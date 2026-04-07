function [Source_Estimate, Parameters] = seal_uSTAR(Data, L, sourceSpaceInfo, varargin)
%SEAL_USTAR Solves the ESI problem using the u-STAR or STAR framework.
%   [Source_Estimate, Parameters] = seal_uSTAR(Data, L, sourceSpaceInfo, 'ParameterName', ParameterValue, ...)
%
%   This function implements the u-STAR framework. It has two modes:
%   1. STAR Mode (Default): If no microstate maps are provided, it runs the
%      core STAR algorithm, which models source activity as S = W*Phi + E,
%      where W are spatial weights, Phi are temporal basis functions (TBFs),
%      and E is a residual source component. All parameters are learned via
%      a Variational Bayes scheme.
%   2. u-STAR Mode: If microstate maps are provided, the data is first
%      segmented into microstate periods. The STAR algorithm is then applied
%      to each segment independently, allowing the TBFs to adapt to the
%      dynamics of each microstate.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data (B).
%       L (double matrix): Nchannels x Nsources lead field matrix.
%       sourceSpaceInfo (struct): Struct with source space geometry. Required:
%           .VertConn (sparse matrix): Nsources x Nsources adjacent matrix.
%           .VertLoc (optional): if LAURA prior is used.
%
%   Optional Name-Value Pair Inputs:
%       'MicrostateMaps' (matrix): Nchannels x Nmaps matrix of microstate topographies.
%                                  Providing this enables u-STAR mode. Default: [].
%       'MicrostateOptions' (struct): Options for microstate segmentation, passed
%                                     to the 'seal_microsmooth' function.
%                                     e.g., .minTime=5, .polarity=0. Default: struct().
%       'InitialTBFs' (matrix): User-provided initial TBFs (K x T). If empty,
%                               TBFs are learned from data. Default: [].
%       'NumTBFs' (integer): Number of TBFs to extract from data via SVD
%                            if 'InitialTBFs' is not provided. Default: 4.
%       'MaxIterations' (integer): Maximum iterations for the core STAR algorithm. Default: 50.
%       'Tolerance' (double): Convergence tolerance for the cost function. Default: 1e-4.
%       'RunAlgorithm' (string):Algorithm type, 'STAR'/'sSTARTS'. Default:'STAR'.
%       'PriorType' (string): Prior type, 'Laplacian'/'LAURA'. Default:'Laplacian'.
%   Outputs:
%       Source_Estimate (double matrix): Nsources x Ntimepoints estimated source activity.
%       Parameters (struct): A struct containing detailed output:
%           .AlgorithmMode ('STAR' or 'u-STAR')
%           .SegmentResults (cell array): Contains detailed results for each
%                                         segment in u-STAR mode.
%           .FinalW (matrix): Spatial weights (for STAR mode).
%           .FinalPhi (matrix): Temporal basis functions (for STAR mode).
%           .FinalE (matrix): Residual source component (for STAR mode).
%
%   Reference: Feng et al. (2022). "u-STAR: a novel framework for
%              spatio-temporal E/MEG source imaging optimized by microstates".


%% --- Input Parsing ---
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);

addParameter(p, 'MicrostateMaps', [], @ismatrix);
addParameter(p, 'MicrostateOptions', struct('minTime', 5, 'polarity', 0), @isstruct);
addParameter(p, 'InitialTBFs', [], @ismatrix);
addParameter(p, 'NumTBFs', 4, @(x) isscalar(x) && x > 0);
addParameter(p, 'MaxIterations', 50, @(x) isscalar(x) && x > 0);
addParameter(p, 'Tolerance', 1e-6, @(x) isscalar(x) && x > 0);
addParameter(p, 'RunAlgorithm', 'STAR', @(x) isstring(x) || ischar(x));
addParameter(p, 'PriorType', 'Laplacian', @(x) isstring(x) || ischar(x));

try
    parse(p, Data, L, sourceSpaceInfo, varargin{:});
catch ME
    disp('Error parsing inputs for seal_uSTAR:');
    rethrow(ME);
end

Parameters.OptionsPassed = p.Results;
RunAlgorithm = p.Results.RunAlgorithm;
PriorType = p.Results.PriorType;
MSmap = p.Results.MicrostateMaps;
if strcmpi(PriorType,'LAURA') 
    fprintf('Construct local regressive prior.\n');
    VertConn = sourceSpaceInfo.VertConn; VertLoc = sourceSpaceInfo.Vertices;
    Nsources_locations = size(VertConn,1);
    dist_mat = inf(size(VertConn));
    for iter = 1:Nsources_locations
        t_dist = sqrt(sum((repmat(VertLoc(iter,:),Nsources_locations,1) - VertLoc).^2,2));
        dist_mat(VertConn(:,iter),iter) = t_dist(VertConn(:,iter));
    end
    D_min = min(min(dist_mat));
    dist_mat = dist_mat./sqrt(D_min);
    dist_mat = -dist_mat.^(-2);
    dist_mat(isinf(dist_mat)) = 0;
    dist_mat = dist_mat + diag((9 ./ sum(dist_mat<0,2)).*sum(-dist_mat,2));
    M_prior = sparse(dist_mat);
else
    % Construct spatial prior matrix M from VertConn (Graph Laplacian)
    fprintf('Construct Graph Laplacian prior.\n');
    VertConn = sourceSpaceInfo.VertConn;
    M_prior = spdiags(sum(VertConn, 2), 0, size(VertConn,1), size(VertConn,2)) - 0.9 * VertConn;
end

if isempty(MSmap)
    %% --- STAR Mode ---
    Parameters.AlgorithmMode = 'STAR';
    fprintf('No microstate maps provided. Running in STAR mode.\n');
    [~,~,Dic] = svd(Data, 'econ');
    Phi = Dic(:, 1:min(p.Results.NumTBFs, size(Dic,2)))';
    opts = p.Results;
    opts.InitialTBFs = Phi;
    % Call the core STAR algorithm
    if strcmpi(RunAlgorithm, 'sSTARTS')
        [Source_Estimate, star_params] = run_sSTARTS(L, Data, M_prior, opts);
    else
        [Source_Estimate, star_params] = run_STAR(L, Data, M_prior, opts);
        Parameters.FinalE = star_params.E;
    end
        
    Parameters.FinalW = star_params.W;
    Parameters.FinalPhi = star_params.Phi;
    
    Parameters.CostHistory = star_params.cost;
    
else
    %% --- u-STAR Mode ---
    Parameters.AlgorithmMode = 'u-STAR';
    fprintf('Microstate maps provided. Running in u-STAR mode.\n');  
    
    [~, nSamp] = size(Data);
    Source_Estimate = zeros(size(L, 2), nSamp);
    
    % 1. Segment data using microstate maps
    mslabel = seal_microsmooth(Data, MSmap, p.Results.MicrostateOptions);
    
    % 2. Find segment boundaries
    tmark = [true, diff(mslabel) ~= 0];
    segment_boundaries = [find(tmark), nSamp + 1];
    
    Parameters.SegmentResults = cell(length(segment_boundaries) - 1, 1);
    fprintf('Found %d microstate segments. Processing each segment...\n', length(segment_boundaries) - 1);
    
    % 3. Process each segment
    for i_seg = 1:length(segment_boundaries) - 1
        start_idx = segment_boundaries(i_seg);
        end_idx = segment_boundaries(i_seg+1) - 1;
        segment_data = Data(:, start_idx:end_idx);
        
        fprintf('  Processing segment %d (samples %d to %d)...\n', i_seg, start_idx, end_idx);
        
        % Automatically determine TBFs for this segment
        [~,~,Dic] = svd(segment_data, 'econ');
        segment_phi = Dic(:, 1:min(p.Results.NumTBFs, size(Dic,2)))';
        
        segment_opts = p.Results;
        segment_opts.InitialTBFs = segment_phi; % Override with segment-specific TBFs

        if strcmpi(RunAlgorithm, 'sSTARTS')
            [S_seg, seg_params] = run_sSTARTS(L, segment_data, M_prior, segment_opts);
        else
            [S_seg, seg_params] = run_STAR(L, segment_data, M_prior, segment_opts);
        end       
        Source_Estimate(:, start_idx:end_idx) = S_seg;
        Parameters.SegmentResults{i_seg} = seg_params;
    end
    
    % 4. Apply final smoothing filter
    Source_Estimate = fir_smooth(1/3*ones(1,3), Source_Estimate);
    fprintf('Finished processing all segments.\n');
end

end

%% --- Core STAR Algorithm Implementation ---
function [S, par] = run_STAR(L, B, M, opts)
% This is the refactored implementation of STAR.m

[nSensor, nSource] = size(L);
nSamp = size(B, 2);

% Initialize parameters from opts
Phi = opts.InitialTBFs;
maxIter = opts.MaxIterations;
epsilon = opts.Tolerance;

if isempty(Phi)
    [~,~,tD] = svd(B, 'econ');
    Phi = tD(:, 1:min(opts.NumTBFs, size(tD,2)))';
end

C_noise = eye(nSensor);
F = L / M; % L * inv(M)
numW = size(Phi, 1);
W = zeros(nSource, numW);
W_t = W;
C_w = cell(numW, 1);
C_phi = eye(numW) * 1e-6;

lambda = ones(nSource, 1);
gamma = ones(nSource, 1);
alpha = ones(numW, 1);

diagCk = zeros(numW, 1);
E = zeros(nSource, nSamp);
active_dipole_indices = 1:nSource;
prune = [1e-6, 1e-1];
cost = 0;
costlist = [];

for iter = 1:maxIter
    % Pruning steps (simplified for clarity)
    active_tbfs = 1./alpha > (max(1./alpha) * prune(2));
    if ~all(active_tbfs) && any(active_tbfs)
        Phi = Phi(active_tbfs, :);
        W = W(:, active_tbfs); W_t = W_t(:, active_tbfs); C_w = C_w(active_tbfs);
        C_phi = C_phi(active_tbfs, active_tbfs); alpha = alpha(active_tbfs);
        diagCk = diagCk(active_tbfs); numW = size(Phi, 1);
    end
    
    active_sources = 1./gamma(active_dipole_indices) > (max(1./gamma(active_dipole_indices)) * prune(1));
    if ~all(active_sources) && any(active_sources)
        active_dipole_indices = active_dipole_indices(active_sources);
    end
    
    % E-Step: Update Spatial Weights W
    FAF = F .* repmat(1./gamma', nSensor, 1) * F.';
    for k = 1:numW
        tIdx = setdiff(1:numW, k);
        x = (Phi(k,:)*Phi(k,:).') + nSamp * C_phi(k,k);
        C_wk = FAF + C_noise / x;
        r_k = (B - L*E) * Phi(k,:).' - L * sum(W(:,tIdx) .* repmat(Phi(k,:)*Phi(tIdx,:).' + nSamp*C_phi(k,tIdx), nSource, 1), 2);
        W_t(:,k) = repmat(1./gamma, 1, nSensor) .* F.' / C_wk * r_k / x;
        diagCk(k) = trace(FAF/C_noise - FAF/C_wk*FAF/C_noise);
    end
    W = M \ W_t;
    
    % E-Step: Update Temporal Basis Phi
    LW = L * W;
    C_p_inv = (LW.'/C_noise*LW + diag(alpha+diagCk));
    C_phi = pinv(C_p_inv);
    Phi = C_phi * LW.' / C_noise * (B - L*E);
    
    % E-Step: Update Residual E
    C_e = F .* repmat(1./lambda', nSensor, 1) * F.' + C_noise;
    E_t = repmat(1./lambda, 1, nSensor) .* F.' / C_e * (B - L*W*Phi);
    E = M \ E_t;
    
    % M-Step: Update Hyperparameters
    mu = zeros(nSource, 1);
    for k = 1:numW
        x = (Phi(k,:)*Phi(k,:).') + nSamp * C_phi(k,k);
        C_wk = FAF + C_noise / x;
        mu = mu + sum((F.'/C_wk).*F.', 2);
    end
    gamma(active_dipole_indices) = sqrt(mu(active_dipole_indices) ./ sum(W_t(active_dipole_indices,:).^2, 2));
    gamma(gamma > 1e16) = 1e16;
    
    alpha = nSamp ./ (diag(C_phi) + diag(Phi*Phi.'));
    
    beta_val = sum(F.'/C_e.*F.', 2);
    lambda(active_dipole_indices) = sqrt(nSamp*beta_val(active_dipole_indices) ./ sum(E_t(active_dipole_indices,:).^2, 2));
    lambda(lambda > min(lambda)*1e12) = min(lambda)*1e12;
    
    res = sum((B - L*W*Phi - L*E).^2, 2) + nSamp*diag(LW*C_phi*LW.');
    C_noise = diag(max(1e-6, res/nSamp));
    
    % Cost Function (simplified for convergence check)
    tempk = 0;
    temp = 0;
    FAF = F.*repmat(1./gamma',nSensor,1)* F.';
    for k = 1:numW
        x = (Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k);
        C_wk = FAF+C_noise/x;
        
        temp = temp + W(:,k).'.*gamma.'*W(:,k);
        tempk = tempk + chol_logdet(C_wk) + log(x)*nSensor;
    end
    cost_old = cost;
    
    cost = -trace((B-L*W*Phi-L*E)'/C_noise*(B-L*W*Phi-L*E))+nSamp*trace(LW.'/C_noise*LW*C_phi)+nSamp*(sum(log(alpha))+log(det(C_phi)))-nSamp*(log(det(C_e)))...
        -temp-trace(E_t.'.*repmat(lambda',nSamp,1)*E_t)-tempk-nSamp*(chol_logdet(C_noise*2*pi));
    
    costlist = [costlist, cost];
    
    if abs((cost - cost_old)/cost) < epsilon && iter > 1
        break;
    end
end

S = E + W*Phi;
par.W = W;
par.Phi = Phi;
par.E = E;
par.cost = costlist;
end

function [S, par] = run_sSTARTS(L, B, M, opts)
% This is the implementation of sSTARTS

% Initialize
[nSensor, nSource] = size(L);
nSamp = size(B, 2);
% Initialize parameters from opts
Phi = opts.InitialTBFs;
maxIter = opts.MaxIterations;
epsilon = opts.Tolerance;
C_noise =eye(nSensor); % noise covariance

if isfield(opts, 'LearnNoise')
    if opts.LearnNoise == 1
        LF = [L eye(nSensor)];
        M = blkdiag(M,eye(nSensor));
        F = LF/M;
        naSource = nSource+nSensor;
    end
else
    LF = L;
    F = LF/M;
    naSource = nSource;
end

numW = size(Phi,1); % number of TBFs
W = ones(naSource,numW);
W_t = W;
C_phi = eye(numW)*1e-6;

gamma = (ones(naSource,1)+1e-3*randn(naSource,1))*trace(F*F.')*trace(Phi*Phi.')/(trace(B*B'));
alpha = ones(numW,1);
diagCk = zeros(numW,1);

dipIdx = 1:naSource;
prune = [1e-6,1e-1];
cost = -1;
costlist = [];

for iter = 1:maxIter
    % check W and TBF
    tbfIdx = abs(1./alpha)>(max(abs(1./alpha))*prune(2));
    if ~prod(tbfIdx)
        Phi = Phi(tbfIdx,:);
        W = W(:,tbfIdx);
        W_t = W_t(:,tbfIdx);
        C_phi = C_phi(tbfIdx,tbfIdx);
        alpha = alpha(tbfIdx);
        diagCk = diagCk(tbfIdx);
        numW = sum(tbfIdx);
    end
    wIdx = abs(1./gamma(dipIdx))>(max(abs(1./gamma(dipIdx)))*prune(1));
    if ~prod(wIdx)
        dipIdx = dipIdx(wIdx);
    end
    % spatial coefficient
    FAF = F.*repmat(1./gamma',nSensor,1)* F.';
    for k = 1:numW
        tIdx = setdiff(1:numW,k);
        x = ((Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k));
        C_wk = FAF+C_noise/x;
        r_k = B*Phi(k,:).'-LF*sum(W(:,tIdx).*repmat(Phi(k,:)*Phi(tIdx,:).'+nSamp*C_phi(k,tIdx),naSource,1),2);
        W_t(:,k) = repmat(1./gamma,1,nSensor).*F.'/C_wk*r_k/x;
        diagCk(k) = trace(FAF/C_noise-FAF/C_wk*FAF/C_noise);% trace(L'L Cwk)
    end
    W = M\W_t;
    % Phi
    LW = LF*W;
    CR = LW.'/C_noise*LW;
    C_phi = inv(CR+diag(alpha+diagCk));
    Phi = C_phi*LW.'/C_noise*(B);
    
    % gamma
    mu = zeros(naSource,1);
    for k = 1:numW
        x = (Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k);
        C_wk = FAF+C_noise/x; % re-calculate C_w{k}
        mu = mu+sum((F.'/C_wk).*F.',2);
    end
    gamma(dipIdx) = sqrt(mu(dipIdx)./sum(W_t(dipIdx,:).^2,2));
    gamma(gamma>1e12) = 1e12;
    % alpha
    alpha = 1./(diag(C_phi)+diag(Phi*Phi')/nSamp);
    % noise
    res = sum((B-LW*Phi).^2,2)+diag(LW*C_phi*LW'*nSamp);
    C_noise = diag(max(1e-6, res/nSamp));
    
    % free energy
    tempk = 0;
    temp = 0;
    tempbw = 0;
    FAF = F.*repmat(1./gamma',nSensor,1)* F.';
    for k = 1:numW
        x = (Phi(k,:)*Phi(k,:).') + nSamp*C_phi(k,k);
        C_wk = FAF + C_noise/x;
        temp = temp + W_t(:,k).'.*gamma.'*W_t(:,k);
        tempk = tempk + chol_logdet(C_wk) + log(x)*nSensor;
        tempbw = tempbw + FAF - FAF/C_wk * FAF;
    end
    cost_old = cost;
    
    cost = -trace((B-LW*Phi)'/C_noise*(B-LW*Phi))-temp-tempk-nSamp*trace(CR*C_phi)...
        +nSamp*(sum(log(alpha))+log(det((C_phi))))-nSamp*chol_logdet((C_noise*2*pi));
    costlist = [costlist,cost];
    if abs((cost-cost_old)/cost)<epsilon
        break;
    end
    
end
par.W = W;
par.Phi = Phi;
par.dipIdx = dipIdx;
par.alpha = alpha;
par.gamma = gamma;
par.noise = C_noise;
par.cost = costlist;
S = W*Phi;
end



function S_out = fir_smooth(b, S_in)
% A simple wrapper for zero-phase FIR filtering
if isempty(b) || all(b==1) || isempty(S_in)
    S_out = S_in;
    return;
end
try
    S_out = filtfilt(b, 1, S_in')';
catch
    warning('seal_uSTAR:FilterFailed', 'Smoothing filter failed. Returning unfiltered result.');
    S_out = S_in;
end
end

function mslabel = seal_microsmooth(Data, MSmap, options)
% modified from eeglab pop_micro_segment
[~,N,T] = size(Data); 
mslabel = nan(T,N);
for tr = 1:T
    mslabel(tr,:) = reject_segments(Data(:,:,tr), MSmap, options);
end

end

function [logDet, err] = chol_logdet(A)
% Calculates the log-determinant of a symmetric positive definite matrix
% using Cholesky factorization to avoid underflow/overflow.
[R, err] = chol(A);
if err ~= 0
    logDet = -inf;
else
    logDet = 2 * sum(log(diag(R)));
end
end

function L = reject_segments(X,A,opts)
% Reject small segments
% If there is a segment of TFs that is smaller than minTime it gets the next
% best label. This starts with segments of length 1 and then iterates up to
% minTime (in samples).

%% Initialisation
[C,N] = size(X);
K = size(A,2);

% Reading settings
if isfield(opts,'minTime')
    minTime = opts.minTime;
else
    minTime = 3;
end
if isfield(opts,'polarity')
    polarity = opts.polarity;
else
    polarity = 0;
end


%% Calculate the Global map dissimilarity for each Microstate
% Normalise EEG and maps (average reference and gfp = 1 for EEG)
X = (X - repmat(mean(X,1), C, 1)) ./ repmat(std(X,1), C, 1);
A = (A - repmat(mean(A,1), C, 1)) ./ repmat(std(A,1), C, 1);

% Global map dissilarity
GMD = nan(K,N);
for k = 1:K
    GMD(k,:) = sqrt(mean( (X - repmat(A(:,k),1,N)).^2 ));
end

% Account for polarity (recommended 0 for spontaneous EEG)
if polarity == 0
    GMDinvpol = nan(K,N);
    for k = 1:K
        GMDinvpol(k,:) = sqrt(mean( (X - repmat(-A(:,k),1,size(X,2))).^2));
    end
    idx = GMDinvpol < GMD;
    GMD(idx) = GMDinvpol(idx);
end

% Sort the GMD to get the labels
[~ ,labels] = sort(GMD,1);


%% reject small maps
for k = 1:minTime
    cruns = k;
    while sum(cruns <= k) > 0
        idx = [];
        [~, runs] = my_RLE(labels(1,:));
        idx(cumsum([1 runs(:,runs>0)])) = 1;
        cruns = runs(:, cumsum(idx(1:find(idx,1,'last')-1)));
        labels(:,cruns<=k) = circshift(labels(:,cruns<=k),-1);
    end
end


%% Export best labels only
L = labels(1,:);
end

function [d,c]=my_RLE(x)
%% RLE
% This function performs Run Length Encoding to a strem of data x.
% [d,c]=rl_enc(x) returns the element values in d and their number of
% apperance in c. All number formats are accepted for the elements of x.
% This function is built by Abdulrahman Ikram Siddiq in Oct-1st-2011 5:15pm.
if nargin~=1
    error('A single 1-D stream must be used as an input')
end

ind=1;
d(ind)=x(1);
c(ind)=1;

for i=2 :length(x)
    if x(i-1)==x(i)
        c(ind)=c(ind)+1;
    else ind=ind+1;
        d(ind)=x(i);
        c(ind)=1;
    end
end
end