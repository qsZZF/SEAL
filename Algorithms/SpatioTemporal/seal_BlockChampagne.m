function [Source_Estimate, Parameters] = seal_BlockChampagne(Data, L, sourceSpaceInfo, varargin)
%SEAL_BLOCKCHAMPAGNE Solves the ESI problem using the Block Champagne SBL algorithm.
%   [Source_Estimate, Parameters] = seal_BlockChampagne(Data, L, sourceSpaceInfo, 'ParameterName', ParameterValue, ...)
%
%   This function implements a block-sparse Bayesian learning algorithm known
%   as Champagne. It promotes sparse and spatially smooth source estimates by
%   grouping neighboring sources and estimating a common hyperparameter for each group.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data (B).
%       L (double matrix): Nchannels x Nsources lead field matrix.
%       sourceSpaceInfo (struct): Struct with source space geometry. Required:
%           .VertConn (sparse matrix): Nsources x Nsources vertex connectivity.
%
%   Optional Name-Value Pair Inputs:
%       'NoiseCovariance' (double matrix): Initial noise covariance matrix. If not
%                                          provided, it's initialized from data.
%                                          Default: [].
%       'NoiseObservationMatrix' (double matrix): A [Nchannels x Nnoise] matrix to
%                                          be used as the leadfield for noise sources
%                                          when 'LearnNoise' is true.
%                                          Default: eye(Nchannels).
%       'AtlasConstraints' (cell array): A cell array where each cell contains a
%                                        vector of source location indices belonging
%                                        to an atlas parcel. The algorithm will
%                                        learn a full symmetric covariance for
%                                        sources within each parcel.
%                                        Example: {[1:10], [11, 15, 22]}. Default: {}.
%       'MaxIterations' (integer): Maximum number of iterations. Default: 100.
%       'Tolerance' (double): Convergence tolerance for the cost function.
%                             Default: 1e-8.
%       'LearnNoise' (logical): If true, augments the leadfield to learn a
%                               noise component. Default: true.
%       'UpdateNoiseCovariance' (logical): If true, iteratively updates the
%                                          noise covariance matrix. Default: true.
%       'NeighborOrder' (integer): Order of neighbors (k) to include in blocks
%                                  (VertConn^k). Default: 1.
%       'BlockCorrelation' (double): Assumed correlation between sources within
%                                    a block. Default: 0.99.
%       'NumOrientations' (integer): Number of orientations per source location.
%                                    Default: 1.
%
%   Outputs:
%       Source_Estimate (double matrix): Nsources x Ntimepoints estimated source activity.
%       Parameters (struct): A struct containing detailed output:
%           .InverseOperator (double matrix): The final operator mapping data to sources.
%           .LearnedGammas (vector): The final hyperparameters for each source group.
%           .LearnedAtlasCovariances (cell): Final covariance matrices for each atlas region.
%           .LearnedNoiseCovariance (matrix): The final estimated noise covariance.
%           .CostHistory (vector): The value of the cost function at each iteration.
%           .NumIterationsPerformed (integer): Number of iterations completed.
%           .OptionsPassed (struct): The user-provided and default options.
%
%   Reference:
%      


%% --- Input Parsing ---
p = inputParser;
p.CaseSensitive = false;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);

addParameter(p, 'NoiseCovariance', [], @(x) ismatrix(x));
addParameter(p, 'NoiseObservationMatrix', [], @(x) ismatrix(x));
addParameter(p, 'AtlasConstraints', {}, @iscell);
addParameter(p, 'MaxIterations', 50, @(x) isscalar(x) && x > 0);
addParameter(p, 'Tolerance', 1e-8, @(x) isscalar(x) && x > 0);
addParameter(p, 'LearnNoise', true, @islogical);
addParameter(p, 'UpdateNoiseCovariance', true, @islogical);
addParameter(p, 'NeighborOrder', 1, @(x) isscalar(x) && x >= 0);
addParameter(p, 'BlockCorrelation', 0.99, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'NumOrientations', 1, @(x) isscalar(x) && ismember(x, [1, 3]));

try
    parse(p, Data, L, sourceSpaceInfo, varargin{:});
catch ME
    disp('Error parsing inputs for seal_BlockChampagne:');
    rethrow(ME);
end

Parameters.OptionsPassed = p.Results;
maxIter = Parameters.OptionsPassed.MaxIterations;
epsilon = Parameters.OptionsPassed.Tolerance;
learnNoise = Parameters.OptionsPassed.LearnNoise;
updateCn = Parameters.OptionsPassed.UpdateNoiseCovariance;
Knb = Parameters.OptionsPassed.NeighborOrder;
Rcorr = Parameters.OptionsPassed.BlockCorrelation;
nd = Parameters.OptionsPassed.NumOrientations;
Cnoise = Parameters.OptionsPassed.NoiseCovariance;
noiseObsMatrix = Parameters.OptionsPassed.NoiseObservationMatrix;
atlasConstraints = Parameters.OptionsPassed.AtlasConstraints;
if nd > 1 && learnNoise
    error('seal_BlockChampagne:IncompatibleParameters', ...
          'Simultaneously using multiple orientations (NumOrientations > 1) and learning noise (LearnNoise=true) is not supported. Please set one of them to its default value.');
end
%% --- Initialization ---
[nSensor, nSources_total] = size(L);
nSamp = size(Data, 2);
nSources_loc = nSources_total / nd;

% Initialize noise covariance if not provided
if isempty(Cnoise)
    Cnoise = 1e-4 * trace(Data*Data'/nSamp) * eye(nSensor);
end
Covy = Data * Data' / nSamp;

MM = sourceSpaceInfo.VertConn;
if Knb > 0
    MM = MM^Knb;
end

% --- Augment Leadfield to Learn Noise ---
if learnNoise
    if ~isempty(noiseObsMatrix)
        if size(noiseObsMatrix, 1) ~= nSensor
            error('seal_BlockChampagne:DimMismatch', 'NoiseObservationMatrix must have the same number of rows (sensors) as the leadfield L.');
        end
        noiseLeadfield = noiseObsMatrix;
    else
        noiseLeadfield = eye(nSensor);
    end
    
    nNoiseSources = size(noiseLeadfield, 2);
    LF = [L, noiseLeadfield];
    M = blkdiag(MM, zeros(nNoiseSources)); % Zeros matrix for noise part
    M = M + speye(size(M));
    naSource = nSources_loc + nNoiseSources;
else
    LF = L;
    M = MM + speye(size(MM));
    naSource = nSources_loc;
end

% --- Construct Block/Group Matrix G ---
group_sizes = full(sum(M > 0, 2));
G = sparse(naSource, sum(group_sizes));
group_indices = cell(naSource, 1);

start_idx = 1;
for i = 1:naSource
    group_indices{i} = start_idx : (start_idx + group_sizes(i) - 1);
    start_idx = start_idx + group_sizes(i);
    
    tidx = sort(unique([i, find(M(i,:)>0)]));
    
    if group_sizes(i) == 1
        G(tidx, group_indices{i}) = 1;
    else
        tempBlk = Rcorr*ones(group_sizes(i)) + eye(group_sizes(i))*(1-Rcorr);
        [U, S_eig] = eig(tempBlk);
        G(tidx, group_indices{i}) = U * sqrt(S_eig);            
    end
end
if nd > 1
   G = kron(G, speye(nd));
end

LFG = LF * G;
nGroups_total = size(LFG, 2);

% --- Initialize hyperparameters (gammas) ---
initial_gamma_val = trace(Covy) / trace(LF*LF.');
Cgamma = (ones(naSource,1) + 1e-1*randn(naSource,1)) * initial_gamma_val;
gamma = [];
for i = 1:naSource
    gamma = [gamma; Cgamma(i) * ones(group_sizes(i)*nd, 1)];
end

% --- Initialize Atlas Structures ---
isAtlas = ~isempty(atlasConstraints);
if isAtlas
    fprintf('  Initializing atlas constraints...\n');
    nAtlas = numel(atlasConstraints);
    AtlasSet = cell(nAtlas, 1);
    for i_atlas = 1:nAtlas
        curVert = atlasConstraints{i_atlas};
        if any(curVert > nSources_loc) || any(curVert < 1)
             error('Atlas constraint %d contains vertex indices out of bounds.', i_atlas);
        end

        als_cls = [];
        tgamma = [];
        for idp = 1:numel(curVert)
            vertex_idx = curVert(idp);
            als_cls = union(als_cls, group_indices{vertex_idx});
            tgamma = [tgamma; Cgamma(vertex_idx) * ones(group_sizes(vertex_idx), 1)];
        end
        
        als.group_indices_unexpanded = als_cls;
        als_cls_expanded = [];
        for j = 1:nd
            als_cls_expanded = [als_cls_expanded, (als_cls-1)*nd + j];
        end
        als.group_indices_expanded = sort(als_cls_expanded);

        gamma_matrix_unexpanded = spdiags(tgamma, 0, numel(tgamma), numel(tgamma));
        als.gamma_matrix = kron(gamma_matrix_unexpanded, speye(nd));
        AtlasSet{i_atlas} = als;
    end
end


%% --- Main Iterative Loop ---
cost = -inf;
costlist = [];
fprintf('Running Block Champagne algorithm...\n');
Sigma_y = LFG .* repmat(gamma',nSensor,1) * LFG'+ Cnoise;
for iter = 1:maxIter
    % --- E-Step: Calculate Posterior Covariance and Mean ---    
    
    [Ug,Sg] = svd(Sigma_y);
    Sg = max(real(diag(Sg)),0);
    invSg = zeros(size(Sg));
    ff = find(Sg>1e-8 * max(Sg));
    invSg(ff) = 1 ./ Sg(ff);
    Sigma_y_inv = Ug * diag(invSg) * Ug';
    
    W = repmat(gamma,1,nSensor) .* LFG'* Sigma_y_inv;

    % --- M-Step: Update Hyperparameters (Gammas) ---
    for i = 1:naSource
        group_idx_expanded = [];
        if nd > 1 && i <=nSource
            for j = 1:nd
                group_idx_expanded = [group_idx_expanded, (group_indices{i}-1)*nd + j];
            end
        else
            group_idx_expanded = group_indices{i};
        end
       
        ss = trace(W(group_idx_expanded,:) * Covy * W(group_idx_expanded,:)');
        z = trace(LFG(:,group_idx_expanded)' * Sigma_y_inv * LFG(:,group_idx_expanded));
        Cgamma(i) = sqrt(ss) / max(sqrt(z),1e-12);
%         Cgamma(i) = (sqrt(z) / max(z, 1e-8)) * sqrt(ss);
    end
    
    % Clamp gamma values to avoid numerical instability
    Cgamma(Cgamma > 1e16) = 1e16;
    Cgamma(Cgamma < 1e-12) = 1e-12;
    gamma = [];
    for i = 1:naSource
        gamma = [gamma; Cgamma(i) * ones(group_sizes(i)*nd, 1)];
    end
    
    % --- M-Step: Update Atlas-based Symmetric Covariance ---
    Sigma_y = LFG .* repmat(gamma',nSensor,1) * LFG';
    if isAtlas
        
        for i_atlas = 1:nAtlas
            als = AtlasSet{i_atlas};
            als_cls = als.group_indices_expanded;
            atlas_gamma_prev = als.gamma_matrix;
            
            W_als = atlas_gamma_prev * LFG(:, als_cls)' * Sigma_y_inv;
            ss = W_als * Covy * W_als';
            z = LFG(:, als_cls)' * Sigma_y_inv * LFG(:, als_cls);
            temp = atlas_gamma_prev - atlas_gamma_prev * z * atlas_gamma_prev;
            cov_symmetric = ss + temp;
            
            gamma_diag_part = spdiags(gamma(als_cls), 0, numel(als_cls), numel(als_cls));
            Sigma_y = Sigma_y - LFG(:, als_cls) * gamma_diag_part * LFG(:, als_cls)';
            Sigma_y = Sigma_y + LFG(:, als_cls) * cov_symmetric * LFG(:, als_cls)';
            
            als.gamma_matrix = cov_symmetric;
            W(als_cls, :) = cov_symmetric * LFG(:, als_cls)' * Sigma_y_inv;
            AtlasSet{i_atlas} = als;
        end
    end
    
    % --- M-Step: Update Noise Covariance ---
    if updateCn
        fw = eye(nSensor) - LFG*W;
        lam = sqrt(sum(fw*Covy.*fw, 2) ./ diag(Sigma_y_inv));
        Cnoise = diag(max(real(lam), 1e-8));
    end
    
    % --- Cost Function (Marginal Log-Likelihood) ---
    cost_old = cost;
    Sigma_y = Sigma_y + Cnoise;
    [logDet_Sigma_y, chol_err] = chol_logdet(Sigma_y);
    
    if chol_err
        cost = -inf;
    else
        cost = -0.5 * trace(Covy * Sigma_y_inv) - 0.5 * logDet_Sigma_y - (nSensor/2)*log(2*pi);
    end
    costlist = [costlist, cost];
    
    if mod(iter,10) == 0
    fprintf('  Iter %d: Cost = %.4f, Change = %.6f\n', iter, cost, abs((cost - cost_old) / cost));
    end
    
    if abs((cost - cost_old) / cost) < epsilon && iter > 1
        fprintf('Convergence reached.\n');
        break;
    end
end



%% --- Final Outputs ---
% Recalculate the final inverse operator with the final parameters
Gamma_final = spdiags(gamma, 0, nGroups_total, nGroups_total);
if isAtlas
    for i_atlas = 1:nAtlas
        als = AtlasSet{i_atlas};
        Gamma_final(als.group_indices_expanded, als.group_indices_expanded) = als.gamma_matrix;
    end
end

Sigma_y = Cnoise + LFG * Gamma_final * LFG';
W_final = Gamma_final * LFG' / Sigma_y;

% Transform back from group space to source space
InverseOperator = G * W_final;

% Calculate final source estimate
Source_Estimate = InverseOperator * Data;
if learnNoise
    Parameters.NoiseWeights = InverseOperator(nSources_total + 1:end,:);
    Source_Estimate = Source_Estimate(1:nSources_total,:);
    InverseOperator = InverseOperator(1:nSources_total,:);
end

% Populate output parameter struct
Parameters.InverseOperator = InverseOperator;
Parameters.LearnedGammas = Cgamma;
if isAtlas
    Parameters.LearnedAtlasCovariances = cellfun(@(c) c.gamma_matrix, AtlasSet, 'UniformOutput', false);
end
Parameters.LearnedNoiseCovariance = Cnoise;
Parameters.CostHistory = costlist;
Parameters.NumIterationsPerformed = iter;

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

