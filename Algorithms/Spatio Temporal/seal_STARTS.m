function [Source_Estimate, Parameters] = seal_STARTS(Data, L, sourceSpaceInfo, varargin)
%SEAL_STARTS Solves the ESI problem using the STARTS framework.
%   [Source_Estimate, Parameters] = seal_STARTS(Data, L, sourceSpaceInfo, 'ParameterName', ParameterValue, ...)
%
%   This function implements the STARTS algorithm using an alternating
%   optimization scheme that aligns its structure with SBL algorithms like
%   Block Champagne. It models source activity as a product of spatial
%   weights (W) and temporal basis functions (Phi) and learns all parameters
%   within a Variational Bayes framework.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data (B).
%       L (double matrix): Nchannels x Nsources lead field matrix.
%       sourceSpaceInfo (struct): Struct with source space geometry. Required:
%           .VertConn (sparse matrix): Nsources x Nsources vertex connectivity.
%
%   Optional Name-Value Pair Inputs:
%       'NoiseCovariance' (double matrix): Initial noise covariance matrix. Default: [].
%       'InitialTBFs' (matrix): User-provided initial TBFs (K x T). If empty,
%                               TBFs are learned from data. Default: [].
%       'NumTBFs' (integer): Number of initial TBFs to extract from data
%                            via SVD if 'InitialTBFs' is not provided. Default: 4.
%       'MaxIterations' (integer): Maximum number of iterations. Default: 100.
%       'Tolerance' (double): Convergence tolerance for the cost function. Default: 1e-6.
%       'UpdateNoiseCovariance' (0, 1, 2): If 1 (diagonal) or 2(full structre), iteratively updates the
%                                          noise covariance matrix. Default: 1.
%       'NeighborOrder' (integer): Order of neighbors (k) for spatial blocks. Default: 1.
%       'BlockCorrelation' (double): Assumed correlation within a block. Default: 0.99.
%       'NumOrientations' (integer): Number of orientations per source location. Default: 1.
%       'PruneGammas' (double): Threshold for pruning spatial groups. Default: 1e-4.
%       'PruneAlphas' (double): Threshold for pruning temporal components. Default: 1e-4.
%
%   Outputs:
%       Source_Estimate (double matrix): Nsources x Ntimepoints estimated source activity.
%       Parameters (struct): A struct containing detailed output:
%           .InverseOperator (matrix): The final spatial inverse operator W.
%           .TemporalBasisFunctions (matrix): The learned temporal basis Phi.
%           .LearnedGammas (vector): Final spatial hyperparameters.
%           .LearnedAlphas (vector): Final temporal hyperparameters.
%           .LearnedNoiseCovariance (matrix): Final noise covariance.
%           .CostHistory (vector): Value of the cost function at each iteration.
%
%   Reference: Feng, F., et al. (2025). "STARTS: A self-adapted spatio-temporal
%              framework for automatic E/MEG source imaging". TMI.


%% --- Input Parsing ---
p = inputParser;
p.CaseSensitive = false;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);

addParameter(p, 'NoiseCovariance', [], @(x) ismatrix(x) && issymmetric(x));
addParameter(p, 'InitialTBFs', [], @(x) ismatrix(x));
addParameter(p, 'NumTBFs', 4, @(x) isscalar(x) && x > 0);
addParameter(p, 'MaxIterations', 100, @(x) isscalar(x) && x > 0);
addParameter(p, 'Tolerance', 1e-6, @(x) isscalar(x) && x > 0);
addParameter(p, 'NeighborOrder', 1, @(x) isscalar(x) && x >= 0);
addParameter(p, 'BlockCorrelation', 0.99, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'NumOrientations', 1, @(x) isscalar(x) && ismember(x, [1, 3]));
addParameter(p, 'PruneGammas', 1e-6, @(x) isscalar(x) && x>0);
addParameter(p, 'PruneAlphas', 1e-1, @(x) isscalar(x) && x>0);
addParameter(p, 'UpdateNoiseCovariance', 1, @(x) ismember(x, [0,1,2]));

try
    parse(p, Data, L, sourceSpaceInfo, varargin{:});
catch ME
    disp('Error parsing inputs for seal_STARTS:');
    rethrow(ME);
end

Parameters.OptionsPassed = p.Results;
maxIter = Parameters.OptionsPassed.MaxIterations;
epsilon = Parameters.OptionsPassed.Tolerance;
updateCn = Parameters.OptionsPassed.UpdateNoiseCovariance;
Knb = Parameters.OptionsPassed.NeighborOrder;
Rcorr = Parameters.OptionsPassed.BlockCorrelation;
nd = Parameters.OptionsPassed.NumOrientations;
Cnoise = Parameters.OptionsPassed.NoiseCovariance;
Phi = Parameters.OptionsPassed.InitialTBFs;
numTBFs = Parameters.OptionsPassed.NumTBFs;
prune_gamma_thr = Parameters.OptionsPassed.PruneGammas;
prune_alpha_thr = Parameters.OptionsPassed.PruneAlphas;

%% --- Initialization ---
[Nchannels, Nsources_total] = size(L);
Ntimepoints = size(Data, 2);
nSources_loc = Nsources_total / nd;

if isempty(Cnoise)
    Cnoise = 1e-4 * trace(Data*Data'/Ntimepoints) * eye(Nchannels);
end

if isempty(Phi)
    [~, ~, tD] = svd(Data, 'econ');
    Phi = tD(:, 1:min(numTBFs, size(tD, 2)))';
else
    numTBFs = size(Phi, 1);
end
numW = numTBFs;

MM = sourceSpaceInfo.VertConn;
if Knb > 0
    MM = MM^Knb;
end


LF = L;
M = MM + speye(size(MM));
naSource = nSources_loc;


% Construct Spatial Grouping Matrix G
group_sizes_full = full(sum(M > 0, 2));
G = sparse(naSource, sum(group_sizes_full));
group_indices = cell(naSource, 1);
start_idx = 1;
for i = 1:naSource
    group_indices{i} = start_idx : (start_idx + group_sizes_full(i) - 1);
    start_idx = start_idx + group_sizes_full(i);
    tidx = find(M(i,:) > 0);
    if group_sizes_full(i) == 1
        G(tidx, group_indices{i}) = 1;
    else
        tempBlk = Rcorr*ones(group_sizes_full(i)) + eye(group_sizes_full(i))*(1-Rcorr);
        [U, S_eig] = eig(tempBlk);
        G(tidx, group_indices{i}) = U * sqrt(S_eig);
    end
end
if nd > 1
    G = kron(G, speye(nd));
    group_sizes_full = group_sizes_full * nd;
%     group_indices = cellfun(@(x) reshape(repmat(x,nd,1) + repmat((1:nd)',1,length(x))',1,[]), group_indices, 'UniformOutput', false);
end
% Initialize active sets
active_sources_mask = true(naSource, 1);
group_sizes_active = group_sizes_full;
active_group_cols = 1:sum(group_sizes_full);
% Initialize hyperparameters
Cgamma = ones(naSource, 1); % Spatial hyperparameters
alpha = ones(numW, 1);      % Temporal hyperparameters
C_phi = eye(numW)*1e-6;
W = zeros(sum(group_sizes_full), numW);
%% --- Main Iterative Loop ---
cost = -inf;
costlist = [];
fprintf('Running STARTS algorithm...\n');

for iter = 1:maxIter
    % --- Pruning Step ---
    if iter > 1
        % Prune temporal basis functions based on alpha
        new_tbfs_mask = 1./alpha > (max(1./alpha) * prune_alpha_thr);
        if ~all(new_tbfs_mask) && sum(new_tbfs_mask) > 0
            Phi = Phi(new_tbfs_mask, :);
            alpha = alpha(new_tbfs_mask);
            numW = size(Phi, 1);
            C_phi = C_phi(new_tbfs_mask, new_tbfs_mask);
            W = W(:,new_tbfs_mask);
%             fprintf('  Pruned TBFs. %d remaining.\n', numW);
        end

        % Prune spatial groups based on Cgamma
        new_sources_mask = 1./Cgamma > (max(1./Cgamma) * prune_gamma_thr);
        if ~all(new_sources_mask) && sum(new_sources_mask) > 0
            active_sources_mask = new_sources_mask;
%             fprintf('  Pruned spatial groups. %d remaining.\n', sum(active_sources_mask));
            Cgamma = Cgamma(active_sources_mask);
            group_sizes_active = group_sizes_active(active_sources_mask);
            group_indices = group_indices(active_sources_mask);
            active_group_cols = cell2mat(group_indices');
            if nd > 1
                active_group_cols = nd * (repmat(active_group_cols, nd, 1)-1) + repmat((1:nd)', 1, length(active_group_cols));
                active_group_cols = active_group_cols(:);
            end

        end
    end

    % Rebuild effective leadfield for active sources

    W = W(active_group_cols,:);
    LFG = LF * G(:, active_group_cols);


    % --- Spatial Update Step (for W) ---
    gamma = repelem(Cgamma, group_sizes_active);
    FAF = LFG .* repmat(1./gamma',Nchannels,1) * LFG';
    diagCk = zeros(numW,1);
    tr_w = zeros(Nchannels);
    for k = 1:numW
        tIdx = setdiff(1:numW,k);
        E_PhiPhiT_k = Phi(k,:) * Phi(k,:)' + Ntimepoints*C_phi(k,k);
        C_wk = FAF + Cnoise / E_PhiPhiT_k;
        Residual_k = Data * Phi(k,:).' - LFG * sum(W(:,tIdx).*repmat(Phi(k,:)*Phi(tIdx,:).'+Ntimepoints*C_phi(k,tIdx),size(LFG,2),1),2);               
        W(:,k) =  repmat(1./gamma,1,Nchannels) .* LFG.' / C_wk * Residual_k / E_PhiPhiT_k;
        diagCk(k) = trace(FAF / Cnoise - FAF / C_wk * FAF / Cnoise);
        tr_w = tr_w + (FAF - FAF / C_wk * FAF) * E_PhiPhiT_k;
    end   

    % --- Temporal Update Step (for Phi) ---
    LW = LFG * W;
    CR = LW.' / Cnoise * LW;
    C_phi = inv(CR + diag(alpha + diagCk));
    Phi = C_phi * LW.'/ Cnoise * Data;  
    
    % --- M-Step: Update Spatial Hyperparameters Cgamma ---
    for i = 1:sum(active_sources_mask)
        mu = 0;
        if i == 1
            cols = 1:sum(group_sizes_active(i)) ;
        else
            cols = sum(group_sizes_active(1:i-1)) + 1 : sum(group_sizes_active(1:i));
        end
        for k = 1:numW
            E_PhiPhiT_k = Phi(k,:) * Phi(k,:)' + Ntimepoints*C_phi(k,k);
            C_wk = FAF + Cnoise / E_PhiPhiT_k;
            mu = mu + trace(LFG(:,cols)' / C_wk * LFG(:,cols));
        end
                
        Cgamma(i) = sqrt( mu / trace(W(cols,:)' * W(cols,:)) );
    end
    % --- M-Step: Update Temporal Hyperparameters alpha ---
    E_PhiPhiT_diag = diag(Phi * Phi') + Ntimepoints * diag(C_phi);
    alpha = Ntimepoints ./ E_PhiPhiT_diag;
    
    % --- Noise Update Step ---
    Residual = Data - LW * Phi;
    if updateCn == 1
    E_ResResT_diag = sum(Residual.^2, 2) + Ntimepoints*diag(LW * C_phi * LW') + diag(tr_w);
    Cnoise = diag(max(1e-8, E_ResResT_diag/Ntimepoints));
    else
        if updateCn == 2
            Cnoise = Residual * Residual' / Ntimepoints + LW * C_phi * LW' + tr_w;
        end
    end

    % --- Cost Function & Convergence ---
    gamma = repelem(Cgamma, group_sizes_active);
    FAF = LFG .* repmat(1./gamma',Nchannels,1) * LFG';
    KL_w = 0;
    for k = 1:numW
        E_PhiPhiT_k = Phi(k,:) * Phi(k,:)' + Ntimepoints*C_phi(k,k);
        C_wk = FAF + Cnoise / E_PhiPhiT_k;
        KL_w = KL_w - 0.5* (W(:,k).'*(gamma.*W(:,k)) + chol_logdet(C_wk) + log(E_PhiPhiT_k)*Nchannels);      
    end  
    KL_phi = -0.5 * (trace(C_phi * CR) - sum(log(alpha)) - chol_logdet(C_phi)) * Ntimepoints;
    
    cost_old = cost;
    cost = -0.5 * trace(Residual' / Cnoise * Residual) + KL_w + KL_phi - 0.5 * Ntimepoints * ( (sum(log(diag(Cnoise))) + size(Cnoise,1)*log(2*pi)) );
    costlist = [costlist, cost];
    
    if mod(iter,10) == 0
    fprintf('  Iter %d: Cost ~ %.4f, Change = %.6f\n', iter, cost, abs((cost - cost_old) / cost));
    end
    % ----- Recover W -----
    temp = zeros(sum(group_sizes_full), numW);
    temp(active_group_cols,:) = W;
    W = temp; 
    clear temp
    
    if abs((cost - cost_old) / cost) < epsilon && iter > 1
        fprintf('Convergence reached.\n');
        break;
    end
end


%% --- Final Outputs ---

InverseOperator = G * W;

% % Map back to full source space
% InverseOperator = zeros(naSource*nd, numW);
% active_indices_expanded = [];
% for j=1:nd
%     active_indices_expanded = [active_indices_expanded, (find(active_sources_mask)-1)*nd + j];
% end
% InverseOperator(sort(active_indices_expanded), :) = InverseOperator_active;

Source_Estimate = InverseOperator * Phi;

Parameters.InverseOperator = InverseOperator;
Parameters.TemporalBasisFunctions = Phi;
Parameters.LearnedGammas = Cgamma;
Parameters.LearnedAlphas = alpha;
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