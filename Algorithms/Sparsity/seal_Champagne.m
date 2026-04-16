function [Source_SBL, Param_SBL] = seal_Champagne(Data, L, varargin)
%SEAL_CHAMPAGNE Computes Sparse Bayesian Learning (Champagne) source imaging.
%   [Source, Param] = seal_Champagne(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Solves the ESI inverse problem using the Champagne algorithm, iteratively 
%   updating the source covariance (gamma) and noise covariance to maximize 
%   the Bayesian model evidence. Highly robust to correlated sources and 
%   effectively eliminates spatial leakage.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data.
%       L (double matrix): Nchannels x Nsources matrix (lead field).
%
%   Optional Name-Value Pair Inputs:
%       'NoiseCovariance' (double matrix): Initial noise covariance. 
%                                          If empty, initialized as scaled identity.
%       'UpdateNoise' (logical): If true, actively learns heteroscedastic noise 
%                                covariance from the data. Default: true.
%       'NumOrientations' (integer): Number of orientations per source. Default: 1.
%       'MaxIterations' (integer): Maximum EM iterations. Default: 100.
%       'Tolerance' (double): Convergence tolerance for the cost function. Default: 1e-8.
%       'PruneThreshold' (double): Threshold to prune inactive dipoles and speed up 
%                                  computation. Default: 1e-6.
%
%   Outputs:
%       Source_SBL (double matrix): Nsources x Ntimepoints estimated source activity.
%       Param_SBL (struct):
%           .InverseOperator (double matrix): The Nsources x Nchannels spatial filter.
%           .Gamma (double vector): Final learned source variances.
%           .NoiseCovariance (double matrix): Final learned noise covariance.
%           .CostFunction (vector): Log-likelihood history per iteration.
%           .ActiveDipoles (vector): Indices of sources that survived pruning.

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;

    defaultNoiseCov = [];
    defaultUpdateNoise = true;
    defaultNumOrientations = 1;
    defaultMaxIter = 100;
    defaultTol = 1e-4;
    defaultPrune = 1e-8;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) isempty(x) || ismatrix(x));
    addParameter(p, 'UpdateNoise', defaultUpdateNoise, @islogical);
    addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'MaxIterations', defaultMaxIter, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'Tolerance', defaultTol, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'PruneThreshold', defaultPrune, @(x) isnumeric(x) && isscalar(x));

    parse(p, Data, L, varargin{:});

    Param_SBL.OptionsPassed = p.Results;
    Data_orig = p.Results.Data;
    L_orig = p.Results.L;
    C_noise = p.Results.NoiseCovariance;
    updateCn = p.Results.UpdateNoise;
    nd = p.Results.NumOrientations;
    maxIter = p.Results.MaxIterations;
    epsilon = p.Results.Tolerance;
    prune_thresh = p.Results.PruneThreshold;

    [Nchannels, Ntimepoints] = size(Data_orig);
    [Nchannels_L, Nsources_total] = size(L_orig);

    if Nchannels ~= Nchannels_L
        error('seal_Champagne:DimensionMismatch', 'Data and L must have same number of channels.');
    end

    Nsources_locations = Nsources_total / nd;

    %% 2. Initialization
    % Data Covariance
    C_y = (Data_orig * Data_orig') / Ntimepoints;

    % Initialize Noise Covariance
    if isempty(C_noise)
        C_noise = eye(Nchannels);
    end

    % Initialize Source Variance (Gamma)
    gamma_loc = (ones(Nsources_locations, 1) + 1e-3 * randn(Nsources_locations, 1)) * ...
                (trace(C_y) / trace(L_orig * L_orig'));

    % Tracking variables
    active_loc_idx = 1:Nsources_locations;
    cost_list = zeros(maxIter, 1);
    cost_old = -inf;

    %% 3. Iterative EM / Fixed-Point Updates
    for iter = 1:maxIter
        
        % A. Pruning Inactive Dipoles
        % Removes dipoles whose variance has collapsed towards zero
        keep_mask = abs(gamma_loc) > (max(abs(gamma_loc)) * prune_thresh);
        
        if sum(~keep_mask) > 0
            active_loc_idx = active_loc_idx(keep_mask);
            gamma_loc = gamma_loc(keep_mask);
        end
        
        % Expand location-based gamma to full orientation space
        if nd > 1
            gamma_full = reshape(repmat(gamma_loc', nd, 1), [], 1);
            active_full_idx = reshape(bsxfun(@plus, (active_loc_idx-1)*nd, (1:nd)'), [], 1);
        else
            gamma_full = gamma_loc;
            active_full_idx = active_loc_idx;
        end
        
        L_active = L_orig(:, active_full_idx);

        % B. Model Covariance & Inversion
        % Gram matrix: C = L * Gamma * L^T + Sigma_noise
        Gamma_mat = bsxfun(@times, L_active, gamma_full');
        Gram = Gamma_mat * L_active' + C_noise;
        Gram = (Gram + Gram') / 2; % Ensure symmetry
        
        % Robust SVD inversion
        [Ug, Sg_mat] = svd(Gram);
        Sg = max(real(diag(Sg_mat)), 0);
        invSg = zeros(Nchannels, 1);
        invSg(Sg > 1e-12) = 1 ./ Sg(Sg > 1e-12);
        invGram = Ug * diag(invSg) * Ug';

        % C. Compute Spatial Filter (W)
        % W = Gamma * L^T * C^-1
        W_active = Gamma_mat' * invGram;

        % D. Update Source Variances
        % ss = diag(W * C_y * W^T) -> variance of estimated source
        % z = diag(L^T * C^-1 * L) -> sensitivity
        ss = sum((W_active * C_y) .* W_active, 2); 
        z = sum((L_active' * invGram) .* L_active', 2);
        
        % Compress back to location space if nd > 1
        if nd > 1
            ss_loc = sum(reshape(ss, nd, []), 1)';
            z_loc = sum(reshape(z, nd, []), 1)';
        else
            ss_loc = ss;
            z_loc = z;
        end
        
        % Gamma update
        gamma_loc = (sqrt(z_loc) ./ max(z_loc, 1e-12)) .* sqrt(ss_loc);
        
        % Clamp values for numerical stability
        gamma_loc(gamma_loc > 1e16) = 1e16;
        gamma_loc(gamma_loc < 1e-16) = 1e-16;

        % E. Update Noise Covariance
        if updateCn
            fw = eye(Nchannels) - L_active * W_active;
            % Heteroscedastic noise update based on model fit
            lam = sqrt(sum((fw * C_y) .* fw, 2) ./ diag(invGram));
            C_noise = diag(max(real(lam), 1e-12));
        end

        % F. Cost Function (Log-Likelihood) evaluation
        cost = -0.5 * sum(sum(C_y .* invGram)) - 0.5 * (sum(log(max(Sg, 1e-12))) + Nchannels * log(2*pi));
        cost_list(iter) = cost;

        % Convergence Check
        if iter > 1 && abs((cost - cost_old) / cost) < epsilon
            cost_list = cost_list(1:iter);
            break;
        end
        cost_old = cost;
    end

    %% 4. Final Operator Assembly
    % Expand the final pruned filter back to full source space dimensions
    W_final = zeros(Nsources_total, Nchannels);
    W_final(active_full_idx, :) = W_active;

    % Reconstruct the full Gamma vector
    Gamma_final = zeros(Nsources_locations, 1);
    Gamma_final(active_loc_idx) = gamma_loc;

    %% 5. Output Assignment
    Param_SBL.InverseOperator = W_final;
    Param_SBL.Gamma = Gamma_final;
    Param_SBL.NoiseCovariance = C_noise;
    Param_SBL.CostFunction = cost_list;
    Param_SBL.ActiveDipoles = active_loc_idx;

    Source_SBL = W_final * Data_orig;

end