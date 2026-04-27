function [Source_MxNE, Param_MxNE] = seal_MxNE(Data, L, varargin)
%SEAL_MXNE Computes Mixed-Norm Estimate (MxNE) for ESI.
%   [Source_MxNE, Param_MxNE] = seal_MxNE(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Implements the L1/L2 Mixed-Norm Estimate (group Lasso) using Block
%   Coordinate Descent and optional FISTA-based debiasing.
%   Inspired by MNE-Python's mixed_norm implementation.
%   This version focuses on n_mxne_iter=1 (standard MxNE).
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data (Y).
%       L (double matrix): Nchannels x Nsources original lead field matrix (G).
%
%   Optional Name-Value Pair Inputs:
%       'Alpha' (double): Regularization parameter (0-100). Higher values promote more sparsity.
%                         Default: 50.
%       'NumOrientations' (integer): Number of orientations per source location. Default: 1.
%       'NoiseCovariance' (double matrix): Noise covariance (C_noise). Default: eye(Nchannels).
%       'DepthBiasComp' (struct): Parameters for depth weighting. Default: struct('exponent', 0.5, 'limit', 10.0).
%       'Weights' (vector/struct): External spatial weights. Default: [].
%       'WeightsMin' (scalar): Minimum threshold for external weights. Default: 0.0.
%       'SolverMaxIter' (integer): Max iterations for BCD solver. Default: 3000.
%       'SolverTolerance' (double): Convergence tolerance for BCD solver. Default: 1e-4.
%       'Debias' (logical): Perform amplitude debiasing. Default: true.
%       'DebiasMaxIter' (integer): Max iterations for FISTA debiasing. Default: 1000.
%       'DebiasTolerance' (double): Tolerance for FISTA debiasing. Default: 1e-7.
%       'PCAWhitening' (logical/integer): Reduce data rank via PCA after whitening. Default: true.
%       'PickOri' (char/string): 'none' or 'vector'. Default: 'none'.
%       'Verbose' (logical): Print progress info. Default: true.
%
%   Outputs:
%       Source_MxNE (double matrix): Estimated source activity.
%       Param_MxNE (struct): Parameters and intermediate results, including:
%           .Converged (logical): True if solver met tolerance.
%           .ConvergedAt (integer): Iteration of convergence (= MaxIter if not).
%           .E_path (vector): Primal objective values during BCD (sampled).

%% Input Parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = false;

defaultAlpha = 50.0;
if size(L,1) > 0, defaultNoiseCov = eye(size(L,1)); else, defaultNoiseCov = []; end
defaultNumOrientations = 1;
defaultDepthBiasComp = struct('exponent', 0.5, 'limit', 10.0);
defaultWeights = [];
defaultWeightsMin = 0.0;
defaultSolverMaxIter = 500;
defaultSolverTol = 1e-4;
defaultDebias = true;
defaultDebiasMaxIter = 1000;
defaultDebiasTol = 1e-7;
defaultPCAWhitening = true;
defaultPickOri = 'none';

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

addParameter(p, 'Alpha', defaultAlpha, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=100);
addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x) && ismember(x,[1,3]));
addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
addParameter(p, 'DepthBiasComp', defaultDepthBiasComp, @(x) isstruct(x) || isempty(x));
addParameter(p, 'Weights', defaultWeights, @(x) isnumeric(x) || isstruct(x) || isempty(x));
addParameter(p, 'WeightsMin', defaultWeightsMin, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'SolverMaxIter', defaultSolverMaxIter, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'SolverTolerance', defaultSolverTol, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'Debias', defaultDebias, @islogical);
addParameter(p, 'DebiasMaxIter', defaultDebiasMaxIter, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'DebiasTolerance', defaultDebiasTol, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'PCAWhitening', defaultPCAWhitening);
addParameter(p, 'PickOri', defaultPickOri, @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, Data, L, varargin{:});
Param_MxNE.OptionsPassed = p.Results;

Data_orig = p.Results.Data; 
L_orig = p.Results.L; 
alpha_user = p.Results.Alpha;
nd = p.Results.NumOrientations; 
C_noise = p.Results.NoiseCovariance;
depth_params = p.Results.DepthBiasComp; 
user_weights_input = p.Results.Weights;
user_weights_min = p.Results.WeightsMin; 
solver_max_iter = p.Results.SolverMaxIter;
solver_tol = p.Results.SolverTolerance; 
perform_debias = p.Results.Debias;
debias_max_iter = p.Results.DebiasMaxIter; 
debias_tol = p.Results.DebiasTolerance;
pca_whitening_opt = p.Results.PCAWhitening; 
pick_ori_opt = lower(char(p.Results.PickOri));
verbose = logical(p.Results.Verbose);
local_check_option(pick_ori_opt, {'none', 'vector'});

[Nchannels, Ntimepoints] = size(Data_orig); [Nchan_L, Nsources_total] = size(L_orig);
if mod(Nsources_total, nd) ~= 0, error('seal_MxNE:NDMismatch', 'L Nsources not divisible by NumOrientations.'); end
nSources_loc = Nsources_total / nd;

if isempty(C_noise), C_noise = eye(Nchannels); Param_MxNE.OptionsPassed.NoiseCovariance = C_noise; end
if Nchannels ~= Nchan_L, error('seal_MxNE:ChannelMismatch', 'Data and L channel mismatch.'); end

% Constants
ALPHA_MAX_SCALE = 0.01;          % MNE-Python convention: alpha range 0-100

% Initialize convergence tracking (so fields always exist)
Param_MxNE.Converged   = false;
Param_MxNE.ConvergedAt = solver_max_iter;
Param_MxNE.E_path      = [];

%% 1. Prepare Gain Matrix (Whitening, PCA, Depth Comp, External Weights)
iW = eye(Nchannels); Data_w = Data_orig; L_w = L_orig; rank_whitener = Nchannels;
if ~isequal(C_noise, eye(Nchannels))
    try
        [U_c, S_c_diag_mat] = eig((C_noise + C_noise')/2); s_c_vec = diag(S_c_diag_mat);
        [s_c_vec_sorted, sort_idx_c] = sort(s_c_vec, 'descend'); U_c_sorted = U_c(:, sort_idx_c);
        rank_eff = sum(s_c_vec_sorted > (Nchannels * eps(max(s_c_vec_sorted))));
        rank_whitener = min(rank_eff, Nchannels);
        
        s_inv_sqrt = zeros(Nchannels,1);
        valid_s_idx = 1:rank_whitener;
        s_inv_sqrt(valid_s_idx) = 1 ./ sqrt(s_c_vec_sorted(valid_s_idx));
        
        iW_dense = diag(s_inv_sqrt) * U_c_sorted';
        iW = iW_dense(1:rank_whitener, :);
        Data_w = iW * Data_orig; L_w = iW * L_orig;
    catch ME_white 
        warning('seal_MxNE:WhiteningFailed', 'Whitening: %s. Using unwhitened.', ME_white.message); 
        iW = eye(Nchannels); Data_w = Data_orig; L_w = L_orig; rank_whitener = Nchannels; 
    end
end
Param_MxNE.WhiteningMatrix = iW;

M_processed = Data_w; Vh_pca = [];
if pca_whitening_opt
    pca_rank = rank_whitener;
    if isnumeric(pca_whitening_opt) && ~islogical(pca_whitening_opt), pca_rank = min(pca_whitening_opt, rank_whitener); end
    if pca_rank < size(M_processed,1)
        if verbose, fprintf('MxNE: Applying PCA to whitened data, reducing to rank %d\n', pca_rank); end
        [U_data, S_data_mat, V_data] = svd(M_processed, 'econ');
        M_processed = U_data(:, 1:pca_rank) * S_data_mat(1:pca_rank, 1:pca_rank);
        Vh_pca = V_data(:, 1:pca_rank)';
    end
end

source_weighting = ones(Nsources_total, 1);
if ~isempty(depth_params) && isfield(depth_params, 'exponent')
    source_norms = zeros(nSources_loc, 1);
    for i = 1:nSources_loc
        idx = (i-1)*nd + (1:nd); 
        source_norms(i) = norm(L_w(:, idx), 'fro'); 
    end
    source_norms(source_norms < eps) = eps;
    depth_weights_loc = source_norms .^ (-depth_params.exponent);
    if isfield(depth_params,'limit') && depth_params.limit > 0
        limit_val = min(depth_weights_loc) * (depth_params.limit^2);
        depth_weights_loc(depth_weights_loc > limit_val) = limit_val;
    end
    if nd == 1 
        source_weighting = depth_weights_loc; 
    else 
        source_weighting = repelem(depth_weights_loc, nd); 
    end
end
if ~isempty(user_weights_input)
    uw = [];
    if isstruct(user_weights_input) && isfield(user_weights_input, 'data')
        uw = abs(max(user_weights_input.data, [], 2));
    elseif isnumeric(user_weights_input) && isvector(user_weights_input) 
        uw = abs(user_weights_input(:)); 
    end
    if ~isempty(uw) && length(uw) == nSources_loc
        if max(uw) > eps 
            uw = uw / max(uw); 
        end
        if nd == 1 
            source_weighting = source_weighting .* uw; 
        else 
            source_weighting = source_weighting .* repelem(uw, nd); 
        end
    else 
        warning('seal_MxNE:UserWeightsDim', 'User weights ignored.'); 
    end
end
Param_MxNE.SourceWeightingUsedForGain = source_weighting;
G_eff = bsxfun(@times, L_w, source_weighting');

mask_active_sources = true(Nsources_total, 1);
if user_weights_min > 0 && ~isempty(user_weights_input)
    mask_active_sources = source_weighting > (user_weights_min * (max(source_weighting)+eps));
    if ~any(mask_active_sources), error('seal_MxNE:NoSourcesAfterWeightMin', 'No sources after WeightsMin.'); end
    G_eff = G_eff(:, mask_active_sources);
    if verbose
        fprintf('MxNE: Reduced source space to %d component(s) due to WeightsMin.\n', sum(mask_active_sources));
    end
end
% Verify the mask respects block structure (otherwise BCD will misalign)
if mod(size(G_eff, 2), nd) ~= 0
    error('seal_MxNE:MaskBlockMismatch', ...
          'Active source mask does not respect orientation-block structure.');
end
Param_MxNE.ActiveSourceMask_Initial = mask_active_sources;

%% 2. Scaling for Solver
GtM = G_eff' * M_processed;
n_pos_eff = size(GtM,1)/nd;
l2inf_norms_sq_per_group = zeros(n_pos_eff, 1);
for i = 1:n_pos_eff 
    idx = (i-1)*nd + (1:nd); 
    l2inf_norms_sq_per_group(i) = sum(GtM(idx,:).^2, 'all'); 
end
alpha_max_val = sqrt(max(l2inf_norms_sq_per_group)) * ALPHA_MAX_SCALE;
if alpha_max_val < eps
    alpha_max_val = 1.0; 
end

G_solver = G_eff / alpha_max_val;
alpha_solver = alpha_user;
Param_MxNE.AlphaSolver = alpha_solver;
Param_MxNE.AlphaMaxScalingFactor = alpha_max_val;

%% 3. Solve MxNE (L1/L2) using Block Coordinate Descent
if verbose, fprintf('MxNE: Solving L1/L2 problem with BCD (n_mxne_iter=1)...\n'); end
[X_active_pca, active_set_solver_mask, E_path_bcd, bcd_converged, bcd_iter] = ...
    local_mxne_bcd_solver(M_processed, G_solver, alpha_solver, nd, ...
        solver_max_iter, solver_tol, max(1,round(solver_max_iter/10)), verbose);

Param_MxNE.E_path      = E_path_bcd;
Param_MxNE.Converged   = bcd_converged;
Param_MxNE.ConvergedAt = bcd_iter;

if sum(active_set_solver_mask) == 0
    warning('seal_MxNE:NoActiveSourcesFound', 'No active sources found by BCD. Alpha might be too high.');
    Source_MxNE = zeros(Nsources_total, Ntimepoints);
    Param_MxNE.ActiveSetIndices_Solver = []; Param_MxNE.ActiveSetIndices_Locations = [];
    return;
end

%% 4. Debias (if requested)
X_debiased_active_pca = X_active_pca;
if perform_debias && sum(active_set_solver_mask) > 0
    if verbose, fprintf('MxNE: Performing debiasing...\n'); end
    G_for_debias = G_solver(:, active_set_solver_mask);
    debiasing_weights = local_compute_bias(M_processed, G_for_debias, X_active_pca, nd, ...
        debias_max_iter, debias_tol, verbose);
    X_debiased_active_pca = X_active_pca .* debiasing_weights;
    Param_MxNE.DebiasingWeights = debiasing_weights;
end

%% 5. Reconstruct full time series if PCA was used
X_debiased_active_full_time = X_debiased_active_pca;
if ~isempty(Vh_pca)
    X_debiased_active_full_time = X_debiased_active_pca * Vh_pca;
end

%% 6. Reapply scaling and map back to original source space
source_weighting_eff = source_weighting(mask_active_sources);
source_weighting_eff_active = source_weighting_eff(active_set_solver_mask);

X_final_active = bsxfun(@times, X_debiased_active_full_time, source_weighting_eff_active);
X_final_active = X_final_active / alpha_max_val;

%% 7. Format Output STC
X_full_components = zeros(Nsources_total, Ntimepoints);
temp_full_masked = zeros(sum(mask_active_sources), Ntimepoints);
temp_full_masked(active_set_solver_mask, :) = X_final_active;
X_full_components(mask_active_sources, :) = temp_full_masked;

if nd > 1 && strcmpi(pick_ori_opt, 'none')
    Source_MxNE_loc = zeros(nSources_loc, Ntimepoints);
    active_comp_indices_final = find(any(X_full_components,2));
    active_locs_final_unique = unique(ceil(active_comp_indices_final/nd));
    Param_MxNE.ActiveSetIndices_Locations = active_locs_final_unique;
    Param_MxNE.ActiveSetIndices_Components = active_comp_indices_final;
    
    for i_loc = 1:nSources_loc
        if ismember(i_loc, active_locs_final_unique)
            idx_comps = (i_loc-1)*nd + (1:nd);
            Source_MxNE_loc(i_loc, :) = sqrt(sum(X_full_components(idx_comps,:).^2, 1));
        end
    end
    Source_MxNE = Source_MxNE_loc;
else
    Source_MxNE = X_full_components;
    Param_MxNE.ActiveSetIndices_Components = find(any(X_full_components,2));
    if nd > 1, Param_MxNE.ActiveSetIndices_Locations = unique(ceil(Param_MxNE.ActiveSetIndices_Components/nd));
    else, Param_MxNE.ActiveSetIndices_Locations = Param_MxNE.ActiveSetIndices_Components; end
end

end


function [X_solved_active, active_set_mask, E_path, converged, last_iter] = local_mxne_bcd_solver(...
    M_data, G_solver, alpha_penalty_val, O_orient, max_iter, tol, dgap_freq, verbose)
% Solves L1/L2 MxNE using Block Coordinate Descent.

    [~, Nt_proc] = size(M_data);
    [~, Ns_proc] = size(G_solver);
    N_pos_proc = Ns_proc / O_orient;

    X = zeros(Ns_proc, Nt_proc);
    R = M_data;

    % Pre-allocate E_path (max possible length)
    E_path = nan(ceil(max_iter / dgap_freq) + 2, 1);
    n_E = 0;
    converged = false;
    last_iter = max_iter;

    % Lipschitz constants for each block. Use eig() on small SPD matrix
    % instead of norm(.,2) (SVD) for speed; mathematically equivalent.
    lc = zeros(N_pos_proc, 1);
    if O_orient == 1
        lc = sum(G_solver.^2, 1)';
    else
        for j_pos = 1:N_pos_proc
            idx = (j_pos-1)*O_orient + (1:O_orient);
            G_j_block = G_solver(:, idx);
            GtG = G_j_block' * G_j_block;
            GtG = (GtG + GtG') * 0.5;            % symmetrize
            lc(j_pos) = max(eig(GtG));
        end
    end
    lc(lc < eps) = eps;
    alpha_scaled_per_block = alpha_penalty_val ./ lc; 

    if verbose
        fprintf('  BCD Solver: MaxIter=%d, Tol=%.1e, Nsources_proc=%d, Npos_proc=%d\n', ...
                max_iter, tol, Ns_proc, N_pos_proc);
    end

    for iter_bcd = 1:max_iter
        last_iter = iter_bcd;
        X_old_iter = X;
        
        for j_pos = 1:N_pos_proc
            idx_block = (j_pos-1)*O_orient + (1:O_orient);
            G_j = G_solver(:, idx_block);
            X_j_current = X(idx_block, :);
            
            % Gradient step 
            X_j_unthres = X_j_current + (1/lc(j_pos)) * (G_j' * R);
            group_norm_L2F = norm(X_j_unthres, 'fro');
            
            % Soft Thresholding (group-level)
            if group_norm_L2F <= alpha_scaled_per_block(j_pos)
                X_j_new = zeros(O_orient, Nt_proc);
            else
                shrink_factor = 1 - alpha_scaled_per_block(j_pos) / group_norm_L2F;
                X_j_new = X_j_unthres * shrink_factor;
            end
            
            if any(X_j_new(:) ~= X_j_current(:))
                R = R - G_j * (X_j_new - X_j_current);
            end
            
            X(idx_block, :) = X_j_new;
        end
        
        % --- Convergence Check via True Duality Gap ---
        if mod(iter_bcd, dgap_freq) == 0 || iter_bcd == max_iter || iter_bcd == 1
            [actual_gap, p_obj, ~] = local_dgap_l21(M_data, G_solver, X, alpha_penalty_val, O_orient);
            
            n_E = n_E + 1;
            E_path(n_E) = p_obj;
            
            if iter_bcd > 1 
                if verbose
                    fprintf('    BCD Iter %4d: Pobj=%.4e, DGap=%.2e (Tol=%.1e), Nactive_locs=%d\n', ...
                        iter_bcd, p_obj, actual_gap, tol, ...
                        sum(any(reshape(local_compute_active_mask(X, O_orient),O_orient,[]),1)));
                end
                if actual_gap < tol
                    if verbose, fprintf('    BCD converged (Duality Gap).\n'); end
                    converged = true;
                    break;
                end
            end
        end
        
        % Fallback absolute change convergence
        if iter_bcd > 1 && norm(X(:) - X_old_iter(:),'fro') / (norm(X_old_iter(:),'fro')+eps) < tol/100
            if verbose, fprintf('    BCD converged (X change) at iter %d.\n', iter_bcd); end
            converged = true;
            break;
        end
    end
    
    if ~converged
        warning('seal_MxNE:BCDSolverMaxIter', 'BCD solver reached max iterations (%d).', max_iter);
    end
    E_path = E_path(1:n_E);     % trim unused

    % Recompute final active mask from final X (avoids stale flags from
    % blocks that oscillated during BCD).
    active_set_mask = local_compute_active_mask(X, O_orient);
    X_solved_active = X(active_set_mask, :);
end

function active_set_mask = local_compute_active_mask(X, O_orient)
% Determine which entries of X belong to active groups.
    [Ns, ~] = size(X);
    if O_orient == 1
        active_set_mask = any(X ~= 0, 2);
    else
        N_pos = Ns / O_orient;
        % A group is active iff any of its components has any non-zero entry
        active_blocks = false(N_pos, 1);
        for j = 1:N_pos
            idx = (j-1)*O_orient + (1:O_orient);
            active_blocks(j) = any(any(X(idx, :) ~= 0));
        end
        active_set_mask = repelem(active_blocks, O_orient);
    end
end

function [gap, p_obj, d_obj] = local_dgap_l21(M_data, G_full, X_full, alpha_val, nd)
% Computes the rigorous Duality Gap for the L1/L2 mixed-norm.

    GX = G_full * X_full;
    R = M_data - GX;
    nR2 = sum(R(:).^2);

    % --- Primal Objective: 0.5 ||M - GX||_F^2 + alpha * sum_g ||X_g||_F ---
    if nd > 1
        X_tensor = reshape(X_full, nd, [], size(X_full, 2));
        X_norms = sum(sqrt(sum(sum(X_tensor.^2, 1), 3)));
    else
        X_norms = sum(sqrt(sum(X_full.^2, 2)));
    end
    p_obj = 0.5 * nR2 + alpha_val * X_norms;

    % --- Dual Objective ---
    Gt_R = G_full' * R;
    if nd > 1
        Gt_R_tensor = reshape(Gt_R, nd, [], size(Gt_R, 2));
        dual_norms = sqrt(sum(sum(Gt_R_tensor.^2, 1), 3)); 
        dual_norm_l2inf = max(dual_norms);
    else
        dual_norm_l2inf = max(sqrt(sum(Gt_R.^2, 2)));
    end

    scaling = min(alpha_val / max(dual_norm_l2inf, eps), 1.0);
    d_obj = (scaling - 0.5 * (scaling^2)) * nR2 + scaling * sum(sum(R .* GX));
    
    gap = p_obj - d_obj;
    if isnan(gap) && ~isnan(p_obj), gap = p_obj; end
end

function D_debias = local_compute_bias(M_data, G_active_comp, X_active_comp, O_orient, max_iter_debias, tol_debias, verbose)
% FISTA debiasing: solves min_D 0.5*||M - G*diag(D)*X||_F^2 s.t. D>=1 (per group avg).

[N_active_comp, Nt_proc] = size(X_active_comp);
N_pos_active = N_active_comp / O_orient;

Dk_comp = ones(N_active_comp, 1);
Yk_comp = ones(N_active_comp, 1);
tk = 1.0;

L_fista_val = 1.1 * local_power_iteration_kron( ...
        G_active_comp, X_active_comp, 1000, 1e-3, 0);
if verbose
    fprintf('  Debiasing (FISTA): MaxIter=%d, Tol=%.1e, L_fista=%.2e\n', ...
            max_iter_debias, tol_debias, L_fista_val);
end

iter_d = 0;
debias_converged = false;
for iter_d = 1:max_iter_debias
    D0_comp = Dk_comp;
    
    M_pred = G_active_comp * bsxfun(@times, X_active_comp, Yk_comp);
    Residual = M_data - M_pred;
    
    Gt_Res = G_active_comp' * Residual;
    Grad_D_comp = -sum(Gt_Res .* X_active_comp, 2);
    
    D_update = Yk_comp - Grad_D_comp / L_fista_val;
    
    if O_orient > 1
        D_reshaped_loc = mean(reshape(D_update, O_orient, N_pos_active), 1)';
        D_reshaped_loc = max(D_reshaped_loc, 1.0);
        Dk_comp = repelem(D_reshaped_loc, O_orient);
        Dk_comp = Dk_comp(:);
    else
        Dk_comp = max(D_update, 1.0);
    end
    
    t_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * tk^2));
    Yk_comp = Dk_comp + ((tk - 1.0) / t_new) * (Dk_comp - D0_comp);
    tk = t_new;
    
    d_diff = max(abs(Dk_comp - D0_comp));
    if d_diff < tol_debias && iter_d > 1
        if verbose, fprintf('    Debiasing converged after %d iterations (Diff=%.1e)\n', iter_d, d_diff); end
        debias_converged = true;
        break;
    end
end
if ~debias_converged
    warning('seal_MxNE:DebiasMaxIter', 'Debiasing did not converge in %d iterations.', max_iter_debias);
end
D_debias = Dk_comp;
end

function local_check_option(value, allowed)
if ~ismember(value, allowed)
    error('seal_MxNE:InvalidOption', 'Invalid value: %s. Allowed: %s', value, strjoin(allowed, ', '));
end
end

function L = local_power_iteration_kron(A, C, max_iter, tol, random_seed)
%LOCAL_POWER_ITERATION_KRON
% Estimate the largest singular value associated with kron(C', A).

    if nargin < 3 || isempty(max_iter), max_iter = 1000; end
    if nargin < 4 || isempty(tol),      tol = 1e-3; end
    if nargin < 5 || isempty(random_seed), random_seed = 0; end

    AS_size = size(C, 1);

    rng_state = rng;
    cleanupObj = onCleanup(@() rng(rng_state)); 
    rng(random_seed, 'twister');

    B = randn(AS_size, AS_size);
    B = B / max(norm(B, 'fro'), eps);

    ATA = A' * A;
    CCT = C * C';

    L0 = inf;
    L = 1.0;

    for iter = 1:max_iter 
        Y = ATA * B * CCT;
        L = norm(Y, 'fro');

        if ~isfinite(L) || L < eps
            L = 1.0;
            break;
        end

        if abs(L - L0) < tol
            break;
        end

        B = Y / L;
        L0 = L;
    end
end