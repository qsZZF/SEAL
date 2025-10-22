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
%       'NumOrientations' (integer): Number of orientations per source location (O). Default: 1.
%       'NoiseCovariance' (double matrix): Noise covariance (C_noise). Default: eye(Nchannels).
%       'DepthBiasComp' (struct): Parameters for depth weighting. Default: struct('exponent', 0.8, 'limit', 10.0).
%       'Weights' (vector/struct): External spatial weights. Default: [].
%       'WeightsMin' (scalar): Minimum threshold for external weights. Default: 0.0.
%       'SolverMaxIter' (integer): Max iterations for BCD solver. Default: 3000.
%       'SolverTolerance' (double): Convergence tolerance for BCD solver. Default: 1e-4.
%       'Debias' (logical): Perform amplitude debiasing. Default: true.
%       'DebiasMaxIter' (integer): Max iterations for FISTA debiasing. Default: 1000.
%       'DebiasTolerance' (double): Tolerance for FISTA debiasing. Default: 1e-7.
%       'PCAWhitening' (logical/integer): Reduce data rank via PCA after whitening. Default: true.
%       'PickOri' (char/string): 'none' or 'vector'. Default: 'none'.
%
%   Outputs:
%       Source_MxNE (double matrix): Estimated source activity.
%       Param_MxNE (struct): Parameters and intermediate results.
%
%   Copyright 2025 Your Name/Institution

%% Input Parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = false;

defaultAlpha = 50.0;
if size(L,1) > 0, defaultNoiseCov = eye(size(L,1)); else, defaultNoiseCov = []; end
defaultNumOrientations = 1;
defaultDepthBiasComp = struct('exponent', 0.8, 'limit', 10.0);
defaultWeights = [];
defaultWeightsMin = 0.0;
defaultSolverMaxIter = 3000;
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
pick_ori_opt = lower(p.Results.PickOri);
local_check_option(pick_ori_opt, {'none', 'vector'});

[Nchan, Nt] = size(Data_orig); [Nchan_L, Nsources_total] = size(L_orig);
S_loc = Nsources_total / nd;

if isempty(C_noise), C_noise = eye(Nchan); Param_MxNE.OptionsPassed.NoiseCovariance = C_noise; end
if Nchan ~= Nchan_L, error('Data and L channel mismatch.'); end
if mod(Nsources_total, nd) ~= 0, error('L Nsources not divisible by NumOrientations.'); end

%% 1. Prepare Gain Matrix (Whitening, PCA, Depth Comp, External Weights)
iW = eye(Nchan); Data_w = Data_orig; L_w = L_orig; rank_whitener = Nchan;
if ~isequal(C_noise, eye(Nchan))
    try
        [U_c, S_c_diag_mat] = eig((C_noise + C_noise')/2); s_c_vec = diag(S_c_diag_mat);
        [s_c_vec_sorted, sort_idx_c] = sort(s_c_vec, 'descend'); U_c_sorted = U_c(:, sort_idx_c);
        rank_eff = sum(s_c_vec_sorted > (Nchan * eps(max(s_c_vec_sorted))));
        rank_whitener = min(rank_eff, Nchan);
        
        s_inv_sqrt = zeros(Nchan,1);
        valid_s_idx = 1:rank_whitener;
        s_inv_sqrt(valid_s_idx) = 1 ./ sqrt(s_c_vec_sorted(valid_s_idx));
        
        iW_dense = diag(s_inv_sqrt) * U_c_sorted';
        iW = iW_dense(1:rank_whitener, :);
        Data_w = iW * Data_orig; L_w = iW * L_orig;
    catch ME_white 
        warning('seal_MxNE:WhiteningFailed', 'Whitening: %s. Using unwhitened.', ME_white.message); 
        iW = eye(Nchan); Data_w = Data_orig; L_w = L_orig; rank_whitener = Nchan; 
    end
end
Param_MxNE.WhiteningMatrix = iW;

M_processed = Data_w; Vh_pca = [];
if pca_whitening_opt
    pca_rank = rank_whitener;
    if isnumeric(pca_whitening_opt), pca_rank = min(pca_whitening_opt, rank_whitener); end
    if pca_rank < size(M_processed,1)
        fprintf('MxNE: Applying PCA to whitened data, reducing to rank %d\n', pca_rank);
        [U_data, S_data_mat, V_data] = svd(M_processed, 'econ'); % V_data is V, V_data_T is V'
        M_processed = U_data(:, 1:pca_rank) * S_data_mat(1:pca_rank, 1:pca_rank);
        Vh_pca = V_data(:, 1:pca_rank)'; % Store V'
    end
end

source_weighting = ones(Nsources_total, 1);
if ~isempty(depth_params) && isfield(depth_params, 'exponent')
    source_norms = zeros(S_loc, 1);
    for i = 1:S_loc
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
    if ~isempty(uw) && length(uw) == S_loc
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
if user_weights_min > 0 && ~isempty(user_weights_input) % Only apply if user_weights were given
    mask_active_sources = source_weighting > (user_weights_min * (max(source_weighting)+eps) );
    if ~any(mask_active_sources), error('seal_MxNE:NoSourcesAfterWeightMin', 'No sources after WeightsMin.'); end
    G_eff = G_eff(:, mask_active_sources);
    fprintf('MxNE: Reduced source space to %d component(s) due to WeightsMin.\n', sum(mask_active_sources));
end
Param_MxNE.ActiveSourceMask_Initial = mask_active_sources; % Mask applied before solver

%% 2. Scaling for Solver
GtM = G_eff' * M_processed;
n_pos_eff = size(GtM,1)/nd;
l2inf_norms_sq_per_group = zeros(n_pos_eff, 1);
for i = 1:n_pos_eff 
    idx = (i-1)*nd + (1:nd); 
    l2inf_norms_sq_per_group(i) = sum(GtM(idx,:).^2, 'all'); 
end % Sum over O and Nt
alpha_max_val = sqrt(max(l2inf_norms_sq_per_group));
alpha_max_val = alpha_max_val * 0.01; % MNE-Python scaling convention
if alpha_max_val < eps
    alpha_max_val = 1.0; 
end

G_solver = G_eff / alpha_max_val;
alpha_solver = alpha_user; % User alpha (0-100) is applied to this scaled G
Param_MxNE.AlphaSolver = alpha_solver;
Param_MxNE.AlphaMaxScalingFactor = alpha_max_val;

%% 3. Solve MxNE (L1/L2) using Block Coordinate Descent
fprintf('MxNE: Solving L1/L2 problem with BCD (n_mxne_iter=1)...\n');
[X_active_pca, active_set_solver_mask, ~] = local_mxne_bcd_solver(...
    M_processed, G_solver, alpha_solver, nd, ...
    solver_max_iter, solver_tol, max(1,round(solver_max_iter/100))); % dgap_freq

if sum(active_set_solver_mask) == 0
    warning('seal_MxNE:NoActiveSourcesFound', 'No active sources found by BCD. Alpha might be too high.');
    Source_MxNE = zeros(Nsources_total, Nt);
    Param_MxNE.ActiveSetIndices_Solver = []; Param_MxNE.ActiveSetIndices_Locations = [];
    return;
end

%% 4. Debias (if requested)
X_debiased_active_pca = X_active_pca; % Initialize with BCD solution
if perform_debias && sum(active_set_solver_mask) > 0
    fprintf('MxNE: Performing debiasing...\n');
    G_for_debias = G_solver(:, active_set_solver_mask);
    debiasing_weights = local_compute_bias(M_processed, G_for_debias, X_active_pca, nd, ...
        debias_max_iter, debias_tol);
    X_debiased_active_pca = X_active_pca .* debiasing_weights; % Element-wise scaling by D_j for each component in group j
    Param_MxNE.DebiasingWeights = debiasing_weights;
end

%% 5. Reconstruct full time series if PCA was used
X_debiased_active_full_time = X_debiased_active_pca;
if ~isempty(Vh_pca)
    X_debiased_active_full_time = X_debiased_active_pca * Vh_pca;
end

%% 6. Reapply scaling and map back to original source space
% X_final_active is for sources corresponding to G_eff (i.e., G_w * diag(source_weighting))
% The solver output X is for G_solver = G_eff / alpha_max_val.
% So, X_for_G_eff = X_debiased_active_full_time / alpha_max_val.
% Since G_eff = G_w * diag(source_weighting), the actual source activity for G_w is
% S_for_Gw = diag(source_weighting_active) * X_for_G_eff.
% MNE-Python's _reapply_source_weighting multiplies X by source_weighting.
% X_final = X_sol * W_source_active / alpha_max

source_weighting_eff = source_weighting(mask_active_sources); % Weights for G_eff
source_weighting_eff_active = source_weighting_eff(active_set_solver_mask);

X_final_active = bsxfun(@times, X_debiased_active_full_time, source_weighting_eff_active);
X_final_active = X_final_active / alpha_max_val;

%% 7. Format Output STC
X_full_components = zeros(Nsources_total, Nt);
temp_full_masked = zeros(sum(mask_active_sources), Nt);
temp_full_masked(active_set_solver_mask, :) = X_final_active;
X_full_components(mask_active_sources, :) = temp_full_masked;

if nd > 1 && strcmpi(pick_ori_opt, 'none')
    Source_MxNE_loc = zeros(S_loc, Nt);
    active_comp_indices_final = find(any(X_full_components,2));
    active_locs_final_unique = unique(ceil(active_comp_indices_final/nd));
    Param_MxNE.ActiveSetIndices_Locations = active_locs_final_unique;
    Param_MxNE.ActiveSetIndices_Components = active_comp_indices_final; % Store component indices too
    
    for i_loc = 1:S_loc
        if ismember(i_loc, active_locs_final_unique)
            idx_comps = (i_loc-1)*nd + (1:nd);
            Source_MxNE_loc(i_loc, :) = sqrt(sum(X_full_components(idx_comps,:).^2, 1));
        end
    end
    Source_MxNE = Source_MxNE_loc;
else % O == 1 or pick_ori == 'vector'
    Source_MxNE = X_full_components;
    Param_MxNE.ActiveSetIndices_Components = find(any(X_full_components,2));
    if nd > 1, Param_MxNE.ActiveSetIndices_Locations = unique(ceil(Param_MxNE.ActiveSetIndices_Components/nd));
    else, Param_MxNE.ActiveSetIndices_Locations = Param_MxNE.ActiveSetIndices_Components; end
end

Param_MxNE.ResolutionMatrixSource = NaN; % Non-linear method
end

%% --- Local Helper Functions for MxNE Solver ---

function [X_solved_active, active_set_mask, E_path] = local_mxne_bcd_solver(...
    M_data, G_solver, alpha_penalty_val, O_orient, max_iter, tol, dgap_freq)
% Solves L1/L2 MxNE using Block Coordinate Descent.
% M_data: Nchan_proc x Nt_proc
% G_solver: Nchan_proc x Nsources_active_masked_scaled
% alpha_penalty_val: scalar regularization (0-100, applied to scaled G_solver)

[N_proc, Nt_proc] = size(M_data);
[~, Ns_proc] = size(G_solver);
N_pos_proc = Ns_proc / O_orient;

X = zeros(Ns_proc, Nt_proc);
R = M_data;
active_set_mask = false(Ns_proc, 1);

E_path = [];
highest_d_obj = -Inf;

% Calculate Lipschitz constants (lc) for each source position block
lc = zeros(N_pos_proc, 1);
if O_orient == 1
    lc = sum(G_solver.^2, 1)';
else
    for j_pos = 1:N_pos_proc
        idx = (j_pos-1)*O_orient + (1:O_orient);
        G_j_block = G_solver(:, idx);
        lc(j_pos) = norm(G_j_block, 'fro')^2;
        if lc(j_pos) < eps, lc(j_pos) = eps; end
    end
end

alpha_scaled_per_block = alpha_penalty_val ./ lc; % alpha_j / L_j (penalty for block j)

fprintf('  BCD Solver: MaxIter=%d, Tol=%.1e, Nsources_proc=%d, Npos_proc=%d\n', max_iter, tol, Ns_proc, N_pos_proc);

for iter_bcd = 1:max_iter
    X_old_iter = X;
    
    for j_pos = 1:N_pos_proc
        idx_block = (j_pos-1)*O_orient + (1:O_orient);
        G_j = G_solver(:, idx_block);
        X_j_current = X(idx_block, :);
        
        if any(active_set_mask(idx_block))
            R = R + G_j * X_j_current;
        end
        
        X_j_unthres = X_j_current + (1/lc(j_pos)) * (G_j' * R);
        
        group_norm_L2F = norm(X_j_unthres, 'fro');
        
        if group_norm_L2F <= alpha_scaled_per_block(j_pos)
            X_j_new = zeros(O_orient, Nt_proc);
            active_set_mask(idx_block) = false;
        else
            shrink_factor = 1 - alpha_scaled_per_block(j_pos) / group_norm_L2F;
            X_j_new = X_j_unthres * shrink_factor;
            active_set_mask(idx_block) = true;
        end
        
        if any(active_set_mask(idx_block)) || any(X_j_current(:) ~= 0)
            R = R - G_j * (X_j_new - X_j_current);
        end
        X(idx_block, :) = X_j_new;
    end
    
    if mod(iter_bcd, dgap_freq) == 0 || iter_bcd == max_iter || iter_bcd == 1
        current_active_comp_indices = find(active_set_mask);
        if isempty(current_active_comp_indices)
            p_obj = 0.5 * sum(M_data(:).^2); d_obj = 0;
        else
            G_curr_active_comps = G_solver(:, current_active_comp_indices);
            X_curr_active_comps = X(current_active_comp_indices, :);
            [~, p_obj, d_obj, ~] = local_dgap_l21(M_data, G_curr_active_comps, X_curr_active_comps, alpha_penalty_val, O_orient, active_set_mask);
        end
        E_path(end+1) = p_obj;
        highest_d_obj = max(highest_d_obj, d_obj);
        actual_gap = p_obj - highest_d_obj;
        
        if iter_bcd > 1 % Only print/check convergence after first full pass
            fprintf('    BCD Iter %4d: Pobj=%.4e, DGap=%.2e (Tol=%.1e), Nactive_locs=%d\n', ...
                iter_bcd, p_obj, actual_gap, tol, sum(any(reshape(active_set_mask,O_orient,[]),1)));
            if actual_gap < tol
                fprintf('    BCD converged (Duality Gap).\n');
                break;
            end
        end
    end
    
    if iter_bcd > 1 && norm(X(:) - X_old_iter(:),'fro') / (norm(X_old_iter(:),'fro')+eps) < tol/100
        fprintf('    BCD converged (X change) at iter %d.\n', iter_bcd);
        break;
    end
end
if iter_bcd == max_iter, warning('seal_MxNE:BCDSolverMaxIter', 'BCD solver reached max iterations.'); end

X_solved_active = X(active_set_mask, :);
end

function [gap, p_obj, d_obj, R_residual] = local_dgap_l21(M_data, G_active_comp, X_active_comp, alpha_penalty_val, O_orient, active_set_mask_full)
% Duality gap for L1/L2 mixed-norm. M, G, X are for active components.
% active_set_mask_full is for the original G_solver space.
if isempty(X_active_comp) || isempty(G_active_comp)
    R_residual = M_data; nR2 = sum(R_residual(:).^2);
    p_obj = 0.5 * nR2; d_obj = 0; gap = p_obj;
    if isnan(gap), gap = p_obj; end
    return;
end

GX = G_active_comp * X_active_comp;
R_residual = M_data - GX;
nR2 = sum(R_residual(:).^2);

penalty_l21 = 0;
N_pos_active = size(X_active_comp,1) / O_orient; % Number of active positions
for i_pos_iter = 1:N_pos_active % Iterate over active positions
    idx_block_in_X_active = (i_pos_iter-1)*O_orient + (1:O_orient);
    penalty_l21 = penalty_l21 + norm(X_active_comp(idx_block_in_X_active,:), 'fro');
end
p_obj = 0.5 * nR2 + alpha_penalty_val * penalty_l21;

Gt_R = G_active_comp' * R_residual;

% Need to consider Gt_R for ALL original positions to find the max group norm
% This requires mapping G_active_comp back to the full G_solver space, or
% calculating G_solver' * R_residual directly.
% For simplicity, if G_active_comp is G_solver(:,active_set_mask_full(active_set_mask_full_indices_in_masked_space))
% this is tricky. Let's assume G_active_comp is G_solver(:, current_active_set_indices)
% The dual norm should be on G_solver' * R_residual
% G_solver_T_R = G_solver' * R_residual (where G_solver is the one used in BCD)
% This part is simplified in MNE-Python's dgap_l21 as it operates on the current active G.
% For a true duality gap, it should use the full G.
% However, for convergence of BCD on active set, using G_active is standard.

dual_group_norms_F = zeros(N_pos_active,1);
for i_pos_iter = 1:N_pos_active
    idx_block_in_Gt_R = (i_pos_iter-1)*O_orient + (1:O_orient);
    dual_group_norms_F(i_pos_iter) = norm(Gt_R(idx_block_in_Gt_R,:), 'fro');
end
dual_norm_l2inf = max(dual_group_norms_F);
if isempty(dual_norm_l2inf), dual_norm_l2inf = 0; end % If no active sources

scaling = alpha_penalty_val / (dual_norm_l2inf + eps);
scaling = min(scaling, 1.0);

d_obj = (scaling - 0.5 * (scaling^2)) * nR2 + scaling * sum(sum(R_residual .* GX));
gap = p_obj - d_obj;
if isnan(gap) && ~isnan(p_obj), gap = p_obj; end
end

function D_debias = local_compute_bias(M_data, G_active_comp, X_active_comp, O_orient, max_iter_debias, tol_debias)
% Solves: min_D 0.5 * || M_data - G_active_comp * diag(D_loc_avg_repeated) * X_active_comp ||_F^2
% s.t. D_loc_avg_j >= 1 (for each source location j)
% D_debias is N_active_comp x 1, but effectively N_active_pos x 1 then repeated.

[N_active_comp, Nt_proc] = size(X_active_comp);
N_pos_active = N_active_comp / O_orient;

Dk_comp = ones(N_active_comp, 1); % Debiasing weights per component
Yk_comp = ones(N_active_comp, 1); % FISTA extrapolation variable
tk = 1.0;

% Heuristic Lipschitz constant (can be refined with power_iteration_kron)
% L_fista = 1.1 * norm(G_active_comp, 'fro')^2 * norm(X_active_comp, 'fro')^2; % Very rough
% For each D_j (location), the gradient involves G_j X_j terms.
% Let's try to estimate L more carefully for the FISTA on D.
% The gradient w.r.t D_s (scalar for location s) involves terms like
% sum_o (G_so X_so_t)' * (G_so X_so_t). Max of these sum of squares.
L_fista_val = 0;
for i_pos = 1:N_pos_active
    idx_block = (i_pos-1)*O_orient + (1:O_orient);
    G_block = G_active_comp(:, idx_block); % Nchan x O
    X_block = X_active_comp(idx_block, :); % O x Nt
    % Effective "regressor" for D_loc(i_pos) is sum_o (G_io X_io_t)
    % This is (G_block .* X_block_reshaped_for_hadamard_with_G)' * (G_block .* X_block_reshaped)
    % Simpler: use sum of squares of G_block_ij * X_block_jk
    % Python: lipschitz_constant = 1.1 * power_iteration_kron(G, X)
    % This is complex to translate directly. Using a simpler heuristic:
    % Max of || G_j X_j ||_F^2 for each block j
    Y_j_eff = G_block * X_block; % Nchan x Nt
    L_fista_val = max(L_fista_val, norm(Y_j_eff, 'fro')^2);
end
if L_fista_val < eps, L_fista_val = 1.0; else, L_fista_val = 1.1 * L_fista_val; end

fprintf('  Debiasing (FISTA): MaxIter=%d, Tol=%.1e, L_fista=%.2e\n', max_iter_debias, tol_debias, L_fista_val);

for iter_d = 1:max_iter_debias
    D0_comp = Dk_comp;
    
    % Yk_comp is the FISTA extrapolation point for D_comp
    % M_pred = G_active_comp * (diag(Yk_comp) * X_active_comp);
    M_pred = G_active_comp * bsxfun(@times, X_active_comp, Yk_comp); % Equivalent if Yk_comp is column
    Residual = M_data - M_pred;
    
    % Gradient of 0.5*||M - G*diag(D)*X||_F^2 w.r.t D_i (component i)
    % Grad_D_i = - tr( (G(:,i) * X(i,:))' * Residual )
    %            = - sum_{t} (G(:,i)' * Residual(:,t)) * X(i,t)
    %            = - sum_{t} GtRes_it * X_it
    Gt_Res = G_active_comp' * Residual; % N_active_comp x Nt
    Grad_D_comp = -sum(Gt_Res .* X_active_comp, 2); % N_active_comp x 1
    
    D_update = Yk_comp - Grad_D_comp / L_fista_val;
    
    % Proximal step: D_j_loc_avg >= 1
    if O_orient > 1
        D_reshaped_loc = mean(reshape(D_update, O_orient, N_pos_active), 1)'; % N_pos_active x 1
        D_reshaped_loc = max(D_reshaped_loc, 1.0);
        Dk_comp = repelem(D_reshaped_loc, O_orient);
    else % O_orient == 1
        Dk_comp = max(D_update, 1.0);
    end
    
    t_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * tk^2));
    Yk_comp = Dk_comp + ((tk - 1.0) / t_new) * (Dk_comp - D0_comp);
    tk = t_new;
    
    d_diff = max(abs(Dk_comp - D0_comp));
    if d_diff < tol_debias && iter_d > 1
        fprintf('    Debiasing converged after %d iterations (Diff=%.1e)\n', iter_d, d_diff);
        break;
    end
end
if iter_d == max_iter_debias, warning('seal_MxNE:DebiasMaxIter', 'Debiasing did not converge.'); end
D_debias = Dk_comp; % N_active_comp x 1 vector of debiasing weights
end

function local_check_option(key, value, allowed)
if ~ismember(value, allowed)
    error('Invalid value for %s. Allowed options are: %s', key, strjoin(allowed, ', '));
end
end
