function [Source_irMxNE, Param_irMxNE] = seal_irMxNE(Data, L, varargin)
%SEAL_IRMXNE Computes Iterative Reweighted Mixed-Norm Estimate (irMxNE).
%   [Source_irMxNE, Param_irMxNE] = seal_irMxNE(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Implements an iterative reweighted L1/L2 Mixed-Norm Estimate, often used
%   to approximate L0.5/L2 solutions, using Block Coordinate Descent.
%   Inspired by MNE-Python's iterative_mixed_norm_solver.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data (Y).
%       L (double matrix): Nchannels x Nsources original lead field matrix (G).
%
%   Optional Name-Value Pair Inputs:
%       'Alpha' (double): Regularization parameter (0-100). Default: 50.
%       'NumIRMXNEIter' (integer): Number of reweighting iterations. Default: 5.
%       'NumOrientations' (integer): Number of orientations per source location. Default: 1.
%       'NoiseCovariance' (double matrix): C_noise. Default: eye(Nchannels).
%       'DepthBiasComp' (struct): Default: struct('exponent', 0.8, 'limit', 10.0).
%       'Weights' (vector/struct): External spatial weights. Default: [].
%       'WeightsMin' (scalar): Default: 0.0.
%       'SolverMaxIter' (integer): Max iterations for inner BCD solver. Default: 1000.
%       'SolverTolerance' (double): Convergence tolerance for inner BCD. Default: 1e-4.
%       'IRMXNETolerance' (double): Tolerance for outer X change. Default: 1e-4.
%       'IRMXNEMinIter' (integer): Minimum outer iterations before convergence check. Default: 2.
%       'Debias' (logical): Perform amplitude debiasing after iterations. Default: true.
%       'DebiasMaxIter' (integer): Default: 1000.
%       'DebiasTolerance' (double): Default: 1e-7.
%       'PCAWhitening' (logical/integer): Default: true.
%       'RankNoiseCov' (integer/[]): Default: [].
%       'PickOri' (char/string): 'none' or 'vector'. Default: 'none'.
%       'Verbose' (logical): Print progress info. Default: true.
%
%   Outputs:
%       Source_irMxNE (double matrix): Estimated source activity.
%       Param_irMxNE (struct): Parameters and intermediate results, including:
%           .Converged (logical): True if outer loop met IRMXNETolerance.
%           .NumIRMXNEIterActual (integer): Number of outer iterations executed.
%           .E_path_per_iter (cell): Inner BCD primal-objective paths per outer iter.

    %% Input Parsing
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = false;

    defaultAlpha = 50.0;
    defaultNumIRMXNEIter = 5;
    if size(L,1) > 0, defaultNoiseCov = eye(size(L,1)); else, defaultNoiseCov = []; end
    defaultNumOrientations = 1;
    defaultDepthBiasComp = struct('exponent', 0.5, 'limit', 10.0); 
    defaultWeights = [];
    defaultWeightsMin = 0.0;
    defaultSolverMaxIter = 1000;
    defaultSolverTol = 1e-4;
    defaultIRMXNETol = 1e-4;
    defaultIRMXNEMinIter = 2;
    defaultDebias = true;
    defaultDebiasMaxIter = 1000;
    defaultDebiasTol = 1e-7;
    defaultPCAWhitening = true;
    defaultRankNoiseCov = [];
    defaultPickOri = 'none';

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

    addParameter(p, 'Alpha', defaultAlpha, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=100);
    addParameter(p, 'NumIRMXNEIter', defaultNumIRMXNEIter, @(x) isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x) && ismember(x,[1,3]));
    addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
    addParameter(p, 'DepthBiasComp', defaultDepthBiasComp, @(x) isstruct(x) || isempty(x));
    addParameter(p, 'Weights', defaultWeights, @(x) isnumeric(x) || isstruct(x) || isempty(x));
    addParameter(p, 'WeightsMin', defaultWeightsMin, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'SolverMaxIter', defaultSolverMaxIter, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'SolverTolerance', defaultSolverTol, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'IRMXNETolerance', defaultIRMXNETol, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'IRMXNEMinIter', defaultIRMXNEMinIter, @(x) isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p, 'Debias', defaultDebias, @islogical);
    addParameter(p, 'DebiasMaxIter', defaultDebiasMaxIter, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'DebiasTolerance', defaultDebiasTol, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'PCAWhitening', defaultPCAWhitening);
    addParameter(p, 'RankNoiseCov', defaultRankNoiseCov, @(x) (isnumeric(x) && isscalar(x) && x>0) || isempty(x));
    addParameter(p, 'PickOri', defaultPickOri, @(x) ischar(x) || isstring(x));
    addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

    parse(p, Data, L, varargin{:});
    Param_irMxNE.OptionsPassed = p.Results;

    % Unpack frequently used parameters
    Y_orig = p.Results.Data; G_orig = p.Results.L; alpha_user = p.Results.Alpha;
    nd = p.Results.NumOrientations; C_noise = p.Results.NoiseCovariance;
    depth_params = p.Results.DepthBiasComp; user_weights_input = p.Results.Weights;
    user_weights_min = p.Results.WeightsMin; solver_max_iter = p.Results.SolverMaxIter;
    solver_tol = p.Results.SolverTolerance; perform_debias = p.Results.Debias;
    debias_max_iter = p.Results.DebiasMaxIter; debias_tol = p.Results.DebiasTolerance;
    pca_whitening_opt = p.Results.PCAWhitening; rank_noise_cov_opt = p.Results.RankNoiseCov;
    pick_ori_opt = lower(char(p.Results.PickOri));
    n_irmxne_iter = p.Results.NumIRMXNEIter;
    irmxne_tol = p.Results.IRMXNETolerance;
    irmxne_min_iter = p.Results.IRMXNEMinIter;
    verbose = logical(p.Results.Verbose);
    local_check_option(pick_ori_opt, {'none', 'vector'});

    [Nchan, Nt] = size(Y_orig); [Nchan_L, Nsources_total] = size(G_orig);
    if mod(Nsources_total, nd) ~= 0, error('seal_irMxNE:NDMismatch', 'L Nsources not divisible by NumOrientations.'); end
    S_loc = Nsources_total / nd;

    if isempty(C_noise), C_noise = eye(Nchan); Param_irMxNE.OptionsPassed.NoiseCovariance = C_noise; end
    if Nchan ~= Nchan_L, error('seal_irMxNE:ChannelMismatch', 'Data and L channel mismatch.'); end

    % Constants
    ALPHA_MAX_SCALE = 0.01;        % MNE-Python convention

    % Initialize tracking
    Param_irMxNE.Converged           = false;
    Param_irMxNE.NumIRMXNEIterActual = 0;
    Param_irMxNE.E_path_per_iter     = {};

    %% 1. Prepare Gain Matrix (Whitening, PCA, Depth Comp, External Weights)
    iW = eye(Nchan); Y_w = Y_orig; G_w = G_orig; rank_whitener = Nchan;
    if ~isequal(C_noise, eye(Nchan))
        try
            [U_c, S_c_diag_mat] = eig((C_noise + C_noise')/2); s_c_vec = diag(S_c_diag_mat);
            [s_c_vec_sorted, sort_idx_c] = sort(s_c_vec, 'descend'); U_c_sorted = U_c(:, sort_idx_c);
            rank_eff = sum(s_c_vec_sorted > (Nchan * eps(max(s_c_vec_sorted))));
            rank_whitener = min(rank_eff, Nchan);
            if ~isempty(rank_noise_cov_opt) && ~islogical(rank_noise_cov_opt)
                rank_whitener = min(rank_noise_cov_opt, rank_whitener); 
            end
            s_inv_sqrt = zeros(Nchan,1); valid_s_idx = 1:rank_whitener;
            s_inv_sqrt(valid_s_idx) = 1 ./ sqrt(s_c_vec_sorted(valid_s_idx));
            iW_dense = diag(s_inv_sqrt) * U_c_sorted'; iW = iW_dense(1:rank_whitener, :);
            Y_w = iW * Y_orig; G_w = iW * G_orig;
        catch ME_white
            warning('seal_irMxNE:WhiteningFailed', 'Whitening: %s. Using unwhitened.', ME_white.message);
            iW=eye(Nchan); Y_w=Y_orig; G_w=G_orig; rank_whitener=Nchan; 
        end
    end
    Param_irMxNE.WhiteningMatrix = iW;

    M_processed = Y_w; Vh_pca = [];
    if pca_whitening_opt
        pca_rank = rank_whitener;
        if isnumeric(pca_whitening_opt) && ~islogical(pca_whitening_opt), pca_rank = min(pca_whitening_opt, rank_whitener); end
        if pca_rank < size(M_processed,1)
            if verbose, fprintf('irMxNE: Applying PCA to whitened data, reducing to rank %d\n', pca_rank); end
            [U_data, S_data_mat, V_data] = svd(M_processed, 'econ');
            M_processed = U_data(:, 1:pca_rank) * S_data_mat(1:pca_rank, 1:pca_rank);
            Vh_pca = V_data(:, 1:pca_rank)';
        end
    end

    source_weighting_init = ones(Nsources_total, 1);
    if ~isempty(depth_params) && isfield(depth_params, 'exponent')
        source_norms = zeros(S_loc, 1);
        for i = 1:S_loc, idx = (i-1)*nd + (1:nd); source_norms(i) = norm(G_w(:, idx), 'fro'); end
        source_norms(source_norms < eps) = eps;
        depth_weights_loc = source_norms .^ (-depth_params.exponent);
        if isfield(depth_params,'limit') && depth_params.limit > 0
            limit_val = min(depth_weights_loc) * (depth_params.limit^2);
            depth_weights_loc(depth_weights_loc > limit_val) = limit_val;
        end
        if nd == 1, source_weighting_init = depth_weights_loc; else, source_weighting_init = repelem(depth_weights_loc, nd); end
    end
    if ~isempty(user_weights_input)
        uw = [];
        if isstruct(user_weights_input) && isfield(user_weights_input, 'data'), uw = abs(max(user_weights_input.data, [], 2));
        elseif isnumeric(user_weights_input) && isvector(user_weights_input), uw = abs(user_weights_input(:)); end
        if ~isempty(uw) && length(uw) == S_loc
            if max(uw) > eps, uw = uw / max(uw); end
            if nd == 1, source_weighting_init = source_weighting_init .* uw; else, source_weighting_init = source_weighting_init .* repelem(uw, nd); end
        else
            warning('seal_irMxNE:UserWeightsDim', 'User weights ignored.'); 
        end
    end
    Param_irMxNE.InitialSourceWeightingForGain = source_weighting_init;
    G_eff_base = bsxfun(@times, G_w, source_weighting_init'); 

    mask_active_sources = true(Nsources_total, 1);
    if user_weights_min > 0 && ~isempty(user_weights_input) 
        mask_active_sources = source_weighting_init > (user_weights_min * (max(source_weighting_init)+eps));
        if ~any(mask_active_sources), error('seal_irMxNE:NoSourcesAfterWeightMin', 'No sources after WeightsMin.'); end
        G_eff_base = G_eff_base(:, mask_active_sources);
        source_weighting_init_masked = source_weighting_init(mask_active_sources);
        if verbose
            fprintf('irMxNE: Reduced source space to %d component(s) due to WeightsMin.\n', sum(mask_active_sources));
        end
    else
        source_weighting_init_masked = source_weighting_init;
    end
    if mod(size(G_eff_base, 2), nd) ~= 0
        error('seal_irMxNE:MaskBlockMismatch', ...
              'Active source mask does not respect orientation-block structure.');
    end
    Param_irMxNE.ActiveSourceMask_Initial = mask_active_sources;
    
    %% 2. Scaling for Solver
    GtM = G_eff_base' * M_processed;
    n_pos_eff = size(GtM,1)/nd;
    l2inf_norms_sq_per_group = zeros(n_pos_eff, 1);
    for i = 1:n_pos_eff
        idx = (i-1)*nd + (1:nd);
        l2inf_norms_sq_per_group(i) = sum(GtM(idx,:).^2, 'all');
    end
    alpha_max_val = sqrt(max(l2inf_norms_sq_per_group)) * ALPHA_MAX_SCALE;
    if alpha_max_val < eps, alpha_max_val = 1.0; end

    G_solver = G_eff_base / alpha_max_val;
    alpha_solver = alpha_user; 
    Param_irMxNE.AlphaSolver = alpha_solver; 
    Param_irMxNE.AlphaMaxScalingFactor = alpha_max_val;

    %% 3. Iterative Reweighted MxNE Solver
    if verbose
        fprintf('irMxNE: Starting iterative reweighted L1/L2 solution (max %d outer iterations)...\n', n_irmxne_iter);
    end
    
    X_active_iter_pca = zeros(size(G_solver,2), size(M_processed,2));
    current_reweighting_weights = ones(size(G_solver,2), 1);
    active_set_bcd = false(size(G_solver,2),1);
    last_outer_iter = 0;

    for k_iter = 1:n_irmxne_iter
        last_outer_iter = k_iter;
        if verbose, fprintf('  irMxNE Iteration %d/%d...\n', k_iter, n_irmxne_iter); end
        X_prev_iter_pca = X_active_iter_pca;
        
        % Apply current reweighting to columns of G
        G_iter_solver = bsxfun(@times, G_solver, current_reweighting_weights');

        [X_solved_at_iter, active_set_at_iter_mask, E_path_bcd, ~, ~] = ...
            local_mxne_bcd_solver(M_processed, G_iter_solver, alpha_solver, nd, ...
                solver_max_iter, solver_tol, max(1,round(solver_max_iter/10)), verbose);
        Param_irMxNE.E_path_per_iter{k_iter} = E_path_bcd;
        
        % Recover the "true" X under the original (unreweighted) variable
        % via the substitution X = diag(w) * X_solved.
        X_active_iter_pca_new = zeros(size(G_solver,2), size(M_processed,2));
        if any(active_set_at_iter_mask)
            X_active_iter_pca_new(active_set_at_iter_mask,:) = ...
                bsxfun(@times, X_solved_at_iter, current_reweighting_weights(active_set_at_iter_mask));
        end
        X_active_iter_pca = X_active_iter_pca_new;
        active_set_bcd = active_set_at_iter_mask;

        if k_iter < n_irmxne_iter
            if ~any(active_set_bcd)
                if verbose
                    fprintf('    No active sources in irMxNE iteration %d. Stopping reweighting.\n', k_iter);
                end
                Param_irMxNE.Converged = true;   % counts as "settled"
                break; 
            end
            
            % Update reweighting based on group L2 norms of current X.
            % This implements the L0.5 majorization-minimization update:
            %   w_g^{new} = 2 * ||X_g||_F^{1/2}
            X_for_gprime = X_active_iter_pca(active_set_bcd, :);
            n_active_pos = size(X_for_gprime, 1) / nd;
            group_L2_norms_sq = local_groups_norm2(X_for_gprime, nd, n_active_pos);
            g_X = sqrt(sqrt(group_L2_norms_sq));         % = ||X_g||_F^{1/2}
            new_reweights_loc = 2.0 * g_X;
            
            temp_current_reweighting_weights = ones(size(G_solver,2), 1);
            if nd > 1
                temp_current_reweighting_weights(active_set_bcd) = repelem(new_reweights_loc, nd);
            else
                temp_current_reweighting_weights(active_set_bcd) = new_reweights_loc;
            end
            current_reweighting_weights = temp_current_reweighting_weights;
        end
        
        % Outer convergence check (only after IRMXNEMinIter)
        if k_iter >= max(2, irmxne_min_iter)
            change_X = norm(X_active_iter_pca(:) - X_prev_iter_pca(:),'fro') / ...
                       (norm(X_prev_iter_pca(:),'fro') + eps);
            if verbose
                fprintf('    irMxNE Iter %d: X change = %.2e (Tol=%.1e)\n', k_iter, change_X, irmxne_tol);
            end
            if change_X < irmxne_tol
                if verbose
                    fprintf('  irMxNE converged (X change) at iteration %d.\n', k_iter);
                end
                Param_irMxNE.Converged = true;
                break;
            end
        end
    end
    
    if ~Param_irMxNE.Converged && last_outer_iter == n_irmxne_iter && n_irmxne_iter > 1
        % Only warn if we genuinely ran out of iterations without converging
        warning('seal_irMxNE:OuterLoopMaxIter', ...
                'irMxNE outer loop reached max iterations (%d) without meeting tol %.1e.', ...
                n_irmxne_iter, irmxne_tol);
    end
    Param_irMxNE.NumIRMXNEIterActual = last_outer_iter;
    
    X_solved_active_pca = X_active_iter_pca(active_set_bcd,:);
    final_active_set_mask_in_G_solver_space = active_set_bcd;

    %% 4. Debias (if requested)
    X_debiased_active_pca = X_solved_active_pca;
    if perform_debias && sum(final_active_set_mask_in_G_solver_space) > 0
        if verbose, fprintf('irMxNE: Performing debiasing...\n'); end
        G_for_debias = G_solver(:, final_active_set_mask_in_G_solver_space);
        debiasing_weights = local_compute_bias(M_processed, G_for_debias, X_solved_active_pca, nd, ...
                                               debias_max_iter, debias_tol, verbose);
        X_debiased_active_pca = X_solved_active_pca .* debiasing_weights;
        Param_irMxNE.DebiasingWeights = debiasing_weights;
    end

    %% 5. Reconstruct full time series if PCA was used
    X_debiased_active_full_time = X_debiased_active_pca;
    if ~isempty(Vh_pca)
        X_debiased_active_full_time = X_debiased_active_pca * Vh_pca;
    end

    %% 6. Reapply scaling and original source weighting
    active_indices_in_masked_space = find(final_active_set_mask_in_G_solver_space);
    sw_init_for_final_active = source_weighting_init_masked(active_indices_in_masked_space);

    X_final_active = bsxfun(@times, X_debiased_active_full_time, sw_init_for_final_active);
    X_final_active = X_final_active / alpha_max_val;

    %% 7. Format Output STC
    X_full_components = zeros(Nsources_total, Nt);
    temp_full_masked = zeros(sum(mask_active_sources), Nt);
    temp_full_masked(final_active_set_mask_in_G_solver_space, :) = X_final_active;
    X_full_components(mask_active_sources, :) = temp_full_masked;

    if nd > 1 && strcmpi(pick_ori_opt, 'none')
        Source_irMxNE_loc = zeros(S_loc, Nt);
        active_comp_indices_final = find(any(X_full_components,2));
        active_locs_final_unique = unique(ceil(active_comp_indices_final/nd));
        Param_irMxNE.ActiveSetIndices_Locations = active_locs_final_unique;
        Param_irMxNE.ActiveSetIndices_Components = active_comp_indices_final;

        for i_loc = 1:S_loc
            if ismember(i_loc, active_locs_final_unique)
                idx_comps = (i_loc-1)*nd + (1:nd);
                Source_irMxNE_loc(i_loc, :) = sqrt(sum(X_full_components(idx_comps,:).^2, 1));
            end
        end
        Source_irMxNE = Source_irMxNE_loc;
    else 
        Source_irMxNE = X_full_components;
        Param_irMxNE.ActiveSetIndices_Components = find(any(X_full_components,2));
        if nd > 1, Param_irMxNE.ActiveSetIndices_Locations = unique(ceil(Param_irMxNE.ActiveSetIndices_Components/nd));
        else, Param_irMxNE.ActiveSetIndices_Locations = Param_irMxNE.ActiveSetIndices_Components; end
    end
    
    Param_irMxNE.ResolutionMatrixSource = NaN; 
end


function [X_solved_active, active_set_mask, E_path, converged, last_iter] = local_mxne_bcd_solver(...
    M_data, G_solver, alpha_penalty_val, O_orient, max_iter, tol, dgap_freq, verbose)
% Solves L1/L2 MxNE using Block Coordinate Descent.

    [~, Nt_proc] = size(M_data);
    [~, Ns_proc] = size(G_solver);
    N_pos_proc = Ns_proc / O_orient;

    X = zeros(Ns_proc, Nt_proc);
    R = M_data;

    E_path = nan(ceil(max_iter / dgap_freq) + 2, 1);
    n_E = 0;
    converged = false;
    last_iter = max_iter;

    % Lipschitz constants. Use eig() on small SPD matrix (faster than SVD).
    lc = zeros(N_pos_proc, 1);
    if O_orient == 1
        lc = sum(G_solver.^2, 1)';
    else
        for j_pos = 1:N_pos_proc
            idx = (j_pos-1)*O_orient + (1:O_orient);
            G_j_block = G_solver(:, idx);
            GtG = G_j_block' * G_j_block;
            GtG = (GtG + GtG') * 0.5;
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
            
            X_j_unthres = X_j_current + (1/lc(j_pos)) * (G_j' * R);
            group_norm_L2F = norm(X_j_unthres, 'fro');
            
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
        
        if mod(iter_bcd, dgap_freq) == 0 || iter_bcd == max_iter || iter_bcd == 1
            [actual_gap, p_obj, ~] = local_dgap_l21(M_data, G_solver, X, alpha_penalty_val, O_orient);
            n_E = n_E + 1;
            E_path(n_E) = p_obj;
            
            if iter_bcd > 1 
                if verbose
                    fprintf('    BCD Iter %4d: Pobj=%.4e, DGap=%.2e (Tol=%.1e), Nactive_locs=%d\n', ...
                        iter_bcd, p_obj, actual_gap, tol, ...
                        sum(any(reshape(local_compute_active_mask(X, O_orient), O_orient,[]),1)));
                end
                if actual_gap < tol
                    if verbose, fprintf('    BCD converged (Duality Gap).\n'); end
                    converged = true;
                    break;
                end
            end
        end
        
        if iter_bcd > 1 && norm(X(:) - X_old_iter(:),'fro') / (norm(X_old_iter(:),'fro')+eps) < tol/100
            if verbose, fprintf('    BCD converged (X change) at iter %d.\n', iter_bcd); end
            converged = true;
            break;
        end
    end
    
    if ~converged
        warning('seal_irMxNE:BCDSolverMaxIter', 'BCD solver reached max iterations (%d).', max_iter);
    end
    E_path = E_path(1:n_E);

    % Recompute active mask from final X
    active_set_mask = local_compute_active_mask(X, O_orient);
    X_solved_active = X(active_set_mask, :);
end


function active_set_mask = local_compute_active_mask(X, O_orient)
    [Ns, ~] = size(X);
    if O_orient == 1
        active_set_mask = any(X ~= 0, 2);
    else
        N_pos = Ns / O_orient;
        active_blocks = false(N_pos, 1);
        for j = 1:N_pos
            idx = (j-1)*O_orient + (1:O_orient);
            active_blocks(j) = any(any(X(idx, :) ~= 0));
        end
        active_set_mask = repelem(active_blocks, O_orient);
    end
end


function [gap, p_obj, d_obj] = local_dgap_l21(M_data, G_full, X_full, alpha_val, nd)
    GX = G_full * X_full;
    R = M_data - GX;
    nR2 = sum(R(:).^2);

    if nd > 1
        X_tensor = reshape(X_full, nd, [], size(X_full, 2));
        X_norms = sum(sqrt(sum(sum(X_tensor.^2, 1), 3)));
    else
        X_norms = sum(sqrt(sum(X_full.^2, 2)));
    end
    p_obj = 0.5 * nR2 + alpha_val * X_norms;

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
% FISTA debiasing: solves min_D 0.5*||M - G*diag(D)*X||_F^2 s.t. D>=1.

[N_active_comp, ~] = size(X_active_comp);
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
    warning('seal_irMxNE:DebiasMaxIter', 'Debiasing did not converge in %d iterations.', max_iter_debias);
end
D_debias = Dk_comp;
end


function local_check_option(value, allowed)
    if ~ismember(value, allowed)
        error('seal_irMxNE:InvalidOption', 'Invalid value: %s. Allowed: %s', value, strjoin(allowed, ', '));
    end
end


function group_L2_norms_sq_val = local_groups_norm2(A_matrix, n_orient_val, n_positions_val)
% Squared Frobenius norm per group (vectorized).
    if n_orient_val == 1
        group_L2_norms_sq_val = sum(A_matrix.^2, 2);
    else
        % Reshape: (n_orient * n_pos) x Nt -> n_orient x n_pos x Nt
        A_tensor = reshape(A_matrix, n_orient_val, n_positions_val, []);
        group_L2_norms_sq_val = squeeze(sum(sum(A_tensor.^2, 1), 3));
        group_L2_norms_sq_val = group_L2_norms_sq_val(:);   % ensure column
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