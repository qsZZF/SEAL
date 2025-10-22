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
%                                  (MNE-Python often uses 1 for MxNE, >1 for irMxNE)
%       'NumOrientations' (integer): O. Default: 1.
%       'NoiseCovariance' (double matrix): C_noise. Default: eye(Nchannels).
%       'DepthBiasComp' (struct): Default: struct('exponent', 0.8, 'limit', 10.0).
%       'Weights' (vector/struct): External spatial weights. Default: [].
%       'WeightsMin' (scalar): Default: 0.0.
%       'SolverMaxIter' (integer): Max iterations for inner BCD solver. Default: 1000.
%       'SolverTolerance' (double): Convergence tolerance for inner BCD. Default: 1e-4.
%       'IRMXNETolerance' (double): Tolerance for convergence of X between reweighting iterations. Default: 1e-4.
%       'Debias' (logical): Perform amplitude debiasing after iterations. Default: true.
%       'DebiasMaxIter' (integer): Default: 1000.
%       'DebiasTolerance' (double): Default: 1e-7.
%       'PCAWhitening' (logical/integer): Default: true.
%       'RankNoiseCov' (integer/[]): Default: [].
%       'PickOri' (char/string): 'none' or 'vector'. Default: 'none'.
%
%   Outputs:
%       Source_irMxNE (double matrix): Estimated source activity.
%       Param_irMxNE (struct): Parameters and intermediate results.
%
%   Copyright 2025 Your Name/Institution

    %% Input Parsing
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = false; % For simplicity, not passing unmatched to BCD for now

    defaultAlpha = 50.0;
    defaultNumIRMXNEIter = 5; % For irMxNE, typically > 1
    if size(L,1) > 0, defaultNoiseCov = eye(size(L,1)); else, defaultNoiseCov = []; end
    defaultNumOrientations = 1;
    defaultDepthBiasComp = struct('exponent', 0.8, 'limit', 10.0); % MNE-Python like
    defaultWeights = [];
    defaultWeightsMin = 0.0;
    defaultSolverMaxIter = 1000; % Max iter for inner BCD loop
    defaultSolverTol = 1e-4;     % Tol for inner BCD loop
    defaultIRMXNETol = 1e-4;     % Tol for outer reweighting loop X change
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
    addParameter(p, 'Debias', defaultDebias, @islogical);
    addParameter(p, 'DebiasMaxIter', defaultDebiasMaxIter, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'DebiasTolerance', defaultDebiasTol, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'PCAWhitening', defaultPCAWhitening);
    addParameter(p, 'RankNoiseCov', defaultRankNoiseCov, @(x) (isnumeric(x) && isscalar(x) && x>0) || isempty(x));
    addParameter(p, 'PickOri', defaultPickOri, @(x) ischar(x) || isstring(x));

    parse(p, Data, L, varargin{:});
    Param_irMxNE.OptionsPassed = p.Results;

    % Unpack frequently used parameters
    Y_orig = p.Results.Data; G_orig = p.Results.L; alpha_user = p.Results.Alpha;
    O = p.Results.NumOrientations; C_noise = p.Results.NoiseCovariance;
    depth_params = p.Results.DepthBiasComp; user_weights_input = p.Results.Weights;
    user_weights_min = p.Results.WeightsMin; solver_max_iter = p.Results.SolverMaxIter;
    solver_tol = p.Results.SolverTolerance; perform_debias = p.Results.Debias;
    debias_max_iter = p.Results.DebiasMaxIter; debias_tol = p.Results.DebiasTolerance;
    pca_whitening_opt = p.Results.PCAWhitening; rank_noise_cov_opt = p.Results.RankNoiseCov;
    pick_ori_opt = lower(p.Results.PickOri);
    n_irmxne_iter = p.Results.NumIRMXNEIter;
    irmxne_tol = p.Results.IRMXNETolerance;
    local_check_option(pick_ori_opt, {'none', 'vector'});


    [Nchan, Nt] = size(Y_orig); [Nchan_L, Nsources_total] = size(G_orig);
    S_loc = Nsources_total / O;

    if isempty(C_noise), C_noise = eye(Nchan); Param_irMxNE.OptionsPassed.NoiseCovariance = C_noise; end
    if Nchan ~= Nchan_L, error('Data and L channel mismatch.'); end
    if mod(Nsources_total, O) ~= 0, error('L Nsources not divisible by NumOrientations.'); end

    %% 1. Prepare Gain Matrix (Whitening, PCA, Depth Comp, External Weights)
    % This part is identical to seal_MxNE and can be a shared helper function eventually
    iW = eye(Nchan); Y_w = Y_orig; G_w = G_orig; rank_whitener = Nchan;
    if ~isequal(C_noise, eye(Nchan))
        try
            [U_c, S_c_diag_mat] = eig((C_noise + C_noise')/2); s_c_vec = diag(S_c_diag_mat);
            [s_c_vec_sorted, sort_idx_c] = sort(s_c_vec, 'descend'); U_c_sorted = U_c(:, sort_idx_c);
            rank_eff = sum(s_c_vec_sorted > (Nchan * eps(max(s_c_vec_sorted))));
            rank_whitener = min(rank_eff, Nchan);
            if ~isempty(rank_noise_cov_opt), rank_whitener = min(rank_noise_cov_opt, rank_whitener); end
            s_inv_sqrt = zeros(Nchan,1); valid_s_idx = 1:rank_whitener;
            s_inv_sqrt(valid_s_idx) = 1 ./ sqrt(s_c_vec_sorted(valid_s_idx));
            iW_dense = diag(s_inv_sqrt) * U_c_sorted'; iW = iW_dense(1:rank_whitener, :);
            Y_w = iW * Y_orig; G_w = iW * G_orig;
        catch ME_white, warning('seal_irMxNE:WhiteningFailed', 'Whitening: %s. Using unwhitened.', ME_white.message); iW=eye(Nchan); Y_w=Y_orig; G_w=G_orig; rank_whitener=Nchan; 
        end
    end
    Param_irMxNE.WhiteningMatrix = iW;

    M_processed = Y_w; Vh_pca = [];
    if pca_whitening_opt
        pca_rank = rank_whitener;
        if isnumeric(pca_whitening_opt), pca_rank = min(pca_whitening_opt, rank_whitener); end
        if pca_rank < size(M_processed,1)
            fprintf('irMxNE: Applying PCA to whitened data, reducing to rank %d\n', pca_rank);
            [U_data, S_data_mat, V_data] = svd(M_processed, 'econ');
            M_processed = U_data(:, 1:pca_rank) * S_data_mat(1:pca_rank, 1:pca_rank);
            Vh_pca = V_data(:, 1:pca_rank)';
        end
    end

    source_weighting_init = ones(Nsources_total, 1); % Initial combined weights for G
    if ~isempty(depth_params) && isfield(depth_params, 'exponent')
        source_norms = zeros(S_loc, 1);
        for i = 1:S_loc, idx = (i-1)*O + (1:O); source_norms(i) = norm(G_w(:, idx), 'fro'); end
        source_norms(source_norms < eps) = eps;
        depth_weights_loc = source_norms .^ (-depth_params.exponent);
        if isfield(depth_params,'limit') && depth_params.limit > 0
            limit_val = min(depth_weights_loc) * (depth_params.limit^2);
            depth_weights_loc(depth_weights_loc > limit_val) = limit_val;
        end
        if O == 1, source_weighting_init = depth_weights_loc; else, source_weighting_init = repelem(depth_weights_loc, O); end
    end
    if ~isempty(user_weights_input)
        uw = [];
        if isstruct(user_weights_input) && isfield(user_weights_input, 'data'), uw = abs(max(user_weights_input.data, [], 2));
        elseif isnumeric(user_weights_input) && isvector(user_weights_input), uw = abs(user_weights_input(:)); end
        if ~isempty(uw) && length(uw) == S_loc
            if max(uw) > eps, uw = uw / max(uw); end
            if O == 1, source_weighting_init = source_weighting_init .* uw; else, source_weighting_init = source_weighting_init .* repelem(uw, O); end
        else, warning('seal_irMxNE:UserWeightsDim', 'User weights ignored.'); 
        end
    end
    Param_irMxNE.InitialSourceWeightingForGain = source_weighting_init;
    G_eff_base = bsxfun(@times, G_w, source_weighting_init'); 

    mask_active_sources = true(Nsources_total, 1);
    if user_weights_min > 0 && ~isempty(user_weights_input) 
        mask_active_sources = source_weighting_init > (user_weights_min * (max(source_weighting_init)+eps) );
        if ~any(mask_active_sources), error('seal_irMxNE:NoSourcesAfterWeightMin', 'No sources after WeightsMin.'); end
        G_eff_base = G_eff_base(:, mask_active_sources);
        source_weighting_init_masked = source_weighting_init(mask_active_sources); % Keep track of these for reapplication
        fprintf('irMxNE: Reduced source space to %d component(s) due to WeightsMin.\n', sum(mask_active_sources));
    else
        source_weighting_init_masked = source_weighting_init;
    end
    Param_irMxNE.ActiveSourceMask_Initial = mask_active_sources;
    
    %% 2. Scaling for Solver (based on initial G_eff_base)
    GtM = G_eff_base' * M_processed;
    n_pos_eff = size(GtM,1)/O;
    l2inf_norms_sq_per_group = zeros(n_pos_eff, 1);
    for i = 1:n_pos_eff, idx = (i-1)*O + (1:O); l2inf_norms_sq_per_group(i) = sum(GtM(idx,:).^2, 'all'); end
    alpha_max_val = sqrt(max(l2inf_norms_sq_per_group));
    alpha_max_val = alpha_max_val * 0.01; 
    if alpha_max_val < eps, alpha_max_val = 1.0; end

    G_scaled_base = G_eff_base / alpha_max_val; % This G is used in the reweighting loop
    alpha_solver = alpha_user; 
    Param_irMxNE.AlphaSolver = alpha_solver; 
    Param_irMxNE.AlphaMaxScalingFactor = alpha_max_val;

    %% 3. Iterative Reweighted MxNE Solver
    fprintf('irMxNE: Starting iterative reweighted L1/L2 solution (%d iterations)...\n', n_irmxne_iter);
    
    X_active_iter_pca = zeros(size(G_scaled_base,2), size(M_processed,2)); % Initialize X for active sources
    current_reweighting_weights = ones(size(G_scaled_base,2), 1); % Initial reweights are 1 (w_k in MNE-Python)
    active_set_bcd = false(size(G_scaled_base,2),1); % Active set from BCD

    for k_iter = 1:n_irmxne_iter
        fprintf('  irMxNE Iteration %d/%d...\n', k_iter, n_irmxne_iter);
        X_prev_iter_pca = X_active_iter_pca; % Store for convergence check
        
        % Apply current reweighting weights to G_scaled_base
        % G_iter = G_scaled_base * diag(current_reweighting_weights)
        % (MNE-Python: G_tmp = G[:, active_set] * weights[np.newaxis, :])
        % Here, weights are applied to columns of G.
        G_iter_solver = bsxfun(@times, G_scaled_base, current_reweighting_weights');

        [X_solved_at_iter, active_set_at_iter_mask, ~] = local_mxne_bcd_solver(...
            M_processed, G_iter_solver, alpha_solver, O, ...
            solver_max_iter, solver_tol, max(1,round(solver_max_iter/100)));
        
        % Update X_active_iter_pca: solution from BCD is for the reweighted problem.
        % X_solved_at_iter is for G_iter_solver.
        % The "true" X for this iteration is X_solved_at_iter .* current_reweighting_weights(active_set_at_iter_mask)
        X_active_iter_pca_new = zeros(size(G_scaled_base,2), size(M_processed,2));
        if any(active_set_at_iter_mask)
            X_active_iter_pca_new(active_set_at_iter_mask,:) = bsxfun(@times, X_solved_at_iter, current_reweighting_weights(active_set_at_iter_mask));
        end
        X_active_iter_pca = X_active_iter_pca_new;
        active_set_bcd = active_set_at_iter_mask; % Update active set from BCD

        if k_iter < n_irmxne_iter % Don't update weights on the last iteration
            if ~any(active_set_bcd)
                fprintf('    No active sources in irMxNE iteration %d. Stopping reweighting.\n', k_iter);
                break; 
            end
            % Update reweighting_weights for next iteration based on current X_active_iter_pca
            % MNE-Python: weights = gprime(X_current_solution_for_active_set)
            % gprime(w) = 2.0 * repeat(g(w), n_orient).ravel()
            % g(w) = sqrt(sqrt(groups_norm2(w, n_orient)))
            % groups_norm2 sums squares over time and orientations for each location
            
            X_for_gprime = X_active_iter_pca(active_set_bcd,:); % Use only active X for gprime
            if isempty(X_for_gprime)
                 new_weights_for_active = ones(sum(active_set_bcd),1); % Should not happen if break above
            else
                group_L2_norms_sq = local_groups_norm2(X_for_gprime, O, size(X_for_gprime,1)/O); % N_active_pos x 1
                g_X = sqrt(sqrt(group_L2_norms_sq)); % N_active_pos x 1
                g_prime_X_loc = 2.0 * g_X; % N_active_pos x 1
                
                % The weights w_k for G_k w_k X_k should be 1/g'(X_k)
                % MNE-Python applies weights = gprime(X) to G, then solves for X_new.
                % Then X_true = X_new * weights. This implies an L0.5 like penalty.
                % Here, we want weights w_k such that G_eff_k * w_k is used in solver.
                % The penalty is alpha * sum || X_k_solved ||_2.
                % For irLasso, w_k = 1 / (|x_k_prev| + epsilon).
                % For L0.5/L2, penalty is sum sqrt(||X_k||_2). Reweighting uses w_k = 1 / (2*sqrt(||X_k_prev||_2) + epsilon)
                % MNE-Python's gprime is 2 * (||X_k_prev||_2)^0.25.
                % If we want penalty sum w_k ||X_k||_2 to approximate sum (||X_k||_2)^0.5, then w_k ~ (||X_k_prev||_2)^-0.5
                
                new_reweights_loc = zeros(size(g_prime_X_loc));
                valid_gprime = g_prime_X_loc > eps;
                new_reweights_loc(valid_gprime) = 1 ./ g_prime_X_loc(valid_gprime);
                new_reweights_loc(~valid_gprime) = 1 / eps; % Max weight if gprime is zero
                
                % Expand to components
                temp_current_reweighting_weights = ones(size(G_scaled_base,2), 1);
                active_indices_in_G_scaled_base = find(active_set_bcd);
                
                map_active_pos_to_all_pos = ceil( (1:sum(active_set_bcd)) / O); % maps from active_set_bcd iteration to original position index within active set
                active_pos_indices_global = unique(ceil(active_indices_in_G_scaled_base/O)); % Global positions that are active

                temp_weights_for_active_pos = zeros(length(active_pos_indices_global),1);
                for i_ap = 1:length(active_pos_indices_global)
                    % find which entry in new_reweights_loc corresponds to this global active pos
                    % This logic is getting complicated. MNE-Python applies weights to G_active.
                    % Let's simplify: weights are for the current active set.
                    if O > 1
                        temp_current_reweighting_weights(active_set_bcd) = repelem(new_reweights_loc, O);
                    else
                        temp_current_reweighting_weights(active_set_bcd) = new_reweights_loc;
                    end
                end
                current_reweighting_weights = temp_current_reweighting_weights;
            end
        end
        
        % Check convergence of X across outer iterations
        if k_iter > 1
            change_X = norm(X_active_iter_pca(:) - X_prev_iter_pca(:),'fro') / (norm(X_prev_iter_pca(:),'fro') + eps);
            fprintf('    irMxNE Iter %d: X change = %.2e (Tol=%.1e)\n', k_iter, change_X, irmxne_tol);
            if change_X < irmxne_tol
                fprintf('  irMxNE converged (X change) at iteration %d.\n', k_iter);
                break;
            end
        end
    end % End irMxNE iterations
    if k_iter == n_irmxne_iter, warning('seal_irMxNE:OuterLoopMaxIter', 'irMxNE reweighting loop reached max iterations.'); end
    Param_irMxNE.NumIRMXNEIterActual = k_iter;
    
    X_solved_active_pca = X_active_iter_pca(active_set_bcd,:); % Final solution from BCD on active set
    final_active_set_mask_in_G_solver_space = active_set_bcd;

    %% 4. Debias (if requested)
    X_debiased_active_pca = X_solved_active_pca;
    if perform_debias && sum(final_active_set_mask_in_G_solver_space) > 0
        fprintf('irMxNE: Performing debiasing...\n');
        G_for_debias = G_solver(:, final_active_set_mask_in_G_solver_space);
        debiasing_weights = local_compute_bias(M_processed, G_for_debias, X_solved_active_pca, O, ...
                                               debias_max_iter, debias_tol);
        X_debiased_active_pca = X_solved_active_pca .* debiasing_weights;
        Param_irMxNE.DebiasingWeights = debiasing_weights;
    end

    %% 5. Reconstruct full time series if PCA was used
    X_debiased_active_full_time = X_debiased_active_pca;
    if ~isempty(Vh_pca)
        X_debiased_active_full_time = X_debiased_active_pca * Vh_pca;
    end

    %% 6. Reapply scaling and original source weighting
    % The X_debiased_active_full_time is for G_solver = G_eff_base / alpha_max_val
    % And G_eff_base = G_w * diag(source_weighting_init_masked)
    % The reweighting weights were applied inside the loop to G_scaled_base.
    % The final X_debiased_active_full_time needs to be scaled by source_weighting_init_masked
    % and then by 1/alpha_max_val to refer to sources of G_w.
    
    active_indices_in_masked_space = find(final_active_set_mask_in_G_solver_space);
    sw_init_for_final_active = source_weighting_init_masked(active_indices_in_masked_space);

    X_final_active = bsxfun(@times, X_debiased_active_full_time, sw_init_for_final_active);
    X_final_active = X_final_active / alpha_max_val;

    %% 7. Format Output STC
    X_full_components = zeros(Nsources_total, Nt);
    temp_full_masked = zeros(sum(mask_active_sources), Nt); % Size of G_eff_base columns
    
    % Map active_indices_in_masked_space (which are indices into G_solver/G_eff_base)
    % to indices in temp_full_masked.
    % active_set_solver_mask was relative to G_solver.
    % final_active_set_mask_in_G_solver_space is this mask.
    temp_full_masked(final_active_set_mask_in_G_solver_space, :) = X_final_active;
    X_full_components(mask_active_sources, :) = temp_full_masked;

    if O > 1 && strcmpi(pick_ori_opt, 'none')
        Source_irMxNE_loc = zeros(S_loc, Nt);
        active_comp_indices_final = find(any(X_full_components,2));
        active_locs_final_unique = unique(ceil(active_comp_indices_final/O));
        Param_irMxNE.ActiveSetIndices_Locations = active_locs_final_unique;
        Param_irMxNE.ActiveSetIndices_Components = active_comp_indices_final;

        for i_loc = 1:S_loc
            if ismember(i_loc, active_locs_final_unique)
                idx_comps = (i_loc-1)*O + (1:O);
                Source_irMxNE_loc(i_loc, :) = sqrt(sum(X_full_components(idx_comps,:).^2, 1));
            end
        end
        Source_irMxNE = Source_irMxNE_loc;
    else 
        Source_irMxNE = X_full_components;
        Param_irMxNE.ActiveSetIndices_Components = find(any(X_full_components,2));
        if O > 1, Param_irMxNE.ActiveSetIndices_Locations = unique(ceil(Param_irMxNE.ActiveSetIndices_Components/O));
        else, Param_irMxNE.ActiveSetIndices_Locations = Param_irMxNE.ActiveSetIndices_Components; end
    end
    
    Param_irMxNE.ResolutionMatrixSource = NaN; 
end

%% --- Local Helper Functions for MxNE Solver (Copied from seal_MxNE, can be shared) ---
function [X_solved_active, active_set_mask, E_path] = local_mxne_bcd_solver(...
    M_data, G_solver, alpha_penalty_val, O_orient, max_iter, tol, dgap_freq)
% (Identical to local_mxne_bcd_solver in seal_MxNE_v1)
    [~, Nt_proc] = size(M_data);
    [~, Ns_proc] = size(G_solver);
    N_pos_proc = Ns_proc / O_orient;

    X = zeros(Ns_proc, Nt_proc); R = M_data; 
    active_set_mask = false(Ns_proc, 1); 
    E_path = []; highest_d_obj = -Inf;

    lc = zeros(N_pos_proc, 1);
    if O_orient == 1, lc = sum(G_solver.^2, 1)'; 
    else
        for j_pos = 1:N_pos_proc
            idx = (j_pos-1)*O_orient + (1:O_orient); G_j_block = G_solver(:, idx);
            lc(j_pos) = norm(G_j_block, 'fro')^2;
            if lc(j_pos) < eps, lc(j_pos) = eps; end 
        end
    end
    alpha_scaled_per_block = alpha_penalty_val ./ lc;

    for iter_bcd = 1:max_iter
        X_old_iter = X; 
        for j_pos = 1:N_pos_proc 
            idx_block = (j_pos-1)*O_orient + (1:O_orient); G_j = G_solver(:, idx_block); 
            X_j_current = X(idx_block, :); 
            if any(active_set_mask(idx_block)), R = R + G_j * X_j_current; end
            X_j_unthres = X_j_current + (1/lc(j_pos)) * (G_j' * R);
            group_norm_L2F = norm(X_j_unthres, 'fro');
            if group_norm_L2F <= alpha_scaled_per_block(j_pos)
                X_j_new = zeros(O_orient, Nt_proc); active_set_mask(idx_block) = false;
            else
                shrink_factor = 1 - alpha_scaled_per_block(j_pos) / group_norm_L2F;
                X_j_new = X_j_unthres * shrink_factor; active_set_mask(idx_block) = true;
            end
            if any(active_set_mask(idx_block)) || any(X_j_current(:) ~= 0), R = R - G_j * (X_j_new - X_j_current); end
            X(idx_block, :) = X_j_new;
        end 
        if mod(iter_bcd, dgap_freq) == 0 || iter_bcd == max_iter || iter_bcd == 1
            current_active_comp_indices = find(active_set_mask);
            if isempty(current_active_comp_indices), p_obj = 0.5*sum(M_data(:).^2); d_obj = 0;
            else
                G_curr_active_comps = G_solver(:, current_active_comp_indices);
                X_curr_active_comps = X(current_active_comp_indices, :);
                [~, p_obj, d_obj, ~] = local_dgap_l21(M_data, G_curr_active_comps, X_curr_active_comps, alpha_penalty_val, O_orient, active_set_mask);
            end
            E_path(end+1) = p_obj; highest_d_obj = max(highest_d_obj, d_obj); actual_gap = p_obj - highest_d_obj;
            if iter_bcd > 1 && actual_gap < tol, break; end
        end
         if iter_bcd > 1 && norm(X(:)-X_old_iter(:),'fro')/(norm(X_old_iter(:),'fro')+eps) < tol/100, break; end
    end 
    if iter_bcd == max_iter, warning('seal_irMxNE:BCDSolverMaxIter', 'BCD solver reached max iterations.'); end
    X_solved_active = X(active_set_mask, :);
end

function [gap, p_obj, d_obj, R_residual] = local_dgap_l21(M_data, G_active_comp, X_active_comp, alpha_penalty_val, O_orient, ~) % active_set_mask_full not used here
% (Identical to local_dgap_l21 in seal_MxNE_v1)
    if isempty(X_active_comp) || isempty(G_active_comp)
        R_residual = M_data; nR2 = sum(R_residual(:).^2); p_obj = 0.5*nR2; d_obj = 0; gap = p_obj;
        if isnan(gap), gap = p_obj; end; return;
    end
    GX = G_active_comp * X_active_comp; R_residual = M_data - GX; nR2 = sum(R_residual(:).^2);
    penalty_l21 = 0; N_pos_active = size(X_active_comp,1)/O_orient;
    for i=1:N_pos_active, idx=(i-1)*O_orient+(1:O_orient); penalty_l21=penalty_l21+norm(X_active_comp(idx,:),'fro'); end
    p_obj = 0.5*nR2 + alpha_penalty_val*penalty_l21;
    Gt_R = G_active_comp'*R_residual; dual_group_norms_F = zeros(N_pos_active,1);
    for i=1:N_pos_active, idx=(i-1)*O_orient+(1:O_orient); dual_group_norms_F(i)=norm(Gt_R(idx,:),'fro'); end
    dual_norm_l2inf = max(dual_group_norms_F); if isempty(dual_norm_l2inf), dual_norm_l2inf=0; end
    scaling = min(1.0, alpha_penalty_val/(dual_norm_l2inf+eps));
    d_obj = (scaling - 0.5*(scaling^2))*nR2 + scaling*sum(sum(R_residual.*GX));
    gap = p_obj - d_obj; if isnan(gap) && ~isnan(p_obj), gap=p_obj; end
end

function D_debias = local_compute_bias(M_data, G_active_comp, X_active_comp, O_orient, max_iter_debias, tol_debias)
% (Identical to local_compute_bias in seal_MxNE_v1)
    [N_active_comp, ~] = size(X_active_comp); N_pos_active = N_active_comp / O_orient;
    Dk_comp = ones(N_active_comp, 1); Yk_comp = ones(N_active_comp, 1); tk = 1.0;
    L_fista_val = 0;
    for i=1:N_pos_active, idx=(i-1)*O_orient+(1:O_orient); G_b=G_active_comp(:,idx); X_b=X_active_comp(idx,:); Y_j=G_b*X_b; L_fista_val=max(L_fista_val,norm(Y_j,'fro')^2); end
    if L_fista_val < eps, L_fista_val=1.0; else, L_fista_val=1.1*L_fista_val; end
    for iter_d = 1:max_iter_debias
        D0_comp = Dk_comp;
        M_pred = G_active_comp * bsxfun(@times, X_active_comp, Yk_comp); Residual = M_data - M_pred; 
        Gt_Res = G_active_comp' * Residual; Grad_D_comp = -sum(Gt_Res .* X_active_comp, 2);
        D_update = Yk_comp - Grad_D_comp / L_fista_val;
        if O_orient > 1
            D_reshaped_loc = mean(reshape(D_update, O_orient, N_pos_active), 1)';
            D_reshaped_loc = max(D_reshaped_loc, 1.0); Dk_comp = repelem(D_reshaped_loc, O_orient);
        else, Dk_comp = max(D_update, 1.0); 
        end
        t_new = 0.5*(1.0+sqrt(1.0+4.0*tk^2)); Yk_comp = Dk_comp + ((tk-1.0)/t_new)*(Dk_comp-D0_comp); tk = t_new;
        if max(abs(Dk_comp - D0_comp)) < tol_debias && iter_d > 1, break; end
    end
    if iter_d == max_iter_debias, warning('seal_irMxNE:DebiasMaxIter', 'Debiasing did not converge.'); end
    D_debias = Dk_comp;
end

function local_check_option(key, value, allowed)
    if ~ismember(value, allowed)
        error('Invalid value for %s. Allowed options are: %s', key, strjoin(allowed, ', '));
    end
end

function group_L2_norms_sq_val = local_groups_norm2(A_matrix, n_orient_val, n_positions_val)
% Compute squared L2 norms of groups (Frobenius norm for each block over time)
% A_matrix is (n_positions_val*n_orient_val) x Ntime
% Output is n_positions_val x 1
    group_L2_norms_sq_val = zeros(n_positions_val, 1);
    for i_pos = 1:n_positions_val
        idx_block = (i_pos-1)*n_orient_val + (1:n_orient_val);
        block_A = A_matrix(idx_block, :); % n_orient_val x Ntime
        group_L2_norms_sq_val(i_pos) = sum(block_A(:).^2); % Squared Frobenius norm of block
    end
end
