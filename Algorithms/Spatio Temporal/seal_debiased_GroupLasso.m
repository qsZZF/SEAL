function [Source_GL, Param_GL] = seal_debiased_GroupLasso(Data, L, varargin)
%SEAL_GROUPLASSO Computes Spatio-Temporal Group LASSO via Dual ADMM.
%   [Source, Param] = seal_GroupLasso(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Translates the advanced python implementation of gl_ADMM_dual_joint.
%   Features block-wise leadfield whitening, dynamic lambda estimation,
%   Spatio-Temporal L2,1 dual projection, and OLS debiasing.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints measured data.
%       L (double matrix): Nchannels x Nsources lead field.
%
%   Optional Name-Value Pair Inputs:
%       'NumOrientations' (integer): Orientations per location. Default: 1.
%       'RegularizationParameter' (double): Base Lambda. Default: 1.0.
%       'DynamicLambda' (logical): If true, scales lambda iteratively based 
%                                  on residual noise. Default: true.
%       'JointEstNum' (integer): Max loops for Dynamic Lambda. Default: 10.
%       'Debias' (logical): Run post-selection OLS bias correction. Default: true.
%       'MaxIterations' (integer): Max inner ADMM iterations. Default: 1000.
%       'Tolerance' (double): ADMM convergence tolerance. Default: 1e-4.

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    
    addParameter(p, 'NumOrientations', 1, @isscalar);
    addParameter(p, 'RegularizationParameter', 0.1, @isscalar);
    addParameter(p, 'DynamicLambda', true, @islogical);
    addParameter(p, 'JointEstNum', 10, @isscalar);
    addParameter(p, 'Debias', true, @islogical);
    addParameter(p, 'MaxIterations', 1000, @isscalar);
    addParameter(p, 'Tolerance', 1e-4, @isscalar);
    addParameter(p, 'MaxIterationsDebias', 300, @isscalar);
    addParameter(p, 'ToleranceDebias', 1e-5, @isscalar);

    parse(p, Data, L, varargin{:});
    opts = p.Results;

    nd = opts.NumOrientations;
    lam_base = opts.RegularizationParameter;
    max_iter = opts.MaxIterations;
    tol = opts.Tolerance;
    max_iter_debias = opts.MaxIterationsDebias;
    tol_debias = opts.ToleranceDebias;
    
    [Nc, T] = size(Data);
    S = size(L, 2) / nd; % Number of anatomical locations

    %% 2. Block-wise Leadfield Whitening 
    % Standardizes the orientations to be isotropic at every location
    fprintf('Whitening leadfield blocks (nd=%d)...\n', nd);
    L_tilde = zeros(size(L));
    sigma_inv_sqrt = zeros(nd, nd, S);
    
    for s = 1:S
        idx = (s-1)*nd + 1 : s*nd;
        Ls = L(:, idx);
        
        % Block covariance
        cov_s = (Ls' * Ls) / Nc; 
        
        % Inverse Square Root via Eigendecomposition
        [V, D] = eig(cov_s);
        d = diag(D);
        d(d < 1e-12) = 1e-12; % Prevent division by zero
        W_s = V * diag(1./sqrt(d)) * V';
        
        L_tilde(:, idx) = Ls * W_s;
        sigma_inv_sqrt(:,:,s) = W_s;
    end

    %% 3. Dynamic Lambda Outer Loop 
    if opts.DynamicLambda
        joint_loops = opts.JointEstNum;
    else
        joint_loops = 1;
    end
    
    Xt = zeros(S*nd, T);
    sigma_hat_pre = inf;
    lam_new = lam_base * (norm(Data, 'fro') / sqrt(Nc*T));
    sigma_list = [];
    
    fprintf('Running Dual ADMM Group LASSO...\n');
    for j = 1:joint_loops
        lam_use = lam_new * sqrt(T); % Scale by time
        
        % --- Inner Dual ADMM Solver ---
        Xt = solve_dual_admm(Xt, L_tilde, Data, lam_use, nd, S, T, max_iter, tol);
        % ------------------------------
        
        % Update residual sigma and Lambda
        sigma_hat = norm(Data - L_tilde * Xt, 'fro') / sqrt(Nc*T);
        lam_new = lam_base * sigma_hat;
        sigma_list = [sigma_list; sigma_hat];
        
        fprintf('  Joint Iter %d: Sigma = %.4f, Lambda = %.4f\n', j, sigma_hat, lam_use);
        
        if abs(sigma_hat - sigma_hat_pre) < 1e-4
            break;
        end
        sigma_hat_pre = sigma_hat;
    end

    %% 4. Un-whiten and Identify Active Set
    X_out = zeros(S*nd, T);
    for s = 1:S
        idx = (s-1)*nd + 1 : s*nd;
        X_out(idx, :) = sigma_inv_sqrt(:,:,s) * Xt(idx, :);
    end

    % Calculate Spatio-Temporal L2 norm per location
    X_out_tensor = reshape(X_out, nd, S, T);
    X_norms_loc = sqrt(sum(sum(X_out_tensor.^2, 1), 3)); % [1 x S x 1]
    
    % Mask non-zero active locations
    active_mask = X_norms_loc(:) > 1e-6;
    active_locs = find(active_mask);
    
    fprintf('Identified %d active source locations.\n', length(active_locs));

    %% 5. Post-Selection Debiasing 
    if opts.Debias && ~isempty(active_locs)
        fprintf('Applying OLS Bias Correction...\n');
        
        % Map active locations back to the 3N columns
        active_cols = zeros(1, length(active_locs)*nd);
        for a = 1:length(active_locs)
            s = active_locs(a);
            active_cols((a-1)*nd + 1 : a*nd) = (s-1)*nd + 1 : s*nd;
        end
        
        G_active_comp = L(:, active_cols);      % Nc x N_active_comp
        X_active_comp = X_out(active_cols, :);     % N_active_comp x T

        D_debias = local_compute_bias( ...
            Data, G_active_comp, X_active_comp, nd, ...
            max_iter_debias, tol_debias);

        X_debias = zeros(S * nd, T);
        X_debias(active_cols, :) = bsxfun(@times, X_active_comp, D_debias);

        Source_GL = X_debias;
    else
        Source_GL = X_out;
        D_debias = [];
        active_cols = [];
    end

    %% 6. Output Assignment
    Param_GL.RawLassoSource = X_out;
    Param_GL.ActiveLocations = active_locs;
    Param_GL.FinalSigma = sigma_list(end);
    Param_GL.SigmaHistory = sigma_list;
    Param_GL.OptionsPassed = opts;
    Param_GL.DebiasWeights = D_debias;
end

function Xt_full = solve_dual_admm(Xt, L_tilde, Y, lam_use, nd, S, T, max_iter, tol)
% Core ADMM Solver for the Dual Spatio-Temporal Group LASSO Problem
    
    Nc = size(Y, 1);
    stepsize = 1 / (S * nd);
    max_stepsize_tol = stepsize * sqrt(Nc*S*nd);
    min_stepsize_tol = stepsize / sqrt(Nc*S*nd);
    % Woodbury identity matrices
    GGT = L_tilde * L_tilde';
    core = (eye(Nc) + stepsize * GGT);
    coreinvG = core \ L_tilde;
    coreinvY = core \ Y;

    Zt = Y - L_tilde * Xt;
    Ut = L_tilde' * Zt;
    sqrtSOT = sqrt(S*nd*T);
    sqrtNT = sqrt(Nc*T);
    
    utmask = 1:S*nd;
    for iter = 1:max_iter
        Ut_old = Ut;
        % Update Z
        grad_Ut = stepsize * Ut;
        grad_Ut(utmask,:) = grad_Ut(utmask,:) - Xt;
        Zt = coreinvY + coreinvG * grad_Ut;
        GTZ = L_tilde' * Zt;

        % Update U (Spatio-Temporal L2,1 Dual Projection)
        Ut = GTZ;
        if any(utmask)
            Ut(utmask, :) = GTZ(utmask, :) + Xt / stepsize;
        end
        % Reshape to pool norm across orientations AND time
        group_norm_U = compute_group_row_norm(Ut, nd, S, T);
        utmask_new = group_norm_U > lam_use;
        % Project onto L2 ball bounded by lambda
        if any(utmask_new)
            Ut(utmask_new, :) = lam_use .* Ut(utmask_new, :) ./ group_norm_U(utmask_new);
        end
        % Update X
        X_new = stepsize * (GTZ(utmask_new, :) - Ut(utmask_new, :));

        overlap_global = find(utmask(:) & utmask_new(:));
        if ~isempty(overlap_global)
            prev_global = find(utmask);
            new_global = find(utmask_new);
            [~, overlap_in_prev] = ismember(overlap_global, prev_global);
            [~, overlap_in_new]  = ismember(overlap_global, new_global);
            X_new(overlap_in_new, :) = Xt(overlap_in_prev, :) + X_new(overlap_in_new, :);
        end

        Xt = X_new;
        utmask = utmask_new;
        if any(utmask)
            LX = L_tilde(:, utmask) * Xt;
        else
            LX = zeros(Nc, T);
        end
        % Check convergence
        loss_p = norm(GTZ - Ut, 'fro');
        loss_d = norm(stepsize * L_tilde * (Ut - Ut_old), 'fro');
        tol_p = tol * sqrtSOT + tol * max(norm(GTZ, 'fro'), norm(Ut, 'fro'));
        tol_d = tol * sqrtNT + tol * norm(LX, 'fro');
        if mod(iter, 200) == 0
            fprintf(' iter = %d, loss_p = %.6g, loss_d = %.6g, stepsize = %.6g, active_source = %d\n', ...
                iter, loss_p, loss_d, stepsize, sum(utmask)/nd);
        end
        if (loss_p < tol_p) && (loss_d < tol_d)
            % fprintf('ADMM converged at iter %d\n', iter);
            break;
        end
        if loss_d / sqrtNT > 10 * loss_p / sqrtSOT
            stepsize = max(stepsize / 2, min_stepsize_tol);
            core = (eye(Nc) + stepsize * GGT);
            coreinvG = core \ L_tilde;
            coreinvY = core \ Y;
        elseif loss_p / sqrtSOT > 10 * loss_d / sqrtNT
            stepsize = min(stepsize * 2, max_stepsize_tol);
            core = (eye(Nc) + stepsize * GGT);
            coreinvG = core \ L_tilde;
            coreinvY = core \ Y;
        end
    end
    Xt_full = zeros(S * nd, T);
    Xt_full(utmask, :) = Xt;
end

function D_debias = local_compute_bias(M_data, G_active_comp, X_active_comp, ...
                                       O_orient, max_iter_debias, tol_debias)
%LOCAL_COMPUTE_BIAS
% Solves the debiasing problem:
%
%   min_D  0.5 * || M_data - G_active_comp * diag(D_comp) * X_active_comp ||_F^2
%
% subject to:
%   the average debiasing weight of each active source location >= 1
%
% Here:
%   - D_debias is defined on active components (N_active_comp x 1)
%   - but the constraint is enforced at the source-location level
%   - for O_orient > 1, one location shares one averaged debias factor,
%     then that factor is repeated across its O_orient components.
%
% Inputs:
%   M_data          : Nc x Nt
%   G_active_comp   : Nc x N_active_comp
%   X_active_comp   : N_active_comp x Nt
%   O_orient        : number of orientations per source location
%   max_iter_debias : maximum FISTA iterations
%   tol_debias      : stopping tolerance
%
% Output:
%   D_debias        : N_active_comp x 1 debiasing weights

    [N_active_comp, ~] = size(X_active_comp); 

    if mod(N_active_comp, O_orient) ~= 0
        error('N_active_comp must be divisible by O_orient.');
    end

    N_pos_active = N_active_comp / O_orient;

    % Initialization
    Dk_comp = ones(N_active_comp, 1);
    Yk_comp = ones(N_active_comp, 1);
    tk = 1.0;

    %% Estimate Lipschitz constant
    % Heuristic but readable/stable form:
    % use the maximum || G_block * X_block ||_F^2 over active locations,
    % then slightly enlarge it.
    L_fista_val = 0;

    for i_pos = 1:N_pos_active
        idx_block = (i_pos - 1) * O_orient + (1:O_orient);

        G_block = G_active_comp(:, idx_block);   % Nc x O
        X_block = X_active_comp(idx_block, :);   % O x Nt

        Y_j_eff = G_block * X_block;             % Nc x Nt
        L_fista_val = max(L_fista_val, norm(Y_j_eff, 'fro')^2);
    end

    if L_fista_val < eps
        L_fista_val = 1.0;
    else
        L_fista_val = 1.1 * L_fista_val;
    end

    fprintf('  Debiasing (FISTA): MaxIter=%d, Tol=%.1e, L_fista=%.2e\n', ...
        max_iter_debias, tol_debias, L_fista_val);

    %% FISTA iterations
    for iter_d = 1:max_iter_debias
        D0_comp = Dk_comp;

        % Prediction under current extrapolated weights
        % M_pred = G_active_comp * diag(Yk_comp) * X_active_comp
        M_pred = G_active_comp * bsxfun(@times, X_active_comp, Yk_comp);

        Residual = M_data - M_pred;

        % Gradient w.r.t. each component-wise debias weight
        % Grad_D_i = - sum_t (G(:,i)' * Residual(:,t)) * X(i,t)
        Gt_Res = G_active_comp' * Residual;                  % N_active_comp x Nt
        Grad_D_comp = -sum(Gt_Res .* X_active_comp, 2);      % N_active_comp x 1

        % Gradient step
        D_update = Yk_comp - Grad_D_comp / L_fista_val;

        % Prox / projection step:
        % average debias factor at each location must be >= 1
        if O_orient > 1
            D_loc_avg = mean(reshape(D_update, O_orient, N_pos_active), 1)';  % N_pos_active x 1
            D_loc_avg = max(D_loc_avg, 1.0);
            Dk_comp = repelem(D_loc_avg, O_orient);
        else
            Dk_comp = max(D_update, 1.0);
        end

        % FISTA extrapolation
        t_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * tk^2));
        Yk_comp = Dk_comp + ((tk - 1.0) / t_new) * (Dk_comp - D0_comp);
        tk = t_new;

        % Convergence check
        d_diff = max(abs(Dk_comp - D0_comp));
        if d_diff < tol_debias && iter_d > 1
            fprintf('    Debiasing converged after %d iterations (Diff=%.1e)\n', ...
                iter_d, d_diff);
            break;
        end
    end

    if iter_d == max_iter_debias
        warning('seal_GroupLasso:DebiasMaxIter', ...
            'Debiasing did not converge within the maximum iterations.');
    end

    D_debias = Dk_comp;
end

function row_norm = compute_group_row_norm(U, nd, S, T)
% Returns a per-row norm vector used for block masking.
% For nd > 1, the same group norm is repeated nd times.

    if nd == 1
        row_norm = vecnorm(U, 2, 2);
    else
        U_tensor = reshape(U, nd, S, T);
        group_norm = sqrt(sum(sum(U_tensor.^2, 1), 3));   % 1 x S x 1
        group_norm = group_norm(:);                       % S x 1
        row_norm = repelem(group_norm, nd, 1);           % (S*nd) x 1
    end
end