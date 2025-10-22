function [Source_DeESI, Param_DeESI] = seal_DeESI(Data, L, varargin)
%SEAL_DEESI Computes De-biased ESI using ADMM for Group LASSO.
%   [Source_DeESI, Param_DeESI] = seal_DeESI(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Translates the Python ADMM_GroupLASSO.py script into the SEAL
%   framework. This function implements a three-stage process:
%   1. Jointly estimates the noise level and a raw source estimate using a
%      wrapper around the core ADMM solver.
%   2. Identifies the active set of sources.
%   3. Performs a de-biasing step on the active set to correct for
%      amplitude shrinkage.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data (Y).
%       L (double matrix): Nchannels x Nsources original lead field matrix (G).
%
%   Optional Name-Value Pair Inputs:
%       'Lambda' (double): Base regularization parameter. Default: 75.
%       'NumOrientations' (integer): O. Default: 3.
%       'NumJointIter' (integer): Max iterations for joint noise/source estimation. Default: 50.
%       'SolverMaxIter' (integer): Max iterations for inner ADMM solver. Default: 1000.
%       'SolverTolerance' (double): Tol for inner ADMM solver. Default: 1e-8.
%       'JointEstTolerance' (double): Tol for joint estimation convergence. Default: 1e-4.
%       'BiasCorrection' (logical): Perform debiasing stage. Default: true.
%       'InitialSourceEstimate' (matrix): X0. Default: zeros.
%
%   Outputs:
%       Source_DeESI (double matrix): Final estimated source activity.
%       Param_DeESI (struct): Parameters and intermediate results.

    %% Input Parsing
    p = inputParser;
    p.CaseSensitive = false;
    
    [Nchan, Nsources_total] = size(L);
    [~, Ntime] = size(Data);

    defaultLambda = 75.0;
    defaultNumOrientations = 3;
    defaultNumJointIter = 50;
    defaultSolverMaxIter = 1000;
    defaultSolverTol = 1e-8;
    defaultJointEstTol = 1e-4;
    defaultBiasCorrection = true;
    defaultInitialSourceEstimate = zeros(Nsources_total, Ntime);

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    
    addParameter(p, 'Lambda', defaultLambda, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NumJointIter', defaultNumJointIter, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'SolverMaxIter', defaultSolverMaxIter, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'SolverTolerance', defaultSolverTol, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'JointEstTolerance', defaultJointEstTol, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'BiasCorrection', defaultBiasCorrection, @islogical);
    addParameter(p, 'InitialSourceEstimate', defaultInitialSourceEstimate, @(x) isnumeric(x) && ismatrix(x));

    parse(p, Data, L, varargin{:});
    Param_DeESI.OptionsPassed = p.Results;
    
    % Unpack parameters
    Y = p.Results.Data;
    G = p.Results.L;
    lam = p.Results.Lambda;
    O = p.Results.NumOrientations;
    joint_est_num = p.Results.NumJointIter;
    max_iter = p.Results.SolverMaxIter;
    tol = p.Results.SolverTolerance;
    joint_tol = p.Results.JointEstTolerance;
    perform_bias_correction = p.Results.BiasCorrection;
    X0 = p.Results.InitialSourceEstimate;

    %% --- Main DeESI Workflow ---
    
    % 1. Create whitening weights for the lead field
    fprintf('DeESI Stage 1: Calculating Lead Field Covariance Weights...\n');
    G_tensor = local_Gmatrix_to_tensor(G, O);
    [sigma_G_sqrt, sigma_G_sqrt_inv] = local_Gtensor_covariance(G_tensor);
    wlist = {sigma_G_sqrt, sigma_G_sqrt_inv};
    Param_DeESI.LeadFieldWeights = wlist;

    % 2. Run Joint Estimation for raw source estimate and data-driven lambda
    fprintf('DeESI Stage 2: Jointly estimating sources and noise level...\n');
    [X_raw, lam_final, sigma_list] = local_gl_ADMM_dual_joint(joint_est_num, X0, G, Y, lam, O, wlist, ...
                                                             'consecutive', joint_tol, tol, max_iter, true);
    Param_DeESI.RawEstimate = X_raw;
    Param_DeESI.FinalLambda = lam_final;
    Param_DeESI.SigmaTrace = sigma_list;
    
    if ~perform_bias_correction
        Source_DeESI = X_raw;
        fprintf('DeESI finished. Bias correction was skipped.\n');
        return;
    end
    
    % 3. Identify Active Set from Raw Estimate
    fprintf('DeESI Stage 3: Identifying active source set...\n');
    [I_hat_locs, ~, Xt_tensor] = local_cal_I_hat(X_raw, G, O, 'consecutive', 0);
    Param_DeESI.ActiveSetLocations = I_hat_locs;
    
    if isempty(I_hat_locs)
        fprintf('DeESI finished. No active sources found. Returning zero solution.\n');
        Source_DeESI = zeros(size(X_raw));
        return;
    end
    
    % 4. Perform Bias Correction on the Active Set
    fprintf('DeESI Stage 4: Performing bias correction on %d active locations...\n', length(I_hat_locs));
    [Xout_debias, ~, significance_list] = local_gl_ADMM_dual_bias_correction(I_hat_locs, Xt_tensor, G_tensor, Y, lam_final, O, wlist, ...
                                                                          'joint', true, 'consecutive', tol, tol, max_iter, true);
    
    Source_DeESI = Xout_debias;
    Param_DeESI.Significance = significance_list;
    fprintf('DeESI finished.\n');
end

%% --- Local Function Translations from ADMM_GroupLASSO.py ---

function [Xt, lam_new, sigma_list] = local_gl_ADMM_dual_joint(joint_est_num, X0, G, Y, lam, O, wlist, block_method, tol, tol_norm, max_iter, varing_rho)
    [N, T] = size(Y);
    lam_new = lam * norm(Y - G * X0) / sqrt(N * T);
    Xt = X0;
    sigma_hat_pre = inf;
    sigma_list = [];
    
    for i = 1:joint_est_num
        lam_use = lam_new * sqrt(T);
        Xt = local_gl_ADMM_dual(Xt, G, Y, lam_use, O, wlist, block_method, tol_norm, tol_norm, max_iter, varing_rho);
        
        sigma_hat = norm(Y - G * Xt, 'fro') / sqrt(N * T);
        lam_new = lam * sigma_hat;
        sigma_list = [sigma_list, sigma_hat];
        
        if abs(sigma_hat - sigma_hat_pre) < tol
            break;
        else
            sigma_hat_pre = sigma_hat;
        end
    end
end

function G_tensor = local_Gmatrix_to_tensor(G, O, ~)
    N = size(G, 1);
    S = size(G, 2) / O;
    G_tensor = reshape(G, N, O, S);
end

function X_tensor = local_Xmatrix_to_tensor(X, O, ~)
    S = size(X, 1) / O;
    T = size(X, 2);
    X_tensor = reshape(X, O, S, T);
    X_tensor = permute(X_tensor, [2, 1, 3]);
end

function [I_hat, G_tensor, Xt_tensor] = local_cal_I_hat(Xt, G, O, block_method, threshold)
    S = size(G, 2) / O;
    G_tensor = local_Gmatrix_to_tensor(G, O, block_method);
    Xt_tensor = local_Xmatrix_to_tensor(Xt, O, block_method);
    
    block_norm_sq = sum(reshape(Xt, O, S, size(Xt, 2)).^2, [1, 3]);
    block_norm = sqrt(block_norm_sq);
    
    I_hat = find(block_norm > threshold); % 1-based indices
end

function [sigma_G_sqrt, sigma_G_sqrt_inv] = local_Gtensor_covariance(G_tensor, depth)
    if nargin < 2, depth = 1; end
    [N, O, S] = size(G_tensor);
    sigma_G = zeros(O, O, S);
    for s = 1:S
        sigma_G(:, :, s) = G_tensor(:, :, s)' * G_tensor(:, :, s) / N;
    end
    
    sigma_G_eigenvalues = zeros(O, S);
    sigma_G_eigenvectors = zeros(O, O, S);
    for s = 1:S
        [V, D] = eig(sigma_G(:, :, s));
        sigma_G_eigenvectors(:, :, s) = V;
        sigma_G_eigenvalues(:, s) = diag(D);
    end
    
    sigma_G_eigenvalues = sigma_G_eigenvalues.^depth;
    sigma_G_eigenvalues_sqrt = sqrt(sigma_G_eigenvalues);
    sigma_G_eigenvalues_sqrt_inv = 1 ./ sigma_G_eigenvalues_sqrt;
    
    sigma_G_sqrt = zeros(O, O, S);
    sigma_G_sqrt_inv = zeros(O, O, S);
    for s = 1:S
        V = sigma_G_eigenvectors(:, :, s);
        D_sqrt = diag(sigma_G_eigenvalues_sqrt(:, s));
        D_sqrt_inv = diag(sigma_G_eigenvalues_sqrt_inv(:, s));
        sigma_G_sqrt(:, :, s) = V * D_sqrt * V';
        sigma_G_sqrt_inv(:, :, s) = V * D_sqrt_inv * V';
    end
end

function [Xout_debias, Xout_debias_rotate, significance_list] = local_gl_ADMM_dual_bias_correction(I_hat_all, Xt_tensor, G_tensor, Y, lam, O, wlist, bias_correction_method, clear_not_select, block_method, tol, tol_norm, max_iter, varing_rho)
    S = size(G_tensor, 3);
    T = size(Y, 2);
    Xout_debias = Xt_tensor;
    I_hat_all_0based = I_hat_all - 1; % Convert to 0-based for setdiff equivalent
    all_indices_0based = 0:(S-1);
    neg_I_hat_0based = setdiff(all_indices_0based, I_hat_all_0based);
    neg_I_hat = neg_I_hat_0based + 1; % Back to 1-based
    
    significance_list = zeros(S, 1);
    
    sigma_G_sqrt = wlist{1};
    sigma_G_sqrt_inv = wlist{2};
    
    G_tilde = G_tensor;
    for s_idx = 1:length(neg_I_hat)
        s = neg_I_hat(s_idx);
        G_tilde(:, :, s) = G_tensor(:, :, s) * sigma_G_sqrt_inv(:, :, s);
    end
    
    if clear_not_select
        Xout_debias(neg_I_hat, :, :) = 0;
        Xt_tilde = Xout_debias;
        for s_idx = 1:length(neg_I_hat)
            s = neg_I_hat(s_idx);
            Xt_tilde(s, :, :) = permute(sigma_G_sqrt(:, :, s) * permute(Xout_debias(s, :, :), [2, 3, 1]), [3, 1, 2]);
        end
    else
        Xout_debias_rotate = Xout_debias;
        Xout_debias_rotate(neg_I_hat, :, :) = 0;
        Xt_tilde = Xt_tensor;
        for s_idx = 1:length(neg_I_hat)
            s = neg_I_hat(s_idx);
            Xt_tilde(s, :, :) = permute(sigma_G_sqrt(:, :, s) * permute(Xt_tensor(s, :, :), [2, 3, 1]), [3, 1, 2]);
        end
    end
    
    if strcmp(bias_correction_method, 'joint')
        [X_hat_debias_idx, X_hat_debias_rotate_idx, debias_flag, effect_mat] = ...
            local_single_bias_correction(I_hat_all, Xt_tilde, G_tilde, Y, lam, O, tol, tol_norm, max_iter, varing_rho);
        Xout_debias(I_hat_all, :, :) = reshape(X_hat_debias_idx, [length(I_hat_all), O, T]);
        % ... further processing for significance ...
    end
    
    Xout_debias = reshape(permute(Xout_debias, [2, 1, 3]), S * O, T);
    Xout_debias_rotate = Xout_debias; % Simplified for now
end

function [X_hat_debias, X_hat_debias_rotate, debias_flag, effect_mat] = local_single_bias_correction(I_hat, Xt_tilde, G_tilde, Y, lam, O, tol, tol_norm, max_iter, varing_rho)
    [N, ~, S] = size(G_tilde);
    T = size(Y, 2);
    debias_flag = false;
    X_hat_debias = [];
    X_hat_debias_rotate = [];
    effect_mat = 0;
    
    if ~isempty(I_hat) && length(I_hat) < N/O
        debias_flag = true;
        all_indices = 1:S;
        neg_I_hat = setdiff(all_indices, I_hat);
        
        G_I = G_tilde(:, :, I_hat);
        G_neg_I = G_tilde(:, :, neg_I_hat);
        X_neg_I = Xt_tilde(neg_I_hat, :, :);
        
        G_I_flat = reshape(G_I, N, []);
        G_neg_I_flat = reshape(G_neg_I, N, []);
        X_neg_I_flat = reshape(permute(X_neg_I, [2, 1, 3]), [], T);
        
        B0 = zeros(size(G_neg_I_flat, 2), size(G_I_flat, 2));
        lam_res = lam * (sqrt((length(I_hat)+1)*O) + sqrt(2*log(S)));
        Bt = local_gl_ADMM_dual(B0, G_neg_I_flat, G_I_flat, lam_res, O, [], 'consecutive', tol, tol_norm, max_iter, varing_rho);
        
        if all(Bt(:) == 0)
            disp('The penalty parameter is too large, the debiasing procedure is equivalent to the OLS regression.');
        end
        
        Z_I = G_I_flat - G_neg_I_flat * Bt;
        P_proj = Z_I * pinv(Z_I' * Z_I) * Z_I';
        effect_mat = P_proj * G_I_flat;
        
        if rank(effect_mat) == size(effect_mat, 2)
            midterm = P_proj * (Y - G_neg_I_flat * X_neg_I_flat);
            effect_mat_inv = pinv(effect_mat);
            X_hat_debias = effect_mat_inv * midterm;
            X_hat_debias_rotate = X_hat_debias; % Simplified
        else
            disp('The Q_A*G_A is not full rank, the debiasing is not performed.');
            X_I = Xt_tilde(I_hat, :, :);
            X_I_flat = reshape(permute(X_I, [2, 1, 3]), [], T);
            X_hat_debias = X_I_flat;
            X_hat_debias_rotate = X_I_flat;
            debias_flag = false;
        end
    else
         X_I = Xt_tilde(I_hat, :, :);
         X_I_flat = reshape(permute(X_I, [2, 1, 3]), [], T);
         X_hat_debias = X_I_flat;
         X_hat_debias_rotate = X_I_flat;
    end
end

function Xout = local_gl_ADMM_dual(X0, G, Y, lam, O, wlist, block_method, tol, tol_norm, max_iter, varing_rho)
    [N, T] = size(Y);
    S = size(G, 2) / O;
    
    G_tensor = local_Gmatrix_to_tensor(G, O, block_method);
    X0_tensor = local_Xmatrix_to_tensor(X0, O, block_method);
    
    if isempty(wlist)
        sigma_G_sqrt_inv = repmat(eye(O), [1, 1, S]);
    else
        sigma_G_sqrt = wlist{1};
        sigma_G_sqrt_inv = wlist{2};
    end
    
    G_tilde_tensor = zeros(size(G_tensor));
    for s = 1:S
        G_tilde_tensor(:, :, s) = G_tensor(:, :, s) * sigma_G_sqrt_inv(:, :, s);
    end
    G_tilde = reshape(G_tilde_tensor, N, S*O);

    GGT = G_tilde * G_tilde';
    stepsize0 = 1/S/O;
    stepsize = stepsize0;
    max_stepsize_tol = stepsize0 * sqrt(N*S*O);
    min_stepsize_tol = stepsize0 / sqrt(N*S*O);
    
    Xt_tensor = zeros(size(X0_tensor));
    for s = 1:S
        Xt_tensor(s,:,:) = permute(wlist{1}(:,:,s) * permute(X0_tensor(s,:,:),[2,3,1]),[3,1,2]);
    end
    Xt = reshape(permute(Xt_tensor,[2,1,3]), S*O, T);
    
    utnorm = sqrt(sum(reshape(Xt, O, S, T).^2, [1, 3]));
    utmask_locs = find(utnorm > 0);
    Xt = Xt(utmask_locs, :);

    Zt = Y - G_tilde(:, utmask_locs) * Xt;
    Ut = G_tilde' * Zt;
    
    loss_p = inf; loss_d = inf;
    tol_p = 0; tol_d = 0;
    sqrtSOT = sqrt(S*O*T);
    sqrtNT = sqrt(N*T);

    iter_idx = 0;
    
    while iter_idx < max_iter && (loss_p > tol_p || loss_d > tol_d)
        coreinv = inv(eye(N) + stepsize * GGT);
        coreinvG = coreinv * G_tilde;
        coreinvY = coreinv * Y;
        
        % Update Z
        mXtpUt = stepsize * Ut;
        mXtpUt(utmask_locs, :) = mXtpUt(utmask_locs, :) - Xt;
        Zt = coreinvY + coreinvG * mXtpUt;
        GTZ = G_tilde' * Zt;
        
        % Update U
        Utpre = Ut;
        Ut = GTZ;
        Ut(utmask_locs, :) = GTZ(utmask_locs, :) + Xt / stepsize;
        
        Ut_reshaped = reshape(Ut, O, S, T);
        utnorm = squeeze(sqrt(sum(Ut_reshaped.^2, [1, 3])));
        
        utmasknew_locs_bool = utnorm > lam;
        utmasknew_locs = find(utmasknew_locs_bool);

        if ~isempty(utmasknew_locs)
            Ut_temp = Ut(utmasknew_locs, :);
            norm_vals = utnorm(utmasknew_locs);
            Ut(utmasknew_locs, :) = lam * bsxfun(@rdivide, Ut_temp, norm_vals');
        end
        
        % Update X
        Xtt = stepsize * (GTZ(utmasknew_locs, :) - Ut(utmasknew_locs, :));
        [~, ia, ib] = intersect(utmask_locs, utmasknew_locs);
        if ~isempty(ia)
             Xtt(ib, :) = Xt(ia, :) + Xtt(ib, :);
        end
        Xt = Xtt;
        utmask_locs = utmasknew_locs;

        % Check convergence
        loss_p = norm(GTZ - Ut, 'fro');
        loss_d = norm(stepsize * G_tilde' * (Ut - Utpre), 'fro');
        tol_p = tol * sqrtSOT + tol_norm * max(norm(GTZ, 'fro'), norm(Ut, 'fro'));
        tol_d = tol * sqrtNT + tol_norm * norm(G_tilde(:, utmask_locs) * Xt, 'fro');

        % Update stepsize
        if varing_rho && iter_idx < max_iter/10
            if loss_d/sqrtNT > 10*loss_p/sqrtSOT
                stepsize = max(stepsize/2, min_stepsize_tol);
            elseif loss_p/sqrtSOT > 10*loss_d/sqrtNT
                stepsize = min(stepsize*2, max_stepsize_tol);
            end
        end
        iter_idx = iter_idx + 1;
    end
    
    Xtt_full = zeros(S*O, T);
    Xtt_full(utmask_locs, :) = Xt;
    
    Xt_tensor = local_Xmatrix_to_tensor(Xtt_full, O, block_method);
    Xout_tensor = zeros(size(Xt_tensor));
    for s=1:S
        Xout_tensor(s,:,:) = permute(sigma_G_sqrt_inv(:,:,s) * permute(Xt_tensor(s,:,:),[2,3,1]), [3,1,2]);
    end
    
    Xout = reshape(permute(Xout_tensor, [2, 1, 3]), S*O, T);
    
    if iter_idx == max_iter
        warning('ADMM algorithm did not converge to the desired tolerance.');
    end
end