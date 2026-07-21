function [Source_GL, Param_GL] = seal_debiased_GroupLasso(Data, L, varargin)
%SEAL_DEBIASED_GROUPLASSO Debiased group Lasso ESI (DeESI).
%   [Source, Param] = seal_debiased_GroupLasso(Data, L, ...)
%
%   Unified SEAL entry point for the debiased spatial-temporal group Lasso
%   proposed by Tong et al. (2025). The file contains both SEAL-compatible
%   parsing/noise prewhitening and the rotated group Lasso + projection-based
%   bias correction implementation.
%
%   Reference:
%       Tong et al. (2025), Debiased Estimation and Inference for
%       Spatial-Temporal EEG/MEG Source Imaging, IEEE TMI.
%       DOI: 10.1109/TMI.2024.3506596.

    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

    addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'RegularizationParameter', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'DynamicLambda', true, @(x) isLogicalLike(x));
    addParameter(p, 'JointEstNum', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Debias', true, @(x) isLogicalLike(x));
    addParameter(p, 'MaxIterations', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Tolerance', 1e-8, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'ToleranceNorm', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'JointTolerance', 1e-4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'VaryingRho', true, @(x) isLogicalLike(x));
    addParameter(p, 'ClearNotSelect', true, @(x) isLogicalLike(x));
    addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'PrewhitenNoise', true, @(x) isLogicalLike(x));
    addParameter(p, 'NoiseObservationMatrix', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x)); %#ok<NVREPL>
    addParameter(p, 'MaxIterationsDebias', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    addParameter(p, 'ToleranceDebias', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    addParameter(p, 'Verbose', true, @(x) isLogicalLike(x));

    try
        parse(p, Data, L, varargin{:});
    catch ME
        disp('Error parsing inputs for seal_debiased_GroupLasso:');
        rethrow(ME);
    end
    opts = p.Results;

    opts.NumOrientations = round(opts.NumOrientations);
    opts.JointEstNum = round(opts.JointEstNum);
    opts.MaxIterations = round(opts.MaxIterations);
    opts.DynamicLambda = parseLogical(opts.DynamicLambda, 'DynamicLambda');
    opts.Debias = parseLogical(opts.Debias, 'Debias');
    opts.VaryingRho = parseLogical(opts.VaryingRho, 'VaryingRho');
    opts.ClearNotSelect = parseLogical(opts.ClearNotSelect, 'ClearNotSelect');
    opts.PrewhitenNoise = parseLogical(opts.PrewhitenNoise, 'PrewhitenNoise');
    opts.Verbose = parseLogical(opts.Verbose, 'Verbose');

    if ~isempty(opts.MaxIterationsDebias)
        opts.MaxIterations = max(opts.MaxIterations, round(opts.MaxIterationsDebias));
    end
    if ~isempty(opts.ToleranceDebias)
        opts.Tolerance = min(opts.Tolerance, opts.ToleranceDebias);
    end

    Data0 = double(opts.Data);
    L0 = double(opts.L);
    [nChan, nTime] = size(Data0);
    if size(L0, 1) ~= nChan
        error('seal_debiased_GroupLasso:DataLeadfieldMismatch', ...
            'Data rows (%d) must match L rows (%d).', nChan, size(L0, 1));
    end
    nSource = size(L0, 2);
    nd = opts.NumOrientations;
    if mod(nSource, nd) ~= 0
        error('seal_debiased_GroupLasso:NumOrientationsMismatch', ...
            'size(L,2) (%d) must be divisible by NumOrientations (%d).', nSource, nd);
    end
    if any(~isfinite(Data0(:))) || any(~isfinite(L0(:)))
        error('seal_debiased_GroupLasso:InvalidInput', ...
            'Data and L must contain finite values.');
    end

    [DataW, LW, Cnoise, noiseInfo] = prepareNoiseModel( ...
        Data0, L0, opts.NoiseCovariance, opts.PrewhitenNoise);
    opts.NoiseCovariance = Cnoise;

    if opts.Verbose
        fprintf('Running DeESI debiased group Lasso (nd=%d, T=%d)...\n', nd, nTime);
    end

    implArgs = { ...
        'NumOrientations', nd, ...
        'RegularizationParameter', opts.RegularizationParameter, ...
        'DynamicLambda', opts.DynamicLambda, ...
        'JointEstNum', opts.JointEstNum, ...
        'Debias', opts.Debias, ...
        'MaxIterations', opts.MaxIterations, ...
        'Tolerance', opts.Tolerance, ...
        'ToleranceNorm', opts.ToleranceNorm, ...
        'JointTolerance', opts.JointTolerance, ...
        'VaryingRho', opts.VaryingRho, ...
        'ClearNotSelect', opts.ClearNotSelect, ...
        'Verbose', opts.Verbose};

    [Source_GL, implParam] = run_deesi_core(DataW, LW, implArgs{:});

    Param_GL = implParam;
    if isfield(Param_GL, 'OptionsPassed')
        Param_GL.ImplementationOptions = Param_GL.OptionsPassed;
    else
        Param_GL.ImplementationOptions = struct();
    end
    Param_GL.Algorithm = 'DeESI';
    Param_GL.Reference = 'Tong et al. 2025, IEEE TMI, DOI:10.1109/TMI.2024.3506596';
    Param_GL.Implementation = 'seal_debiased_GroupLasso';
    Param_GL.NoiseWhitening = noiseInfo;
    Param_GL.OptionsPassed = opts;

    if isfield(implParam, 'RawLassoSource')
        Param_GL.RawRotatedGroupLassoSource = implParam.RawLassoSource;
    end
    if isfield(implParam, 'RotatedDebias')
        Param_GL.StandardizedSource = implParam.RotatedDebias;
    else
        Param_GL.StandardizedSource = [];
    end
    if isfield(implParam, 'Significance')
        Param_GL.SignificanceStatistic = implParam.Significance;
    else
        Param_GL.SignificanceStatistic = [];
    end
end

function tf = isLogicalLike(x)
    tf = islogical(x) || (isnumeric(x) && isscalar(x)) || ischar(x) || isstring(x);
end

function tf = parseLogical(x, name)
    if islogical(x)
        tf = x;
    elseif isnumeric(x) && isscalar(x)
        tf = x ~= 0;
    elseif ischar(x) || isstring(x)
        x = lower(strtrim(char(x)));
        if any(strcmp(x, {'true', '1', 'yes', 'on'}))
            tf = true;
        elseif any(strcmp(x, {'false', '0', 'no', 'off'}))
            tf = false;
        else
            error('seal_debiased_GroupLasso:InvalidLogical', ...
                '%s must be logical-compatible.', name);
        end
    else
        error('seal_debiased_GroupLasso:InvalidLogical', ...
            '%s must be logical-compatible.', name);
    end
end

function [DataW, LW, Cnoise, noiseInfo] = prepareNoiseModel(Data, L, Cnoise, prewhiten)
    nChan = size(Data, 1);
    noiseInfo = struct('Applied', false, 'Method', 'none', 'Regularization', 0);

    if isempty(Cnoise)
        Cnoise = eye(nChan);
        DataW = Data;
        LW = L;
        return;
    end
    if size(Cnoise, 1) ~= nChan || size(Cnoise, 2) ~= nChan
        error('seal_debiased_GroupLasso:NoiseCovarianceSize', ...
            'NoiseCovariance must be %d x %d.', nChan, nChan);
    end

    Cnoise = double((Cnoise + Cnoise') / 2);
    if ~prewhiten
        DataW = Data;
        LW = L;
        return;
    end

    [V, D] = eig(Cnoise);
    d = real(diag(D));
    floorVal = max(max(d) * 1e-12, eps);
    d(d < floorVal) = floorVal;
    W = diag(1 ./ sqrt(d)) * V';
    DataW = W * Data;
    LW = W * L;
    noiseInfo = struct('Applied', true, 'Method', 'eig', 'Regularization', floorVal);
end

function [Source_GL, Param_GL] = run_deesi_core(Data, L, varargin)
%RUN_DEESI_CORE Core rotated group Lasso and DeESI bias correction.

    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L',    @(x) isnumeric(x) && ismatrix(x));

    addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'RegularizationParameter', 1.0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'DynamicLambda', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'JointEstNum', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Debias', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'MaxIterations', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Tolerance', 1e-8, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'ToleranceNorm', 1e-3, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'JointTolerance', 1e-4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'VaryingRho', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'ClearNotSelect', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

    parse(p, Data, L, varargin{:});
    opts = p.Results;

    Y = Data;
    G = L;

    [N, T] = size(Y);
    O = opts.NumOrientations;
    S = size(G, 2) / O;

    if abs(round(S) - S) > 0
        error('seal_debiased_GroupLasso:NumOrientationsMismatch', ...
            'size(L,2) must be divisible by NumOrientations.');
    end
    S = round(S);

    [sigma_G_sqrt, sigma_G_sqrt_inv, G_tensor] = local_gtensor_covariance_from_matrix(G, O);
    wlist.sigma_G_sqrt = sigma_G_sqrt;
    wlist.sigma_G_sqrt_inv = sigma_G_sqrt_inv;

    X0 = zeros(S * O, T);

    if opts.DynamicLambda
        joint_loops = opts.JointEstNum;
    else
        joint_loops = 1;
    end

    lam_base = opts.RegularizationParameter;
    lam_new = lam_base * norm(Y - G * X0, 'fro') / sqrt(N * T);

    Xt = X0;
    sigma_hat_pre = inf;
    sigma_list = [];

    for j = 1:joint_loops
        lam_use = lam_new * sqrt(T);

        Xt = gl_ADMM_dual_matlab( ...
            Xt, G, Y, lam_use, O, wlist, ...
            opts.Tolerance, opts.ToleranceNorm, opts.MaxIterations, opts.VaryingRho, opts.Verbose);

        sigma_hat = norm(Y - G * Xt, 'fro') / sqrt(N * T);
        lam_new = lam_base * sigma_hat;
        sigma_list(end+1, 1) = sigma_hat; %#ok<AGROW>

        if opts.Verbose
            fprintf('Joint iter %d: sigma = %.6g, lam_use = %.6g\n', j, sigma_hat, lam_use);
        end

        if abs(sigma_hat - sigma_hat_pre) < opts.JointTolerance
            break;
        end
        sigma_hat_pre = sigma_hat;
    end

    [I_hat, ~, Xt_tensor] = local_cal_I_hat(Xt, G, O, 0);

    debias_info = struct( ...
        'Requested', logical(opts.Debias), ...
        'Mode', 'none', ...
        'Applied', false, ...
        'JointFeasible', false, ...
        'ActiveLocationCount', numel(I_hat), ...
        'MaxJointLocationCount', max(0, ceil(N / O) - 1), ...
        'CorrectedLocationCount', 0, ...
        'Reason', 'Debias not requested or no active locations.');

    if opts.Debias && ~isempty(I_hat)
        [Xout_debias, Xout_debias_rotate, significance_list, debias_info] = ...
            gl_ADMM_dual_bias_correction_matlab( ...
                I_hat, Xt_tensor, G_tensor, Y, lam_use, O, wlist, ...
                true, opts.ClearNotSelect, ...
                opts.Tolerance, opts.ToleranceNorm, opts.MaxIterations, opts.VaryingRho, opts.Verbose);
        Source_GL = Xout_debias;
    else
        Source_GL = Xt;
        Xout_debias_rotate = Xt;
        significance_list = zeros(S, 1);
    end

    Param_GL.RawLassoSource = Xt;
    Param_GL.ActiveLocations = I_hat;
    Param_GL.ActiveLocationCount = numel(I_hat);
    Param_GL.FinalSigma = sigma_list(end);
    Param_GL.SigmaHistory = sigma_list;
    Param_GL.FinalLambda = lam_new;
    Param_GL.RotatedDebias = Xout_debias_rotate;
    Param_GL.Significance = significance_list;
    Param_GL.DebiasingInfo = debias_info;
    Param_GL.OptionsPassed = opts;
end

function Xout = gl_ADMM_dual_matlab(X0, G, Y, lam, O, wlist, tol, tol_norm, max_iter, varying_rho, verbose)
%GL_ADMM_DUAL_MATLAB Core ADMM solver for the rotated group Lasso problem.

    [N, SO] = size(G);
    S = SO / O;
    T = size(Y, 2);

    G_tensor = local_Gmatrix_to_tensor(G, O);
    X0_tensor = local_Xmatrix_to_tensor(X0, O);

    sigma_G_sqrt = wlist.sigma_G_sqrt;
    sigma_G_sqrt_inv = wlist.sigma_G_sqrt_inv;

    G_tilde_tensor = zeros(N, O, S);
    for s = 1:S
        G_tilde_tensor(:, :, s) = G_tensor(:, :, s) * sigma_G_sqrt_inv(:, :, s);
    end
    G_tilde = reshape(G_tilde_tensor, [N, S * O]);

    GGT = G_tilde * G_tilde';

    stepsize0 = 1 / (S * O);
    max_stepsize_tol = stepsize0 * sqrt(N * S * O);
    min_stepsize_tol = stepsize0 / sqrt(N * S * O);
    stepsize = stepsize0;

    coreinv = inv(eye(N) + stepsize * GGT);
    coreinvG = coreinv * G_tilde;
    coreinvY = coreinv * Y;

    fastcoeff = 1;

    Xt_tensor = zeros(S, O, T);
    for s = 1:S
        Xt_tensor(s, :, :) = sigma_G_sqrt(:, :, s) * squeeze(X0_tensor(s, :, :));
    end
    Xt_full = reshape(Xt_tensor, [S * O, T]);

    utnorm = local_row_or_group_norm(Xt_full, O, S, T);
    utmask = (utnorm > 0);

    Xt = Xt_full(utmask, :);
    Zt = Y - G_tilde(:, utmask) * Xt;
    Ut = G_tilde' * Zt;

    loss_p = inf;
    loss_d = inf;

    sqrtSOT = sqrt(S * O * T);
    sqrtNT = sqrt(N * T);
    tolsqrtSOT = tol * sqrtSOT;
    tolsqrtNT = tol * sqrtNT;

    iter_idx = 0;
    if verbose
        fprintf('iter = %d, loss_p = %.6g, loss_d = %.6g, lam = %.6g, stepsize = %.6g\n', ...
            iter_idx, loss_p, loss_d, lam, stepsize);
    end

    while iter_idx < max_iter
        mXtpUt = stepsize * Ut;
        mXtpUt(utmask, :) = mXtpUt(utmask, :) - Xt;
        Zt = coreinvY + coreinvG * mXtpUt;

        GTZ = G_tilde' * Zt;

        Utpre = Ut;
        Ut = GTZ;
        Ut(utmask, :) = GTZ(utmask, :) + Xt / stepsize;

        utnorm_new = local_row_or_group_norm(Ut, O, S, T);
        utmasknew = (utnorm_new > lam);

        idx_new = find(utmasknew);
        if ~isempty(idx_new)
            Ut(idx_new, :) = lam .* Ut(idx_new, :) ./ utnorm_new(idx_new);
        end

        Xtt = fastcoeff * stepsize * (GTZ(utmasknew, :) - Ut(utmasknew, :));

        utmask_and = utmask & utmasknew;
        idx_and_global = find(utmask_and);
        if ~isempty(idx_and_global)
            idx_old_local = local_globalmask_to_localmask(utmask, idx_and_global);
            idx_new_local = local_globalmask_to_localmask(utmasknew, idx_and_global);
            Xtt(idx_new_local, :) = Xt(idx_old_local, :) + Xtt(idx_new_local, :);
        end

        Xt = Xtt;
        utmask = utmasknew;

        GX = G_tilde(:, utmask) * Xt;

        loss_p = norm(GTZ - Ut, 'fro');
        loss_d = norm(stepsize * G_tilde * (Ut - Utpre), 'fro');

        tol_p = tolsqrtSOT + tol_norm * max(norm(GTZ, 'fro'), norm(Ut, 'fro'));
        tol_d = tolsqrtNT + tol_norm * norm(GX, 'fro');

        if loss_p <= tol_p && loss_d <= tol_d
            break;
        end

        if varying_rho && iter_idx < max_iter / 10
            if loss_d / sqrtNT > 10 * loss_p / sqrtSOT
                stepsize = max(stepsize / 2, min_stepsize_tol);
                coreinv = inv(eye(N) + stepsize * GGT);
                coreinvG = coreinv * G_tilde;
                coreinvY = coreinv * Y;
            elseif loss_p / sqrtSOT > 10 * loss_d / sqrtNT
                stepsize = min(stepsize * 2, max_stepsize_tol);
                coreinv = inv(eye(N) + stepsize * GGT);
                coreinvG = coreinv * G_tilde;
                coreinvY = coreinv * Y;
            end
        end

        iter_idx = iter_idx + 1;
        if verbose && mod(iter_idx, 200) == 0
            fprintf('iter = %d, loss_p = %.6g, loss_d = %.6g, lam = %.6g, stepsize = %.6g, active_source = %d\n', ...
                iter_idx, loss_p, loss_d, lam, stepsize, sum(utmask) / O);
        end
    end

    Xtt_full = zeros(S * O, T);
    Xtt_full(utmask, :) = Xt;

    Xt_tensor2 = local_Xmatrix_to_tensor(Xtt_full, O);
    Xout_tensor = zeros(S, O, T);
    for s = 1:S
        Xout_tensor(s, :, :) = sigma_G_sqrt_inv(:, :, s) * squeeze(Xt_tensor2(s, :, :));
    end
    Xout = reshape(Xout_tensor, [S * O, T]);

    if iter_idx == max_iter
        warning('seal_debiased_GroupLasso:ADMMMaxIter', ...
            ['ADMM did not converge within %d iterations. Final primal %.3g (tol %.3g), ' ...
             'dual %.3g (tol %.3g). Consider increasing MaxIterations or relaxing ToleranceNorm.'], ...
            max_iter, loss_p, tol_p, loss_d, tol_d);
    end
end

function [Xout_debias, Xout_debias_rotate, significance_list, debias_info] = ...
    gl_ADMM_dual_bias_correction_matlab(I_hat_all, Xt_tensor, G_tensor, Y, lam, O, wlist, ...
                                        bias_correction_joint, clear_not_select, ...
                                        tol, tol_norm, max_iter, varying_rho, verbose)
%GL_ADMM_DUAL_BIAS_CORRECTION_MATLAB Projection-based DeESI correction.

    [N, ~, S] = size(G_tensor);
    T = size(Y, 2);
    active_count = numel(I_hat_all);
    max_joint_count = max(0, ceil(N / O) - 1);
    joint_feasible = active_count > 0 && active_count < N / O;

    Xout_debias = Xt_tensor;
    significance_list = zeros(S, 1);
    debias_info = struct( ...
        'Requested', true, ...
        'Mode', 'none', ...
        'Applied', false, ...
        'JointFeasible', joint_feasible, ...
        'ActiveLocationCount', active_count, ...
        'MaxJointLocationCount', max_joint_count, ...
        'CorrectedLocationCount', 0, ...
        'Reason', '');

    neg_I_hat = setdiff(1:S, I_hat_all(:)');

    sigma_G_sqrt = wlist.sigma_G_sqrt;
    sigma_G_sqrt_inv = wlist.sigma_G_sqrt_inv;

    G_tilde = G_tensor;
    for s = neg_I_hat
        G_tilde(:, :, s) = G_tensor(:, :, s) * sigma_G_sqrt_inv(:, :, s);
    end

    if clear_not_select
        Xout_debias(neg_I_hat, :, :) = 0;
        Xout_debias_rotate = Xout_debias;

        Xt_tilde = Xout_debias;
        for s = neg_I_hat
            Xt_tilde(s, :, :) = sigma_G_sqrt(:, :, s) * squeeze(Xout_debias(s, :, :));
        end
    else
        Xout_debias_rotate = Xout_debias;
        Xout_debias_rotate(neg_I_hat, :, :) = 0;

        Xt_tilde = Xt_tensor;
        for s = neg_I_hat
            Xt_tilde(s, :, :) = sigma_G_sqrt(:, :, s) * squeeze(Xt_tensor(s, :, :));
        end
    end

    if bias_correction_joint && joint_feasible
        I_hat = I_hat_all(:)';
        [X_hat_debias_idx, X_hat_debias_rotate_idx, debias_flag_idx, effect_mat] = ...
            local_gl_ADMM_dual_bias_correction(I_hat, Xt_tilde, G_tilde, Y, lam, O, ...
                                               tol, tol_norm, max_iter, varying_rho, verbose);

        if ~isempty(I_hat)
            Xout_debias(I_hat, :, :) = reshape(X_hat_debias_idx, [length(I_hat), O, T]);
            Xout_debias_rotate(I_hat, :, :) = reshape(X_hat_debias_rotate_idx, [length(I_hat), O, T]);
        end

        if debias_flag_idx
            debias_info.Mode = 'joint';
            debias_info.Applied = true;
            debias_info.CorrectedLocationCount = length(I_hat);
            effect_mat_tensor = local_Gmatrix_to_tensor(effect_mat, O);
            for k = 1:length(I_hat)
                s = I_hat(k);
                significance_list(s) = norm(effect_mat_tensor(:, :, k) * squeeze(Xout_debias(s, :, :)), 'fro');
            end
        else
            debias_info.Mode = 'joint';
            debias_info.Reason = 'Projected active lead field was rank deficient.';
        end
    elseif bias_correction_joint
        debias_info.Mode = 'skipped';
        debias_info.Reason = sprintf( ...
            'Joint debias requires active locations < N/O; active=%d, maximum=%d.', ...
            active_count, max_joint_count);
        if verbose
            fprintf('Skipping joint debias: active locations = %d, maximum = %d for N = %d, O = %d.\n', ...
                active_count, max_joint_count, N, O);
        end
    else
        debias_info.Mode = 'single';
        for ii = 1:length(I_hat_all)
            I_hat = I_hat_all(ii);
            [X_hat_debias_idx, X_hat_debias_rotate_idx, debias_flag_idx, effect_mat] = ...
                local_gl_ADMM_dual_bias_correction(I_hat, Xt_tilde, G_tilde, Y, lam, O, ...
                                                   tol, tol_norm, max_iter, varying_rho, verbose);

            Xout_debias(I_hat, :, :) = reshape(X_hat_debias_idx, [1, O, T]);
            Xout_debias_rotate(I_hat, :, :) = reshape(X_hat_debias_rotate_idx, [1, O, T]);
            if debias_flag_idx
                debias_info.Applied = true;
                debias_info.CorrectedLocationCount = debias_info.CorrectedLocationCount + 1;
                significance_list(I_hat) = norm(effect_mat * X_hat_debias_idx, 'fro');
            end
        end
        if ~debias_info.Applied
            debias_info.Reason = 'No single-source debias correction satisfied the rank condition.';
        end
    end

    Xout_debias = reshape(Xout_debias, [S * O, T]);
    Xout_debias_rotate = reshape(Xout_debias_rotate, [S * O, T]);
end

function [X_hat_debias, X_hat_debias_rotate, debias_flag, effect_mat] = ...
    local_gl_ADMM_dual_bias_correction(I_hat, Xt_tilde, G_tilde, Y, lam, O, ...
                                       tol, tol_norm, max_iter, varying_rho, verbose)
%LOCAL_GL_ADMM_DUAL_BIAS_CORRECTION Correct active groups using Q_A.

    [N, ~, S] = size(G_tilde);
    T = size(Y, 2);

    I_hat = I_hat(:)';
    A = length(I_hat);

    if A > 0 && A < N / O
        debias_flag = true;
        neg_I_hat = setdiff(1:S, I_hat);

        G_I = G_tilde(:, :, I_hat);
        G_neg_I = G_tilde(:, :, neg_I_hat);
        X_neg_I = Xt_tilde(neg_I_hat, :, :);

        G_neg_I_flat = reshape(G_neg_I, [N, (S - A) * O]);
        G_I_flat = reshape(G_I, [N, A * O]);
        X_neg_I_flat = reshape(X_neg_I, [(S - A) * O, T]);

        B0 = zeros((S - A) * O, A * O);

        lam_res = lam * (sqrt((A + 1) * O) + sqrt(2 * log(S)));

        Bt = gl_ADMM_dual_matlab( ...
            B0, G_neg_I_flat, G_I_flat, lam_res, O, local_identity_wlist(S - A, O), ...
            tol, tol_norm, max_iter, varying_rho, verbose);

        if verbose && all(Bt(:) == 0)
            fprintf('Penalty too large in debias stage; close to OLS correction.\n');
        end

        Z_I = G_I_flat - G_neg_I_flat * Bt;
        P_proj = Z_I / (Z_I' * Z_I) * Z_I';
        effect_mat = P_proj * G_I_flat;

        if rank(effect_mat) == size(effect_mat, 2)
            midterm = P_proj * (Y - G_neg_I_flat * X_neg_I_flat);
            effect_mat_inv = pinv(effect_mat);
            X_hat_debias = effect_mat_inv * midterm;

            debias_var = effect_mat_inv * P_proj * effect_mat_inv';
            d = diag(debias_var);
            d(d <= 0) = eps;
            debias_var_inv_square = diag(1 ./ sqrt(d));
            X_hat_debias_rotate = debias_var_inv_square * X_hat_debias;
        else
            if verbose
                fprintf('Q_A * G_A is not full rank; skip debias.\n');
            end
            X_I = Xt_tilde(I_hat, :, :);
            X_hat_debias = reshape(X_I, [A * O, T]);
            X_hat_debias_rotate = X_hat_debias;
            debias_flag = false;
            effect_mat = 0;
        end
    else
        debias_flag = false;
        X_I = Xt_tilde(I_hat, :, :);
        X_hat_debias = reshape(X_I, [A * O, T]);
        X_hat_debias_rotate = X_hat_debias;
        effect_mat = 0;
    end
end

function [sigma_G_sqrt, sigma_G_sqrt_inv, G_tensor] = local_gtensor_covariance_from_matrix(G, O)
    G_tensor = local_Gmatrix_to_tensor(G, O);
    [N, ~, S] = size(G_tensor);

    sigma_G_sqrt = zeros(O, O, S);
    sigma_G_sqrt_inv = zeros(O, O, S);

    for s = 1:S
        Ls = G_tensor(:, :, s);
        sigma_G = (Ls' * Ls) / N;

        [V, D] = eig(sigma_G);
        d = diag(D);
        d(d < 1e-12) = 1e-12;

        sigma_G_sqrt(:, :, s) = V * diag(sqrt(d)) * V';
        sigma_G_sqrt_inv(:, :, s) = V * diag(1 ./ sqrt(d)) * V';
    end
end

function G_tensor = local_Gmatrix_to_tensor(G, O)
    [N, SO] = size(G);
    S = SO / O;
    G_tensor = zeros(N, O, S);
    for s = 1:S
        idx = (s - 1) * O + (1:O);
        G_tensor(:, :, s) = G(:, idx);
    end
end

function X_tensor = local_Xmatrix_to_tensor(X, O)
    [SO, T] = size(X);
    S = SO / O;
    X_tensor = zeros(S, O, T);
    for s = 1:S
        idx = (s - 1) * O + (1:O);
        X_tensor(s, :, :) = X(idx, :);
    end
end

function [I_hat, G_tensor, Xt_tensor] = local_cal_I_hat(Xt, G, O, threshold)
    G_tensor = local_Gmatrix_to_tensor(G, O);
    Xt_tensor = local_Xmatrix_to_tensor(Xt, O);

    [S, ~, ~] = size(Xt_tensor);
    block_norm = zeros(S, 1);
    for s = 1:S
        block_norm(s) = norm(squeeze(Xt_tensor(s, :, :)), 'fro');
    end
    I_hat = find(block_norm > threshold).';
end

function utnorm = local_row_or_group_norm(U, O, S, T)
    if O == 1
        utnorm = vecnorm(U, 2, 2);
    else
        U3 = reshape(U, [O, S, T]);
        gnorm = squeeze(sqrt(sum(sum(U3.^2, 1), 3)));
        utnorm = repelem(gnorm, O, 1);
        utnorm = utnorm(:);
    end
end

function idx_local = local_globalmask_to_localmask(mask_global, idx_global)
    global_pos = find(mask_global);
    [~, idx_local] = ismember(idx_global, global_pos);
end

function wlist = local_identity_wlist(S, O)
    sigma_G_sqrt = zeros(O, O, S);
    sigma_G_sqrt_inv = zeros(O, O, S);
    I = eye(O);
    for s = 1:S
        sigma_G_sqrt(:, :, s) = I;
        sigma_G_sqrt_inv(:, :, s) = I;
    end
    wlist.sigma_G_sqrt = sigma_G_sqrt;
    wlist.sigma_G_sqrt_inv = sigma_G_sqrt_inv;
end
