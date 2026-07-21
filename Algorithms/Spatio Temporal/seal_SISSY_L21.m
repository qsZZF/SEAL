function [S_final, Param_SISSY] = seal_SISSY_L21(Data, L, sourceSpaceInfo, varargin)
% SEAL_SISSY_L21 L1,2-SISSY source imaging solved with ADMM.
%
% Solves:
%   min_S 0.5*||X - L*S||_F^2 + lambda*(||V*S||_1,2 + alpha*||S||_1,2)
%
% The original SISSY paper uses fixed-orientation cortical sources. This
% implementation also supports grouped multi-orientation sources as a SEAL
% extension when NumOrientations > 1.
%
% Required inputs:
%   Data            Nchannels x Ntimepoints data matrix.
%   L               Nchannels x Nsources leadfield matrix.
%   sourceSpaceInfo Struct with source-space geometry and VertConn.
%
% Name-value options:
%   lambda             Base regularization parameter. Default: 1.0.
%   DynamicLambda      If true, scale lambda by residual RMS noise. Default: true.
%   JointEstNum        Outer residual-noise update count. Default: 10.
%   alpha              Source sparsity weight. Default: 0.07.
%   rho                ADMM penalty parameter. Default: 1.0.
%   NumOrientations    Orientations per source location. Default: 1.
%   max_iter           Maximum ADMM iterations per outer loop. Default: 200.
%   tol                ADMM convergence tolerance. Default: 1e-4.
%   NoiseCovariance    Optional channel noise covariance for direct calls.
%   PrewhitenNoise     Whiten Data/L when NoiseCovariance is non-identity. Default: true.
%   NormalizeLeadfield Normalize leadfield blocks and restore source units. Default: true.

    %% Input parsing and validation
    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'sourceSpaceInfo', @isstruct);

    addParameter(p, 'lambda', 1.0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'DynamicLambda', true, @isLogicalLike);
    addParameter(p, 'JointEstNum', 10, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'alpha', 0.07, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'rho', 1.0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'max_iter', 200, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'tol', 1e-4, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'PrewhitenNoise', true, @isLogicalLike);
    addParameter(p, 'NormalizeLeadfield', true, @isLogicalLike);

    parse(p, Data, L, sourceSpaceInfo, varargin{:});

    Data = double(Data);
    L = double(L);

    lam_base = double(p.Results.lambda);
    dynamic_lambda = localToLogical(p.Results.DynamicLambda);
    joint_loops = max(1, round(double(p.Results.JointEstNum)));
    alpha = double(p.Results.alpha);
    rho = double(p.Results.rho);
    nd = round(double(p.Results.NumOrientations));
    max_iter = max(1, round(double(p.Results.max_iter)));
    tol = double(p.Results.tol);
    prewhiten_noise = localToLogical(p.Results.PrewhitenNoise);
    normalize_leadfield = localToLogical(p.Results.NormalizeLeadfield);

    validateScalar(lam_base, 'lambda', 0, inf);
    validateScalar(alpha, 'alpha', 0, 1);
    validateScalar(rho, 'rho', eps, inf);
    validateScalar(tol, 'tol', eps, inf);
    if nd < 1
        error('seal_SISSY_L21:InvalidNumOrientations', ...
            'NumOrientations must be a positive integer.');
    end
    if ~all(isfinite(Data(:))) || ~all(isfinite(L(:)))
        error('seal_SISSY_L21:InvalidInput', ...
            'Data and L must contain only finite values.');
    end

    [Nchannels, Ntimepoints] = size(Data);
    if size(L, 1) ~= Nchannels
        error('seal_SISSY_L21:DimensionMismatch', ...
            'Data channels (%d) and leadfield rows (%d) must match.', ...
            Nchannels, size(L, 1));
    end

    Nsources_total = size(L, 2);
    if mod(Nsources_total, nd) ~= 0
        error('seal_SISSY_L21:OrientationMismatch', ...
            'Leadfield columns (%d) must be divisible by NumOrientations (%d).', ...
            Nsources_total, nd);
    end
    Nsources = Nsources_total / nd;

    if ~isfield(sourceSpaceInfo, 'VertConn') || isempty(sourceSpaceInfo.VertConn)
        error('seal_SISSY_L21:MissingVertConn', ...
            'sourceSpaceInfo.VertConn is required for SISSY spatial variation regularization.');
    end
    [V_loc, graphInfo] = VariationEdge(sourceSpaceInfo.VertConn, Nsources);

    if nd > 1
        V = kron(V_loc, speye(nd));
    else
        V = V_loc;
    end
    [E_total, ~] = size(V);

    [Data, L, whiteningInfo] = applyNoiseWhitening( ...
        Data, L, p.Results.NoiseCovariance, prewhiten_noise);

    [L, leadfieldScale] = normalizeLeadfieldBlocks(L, Nsources, nd, normalize_leadfield);

    %% Initialization and Woodbury precomputation
    fprintf('Initializing SISSY ADMM (locations=%d, nd=%d, edges=%d)...\n', ...
        Nsources, nd, graphInfo.NumEdges);

    S = zeros(Nsources_total, Ntimepoints);
    Y = zeros(E_total, Ntimepoints);
    Z = zeros(Nsources_total, Ntimepoints);
    U = zeros(E_total, Ntimepoints);
    W = zeros(Nsources_total, Ntimepoints);

    Q = rho * (V' * V + speye(Nsources_total));
    dA = decomposition(Q, 'chol');

    Q_inv_L = dA \ L';
    C_core = eye(Nchannels) + L * Q_inv_L;
    LTX = L' * Data;

    history.primal_res_norm = [];
    history.dual_res_norm = [];
    history.sigma_list = [];
    history.lambda_list = [];

    %% Dynamic lambda outer loop
    if ~dynamic_lambda
        joint_loops = 1;
        lam_new = lam_base;
        lambdaStrategy = 'fixed';
    else
        lam_new = lam_base * (norm(Data, 'fro') / sqrt(Nchannels * Ntimepoints));
        lambdaStrategy = 'residual_noise_rms';
    end

    sigma_hat_pre = inf;
    total_iter_count = 0;

    fprintf('Starting SISSY optimization (lambda strategy: %s)...\n', lambdaStrategy);
    for j = 1:joint_loops
        lam_use = lam_new;

        %% Inner ADMM loop
        for k = 1:max_iter
            Y_old = Y;
            Z_old = Z;

            % S-update: solve (L'*L + rho*(V'*V + I))*S = rhs by Woodbury.
            S_rhs = LTX + rho * V' * (Y - U) + rho * (Z - W);
            Q_inv_rhs = dA \ S_rhs;
            S = Q_inv_rhs - Q_inv_L * (C_core \ (L * Q_inv_rhs));

            Y = prox_L21_group(V * S + U, lam_use / rho, nd);
            Z = prox_L21_group(S + W, lam_use * alpha / rho, nd);

            U = U + V * S - Y;
            W = W + S - Z;

            r_primal_Y = V * S - Y;
            r_primal_Z = S - Z;
            prim_res = sqrt(norm(r_primal_Y, 'fro')^2 + norm(r_primal_Z, 'fro')^2);

            s_dual_Y = -rho * V' * (Y - Y_old);
            s_dual_Z = -rho * (Z - Z_old);
            dual_res = sqrt(norm(s_dual_Y, 'fro')^2 + norm(s_dual_Z, 'fro')^2);

            eps_pri = sqrt((E_total + Nsources_total) * Ntimepoints) * tol + ...
                tol * max( ...
                    sqrt(norm(V * S, 'fro')^2 + norm(S, 'fro')^2), ...
                    sqrt(norm(Y, 'fro')^2 + norm(Z, 'fro')^2));
            eps_dual = sqrt(Nsources_total * Ntimepoints) * tol + ...
                tol * sqrt(norm(rho * V' * U, 'fro')^2 + norm(rho * W, 'fro')^2);

            history.primal_res_norm = [history.primal_res_norm; prim_res]; %#ok<AGROW>
            history.dual_res_norm = [history.dual_res_norm; dual_res]; %#ok<AGROW>

            if (prim_res < eps_pri) && (dual_res < eps_dual)
                break;
            end
        end
        total_iter_count = total_iter_count + k;

        sigma_hat = norm(Data - L * S, 'fro') / sqrt(Nchannels * Ntimepoints);
        history.sigma_list = [history.sigma_list; sigma_hat]; %#ok<AGROW>
        history.lambda_list = [history.lambda_list; lam_use]; %#ok<AGROW>

        if dynamic_lambda
            lam_new = lam_base * sigma_hat;
        end

        fprintf('  Outer Loop %d | ADMM Iters: %d | Sigma: %.4g | Lambda: %.4g\n', ...
            j, k, sigma_hat, lam_use);

        if dynamic_lambda && abs(sigma_hat - sigma_hat_pre) < 1e-4
            break;
        end
        sigma_hat_pre = sigma_hat;
    end

    fprintf('Optimization finished. Total ADMM iterations: %d\n', total_iter_count);

    S_final = restoreLeadfieldScale(S, leadfieldScale, nd);

    Param_SISSY = struct();
    Param_SISSY.OptionsPassed = p.Results;
    Param_SISSY.OptionsPassed.DynamicLambda = dynamic_lambda;
    Param_SISSY.OptionsPassed.PrewhitenNoise = prewhiten_noise;
    Param_SISSY.OptionsPassed.NormalizeLeadfield = normalize_leadfield;
    Param_SISSY.history = history;
    Param_SISSY.FinalSigma = history.sigma_list(end);
    Param_SISSY.FinalLambda = history.lambda_list(end);
    Param_SISSY.iterations_run = total_iter_count;
    Param_SISSY.LambdaStrategy = lambdaStrategy;
    Param_SISSY.NoiseWhitening = whiteningInfo;
    Param_SISSY.LeadfieldNormalization = struct( ...
        'Applied', normalize_leadfield, ...
        'Scale', leadfieldScale, ...
        'ReturnedOriginalUnits', normalize_leadfield);
    Param_SISSY.Graph = graphInfo;
end

%% Helper functions

function Y_out = prox_L21_group(Y_in, beta, nd)
% Proximal operator for grouped L1,2 norm over orientations and time.
    if isempty(Y_in) || beta <= 0
        Y_out = Y_in;
        return;
    end

    [N_total, T] = size(Y_in);
    if mod(N_total, nd) ~= 0
        error('seal_SISSY_L21:GroupSizeMismatch', ...
            'Input rows (%d) must be divisible by group size (%d).', N_total, nd);
    end

    N_locs = N_total / nd;
    Y_tensor = reshape(Y_in, nd, N_locs, T);
    loc_norms = sqrt(sum(sum(Y_tensor .^ 2, 1), 3));
    shrinkage = max(1 - beta ./ max(loc_norms, eps), 0);
    Y_out_tensor = bsxfun(@times, Y_tensor, shrinkage);
    Y_out = reshape(Y_out_tensor, N_total, T);
end

function [V, info] = VariationEdge(VertConn, nSources)
% Build an incidence matrix from a sanitized undirected cortical graph.
    if size(VertConn, 1) ~= size(VertConn, 2)
        error('seal_SISSY_L21:InvalidVertConn', ...
            'VertConn must be a square adjacency matrix.');
    end
    if size(VertConn, 1) ~= nSources
        error('seal_SISSY_L21:VertConnSizeMismatch', ...
            'VertConn size (%d) must match source locations (%d).', ...
            size(VertConn, 1), nSources);
    end

    A_in = sparse(VertConn);
    inputEdges = nnz(A_in);
    A = spones(A_in + A_in');
    A = A - spdiags(diag(A), 0, nSources, nSources);
    A = spones(A);

    [i, j] = find(tril(A, -1));
    nEdges = numel(i);
    V = sparse(nEdges, nSources);
    if nEdges > 0
        edgeIdx = (1:nEdges)';
        V = sparse( ...
            [edgeIdx; edgeIdx], ...
            [i; j], ...
            [ones(nEdges, 1); -ones(nEdges, 1)], ...
            nEdges, nSources);
    end

    info = struct();
    info.NumLocations = nSources;
    info.NumEdges = nEdges;
    info.InputNonzeros = inputEdges;
    info.SanitizedNonzeros = nnz(A);
    info.Symmetrized = ~isequal(spones(A_in), A);
end

function [DataOut, LOut, info] = applyNoiseWhitening(Data, L, NoiseCovariance, prewhitenNoise)
    DataOut = Data;
    LOut = L;
    info = struct( ...
        'Applied', false, ...
        'Reason', 'not_requested', ...
        'Rank', [], ...
        'ConditionNumber', NaN);

    if ~prewhitenNoise
        return;
    end

    if isempty(NoiseCovariance)
        info.Reason = 'empty_noise_covariance';
        return;
    end

    nChannels = size(Data, 1);
    if size(NoiseCovariance, 1) ~= size(NoiseCovariance, 2)
        error('seal_SISSY_L21:InvalidNoiseCovariance', ...
            'NoiseCovariance must be square.');
    end
    if size(NoiseCovariance, 1) ~= nChannels
        error('seal_SISSY_L21:NoiseCovarianceMismatch', ...
            'NoiseCovariance size (%d) must match data channels (%d).', ...
            size(NoiseCovariance, 1), nChannels);
    end

    NoiseCovariance = full((double(NoiseCovariance) + double(NoiseCovariance)') / 2);
    if isIdentityLike(NoiseCovariance)
        info.Reason = 'identity_noise_covariance';
        return;
    end

    [U, S] = svd(NoiseCovariance, 'econ');
    singularValues = diag(S);
    if isempty(singularValues) || ~isfinite(max(singularValues)) || max(singularValues) <= 0
        error('seal_SISSY_L21:InvalidNoiseCovariance', ...
            'NoiseCovariance must have at least one positive finite singular value.');
    end

    sMax = max(singularValues);
    rankTol = max(size(NoiseCovariance)) * eps(sMax);
    rankNoise = sum(singularValues > rankTol);
    lowerBound = sMax * 1e-12;
    singularValues(singularValues < lowerBound) = lowerBound;

    Wnoise = diag(1 ./ sqrt(singularValues)) * U';
    DataOut = Wnoise * Data;
    LOut = Wnoise * L;

    info.Applied = true;
    info.Reason = 'noise_covariance';
    info.Rank = rankNoise;
    info.ConditionNumber = sMax / min(singularValues);
end

function [Lout, scale] = normalizeLeadfieldBlocks(Lin, nSources, nd, doNormalize)
    Lout = Lin;
    scale = ones(nSources, 1);

    if ~doNormalize
        return;
    end

    for s = 1:nSources
        idx = (s - 1) * nd + 1 : s * nd;
        blockNorm = norm(Lin(:, idx), 'fro');
        if ~isfinite(blockNorm) || blockNorm <= 0
            error('seal_SISSY_L21:ZeroLeadfieldBlock', ...
                'Leadfield block %d has zero or invalid norm.', s);
        end
        scale(s) = blockNorm;
        Lout(:, idx) = Lin(:, idx) ./ blockNorm;
    end
end

function S_out = restoreLeadfieldScale(S_in, scale, nd)
    S_out = S_in;
    if isempty(scale) || all(abs(scale - 1) < eps)
        return;
    end
    for s = 1:numel(scale)
        idx = (s - 1) * nd + 1 : s * nd;
        S_out(idx, :) = S_in(idx, :) ./ scale(s);
    end
end

function validateScalar(value, name, lowerBound, upperBound)
    if ~isscalar(value) || ~isfinite(value) || value < lowerBound || value > upperBound
        error('seal_SISSY_L21:InvalidParameter', ...
            '%s must be a finite scalar in [%g, %g].', name, lowerBound, upperBound);
    end
end

function tf = isLogicalLike(value)
    tf = (islogical(value) && isscalar(value)) || ...
        (isnumeric(value) && isscalar(value)) || ...
        ischar(value) || ...
        (isstring(value) && isscalar(value));
end

function tf = localToLogical(value)
    if islogical(value)
        tf = value;
    elseif isnumeric(value)
        tf = value ~= 0;
    else
        value = lower(strtrim(char(string(value))));
        tf = any(strcmp(value, {'true', '1', 'yes', 'on'}));
    end
    tf = logical(tf(1));
end

function tf = isIdentityLike(A)
    n = size(A, 1);
    delta = A - speye(n);
    tf = norm(delta, 'fro') <= 1e-12 * max(1, norm(A, 'fro'));
end
