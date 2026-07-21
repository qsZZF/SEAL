function [S_final, Param_dMAP] = seal_dMAP(X, L, VertConn, varargin)
%SEAL_DMAP Low-memory approximate dynamic MAP ESI solver.
%   [S, Param] = seal_dMAP(X, L, VertConn, ...)
%   [S, Param] = seal_dMAP(X, L, sourceSpaceInfo, ...)
%
%   This implementation follows the dMAP state-space idea with a diagonal
%   mean-field covariance approximation. It is intentionally lighter than a
%   full Kalman/RTS EM implementation: full source covariance matrices are
%   never formed.
%
% REQUIRED INPUTS
%   X        - [N_chan x N_time] sensor time series.
%   L        - [N_chan x N_src] lead field matrix, where
%              N_src = NumOrientations * N_loc.
%   VertConn - [N_loc x N_loc] cortical adjacency, or a struct containing
%              .VertConn. The graph is defined on physical source locations,
%              not orientation-expanded source rows.
%
% OPTIONAL NAME-VALUE PARAMETERS
%   'NoiseCovariance'      - [N_chan x N_chan] sensor noise covariance R.
%                            Default: eye(N_chan).
%   'InitialSourceCovDiag' - [N_src x 1] or [N_loc x 1] initial diagonal Q.
%                            Default: 1e-3 * ones(N_src,1).
%   'NumOrientations'      - orientations per location. Default: 1.
%   'lambda_F'             - transition decay in (0,1]. Default: 0.95.
%   'alpha_F'              - self/neighbor split in [0,1]. Default: 0.51.
%   'max_iter'             - maximum EM iterations. Default: 15.
%   'tol'                  - relative log-posterior tolerance. Default: 1e-5.
%   'q_min'                - lower bound for diagonal source covariance.
%                            Default: 1e-12.
%   'InnovationJitter'     - relative diagonal loading for innovation S.
%                            Default: 1e-10.
%
% OUTPUTS
%   S_final    - [N_src x N_time] smoothed source time courses.
%   Param_dMAP - struct with learned covariance diagonal and diagnostics.

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;

    [N_chan0, ~] = size(X);
    [~, N_src0] = size(L);

    addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'VertConn', @local_is_valid_graph_input);
    addParameter(p, 'NoiseCovariance', eye(N_chan0), @(x) isempty(x) || ismatrix(x));
    addParameter(p, 'InitialSourceCovDiag', 1e-3 * ones(N_src0, 1), @isvector);
    addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'lambda_F', 0.95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    addParameter(p, 'alpha_F', 0.51, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    addParameter(p, 'max_iter', 15, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'tol', 1e-5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'q_min', 1e-12, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'InnovationJitter', 1e-10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, X, L, VertConn, varargin{:});

    max_iter = round(p.Results.max_iter);
    tol = p.Results.tol;
    nd = round(p.Results.NumOrientations);
    q_min = p.Results.q_min;
    innovation_jitter = p.Results.InnovationJitter;

    %% 2. Dimension Checks and Initialization
    [N_chan, N_time] = size(X);
    [N_chanL, N_src] = size(L);
    if N_chan ~= N_chanL
        error('seal_dMAP:DimensionMismatch', 'X and L channel dimensions do not match.');
    end
    if mod(N_src, nd) ~= 0
        error('seal_dMAP:OrientationMismatch', 'size(L,2) must be divisible by NumOrientations.');
    end

    G_adj = local_extract_vertconn(p.Results.VertConn);
    N_loc = N_src / nd;
    if ~isequal(size(G_adj), [N_loc, N_loc])
        error('seal_dMAP:VertConnSize', ...
            'VertConn must be N_loc x N_loc, where N_loc = size(L,2)/NumOrientations.');
    end

    R = p.Results.NoiseCovariance;
    if isempty(R)
        R = eye(N_chan);
    end
    if ~isequal(size(R), [N_chan, N_chan])
        error('seal_dMAP:NoiseCovarianceSize', 'NoiseCovariance must be N_chan x N_chan.');
    end
    R = double((R + R') / 2);

    fprintf('dMAP: Constructing state transition matrix F...\n');
    F_loc = local_construct_F(G_adj, p.Results.lambda_F, p.Results.alpha_F);
    if nd > 1
        F = kron(F_loc, speye(nd));
    else
        F = F_loc;
    end
    F_sq = F .^ 2;

    q0 = double(p.Results.InitialSourceCovDiag(:));
    if numel(q0) == N_loc
        theta_loc = q0;
    elseif numel(q0) == N_src
        theta_loc = mean(reshape(q0, nd, []), 1)';
    else
        error('seal_dMAP:InitialSourceCovDiagSize', ...
            'InitialSourceCovDiag must have length N_src or N_loc.');
    end
    theta_loc = max(theta_loc, q_min);
    theta = repelem(theta_loc, nd);
    x0 = zeros(N_src, 1);

    logPost_hist = zeros(max_iter, 1);
    fprintf('Starting Diagonal Mean-Field dMAP (nd=%d)...\n', nd);

    %% 3. Expectation-Maximization Loop
    for iter = 1:max_iter
        % A. E-step: diagonal steady-state Riccati approximation.
        p_pred = theta;
        for k = 1:100
            [~, p_filt, ~, ~, ~] = local_filter_step(L, R, p_pred, q_min, innovation_jitter);
            p_pred_new = max(real(F_sq * p_filt + theta), q_min);
            rel_change = norm(p_pred_new - p_pred) / max(norm(p_pred), eps);
            p_pred = p_pred_new;
            if rel_change < 1e-5
                break;
            end
        end
        [K, p_filt, S, Ls, logdetS] = local_filter_step(L, R, p_pred, q_min, innovation_jitter);
        p_pred_safe = max(p_pred, q_min);

        % B. Backward diagonal Riccati approximation.
        % RTS mean gain under diagonal P: J = diag(p_filt) * F' * diag(1 ./ p_pred).
        J_sq = spdiags(p_filt .^ 2, 0, N_src, N_src) * ...
               (F' .^ 2) * ...
               spdiags(1 ./ (p_pred_safe .^ 2), 0, N_src, N_src);
        p_smooth = p_filt;
        for k = 1:100
            p_smooth_new = p_filt + J_sq * (p_smooth - p_pred_safe);
            p_smooth_new = max(real(p_smooth_new), q_min);
            rel_change = norm(p_smooth_new - p_smooth) / max(norm(p_smooth), eps);
            p_smooth = p_smooth_new;
            if rel_change < 1e-5
                break;
            end
        end

        diag_cross = p_smooth .* (F_sq * p_filt) ./ p_pred_safe;
        diag_cov_err = p_smooth - 2 * diag_cross + F_sq * p_smooth;
        diag_cov_err = max(real(diag_cov_err), 0);

        % C. Forward data pass with steady-state gain.
        x_filt = zeros(N_src, N_time);
        x_k = x0;
        logL = 0;
        for t = 1:N_time
            v_k = X(:, t) - L * x_k;
            x_filt_k = x_k + K * v_k;
            x_filt(:, t) = x_filt_k;
            x_k = F * x_filt_k;

            y = Ls \ v_k;
            quad = y' * y;
            logL = logL - 0.5 * (logdetS + quad + N_chan * log(2*pi));
        end

        % D. Backward mean smoothing.
        x_smooth = zeros(N_src, N_time);
        x_smooth(:, N_time) = x_filt(:, N_time);
        for t = N_time-1:-1:1
            dx = x_smooth(:, t+1) - F * x_filt(:, t);
            x_smooth(:, t) = x_filt(:, t) + p_filt .* (F' * (dx ./ p_pred_safe));
        end

        % E. M-step: update diagonal Q from smoothed process increments.
        if N_time > 1
            x_t = x_smooth(:, 2:end);
            x_tm1 = x_smooth(:, 1:end-1);
            data_err_pow = sum((x_t - F * x_tm1) .^ 2, 2);
            n_transitions = N_time - 1;
            total_err_pow = n_transitions * diag_cov_err + data_err_pow;
        else
            n_transitions = 1;
            total_err_pow = theta;
        end

        if nd > 1
            loc_pow = sum(reshape(total_err_pow, nd, []), 1)';
            theta_loc = loc_pow / (nd * n_transitions);
            theta = repelem(theta_loc, nd);
        else
            theta = total_err_pow / n_transitions;
        end
        theta = max(real(theta), q_min);

        % F. Convergence check.
        theta_safe = max(theta, q_min);
        log_prior_Q = -0.5 * sum(log(theta_safe) + 1 ./ theta_safe);
        logPost = logL + log_prior_Q;
        logPost_hist(iter) = logPost;

        if iter > 1
            change = abs(logPost_hist(iter) - logPost_hist(iter-1)) / max(abs(logPost_hist(iter-1)), 1);
            fprintf('  Iter %d: LogPost = %.4e, Change = %.4e\n', iter, logPost, change);
            if change < tol
                break;
            end
        else
            fprintf('  Iter %d: LogPost = %.4e\n', iter, logPost);
        end
    end

    S_final = x_smooth;
    Param_dMAP.FinalSourceCovDiag = theta;
    Param_dMAP.LogPosteriorHistory = logPost_hist(1:iter);
    Param_dMAP.NumIterationsPerformed = iter;
    Param_dMAP.Approximation = 'diagonal_mean_field_steady_state_rts';
    Param_dMAP.InnovationCovariance = S;
end

function tf = local_is_valid_graph_input(x)
    tf = isstruct(x) || ((isnumeric(x) || islogical(x) || issparse(x)) && ismatrix(x));
end

function VertConn = local_extract_vertconn(inputArg)
    if isstruct(inputArg)
        if isfield(inputArg, 'VertConn') && ~isempty(inputArg.VertConn)
            inputArg = inputArg.VertConn;
        else
            error('seal_dMAP:MissingVertConn', ...
                'The third input must be VertConn or a struct containing .VertConn.');
        end
    end
    if isempty(inputArg) || ~ismatrix(inputArg)
        error('seal_dMAP:InvalidVertConn', 'VertConn must be a non-empty matrix.');
    end
    VertConn = sparse(inputArg);
    VertConn = spones(VertConn);
    VertConn = VertConn - spdiags(diag(VertConn), 0, size(VertConn, 1), size(VertConn, 2));
    VertConn = spones(VertConn | VertConn');
end

function F = local_construct_F(G_adj, lambda, alpha_n)
    [N_loc, ~] = size(G_adj);
    G_adj = spones(sparse(G_adj));
    F = spalloc(N_loc, N_loc, max(nnz(G_adj) + N_loc, N_loc));
    [row, col] = find(G_adj);
    for n = 1:N_loc
        neighbors = col(row == n);
        F(n, n) = lambda * alpha_n;
        if ~isempty(neighbors)
            neighbor_weight = lambda * (1 - alpha_n) / numel(neighbors);
            F(n, neighbors) = neighbor_weight;
        end
    end
end

function [K, p_filt, S, Ls, logdetS] = local_filter_step(L, R, p_pred, q_min, innovation_jitter)
    U = L' .* p_pred;
    S = L * U + R;
    [Ls, S, logdetS] = local_chol_spd(S, innovation_jitter);
    K = (Ls' \ (Ls \ U'))';
    p_filt = p_pred - p_pred .* sum(K .* L', 2);
    p_filt = max(real(p_filt), q_min);
end

function [Ls, S, logdetS] = local_chol_spd(S, innovation_jitter)
    n = size(S, 1);
    S = double((S + S') / 2);
    scale = trace(S) / max(n, 1);
    if ~isfinite(scale) || scale <= 0
        scale = 1;
    end
    loading = max(innovation_jitter * scale, eps(scale));
    eyeN = eye(n);
    for attempt = 1:8
        [Ls, flag] = chol(S, 'lower');
        if flag == 0 && all(diag(Ls) > 0)
            logdetS = 2 * sum(log(diag(Ls)));
            return;
        end
        S = S + loading * eyeN;
        loading = loading * 10;
    end
    error('seal_dMAP:InnovationCovarianceNotPD', ...
        'Innovation covariance is not positive definite after diagonal loading.');
end
