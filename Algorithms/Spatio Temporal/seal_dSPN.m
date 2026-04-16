function [S_final, Param_dSPN] = seal_dSPN(X, L, VertConn, varargin)
% SEAL_DSPN
% Low-memory approximate dynamic sparse ESI using a steady-state Kalman-like filter.
%
% Main features:
%   1) low-memory steady-state gain approximation
%   2) graph-regularized state transition
%   3) Q update based on process increments, not correction terms
%
% Notes:
%   - This is still an approximation to Kalman/EM, not a full covariance RTS-EM method.
%   - Memory footprint stays O(N_src * T + N_src * N_chan).

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;

    [N_chan0, ~] = size(X);
    [N_chanL0, N_src0] = size(L);

    addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'VertConn', @(x) isnumeric(x) || islogical(x));

    addParameter(p, 'NoiseCovariance', eye(N_chan0), @ismatrix);
    addParameter(p, 'InitialSourceCov', 1e-3 * ones(N_src0, 1), @isvector);
    addParameter(p, 'NumOrientations', 1, @isscalar);
    addParameter(p, 'max_iter', 30, @isscalar);
    addParameter(p, 'tol', 1e-5, @isscalar);
    addParameter(p, 'rho', 0.99, @isscalar);
    addParameter(p, 'q_min', 1e-10, @isscalar);
    addParameter(p, 'sigma_min', 1e-6, @isscalar);
    addParameter(p, 'riccati_iter', 200, @isscalar);
    addParameter(p, 'riccati_tol', 1e-5, @isscalar);
    addParameter(p, 'smooth_pass', true, @(x) islogical(x) || isnumeric(x));

    parse(p, X, L, VertConn, varargin{:});

    R = p.Results.NoiseCovariance;
    q0 = p.Results.InitialSourceCov(:);
    nd = p.Results.NumOrientations;
    max_iter = p.Results.max_iter;
    tol = p.Results.tol;
    rho = p.Results.rho;
    q_min = p.Results.q_min;
    sigma_min = p.Results.sigma_min;
    riccati_iter = p.Results.riccati_iter;
    riccati_tol = p.Results.riccati_tol;
    do_smooth = logical(p.Results.smooth_pass);

    %% 2. Dimension Checks
    [N_chan, N_time] = size(X);
    [N_chanL, N_src] = size(L);

    if N_chan ~= N_chanL
        error('seal_dSPN:DimensionMismatch', 'X and L channel dimensions do not match.');
    end
    if mod(N_src, nd) ~= 0
        error('seal_dSPN:OrientationMismatch', 'size(L,2) must be divisible by NumOrientations.');
    end
    if length(q0) ~= N_src
        error('seal_dSPN:InitialSourceCovSize', ...
            'InitialSourceCov must have length equal to size(L,2).');
    end
    if ~isequal(size(R), [N_chan, N_chan])
        error('seal_dSPN:NoiseCovarianceSize', ...
            'NoiseCovariance must be N_chan x N_chan.');
    end

    N_loc = N_src / nd;
    if ~isequal(size(VertConn), [N_loc, N_loc])
        error('seal_dSPN:VertConnSize', ...
            'VertConn must be N_loc x N_loc, where N_loc = size(L,2)/NumOrientations.');
    end

    %% 3. Graph Normalization with Zero-Degree Protection
    deg = sum(VertConn, 2);
    deg_safe = deg;
    deg_safe(deg_safe <= 0) = 1;   % avoid divide-by-zero for isolated nodes
    Fg_loc = bsxfun(@rdivide, VertConn, deg_safe) * 0.99;

    if nd > 1
        Fg = kron(Fg_loc, eye(nd));
    else
        Fg = Fg_loc;
    end

    %% 4. Initialization
    q_loc = mean(reshape(q0, nd, []), 1)';
    q = repelem(q_loc, nd);
    q = max(q, q_min);

    R = (R + R') / 2;
    R = R + sigma_min * eye(N_chan);

    % Initial state via regularized pseudoinverse
    lambda_init = trace(L * L') / max(N_chan * 100, 1);
    x0 = L' / (L * L' + lambda_init * eye(N_chan)) * X(:, 1);

    logL_hist = zeros(max_iter, 1);

    fprintf('Starting low-memory dSPN (nd=%d)...\n', nd);

    %% 5. EM-like Outer Loop
    for iter = 1:max_iter

        % ------------------------------------------------------------
        % A. Build transition matrix
        % ------------------------------------------------------------
        % Same design spirit as your original approximate SSM:
        %   F = rho*I + (1-rho)*Fg
        %
        % If you want a tunable scalar A again, it can be reintroduced later,
        % but for stability and low-memory consistency this form is cleaner.
        F = rho * eye(N_src) + (1 - rho) * Fg;

        % ------------------------------------------------------------
        % B. Steady-state low-rank Riccati iteration
        % ------------------------------------------------------------
        % U approximates P * L'
        U = L' .* q;   % N_src x N_chan

        for k = 1:riccati_iter
            S = L * U + R;
            S = (S + S') / 2;

            % stable solve through Cholesky
            [Ls, flag] = chol(S, 'lower');
            if flag ~= 0
                S = S + 1e-8 * eye(N_chan);
                [Ls, flag] = chol(S, 'lower');
                if flag ~= 0
                    error('seal_dSPN:InnovationCovarianceNotPD', ...
                        'Innovation covariance is not positive definite.');
                end
            end

            invS_R = Ls' \ (Ls \ R);
            U_new = rho^2 * U * invS_R + L' .* q;

            relU = norm(U_new - U, 'fro') / max(norm(U, 'fro'), eps);
            U = U_new;
            if relU < riccati_tol
                break;
            end
        end

        % final steady-state quantities
        S = L * U + R;
        S = (S + S') / 2;
        [Ls, flag] = chol(S, 'lower');
        if flag ~= 0
            S = S + 1e-8 * eye(N_chan);
            [Ls, flag] = chol(S, 'lower');
            if flag ~= 0
                error('seal_dSPN:InnovationCovarianceNotPD_Final', ...
                    'Final innovation covariance is not positive definite.');
            end
        end

        % K = U / S, solved stably
        K = (Ls' \ (Ls \ U'))';   % N_src x N_chan

        logdetS = 2 * sum(log(diag(Ls)));

        % ------------------------------------------------------------
        % C. Forward filtering
        % ------------------------------------------------------------
        x_pred = zeros(N_src, N_time);
        x_filt = zeros(N_src, N_time);
        innov = zeros(N_chan, N_time);

        x_prev = x0;
        logL = 0;

        for t = 1:N_time
            x_pred(:, t) = F * x_prev;
            innov(:, t) = X(:, t) - L * x_pred(:, t);
            x_filt(:, t) = x_pred(:, t) + K * innov(:, t);
            x_prev = x_filt(:, t);

            quad = innov(:, t)' * (Ls' \ (Ls \ innov(:, t)));
            logL = logL - 0.5 * (logdetS + quad + N_chan * log(2*pi));
        end

        % ------------------------------------------------------------
        % D. Low-memory backward mean smoothing
        % ------------------------------------------------------------
        if do_smooth
            x_smooth = zeros(N_src, N_time);
            x_smooth(:, N_time) = x_filt(:, N_time);

            % low-memory steady-state smoother gain approximation
            J = F - K * (L * F);

            for t = N_time-1:-1:1
                dx = x_smooth(:, t+1) - x_pred(:, t+1);
                x_smooth(:, t) = x_filt(:, t) + J * dx;
            end
        else
            x_smooth = x_filt;
        end

        logL_hist(iter) = logL;

        % ------------------------------------------------------------
        % E. Update Q from state increments (crucial fix)
        % ------------------------------------------------------------
        % Approximate process noise:
        %   w_t = x_t - F x_{t-1}
        %
        % This is much more consistent with the state model than using
        % correction terms based on innovation rhs.
        if N_time > 1
            Wproc = x_smooth(:, 2:end) - F * x_smooth(:, 1:end-1);
            q_update = mean(Wproc.^2, 2);
        else
            q_update = q;
        end

        % enforce group-sharing across orientations
        if nd > 1
            q_loc = mean(reshape(q_update, nd, []), 1)';
            q_loc = max(q_loc, q_min);
            q = repelem(q_loc, nd);
        else
            q = max(q_update, q_min);
        end

        % ------------------------------------------------------------
        % F. Update Sigma from smoothed measurement residual
        % ------------------------------------------------------------
        residual = X - L * x_smooth;
        sigma_diag = mean(residual.^2, 2);
        sigma_diag = max(sigma_diag, sigma_min);
        R = diag(sigma_diag);

        % ------------------------------------------------------------
        % G. Convergence
        % ------------------------------------------------------------
        if iter > 1
            rel_change = abs(logL_hist(iter) - logL_hist(iter-1)) / max(abs(logL_hist(iter-1)), 1);
            fprintf('  Iter %d: LogL = %.4e, Change = %.4e\n', iter, logL_hist(iter), rel_change);
            if rel_change < tol
                break;
            end
        else
            fprintf('  Iter %d: LogL = %.4e\n', iter, logL_hist(iter));
        end
    end

    %% 6. Output
    S_final = x_smooth;
    Param_dSPN.FinalSourceCov = spdiags(q, 0, N_src, N_src);
    Param_dSPN.FinalNoiseCov = R;
    Param_dSPN.KalmanGain = K;
    Param_dSPN.InnovationCov = S;
    Param_dSPN.LogLikelihoodHistory = logL_hist(1:iter);
end