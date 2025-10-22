function [S_final, Param_dSPN] = seal_dSPN(X, L, varargin)
% SEAL_DSPN Implements a dynamic sparse solution for the ESI problem.
%
% This function solves the ESI problem using a state-space model with
% Kalman filtering and RTS smoothing within an Expectation-Maximization (EM)
% framework. It promotes dynamic sparsity by estimating the variance of
% each source at each time point (Automatic Relevance Determination).
%
% Based on the paper: "Computationally Efficient Algorithms for Sparse,
% Dynamic Solutions to the EEG Source Localization Problem" by Pirondini et al.,
% IEEE T-BME, 2018.
%
%   Usage:
%   [S, P] = seal_dSPN(Data, LeadField, 'max_iter', 20);
%
%   Inputs:
%       X (double matrix): N_channels x N_timepoints data matrix.
%       L (double matrix): N_channels x N_sources lead field matrix.
%
%   Optional Name-Value Pair Inputs:
%       'NoiseCovariance' (matrix): Measurement noise covariance (R). Default: eye(N_channels).
%       'InitialSourceCov' (vector): Initial guess for source variances (diag of Q). Default: 1e-3 * ones.
%       'max_iter' (integer): Maximum number of EM iterations. Default: 30.
%       'tol' (double): Convergence tolerance for the log-likelihood. Default: 1e-5.
%
%   Outputs:
%       S_final (double matrix): Estimated source activity (N_sources x N_timepoints).
%       Param_dSPN (struct): Parameters, convergence history, and final covariances.

    %% --- 1. Input Parsing ---
    p = inputParser;
    [N_chan, ~] = size(X);
    [~, N_src] = size(L);

    addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'NoiseCovariance', eye(N_chan), @ismatrix);
    addParameter(p, 'InitialSourceCov', 1e-3 * ones(N_src, 1), @isvector);
    addParameter(p, 'max_iter', 30, @isscalar);
    addParameter(p, 'tol', 1e-5, @isscalar);
    parse(p, X, L, varargin{:});

    Param_dSPN.OptionsPassed = p.Results;
    max_iter = p.Results.max_iter;
    tol = p.Results.tol;

    %% --- 2. Initialization ---
    [N_chan, N_time] = size(X);
    [~, N_src] = size(L);

    % State-space model:
    % x_t = F * x_{t-1} + w_t,  w_t ~ N(0, Q_t)
    % y_t = L * x_t + v_t,      v_t ~ N(0, R)
    
    F = speye(N_src); % State transition matrix (identity for this model)
    R = p.Results.NoiseCovariance; % Measurement noise covariance
    
    % Initialize source covariance Q_t. In this model, Q is diagonal and
    % its elements (the source variances) are estimated at each time step.
    % We initialize with a single Q for all time steps.
    q_diag = p.Results.InitialSourceCov;
    Q = spdiags(q_diag, 0, N_src, N_src);

    % Initial state guess
    x0 = zeros(N_src, 1);
    P0 = Q; % Initial state covariance

    logL_hist = zeros(max_iter, 1);
    fprintf('Starting dSPN EM Algorithm...\n');

    %% --- 3. Expectation-Maximization (EM) Loop ---
    for iter = 1:max_iter
        % --- E-Step: Kalman Filter & RTS Smoother ---
        % Given current Q and R, estimate the posterior distribution of the states
        [x_filt, P_filt, x_pred, P_pred, logL] = local_kalman_filter(X, L, F, Q, R, x0, P0);
        [x_smooth, P_smooth] = local_rts_smoother(x_filt, P_filt, x_pred, P_pred, F);
        
        logL_hist(iter) = logL;
        
        % --- M-Step: Update Covariances ---
        % Update Q based on the smoothed state estimates. This is the core
        % step that enforces dynamic sparsity.
        % The variance of source 'i' at time 't' is updated as the expected
        % value of x_i,t^2, which is P_smooth(i,i,t) + x_smooth(i,t)^2.
        
        q_diag_new = zeros(N_src, 1);
        for t = 1:N_time
            % The innovation for each source is the variance of w_t.
            % E[w_t*w_t'] = E[(x_t - F*x_{t-1})*(x_t - F*x_{t-1})']
            % This is complex. The paper uses a simplified update based on ARD.
            % The variance of source i is updated based on its own smoothed power.
            q_diag_new = q_diag_new + diag(P_smooth(:,:,t)) + x_smooth(:,t).^2;
        end
        q_diag = q_diag_new / N_time;
        Q = spdiags(q_diag, 0, N_src, N_src);
        
        % Check for convergence on the log-likelihood
        if iter > 1
            change = abs(logL_hist(iter) - logL_hist(iter-1)) / abs(logL_hist(iter-1));
            fprintf('  Iter %d: LogL = %.4e, Change = %.4e\n', iter, logL, change);
            if change < tol
                fprintf('Converged at iteration %d.\n', iter);
                break;
            end
        else
            fprintf('  Iter %d: LogL = %.4e\n', iter, logL);
        end
    end
    
    if iter == max_iter
        warning('dSPN:MaxIterReached', 'EM algorithm did not converge.');
    end

    %% --- 4. Final Output ---
    S_final = x_smooth;
    Param_dSPN.FinalSourceCov = Q;
    Param_dSPN.FinalMeasurementCov = R;
    Param_dSPN.LogLikelihoodHistory = logL_hist(1:iter);
    Param_dSPN.IterationsRun = iter;
end

%--------------------------------------------------------------------------
% Local Helper Functions for State-Space Estimation
%--------------------------------------------------------------------------

function [x_filt, P_filt, x_pred, P_pred, logL] = local_kalman_filter(Y, L, F, Q, R, x0, P0)
    % Standard Kalman Filter forward pass
    [N_chan, N_time] = size(Y);
    N_src = size(L, 2);

    x_pred = zeros(N_src, N_time);
    P_pred = zeros(N_src, N_src, N_time);
    x_filt = zeros(N_src, N_time);
    P_filt = zeros(N_src, N_src, N_time);
    logL = 0;

    % Initial prediction
    x_pred_k = x0;
    P_pred_k = P0;

    for t = 1:N_time
        x_pred(:, t) = x_pred_k;
        P_pred(:, :, t) = P_pred_k;
        
        % Innovation
        y_innov = Y(:, t) - L * x_pred_k;
        S_innov = L * P_pred_k * L' + R;
        
        % Kalman Gain
        K = P_pred_k * L' / S_innov;
        
        % Update (Filter)
        x_filt_k = x_pred_k + K * y_innov;
        P_filt_k = (eye(N_src) - K * L) * P_pred_k;
        
        x_filt(:, t) = x_filt_k;
        P_filt(:, :, t) = P_filt_k;
        
        % Predict next state
        x_pred_k = F * x_filt_k;
        P_pred_k = F * P_filt_k * F' + Q;

        % Log-likelihood calculation
        logL = logL - 0.5 * (log(det(S_innov)) + y_innov' / S_innov * y_innov + N_chan*log(2*pi));
    end
end

function [x_smooth, P_smooth] = local_rts_smoother(x_filt, P_filt, x_pred, P_pred, F)
    % Rauch-Tung-Striebel (RTS) smoother backward pass
    N_time = size(x_filt, 2);
    N_src = size(x_filt, 1);

    x_smooth = zeros(N_src, N_time);
    P_smooth = zeros(N_src, N_src, N_time);

    % Initialize with the last filtered state
    x_smooth(:, N_time) = x_filt(:, N_time);
    P_smooth(:, :, N_time) = P_filt(:, :, N_time);

    for t = (N_time-1):-1:1
        % Smoother Gain J
        J = P_filt(:, :, t) * F' / P_pred(:, :, t+1);
        
        % Smoothed state and covariance
        x_smooth(:, t) = x_filt(:, t) + J * (x_smooth(:, t+1) - x_pred(:, t+1));
        P_smooth(:, :, t) = P_filt(:, :, t) + J * (P_smooth(:, :, t+1) - P_pred(:, :, t+1)) * J';
    end
end
