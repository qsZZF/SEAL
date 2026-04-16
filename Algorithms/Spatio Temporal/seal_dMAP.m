function [S_final, Param_dMAP] = seal_dMAP(X, L, G_adj, varargin)
% SEAL_DMAP Spatiotemporal dynamic ESI using a Diagonal Mean-Field Approximation.
%
%   Features a perfectly memory-optimized Diagonal Riccati tracker that preserves 
%   the lateral spatial spreading dynamics of the cortical graph (F). Memory is O(N).
% REQUIRED INPUTS
%   X       - [N_chan x N_time]   Sensor time series (EEG/MEG).
%   L       - [N_chan x N_src]    Lead field matrix. N_src = nd * N_loc,
%                                 where nd is the number of orientations
%                                 per cortical location.
%   G_adj   - [N_loc x N_loc]     Sparse cortical adjacency (undirected).
%                                 Defined on physical locations, not on
%                                 orientation-expanded sources.
%
% OPTIONAL NAME-VALUE PARAMETERS
%   'NoiseCovariance'      - [N_chan x N_chan]  Sensor noise covariance R.
%                                               Default: eye(N_chan).
%   'InitialSourceCovDiag' - [N_src x 1]        Initial diagonal of the
%                                               state noise covariance.
%                                               Default: 1e-3 * ones(N_src,1).
%   'NumOrientations'      - scalar             Dipole orientations per
%                                               location (1 = fixed,
%                                               3 = free). Default: 1.
%   'lambda_F'             - scalar in (0,1]    Temporal decay factor in F.
%                                               Default: 0.95.
%   'alpha_F'              - scalar in [0,1]    Self- vs. neighbor-weight
%                                               split in F (alpha = self,
%                                               1-alpha distributed to
%                                               neighbors). Default: 0.51.
%   'max_iter'             - scalar             Maximum EM iterations.
%                                               Default: 15.
%   'tol'                  - scalar             Relative tolerance on the
%                                               log-posterior for EM
%                                               convergence. Default: 1e-5.
%
% OUTPUTS
%   S_final    - [N_src x N_time]  Smoothed source time courses (posterior
%                                  mean from the RTS-style backward pass).
%   Param_dMAP - struct with fields:
%       .FinalSourceCovDiag    - [N_src x 1]  Learned diagonal of Q.
%       .LogPosteriorHistory   - [iter x 1]   Log-posterior trace over EM
%                                             iterations (truncated to the
%                                             actual number of iterations
%                                             performed).
    %% 1. Input Parsing
    p = inputParser;
    [N_chan, ~] = size(X);
    [~, N_src] = size(L);

    addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'G_adj', @(x) issparse(x) && ismatrix(x));
    addParameter(p, 'NoiseCovariance', eye(N_chan), @ismatrix);
    addParameter(p, 'InitialSourceCovDiag', 1e-3 * ones(N_src, 1), @isvector);
    addParameter(p, 'NumOrientations', 1, @isscalar);
    addParameter(p, 'lambda_F', 0.95, @isscalar);
    addParameter(p, 'alpha_F', 0.51, @isscalar);
    addParameter(p, 'max_iter', 15, @isscalar);
    addParameter(p, 'tol', 1e-5, @isscalar);
    parse(p, X, L, G_adj, varargin{:});

    max_iter = p.Results.max_iter;
    tol = p.Results.tol;
    nd = p.Results.NumOrientations;

    %% 2. Initialization
    [~, N_time] = size(X);
    R = p.Results.NoiseCovariance;
    
    fprintf('dMAP: Constructing state transition matrix F...\n');
    F_loc = local_construct_F(G_adj, p.Results.lambda_F, p.Results.alpha_F);
    if nd > 1
        F = kron(F_loc, speye(nd));
    else
        F = F_loc;
    end
    F_sq = F .^ 2; % For rapid diagonal variance propagation
    
    theta_loc = mean(reshape(p.Results.InitialSourceCovDiag, nd, []), 1)';
    theta = repelem(theta_loc, nd);
    x0 = zeros(N_src, 1); 
    
    logPost_hist = zeros(max_iter, 1);
    fprintf('Starting Diagonal Mean-Field dMAP (nd=%d)...\n', nd);

    %% 3. Expectation-Maximization (EM) Loop
    for iter = 1:max_iter
        
        % --- A. E-Step: Diagonal Steady-State Filter ---
        
        p_pred = theta;
        % 1. Forward Diagonal Riccati
        for i = 1:100
            S = L * (L' .* p_pred) + R; % O(N) memory
            invS = inv(S);
            K = (L' .* p_pred) * invS; 
            
            p_filt = p_pred - sum(K .* (L .* p_pred')', 2);
            p_pred_new = F_sq * p_filt + theta;
            
            if norm(p_pred_new - p_pred) / max(norm(p_pred), eps) < 1e-5, break; end
            p_pred = p_pred_new;
        end
        
        % 2. Backward Diagonal Riccati
        J_sq = (F' .^ 2) .* (p_filt ./ p_pred);
        p_smooth = p_filt;
        for i = 1:100
            p_smooth_new = p_filt + J_sq * (p_smooth - p_pred);
            if norm(p_smooth_new - p_smooth) / max(norm(p_smooth), eps) < 1e-5, break; end
            p_smooth = p_smooth_new;
        end
        
        % 3. Extract Covariance Error Term Analytically
        diag_cross = p_smooth .* (F_sq' * p_filt) ./ p_pred;
        diag_cov_err = p_smooth - 2 * diag_cross + F_sq * p_smooth;
        
        % 4. Fast Forward Data Pass
        x_filt = zeros(N_src, N_time);
        x_k = x0;
        logL = 0;
        for t = 1:N_time
            v_k = X(:, t) - L * x_k;
            x_filt_k = x_k + K * v_k;
            x_filt(:, t) = x_filt_k;
            x_k = F * x_filt_k; 
            logL = logL - 0.5 * (log(det(S)) + v_k' * invS * v_k + N_chan*log(2*pi));
        end
        
        % 5. Fast Backward Data Pass
        x_smooth = zeros(N_src, N_time);
        x_smooth(:, N_time) = x_filt(:, N_time);
        for t = N_time-1:-1:1
            dx = x_smooth(:, t+1) - F * x_filt(:, t);
            % Diagonal J approximation
            x_smooth(:, t) = x_filt(:, t) + (p_filt ./ p_pred) .* (F' * dx);
        end
        
        % --- B. M-Step: Update State Noise Covariance Q ---
        
        x_t = x_smooth(:, 2:end);
        x_tm1 = x_smooth(:, 1:end-1);
        data_err_pow = sum((x_t - F * x_tm1).^2, 2);
        
        total_err_pow = (N_time-1) * diag_cov_err + data_err_pow;
        
        % Group Sparsity Pooling
        if nd > 1
            loc_pow = sum(reshape(total_err_pow, nd, []), 1)';
            theta_loc = loc_pow / (nd * N_time);
            theta = repelem(theta_loc, nd);
        else
            theta = total_err_pow / N_time;
        end
        
        % Convergence Check
        log_prior_Q = -0.5 * sum(log(theta) + 1./(theta + eps));
        logPost = logL + log_prior_Q;
        logPost_hist(iter) = logPost;

        if iter > 1
            change = abs(logPost_hist(iter) - logPost_hist(iter-1)) / abs(logPost_hist(iter-1));
            fprintf('  Iter %d: LogPost = %.4e, Change = %.4e\n', iter, logPost, change);
            if change < tol, break; end
        else
            fprintf('  Iter %d: LogPost = %.4e\n', iter, logPost);
        end
    end
    
    S_final = x_smooth;
    Param_dMAP.FinalSourceCovDiag = theta;
    Param_dMAP.LogPosteriorHistory = logPost_hist(1:iter);
end

function F = local_construct_F(G_adj, lambda, alpha_n)
    [N_src, ~] = size(G_adj);
    F = spalloc(N_src, N_src, N_src * 7); 
    [row, col] = find(G_adj);
    for n = 1:N_src
        neighbors = col(row == n);
        F(n, n) = lambda * alpha_n;
        if ~isempty(neighbors)
            d_ni = 1 / length(neighbors); 
            for i_idx = 1:length(neighbors)
                F(n, neighbors(i_idx)) = lambda * (1 - alpha_n) * d_ni;
            end
        end
    end
end