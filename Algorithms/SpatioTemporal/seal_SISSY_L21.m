function [S_final, Param_SISSY] = seal_SISSY_L21(X, G, V, varargin)
% SEAL_SISSY_L21 Implements the L1,2-SISSY algorithm using ADMM.
%
% This function solves the optimization problem:
%  min_S 0.5*||X - GS||_F^2 + lambda * (||VS||_1,2 + alpha*||S||_1,2)
%
% Based on the paper: "SISSY: An efficient and automatic algorithm for the
% analysis of EEG sources based on structured sparsity" by H. Becker et al.,
% NeuroImage, 2017. 
%
%   Usage:
%   [S_final, Param_SISSY] = seal_SISSY_L21(X, G, V, 'lambda', 10, 'alpha', 0.1);
%
%   Inputs:
%       X (double matrix): N_channels x N_timepoints data matrix.
%       G (double matrix): N_channels x N_sources lead field matrix.
%       V (sparse matrix): N_edges x N_sources matrix computing the spatial
%                          gradient on the cortical mesh. 
%
%   Optional Name-Value Pair Inputs:
%       'lambda' (double): Overall regularization parameter. Default: 1.0.
%       'alpha' (double): Balances between source and gradient sparsity. Default: 0.1.
%       'rho' (double): ADMM penalty parameter. Default: 1.0.
%       'max_iter' (integer): Maximum number of ADMM iterations. Default: 200.
%       'tol' (double): Convergence tolerance for ADMM. Default: 1e-4.
%
%   Outputs:
%       S_final (double matrix): Estimated source activity (N_sources x N_timepoints).
%       Param_SISSY (struct): Parameters and convergence history.

    %% --- 1. Input Parsing ---
    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'G', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'V', @(x) isnumeric(x) && (ismatrix(x) || issparse(x)));

    addParameter(p, 'lambda', 1.0, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'alpha', 0.1, @(x) isnumeric(x) && isscalar(x) && x>=0);
    addParameter(p, 'rho', 1.0, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'max_iter', 200, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'tol', 1e-4, @(x) isnumeric(x) && isscalar(x) && x>0);

    parse(p, X, G, V, varargin{:});

    lambda = p.Results.lambda;
    alpha = p.Results.alpha;
    rho = p.Results.rho;
    max_iter = p.Results.max_iter;
    tol = p.Results.tol;

    Param_SISSY.OptionsPassed = p.Results;

    %% --- 2. Initialization ---
    [~, D] = size(G);
    [E, T] = size(V);

    % Initialize primal and dual variables
    S = zeros(D, T);
    Y = zeros(E, T); % Corresponds to VS
    Z = zeros(D, T); % Corresponds to S
    U = zeros(E, T); % Lagrange multiplier for Y
    W = zeros(D, T); % Lagrange multiplier for Z

    % Pre-compute matrix inverse for the S-update step
    % P = G'G + rho*(V'V + I)
    % We use the matrix inversion lemma for efficiency if N < D
    fprintf('Initializing ADMM...\n');
    I_D = speye(D);
    P_inv = (G'*G + rho*(V'*V + I_D)) \ I_D;
    
    % Pre-compute terms for the S-update
    GTX = G' * X;

    history.primal_res_norm = zeros(max_iter, 1);
    history.dual_res_norm = zeros(max_iter, 1);
    history.eps_pri = zeros(max_iter, 1);
    history.eps_dual = zeros(max_iter, 1);

    %% --- 3. ADMM Iterations ---
    fprintf('Starting ADMM iterations...\n');
    for k = 1:max_iter
        % Store previous values for convergence check
        Y_old = Y;
        Z_old = Z;

        % --- S-update ---
        % Solves: argmin_S 0.5||X-GS||^2 + 0.5*rho*||VS - (Y-U)||^2 + 0.5*rho*||S - (Z-W)||^2
        % This is a standard least-squares problem with the solution:
        % S = (G'G + rho*V'V + rho*I)^-1 * (G'X + rho*V'(Y-U) + rho*(Z-W))
        % S = P_inv * (GTX + rho * V' * (Y - U) + rho * (Z - W));
        
        S_update_rhs = GTX + rho * V' * (Y - U) + rho * (Z - W);
        S = P_inv * S_update_rhs;
        
        % --- Y-update (for gradient) ---
        % Solves: argmin_Y lambda*||Y||_1,2 + 0.5*rho*||Y - (VS+U)||^2
        Y = prox_L21(V*S + U, lambda/rho);

        % --- Z-update (for sources) ---
        % Solves: argmin_Z lambda*alpha*||Z||_1,2 + 0.5*rho*||Z - (S+W)||^2
        Z = prox_L21(S + W, lambda*alpha/rho);

        % --- U and W update (dual variables) ---
        U = U + V*S - Y;
        W = W + S - Z;

        % --- Convergence Check ---
        % Following Boyd et al. (2010), Section 3.3.1
        % Primal residual
        r_primal_Y = V*S - Y;
        r_primal_Z = S - Z;
        history.primal_res_norm(k) = sqrt(norm(r_primal_Y, 'fro')^2 + norm(r_primal_Z, 'fro')^2);
        
        % Dual residual
        s_dual_Y = -rho * V' * (Y - Y_old);
        s_dual_Z = -rho * (Z - Z_old);
        history.dual_res_norm(k) = sqrt(norm(s_dual_Y, 'fro')^2 + norm(s_dual_Z, 'fro')^2);
        
        % Tolerances
        history.eps_pri(k) = sqrt(2*D*T)*tol + tol * max(sqrt(norm(V*S,'fro')^2 + norm(S,'fro')^2), sqrt(norm(Y,'fro')^2 + norm(Z,'fro')^2));
        history.eps_dual(k) = sqrt(D*T)*tol + tol * sqrt(norm(rho*V'*U,'fro')^2 + norm(rho*W,'fro')^2);
        
        if mod(k,10) == 0
            fprintf('Iter %d: Primal Res = %.4e (tol %.4e), Dual Res = %.4e (tol %.4e)\n', ...
                k, history.primal_res_norm(k), history.eps_pri(k), ...
                history.dual_res_norm(k), history.eps_dual(k));
        end

        if (history.primal_res_norm(k) < history.eps_pri(k)) && (history.dual_res_norm(k) < history.eps_dual(k))
            fprintf('Converged at iteration %d.\n', k);
            break;
        end
    end

    if k == max_iter
        warning('SISSY:MaxIterReached', 'ADMM did not converge within the maximum number of iterations.');
    end
    
    S_final = Z; % The final estimate for S is the output of the proximal operator for the source term
    Param_SISSY.history = history;
    Param_SISSY.iterations_run = k;

end

function Y_out = prox_L21(Y_in, beta)
% Proximal operator for the L2,1-norm (group-wise soft thresholding).
% Solves: argmin_X beta*||X||_1,2 + 0.5*||X - Y_in||_F^2
% ||X||_1,2 is the sum of the L2-norms of the rows of X.
    
    Y_out = zeros(size(Y_in));
    
    % Calculate L2 norm for each row (across time)
    row_norms = sqrt(sum(Y_in.^2, 2));
    
    % Apply soft thresholding
    shrinkage_factor = 1 - beta ./ row_norms;
    shrinkage_factor(shrinkage_factor < 0) = 0; % Ensure factor is non-negative
    
    % Apply shrinkage to rows that are not zero
    active_rows = row_norms > 0;
    Y_out(active_rows, :) = Y_in(active_rows, :) .* shrinkage_factor(active_rows);
end