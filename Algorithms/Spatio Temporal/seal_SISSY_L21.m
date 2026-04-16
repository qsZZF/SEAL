function [S_final, Param_SISSY] = seal_SISSY_L21(Data, L, sourceSpaceInfo, varargin)
% SEAL_SISSY_L21 Implements the L1,2-SISSY algorithm using ADMM (Dynamic Lambda).
%
% Solves: min_S 0.5*||X - LS||_F^2 + lambda * (||VS||_1,2 + alpha*||S||_1,2)
%
%   Inputs:
%       Data (double matrix): N_channels x N_timepoints data matrix.
%       L (double matrix): N_channels x N_sources lead field matrix.
%       sourceSpaceInfo (struct): Struct with source space geometry. Required:
%           .VertConn (sparse matrix): Nlocations x Nlocations vertex connectivity. 
%
%   Optional Name-Value Pair Inputs:
%       'lambda' (double): Base regularization parameter. Default: 1.0.
%       'DynamicLambda' (logical): If true, scales lambda iteratively based 
%                                  on the residual noise variance. Default: true.
%       'JointEstNum' (integer): Max outer loops for Dynamic Lambda. Default: 10.
%       'alpha' (double): Balances source and gradient sparsity. Default: 0.51.
%       'rho' (double): ADMM penalty parameter. Default: 1.0.
%       'NumOrientations' (integer): Orientations per location. Default: 1.
%       'max_iter' (integer): Maximum ADMM iterations per outer loop. Default: 200.
%       'tol' (double): Convergence tolerance for ADMM. Default: 1e-4.

    %% --- 1. Input Parsing ---
    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'sourceSpaceInfo', @isstruct);

    addParameter(p, 'lambda', 1, @isscalar);
    addParameter(p, 'DynamicLambda', true, @islogical);
    addParameter(p, 'JointEstNum', 10, @isscalar);
    addParameter(p, 'alpha', 0.5, @isscalar);
    addParameter(p, 'rho', 1.0, @isscalar);
    addParameter(p, 'NumOrientations', 1, @isscalar);
    addParameter(p, 'max_iter', 200, @isscalar);
    addParameter(p, 'tol', 1e-4, @isscalar);

    parse(p, Data, L, sourceSpaceInfo, varargin{:});
    
    lam_base = p.Results.lambda;
    dynamic_lambda = p.Results.DynamicLambda;
    joint_loops = p.Results.JointEstNum;
    alpha = p.Results.alpha;
    rho = p.Results.rho;
    nd = p.Results.NumOrientations;
    max_iter = p.Results.max_iter;
    tol = p.Results.tol;
    Param_SISSY.OptionsPassed = p.Results;
    [Nchannels, Ntimepoints] = size(Data);

    Nsources_total = size(L, 2);
    V_loc = VariationEdge(sourceSpaceInfo.VertConn);
    
    if nd > 1
        V = kron(V_loc, speye(nd));
    else
        V = V_loc;
    end
    [E_total, ~] = size(V);
    
    
    Nsources = Nsources_total / nd;
%     sigma_inv_sqrt = zeros(nd, nd, Nsources);
%     for s = 1:Nsources
%         idx = (s-1)*nd + 1 : s*nd;
%         Ls = L(:, idx);
%         
%         cov_s = (Ls' * Ls) / Nchannels; 
%         
%         [U_L, D] = eig(cov_s);
%         d = diag(D);
%         d(d < 1e-12) = 1e-12; % Prevent division by zero
%         W_s = U_L * diag(1./sqrt(d)) * U_L';
%         
%         L(:, idx) = Ls * W_s;
%         sigma_inv_sqrt(:,:,s) = W_s;
%     end
    for s = 1:Nsources
        idx = (s-1)*nd + 1 : s*nd;
        Ls = L(:, idx);
        L(:, idx) = L(:, idx) ./ norm(Ls', 'fro');
    end
    %% --- 3. Initialization & Woodbury Precomputation ---
    fprintf('Initializing SISSY ADMM (nd=%d)...\n', nd);
    
    Y = zeros(E_total, Ntimepoints); 
    Z = zeros(Nsources_total, Ntimepoints); 
    U = zeros(E_total, Ntimepoints); 
    W = zeros(Nsources_total, Ntimepoints); 

    % --- Woodbury Identity Optimization ---
    % Computed ONLY ONCE. Entirely independent of Lambda.
    Q = rho * (V'*V + speye(Nsources_total));
    dA = decomposition(Q, 'chol');
    
    Q_inv_L = dA \ L';
    C_core = (eye(Nchannels) + L * Q_inv_L);
    LTX = L' * Data;

    history.primal_res_norm = [];
    history.dual_res_norm = [];
    history.sigma_list = [];

    %% --- 4. Dynamic Lambda Outer Loop ---
    if ~dynamic_lambda
        joint_loops = 1;
    end
    
    % Initial lambda estimate assuming S=0 (Noise = Data)
    lam_new = lam_base * (norm(Data, 'fro') / sqrt(Nchannels * Ntimepoints));
    sigma_hat_pre = inf;
    total_iter_count = 0;

    fprintf('Starting Dynamic Lambda SISSY Optimization...\n');
    for j = 1:joint_loops
        
        lam_use = lam_new;
        
        % --- 5. Inner ADMM Loop (Warm Started) ---
        for k = 1:max_iter
            Y_old = Y;
            Z_old = Z;

            % S-update (Woodbury Fast Solve)
            S_rhs = LTX + rho * V' * (Y - U) + rho * (Z - W);
            Q_inv_rhs = dA \ S_rhs;
            S = Q_inv_rhs - Q_inv_L * (C_core \ (L * Q_inv_rhs));
            
            % Y-update (Gradient Sparsity) using dynamic lam_use
            Y = prox_L21_group(V*S + U, lam_use/rho, nd);

            % Z-update (Source Sparsity) using dynamic lam_use
            Z = prox_L21_group(S + W, lam_use*alpha/rho, nd);

            % Dual Updates
            U = U + V*S - Y;
            W = W + S - Z;

            % Convergence Check for Inner Loop
            r_primal_Y = V*S - Y;
            r_primal_Z = S - Z;
            prim_res = sqrt(norm(r_primal_Y, 'fro')^2 + norm(r_primal_Z, 'fro')^2);
            
            s_dual_Y = -rho * V' * (Y - Y_old);
            s_dual_Z = -rho * (Z - Z_old);
            dual_res = sqrt(norm(s_dual_Y, 'fro')^2 + norm(s_dual_Z, 'fro')^2);
            
            eps_pri = sqrt(2*Nsources_total*Ntimepoints)*tol + tol * max(sqrt(norm(V*S,'fro')^2 + norm(S,'fro')^2), sqrt(norm(Y,'fro')^2 + norm(Z,'fro')^2));
            eps_dual = sqrt(Nsources_total*Ntimepoints)*tol + tol * sqrt(norm(rho*V'*U,'fro')^2 + norm(rho*W,'fro')^2);
            
            history.primal_res_norm = [history.primal_res_norm; prim_res];
            history.dual_res_norm = [history.dual_res_norm; dual_res];

            if (prim_res < eps_pri) && (dual_res < eps_dual)
                break;
            end
        end
        total_iter_count = total_iter_count + k;

        % --- 6. Update Residual Noise and Lambda ---
        sigma_hat = norm(Data - L * S, 'fro') / sqrt(Nchannels * Ntimepoints);
        lam_new = lam_base * sigma_hat;
        history.sigma_list = [history.sigma_list; sigma_hat];
        
        fprintf('  Outer Loop %d | ADMM Iters: %d | Sigma: %.4f | Lambda: %.4f\n', j, k, sigma_hat, lam_use);

        % Convergence check for Dynamic Lambda
        if abs(sigma_hat - sigma_hat_pre) < 1e-4
            break;
        end
        sigma_hat_pre = sigma_hat;
    end

    fprintf('Optimization Finished. Total ADMM iterations: %d\n', total_iter_count);
%     for s = 1:Nsources
%         idx = (s-1)*nd + 1 : s*nd;
%         S(idx, :) = sigma_inv_sqrt(:,:,s) * S(idx, :);
%     end
    S_final = S; 
    Param_SISSY.history = history;
    Param_SISSY.FinalSigma = history.sigma_list(end);
    Param_SISSY.FinalLambda = lam_use;
    Param_SISSY.iterations_run = total_iter_count;

end

%% --- Helper Functions ---

function Y_out = prox_L21_group(Y_in, beta, nd)
% Spatio-Temporal Group Proximal Operator
    [N_total, T] = size(Y_in);
    N_locs = N_total / nd;
    
    Y_tensor = reshape(Y_in, nd, N_locs, T);
    loc_norms = sqrt(sum(sum(Y_tensor.^2, 1), 3)); 
    
    shrinkage = max(1 - beta ./ max(loc_norms, eps), 0);
    Y_out_tensor = bsxfun(@times, Y_tensor, shrinkage);
    
    Y_out = reshape(Y_out_tensor, N_total, T);
end

function V = VariationEdge(VertConn)
% VariationEdge computes the incidence matrix V for the cortical graph.
    Nsource = size(VertConn, 1);
    Nedge = numel(find(VertConn(:) ~= 0)) / 2; 
    
    V = sparse(Nedge, Nsource);
    edge_idx = 0;
    
    for i = 1:Nsource
        idx = find(VertConn(i, :) ~= 0);
        idx = idx(idx < i); 
        
        for j = 1:numel(idx)
            edge_idx = edge_idx + 1;
            V(edge_idx, i) = 1;
            V(edge_idx, idx(j)) = -1;
        end
    end
end