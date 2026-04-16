function [Source_Lq, Param_Lq] = seal_Lq_ADMM(Data, L, varargin)
%SEAL_LQ_ADMM Computes Lq-regularized source imaging using Spatio-Temporal ADMM.
%   [Source, Param] = seal_Lq_ADMM(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Solves the sparse inverse problem using Alternating Direction Method of 
%   Multipliers (ADMM). Uses a mixed L2,q norm across orientations AND time 
%   (or Temporal Basis Functions) to enforce spatio-temporal group sparsity.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints measured data.
%       L (double matrix): Nchannels x Nsources lead field.
%
%   Optional Name-Value Pair Inputs:
%       'Phi' (double matrix): K x Ntimepoints Temporal Basis Functions. If 
%                              provided, the data is projected onto these bases 
%                              before solving the L2,q penalty. Default: [].
%       'DataFidelity' (char): 'L2' (Standard Least Squares) or 
%                              'L1' (Robust Absolute Deviations). Default: 'L2'.
%       'q_norm' (double): The q parameter for Lq norm (0 <= q <= 1). Default: 0.5.
%       'RegularizationParameter' (double): Lambda. Default: 0.1.
%       'ADMM_Rho' (double): Penalty parameter rho. Default: 1.0 (L2) or 100 (L1).
%       'NumOrientations' (integer): Orientations per location. Default: 1.
%       'MaxIterations' (integer): Maximum ADMM iterations. Default: 500.
%       'Tolerance' (double): Convergence tolerance. Default: 1e-5.

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    
    addParameter(p, 'Phi', [], @(x) isempty(x) || ismatrix(x));
    addParameter(p, 'DataFidelity', 'L2', @(x) ischar(x) || isstring(x));
    addParameter(p, 'q_norm', 0.5, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    addParameter(p, 'RegularizationParameter', 0.1, @isscalar);
    addParameter(p, 'ADMM_Rho', [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'NumOrientations', 1, @isscalar);
    addParameter(p, 'MaxIterations', 500, @isscalar);
    addParameter(p, 'Tolerance', 1e-5, @isscalar);

    parse(p, Data, L, varargin{:});
    opts = p.Results;

    fidelity = upper(opts.DataFidelity);
    q = opts.q_norm;
    lamda = opts.RegularizationParameter;
    nd = opts.NumOrientations;
    Phi = opts.Phi;
    MAX_ITER = opts.MaxIterations;
    TOL = opts.Tolerance;
    
    [Nchannels, Ntimepoints] = size(Data);  

    nSources_total = size(L, 2);
    nSources_loc = nSources_total / nd;

    %% 3. Temporal Basis Projection (TBF)
    if ~isempty(Phi)
        % Project Data: Y_tilde = Y * pinv(Phi)
        Y_tilde = Data * pinv(Phi);
        K_dim = size(Phi, 1);
        fprintf('Projected data onto %d Temporal Basis Functions.\n', K_dim);
    else
        % Use raw time points
        Y_tilde = Data;
        K_dim = Ntimepoints;
    end

    %% 4. L2 Data Fidelity Solver (Least Squares + Spatio-Temporal L2,q)
    if strcmp(fidelity, 'L2')
        if isempty(opts.ADMM_Rho), rho = 1; else, rho = opts.ADMM_Rho; end
        
        fprintf('Running Spatio-Temporal L2,q-L2 ADMM (q=%g, lambda=%g, nd=%d, K=%d)...\n', q, lamda, nd, K_dim);
        
        % Woodbury Identity Optimization (Channel Space Inversion)
        inv_core = inv((rho/2) * eye(Nchannels) + L * L');
        
        % Initialize Matrix States (No time loop needed)
        Lty = L' * Y_tilde;
        X = zeros(nSources_total, K_dim);
        Z = zeros(nSources_total, K_dim);
        U = zeros(nSources_total, K_dim);
        
        for iter = 1:MAX_ITER
            Zm1 = Z;
            
            % X-update (Woodbury Direct Solve on all bases simultaneously)
            RHS = 2 * Lty + rho * (Z - U);
            X = (1/rho) * RHS - (1/rho) * L' * (inv_core * (L * RHS));
            
            % Z-update (Spatio-Temporal L2,q Group Shrinkage)
            Z = spatio_temporal_shrinkage(X + U, q, lamda, rho, nd, nSources_loc, K_dim);
            
            % U-update (Dual ascent)
            U = U + (X - Z);
            
            % Convergence Check
            if (norm(rho * (Z - Zm1), 'fro') < sqrt(nSources_total * K_dim) * TOL) && ...
               (norm(X - Z, 'fro') < sqrt(nSources_total * K_dim) * TOL)
                Param_Lq.ConvergedAt = iter;
                break;
            end
        end
        W_est = Z;

    %% 5. L1 Data Fidelity Solver (Robust Absolute Deviation + Spatio-Temporal L2,q)
    elseif strcmp(fidelity, 'L1')
        if isempty(opts.ADMM_Rho), rho = 100; else, rho = opts.ADMM_Rho; end
        
        fprintf('Running Spatio-Temporal L2,q-L1 ADMM (q=%g, lambda=%g, nd=%d, K=%d)...\n', q, lamda, nd, K_dim);
        
        % Spectral Normalization
        scale_factor = 1.01 * norm(L, 2);
        Ls = L / scale_factor;
        Ys = Y_tilde / scale_factor;
        
        ep = 1e-3;
        tao1 = 0.99; 
        tao2 = 2 * ep;
        rho_cov = 3.2 / ep / lamda;
        
        % Initialize Matrix States
        X = zeros(nSources_total, K_dim);
        W_dual = zeros(Nchannels, K_dim);
        V = zeros(Nchannels, K_dim);
        current_rho = rho;
        
        for iter = 1:MAX_ITER
            if (current_rho < rho_cov && q < 1)
                current_rho = current_rho * 1.01;
            end
            Xm1 = X;
            
            % X-step (Proximal Gradient)
            grad_Z = X - tao1 * (Ls' * (Ls * X - Ys - V - W_dual / current_rho));
            X = spatio_temporal_shrinkage(grad_Z, q, tao1, current_rho, nd, nSources_loc, K_dim);
            
            % V-step (Huber-like smoothing)
            d2 = V ./ sqrt(V.^2 + ep^2);
            V = 1 / (1/tao2 + current_rho * lamda) * ...
                (1/tao2 * V - d2 + current_rho * lamda * (Ls * X - Ys - W_dual / current_rho));
            
            % W-step (Dual ascent)
            LsX = Ls * X;
            W_dual = W_dual - current_rho * (LsX - Ys - V);
            
            % Convergence Check
            if (norm(X - Xm1, 'fro') < sqrt(nSources_total * K_dim) * TOL) && ...
               (norm(LsX - Ys - V, 'fro') < sqrt(Nchannels * K_dim) * TOL)
                Param_Lq.ConvergedAt = iter;
                break;
            end
        end
        W_est = X;
    end

    %% 6. Final State Assembly
    if ~isempty(Phi)
        Source_Lq = W_est * Phi;
    else
        Source_Lq = W_est;
    end
    
    Param_Lq.BasisCoefficients = W_est;
    Param_Lq.OptionsPassed = opts;
end

%% --- Helper Functions ---

function [Z_new] = spatio_temporal_shrinkage(Z, q, lam, rho_penalty, nd, nLoc, K_dim)
% Applies L2,q group shrinkage across both Orientation (nd) and Time/Basis (K_dim)
    
    % Reshape to [Orientations x Locations x Time/Basis]
    Z_tensor = reshape(Z, nd, nLoc, K_dim);
    
    % Calculate L2 norm (energy) across BOTH orientations and time for each location
    % Output size: [1 x nLoc x 1]
    Z_norms_sq = sum(sum(Z_tensor.^2, 1), 3);
    Z_norms = sqrt(Z_norms_sq(:)); 
    
    % Apply 1D Lq shrinkage to the scalar magnitude of each location
    X_norms = scalar_shrinkage_Lq(Z_norms, q, lam, rho_penalty);
    
    % Compute the isotropic scaling factor for the surviving locations
    scale = X_norms ./ max(Z_norms, eps);
    
    % Apply the scale uniformly across the orientation and time axes
    scale_tensor = reshape(scale, 1, nLoc, 1);
    Z_tensor_scaled = bsxfun(@times, Z_tensor, scale_tensor);
    
    % Flatten back to [Nsources x K_dim]
    Z_new = reshape(Z_tensor_scaled, nLoc * nd, K_dim);
end

function [x] = scalar_shrinkage_Lq(b, q, lam, rho_penalty)
% Evaluates the 1D proximal operator for the Lq norm.
    if (q == 0) % L0 Hard Thresholding
        x = b;
        i1 = abs(b) <= sqrt(2 * lam / rho_penalty);
        x(i1) = 0;
        
    elseif (q < 1 && q > 0) % Fractional Lq
        max_iter = 10;
        ABSTOL   = 1e-4;
        x    = zeros(length(b), 1);
        
        beta = (rho_penalty / (lam * q * (1 - q)))^(1 / (q - 2));
        f1   = lam * q * beta^(q - 1) + rho_penalty * beta - rho_penalty * abs(b);
        i0   = find(f1 < 0);
        
        if ~isempty(i0) 
            b_u = abs(b(i0));
            x_u = b_u;
            for i = 1:max_iter              
                deta_x = (lam * q * x_u.^(q-1) + rho_penalty * x_u - rho_penalty * b_u) ./ ...
                         (lam * q * (q-1) * x_u.^(q-2) + rho_penalty);
                x_u = x_u - deta_x;
                if (norm(deta_x) < sqrt(length(x_u)) * ABSTOL)
                    break;
                end
            end
            x_u = x_u .* sign(b(i0));
            
            % Check objective value to ensure global minimum
            obj_zero = rho_penalty/2 * b_u.^2;
            obj_val = lam * abs(x_u).^q + rho_penalty/2 * (x_u - b(i0)).^2;
            
            i1 = find(obj_zero - obj_val < 0);
            x_u(i1) = 0;
            x(i0) = x_u;
        end 
        
    elseif (q == 1) % L1 Soft Thresholding
        x = sign(b) .* max(abs(b) - lam / rho_penalty, 0);
    else
        error('q must be between 0 and 1');
    end
end