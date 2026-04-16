function [Source_BESTIES, Param_BESTIES] = seal_BESTIES(Data, L, Phi, VertConn, varargin)
%SEAL_BESTIES Computes source imaging using MRF and Temporal Basis Functions.
%   [Source, Param] = seal_BESTIES(Data, L, Phi, VertConn, 'ParameterName', ParameterValue, ...)
%
%   Formulates the inverse problem using a hierarchical Bayesian model with
%   multi-orientation support. Enforces spatial smoothness using an MRF 
%   (VertConn) and temporal continuity via Temporal Basis Functions (Phi). 
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints measured data.
%       L (double matrix): Nchannels x Nsources lead field.
%       Phi (double matrix): K x Ntimepoints temporal basis functions. 
%                            Can be empty [] for spatial-only MRF.
%       VertConn (logical/double matrix): Nlocations x Nlocations adjacency matrix.
%
%   Optional Name-Value Pair Inputs:
%       'NumOrientations' (integer): Orientations per location. Default: 1.
%       'OrientationTransform' (matrix): Matrix B to rotate/constrain orientations.
%       'NoiseCovariance' (double matrix): Sensor noise covariance. Default: eye.
%       'MaxIterations' (integer): Maximum VB iterations. Default: 200.
%       'Tolerance' (double): Convergence tolerance. Default: 1e-6.
%       'LearnTau' (logical): Dynamically learn MRF spatial parameter. Default: true.
%       'InitialTau' (double): Initial tau (0 < tau < 1). Default: 0.99.
%       'ComputeFreeEnergy' (logical): Use the original ELBO / Free Energy 
%                                      for the cost function. Default: true.

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'Phi', @(x) isnumeric(x) || isempty(x));
    addRequired(p, 'VertConn', @(x) isnumeric(x) || islogical(x));

    addParameter(p, 'NumOrientations', 1, @isscalar);
    addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || ismatrix(x));
    addParameter(p, 'MaxIterations', 100, @isscalar);
    addParameter(p, 'Tolerance', 1e-3, @isscalar);
    addParameter(p, 'LearnTau', true, @islogical);
    addParameter(p, 'InitialTau', 0.99, @isscalar);
    addParameter(p, 'ComputeFreeEnergy', true, @islogical);

    parse(p, Data, L, Phi, VertConn, varargin{:});
    opts = p.Results;

    Data_orig = opts.Data;
    L_orig = opts.L;
    Phi_basis = opts.Phi;
    Adj_loc = opts.VertConn;
    nd = opts.NumOrientations;    
    [nChan, nSnap] = size(Data_orig);
    
 
    nSource = size(L_orig, 2);
    
    if isempty(Phi_basis), K = 0; else, K = size(Phi_basis, 1); end
    if isempty(opts.NoiseCovariance), C_noise = eye(nChan); else, C_noise = opts.NoiseCovariance; end

    %% 3. Initialization
    Cost = zeros(opts.MaxIterations, 1);
    Cost_old = 0;
    S_old = zeros(nSource, nSnap);
    
    W = zeros(nSource, K); W1 = W;
    V = zeros(nSource, nSnap); V1 = V;
    keeplist = (1:K)';
    c = ones(K, 1);
    
    % Safe Graph Normalization (1D Location Space)
    NVertConn = sum(Adj_loc, 2);
    NVertConn(NVertConn == 0) = 1; 
    N_loc = bsxfun(@rdivide, Adj_loc, NVertConn);
    
    % Expand Graph to 3D Orientation Space
    if nd > 1
        N = kron(N_loc, eye(nd));
    else
        N = N_loc;
    end
    
    tau = opts.InitialTau;
    M = speye(nSource) - tau * N;
    LL = L_orig / M; 
    
    % Initialize Precisions
    tr_phi = 0; if K > 0, tr_phi = trace(Phi_basis' * Phi_basis); end
    init_prec = 1.0 * trace(LL * LL') * (tr_phi + nSnap) / (nSnap * trace(C_noise));
    gamma_prec = init_prec * ones(nSource, 1);
    alpha_prec = init_prec * ones(nSource, 1);
    
    c = c .* norm(alpha_prec);
    alpha_prec = alpha_prec ./ norm(alpha_prec);

    %% 4. Variational Bayes Loop
    fprintf('Running BESTIES (nd=%d)...\n', nd);
    for iter = 1:opts.MaxIterations
        
        if ~isempty(Phi_basis)
            index = find(abs(1./c) > max(abs(1./c)) * 1e-3);
            c = c(index); keeplist = keeplist(index);
            Phi_basis = Phi_basis(index, :);
            W = W(:, index); W1 = W1(:, index);
            K = numel(index); C2 = cell(numel(c), 1);
        end

        if opts.LearnTau   
           M = speye(nSource) - tau * N;
           LL = L_orig / M;
        end    
        
        % W Step
        if ~isempty(Phi_basis)
            TrW = zeros(nSource, 1);
            for k = numel(c):-1:1
                C2{k} = LL .* repmat(1./(c(k)*alpha_prec)', nChan, 1) * LL' + C_noise / (Phi_basis(k,:) * Phi_basis(k,:)');
                residual = (Data_orig - L_orig * V) * Phi_basis(k,:)' - ...
                    L_orig * sum((W(:, setdiff(1:K, k)) .* repmat(Phi_basis(k,:) * Phi_basis(setdiff(1:K, k), :)', nSource, 1)), 2);
                Kernel = repmat(1./(c(k)*alpha_prec), 1, nChan) .* LL' / C2{k};
                W1(:, k) = Kernel * residual / (Phi_basis(k,:) * Phi_basis(k,:)');
                if opts.LearnTau
                    G = M \ Kernel * LL;
                    TrW = TrW + sum(N .* G', 2);
                end
            end           
            W = M \ W1;
        end
        
        % V Step
        C1 = LL .* repmat(1./gamma_prec', nChan, 1) * LL' + C_noise;
        if ~isempty(Phi_basis)
            V1 = repmat(1./gamma_prec, 1, nChan) .* LL' / C1 * (Data_orig - L_orig * W * Phi_basis);
        else 
            V1 = repmat(1./gamma_prec, 1, nChan) .* LL' / C1 * Data_orig;      
        end
        E = M \ (repmat(1./gamma_prec, 1, nChan) .* LL' / C1 * LL);
        V = M \ V1;
        
        % Gamma Step (Grouped by Location)
        rou = diag(LL' / C1 * LL);
        if nd > 1
            rou_loc = sum(reshape(rou, nd, []), 1)';
            V1_sq_loc = sum(reshape(sum(V1.^2, 2), nd, []), 1)';
            gamma_loc = sqrt(nSnap * rou_loc ./ max(V1_sq_loc, 1e-12));
            gamma_prec = reshape(repmat(gamma_loc', nd, 1), [], 1);
        else
            gamma_prec = sqrt(nSnap * rou ./ max(sum(V1.^2, 2), 1e-12));
        end
%         gamma_prec(gamma_prec > min(gamma_prec)*1e12) = min(gamma_prec)*1e12;

        % Alpha Step (Grouped by Location)
        if ~isempty(Phi_basis)
            kesi = zeros(nSource, 1);
            for k = 1:numel(c)
                kesi = kesi + diag(LL' / C2{k} * LL) / c(k);
            end
            if nd > 1
                kesi_loc = sum(reshape(kesi, nd, []), 1)';
                W1_sq_loc = sum(reshape(sum(W1 .* repmat(c', nSource, 1) .* W1, 2), nd, []), 1)';
                alpha_loc = sqrt(kesi_loc ./ max(W1_sq_loc, 1e-12));
                alpha_prec = reshape(repmat(alpha_loc', nd, 1), [], 1);
            else
                alpha_prec = sqrt(kesi ./ max(sum(W1 .* repmat(c', nSource, 1) .* W1, 2), 1e-12));
            end
%             alpha_prec(alpha_prec > min(alpha_prec)*1e12) = min(alpha_prec)*1e12;
            
            % C Step
            alpha_prec = alpha_prec ./ norm(alpha_prec);
            for k = 1:numel(c)
                omega = diag(LL' / C2{k} * LL)' * (1./alpha_prec);
                c(k) = sqrt(omega / max(W1(:, k)' .* alpha_prec' * W1(:, k), 1e-12));
            end
        end
        
        % Tau Step
        if opts.LearnTau
             if ~isempty(Phi_basis), W_mean = N * W; end
             V_mean = N * V;
             TrV = nSnap * sum(E' .* N, 2);
             
             if ~isempty(Phi_basis)
                 Num = sum(W .* repmat(c', nSource, 1) .* W_mean .* repmat(alpha_prec, 1, numel(c)), 2) - TrW + ...
                       sum(V .* V_mean .* repmat(gamma_prec, 1, nSnap), 2) - TrV;
                 Det = sum(repmat(c', nSource, 1) .* (W_mean.^2) .* repmat(alpha_prec, 1, numel(c)), 2) + ...
                       sum((V_mean.^2) .* repmat(gamma_prec, 1, nSnap), 2);
             else
                 Num = sum(V .* V_mean .* repmat(gamma_prec, 1, nSnap), 2) - TrV;
                 Det = sum((V_mean.^2) .* repmat(gamma_prec, 1, nSnap), 2);
             end
             tau = sum(Num) / max(sum(Det), 1e-12);
             tau = max(min(tau, 0.99), 0.05); 
        end 
        
        if ~isempty(Phi_basis), S_est = W * Phi_basis + V; else, S_est = V; end
        
        %% J. Convergence Check (Original ELBO / Free Energy)
        if ~opts.ComputeFreeEnergy
            MSE = norm(S_est - S_old, 'fro')^2 / max(norm(S_est, 'fro')^2, eps);
            S_old = S_est;
            Cost(iter) = MSE;
        else
            % Limit precisions to prevent log(0) or Inf inside logdet
            TOL_FE = 1e-14;
            alpha_fe = max(min(alpha_prec, 1/TOL_FE), TOL_FE);
            gamma_fe = max(min(gamma_prec, 1/TOL_FE), TOL_FE);
            
            temp = 0; logCWprior = 0; logCW = 0;
            if ~isempty(Phi_basis)
                for k = 1:K
                    temp = temp + c(k) * W1(:, k)' .* alpha_fe' * W1(:, k);
                    logCWprior = logCWprior + sum(log(c(k) * alpha_fe));
                    
                    x_phi = Phi_basis(k,:) * Phi_basis(k,:)';
                    LALT = LL .* repmat(1./(c(k)*alpha_fe)', nChan, 1) * LL';
                    
                    % Original Identity: |A + U W V'| = |W^(-1) + V'*A^(-1)*U| * |W| * |A|
                    logCW = logCW + robust_logdet(eye(nChan)/x_phi + LALT) + nChan*log(x_phi) + sum(log(c(k)*alpha_fe));
                end
            end
            
            logCV = robust_logdet(eye(nChan) + LL .* repmat(1./gamma_fe', nChan, 1) * LL') + sum(log(gamma_fe));
            logCVprior = sum(log(gamma_fe));

            % Final Original ELBO Formulation
            FreeEnergy = -norm(Data_orig - L_orig * S_est, 'fro')^2 ...
                         - trace(V1' .* repmat(gamma_fe', nSnap, 1) * V1) ...
                         - temp ...
                         + (logCWprior - logCW) ...
                         + nSnap * (logCVprior - logCV);
                     
            Cost(iter) = FreeEnergy;
            
            if iter > 1
                MSE = abs(FreeEnergy - Cost_old) / abs(FreeEnergy);
            else
                MSE = inf;
            end
            Cost_old = FreeEnergy;
            if mod(iter, 10) == 0
                fprintf('BESTIES iteration: %d, tau: %.4f, Diff: %g\n', iter, tau, MSE);
            end
        end
        
        if abs(MSE) < opts.Tolerance
            fprintf('BESTIES Converged at iteration: %d, tau: %.4f, Diff: %g\n', iter, tau, MSE);
            Cost = Cost(1:iter);
            break;
        end
    end
    
    Source_BESTIES = S_est;
    Param_BESTIES.BasisCoefficients_W = W;
    Param_BESTIES.SpatialPrecision_alpha = alpha_prec;
    Param_BESTIES.SourceErrorPrecision_gamma = gamma_prec;
    Param_BESTIES.TemporalPrecision_c = c;
    Param_BESTIES.CostFunction = Cost;
    Param_BESTIES.OptionsPassed = opts;
end

%% --- Helper Function ---
function ld = robust_logdet(C)
% Computes log(det(C)) robustly to safely replace spm_logdet
    try
        R = chol(C);
        ld = 2 * sum(log(diag(R)));
    catch
        % Fallback using SVD if matrix is not perfectly positive definite
        C_sym = (C + C') / 2;
        [~, S, ~] = svd(C_sym);
        s_diag = diag(S);
        ld = sum(log(max(s_diag, 1e-12)));
    end
end