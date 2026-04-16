function [Source_Estimate, Parameters] = seal_FASTIRES(Data, L, sourceSpaceInfo, varargin)
%SEAL_FASTIRES Solves the ESI problem using the FAST-IRES algorithm.
% min alpha*||W_1 J||_1 + ||W_2 VJ||_1
%   s.t. trace[(b-LJA)'Cnoise^{-1}(b-LJA)]<= beta^2, A is temporal basis
%   function, J is spatial coefficients
%   [Source_Estimate, Parameters] = seal_FASTIRES(Data, L, sourceSpaceInfo, 'ParameterName', ParameterValue, ...)
%
%   Implements the FAST-IRES algorithm (modified from the code by Sohrabpour et al., 2017)
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data.
%       L (double matrix): Nchannels x Nsources lead field matrix.
%       sourceSpaceInfo (struct): Struct with source space geometry. Required:
%           .EdgeOperator (sparse matrix): Nedges x Nsources_loc edge operator (V).
%
%   Optional Name-Value Pair Inputs:
%       'NoiseCovariance' (double matrix): Initial noise covariance matrix. Default: [].
%       'InitialTBFs' (matrix): Initial temporal basis functions (K x T). Default: [].
%       'NumTBFs' (integer): Number of TBFs to extract via SVD if InitialTBFs is empty. Default: 4.
%       'Alpha' (double or vector): Regularization parameter(s) for source sparsity. Default: 0.08.
%       'Lambda' (double): ADMM balancing parameter. Default: computed from data.
%       'BettaTilda' (double): Noise power constraint. Default: computed from SNR.
%       'Epsilon' (double): Reweighting stabilization parameter. Default: 1e-3.
%       'MaxIterations' (integer): Maximum reweighting iterations. Default: 10.
%       'MaxInnerIterationsX' (integer): Maximum iterations for x-update. Default: 20.
%       'MaxInnerIterationsY' (integer): Maximum iterations for y-update. Default: 20.
%       'MaxInnerIterationsProj' (integer): Maximum iterations for projection step. Default: 20.
%       'NumOrientations' (integer): Number of dipole orientations (1 or 3). Default: 3.
%       'Tolerance' (double): Convergence tolerance for inner loops. Default: 1e-4.
%
%   Outputs:
%       Source_Estimate (double matrix): Nsources x Ntimepoints estimated source activity.
%       Parameters (struct): Contains:
%           .InverseOperator (matrix): Spatial mapping matrix.
%           .EdgeEstimate (matrix): Edge-based sparse solution (y).
%           .LearnedAlpha (vector): Final regularization parameters.
%           .LearnedBettaTilda (double): Final noise power constraint.
%           .LearnedLambda (double): Final ADMM balancing parameter.
%           .CostHistory (vector): Cost function values per iteration.
%           .NumIterationsPerformed (integer): Total reweighting iterations.
%
%   Reference: Sohrabpour et al., "Noninvasive Electromagnetic Source Imaging of Spatio-temporally Distributed Epileptogenic Brain Sources", 2020.

%% --- Input Parsing ---
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;
addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);

addParameter(p, 'NoiseCovariance', [], @(x) ismatrix(x) && issymmetric(x));
addParameter(p, 'InitialTBFs', [], @(x) ismatrix(x));
addParameter(p, 'NumTBFs', 4, @(x) isscalar(x) && x > 0);
addParameter(p, 'Alpha', 0.08, @(x) isvector(x) && all(x > 0));
addParameter(p, 'Lambda', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(p, 'BettaTilda', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(p, 'Epsilon', 1e-3, @(x) isscalar(x) && x > 0);
addParameter(p, 'MaxIterations', 10, @(x) isscalar(x) && x > 0);
addParameter(p, 'MaxInnerIterationsX', 20, @(x) isscalar(x) && x > 0);
addParameter(p, 'MaxInnerIterationsY', 20, @(x) isscalar(x) && x > 0);
addParameter(p, 'MaxInnerIterationsProj', 20, @(x) isscalar(x) && x > 0);
addParameter(p, 'NumOrientations', 3, @(x) isscalar(x) && ismember(x, [1, 3]));
addParameter(p, 'Tolerance', 1e-4, @(x) isscalar(x) && x > 0);

try
    parse(p, Data, L, sourceSpaceInfo, varargin{:});
catch ME
    disp('Error parsing inputs for seal_FASTIRES:');
    rethrow(ME);
end

Parameters.OptionsPassed = p.Results;
maxIter = Parameters.OptionsPassed.MaxIterations;
numItX = Parameters.OptionsPassed.MaxInnerIterationsX;
numItY = Parameters.OptionsPassed.MaxInnerIterationsY;
numItProj = Parameters.OptionsPassed.MaxInnerIterationsProj;
epsilon = Parameters.OptionsPassed.Tolerance;
nd = Parameters.OptionsPassed.NumOrientations;
alpha = Parameters.OptionsPassed.Alpha;
lambda = Parameters.OptionsPassed.Lambda;
betta_tilda = Parameters.OptionsPassed.BettaTilda;
eps_reweight = Parameters.OptionsPassed.Epsilon;
Cnoise = Parameters.OptionsPassed.NoiseCovariance;
Phi = Parameters.OptionsPassed.InitialTBFs;
numTBFs = Parameters.OptionsPassed.NumTBFs;

%% --- Initialization ---
[nSensor, nSources_total] = size(L);
nSamp = size(Data, 2);
nSources_loc = nSources_total / nd;

if mod(nSources_total, nd) ~= 0
    error('Lead field columns (%d) must be divisible by NumOrientations (%d).', nSources_total, nd);
end

% Initialize noise covariance
if isempty(Cnoise)
    % Robust diagonal noise estimate from MAD
    sigma_hat = median(abs(Data - median(Data,2)), 2) / 0.6745;
    sigma_hat = max(sigma_hat, 1e-8);
    Cnoise = diag(sigma_hat.^2);
end

Cnoise = (Cnoise + Cnoise') / 2;
Sigma_inv_half = chol(inv(Cnoise + 1e-8 * eye(nSensor)), 'lower');

% Initialize TBFs
if isempty(Phi)
    [~, ~, tD] = svd(Sigma_inv_half * Data, 'econ');
    Phi = tD(:, 1:min(numTBFs, size(tD, 2)))';
end
numTBFs = size(Phi, 1);

% Compute normalized data and lead field
Data_norm = Sigma_inv_half * Data;
L_norm = Sigma_inv_half * L;

% Edge operator
V = sourceSpaceInfo.EdgeOperator;
if size(V, 2) ~= nSources_total
    error('EdgeOperator must have %d columns for %d source locations.', nSources_loc, nSources_loc);
end
% V = kron(V, speye(nd));
nEdges = size(V, 1);

% Precompute matrices
C_t = Phi * Phi' + 1e-8 * eye(numTBFs);
TBF_norm = Phi.'*(C_t\Phi);
T_norm = pinv(C_t) * Phi;
W_v = V' * V;
L_v = 1.01 * svds(W_v, 1);
Abbas = L_norm * L_norm';
[U, D, ~] = svd(Abbas);
if any(diag(D) < -1e-10)
    error('Lead field covariance must be positive semi-definite.');
end
K_U = U' * L_norm;

% Noise power (betta_tilda)
% if isempty(betta_tilda)
%     Noise = Data_norm(:, 1:floor(nSamp/2));
%     SNR = norm(Data_norm, 'fro')^2 / norm(Noise, 'fro')^2;
%     power = sum(Noise.^2, 1);
%     [X, B] = hist(power, 100);
%     Sum_X = cumsum(X) / sum(X);
%     ind_90 = find(Sum_X > 0.90, 1);
%     beta = icdf('chi2', 0.95, nSensor-1) * numTBFs;
%     CF_min = beta / (norm(Data_norm, 'fro')^2 / nSources_total);
%     correction_factor = max(beta / B(ind_90), CF_min);
%     betta_tilda = norm(Data_norm, 'fro')^2 / SNR / correction_factor;
% end
if isempty(betta_tilda)
    dof_eff = nSensor * nSamp;
    conf_level = 0.99;
    betta_tilda = chi2inv(conf_level, dof_eff);
end
% Lambda
if isempty(lambda)
    lambda = 1.01 * max(1, svds(L_norm, 1)^2 * sqrt(nSensor) / max(norms(L_norm' * Data_norm)));
end

% Initialize solution
D_W = 1 ./ (norms(L_norm) + 1e-20);
x_0 = (L_norm' * pinv(L_norm * L_norm' + mean(diag(D)) * eye(nSensor)) * (Data_norm * Phi'/C_t)) .* repmat(D_W', 1, numTBFs);
W = ones(nSources_total, numTBFs);
W_d = ones(nEdges, numTBFs);

% Alpha selection
if isscalar(alpha)
    alpha_vec = alpha;
else
    alpha_vec = alpha;
end
num_alpha = numel(alpha_vec);

%% --- Main Iterative Loop ---

J_sol = zeros(nSources_total, numTBFs, maxIter, num_alpha);
Y_sol = zeros(nEdges, numTBFs, maxIter, num_alpha);
Parameters.LearnedAlpha = zeros(num_alpha, 1);

fprintf('Running FAST-IRES algorithm...\n');
for i_alpha = 1:num_alpha
    costlist = [];
    alpha = alpha_vec(i_alpha);
    x_0_iter = x_0;
    W_iter = W;
    W_d_iter = W_d;
    weight_it = 0;
    stop_itr = false;

    while ~stop_itr
        weight_it = weight_it + 1;
        if weight_it > maxIter
            stop_itr = true;
        end

        % FISTA-ADMM Core
        y_cond = true;
        y_it = 0;
        x_tild = x_0_iter;
        x_old = x_0_iter;
        y = V * x_tild;
        u = zeros(size(y));
        e_rel = 1e-4;
        e_abs = 1e-8;
        tau = 2;
        mu = 10;

        while y_cond
            y_it = y_it + 1;
            if y_it > numItY
                y_cond = false;
            end

            t_old = 1;
            x_it = 0;
            x_cond = true;

            while x_cond
                x_it = x_it + 1;
                if x_it > numItX
                    x_cond = false;
                end

                % x-update with soft thresholding
                Alpha = diag(sum(abs(x_0_iter)) / max(sum(abs(x_0_iter))));
                x_new_tran = wthresh(x_tild - (W_v * x_tild - V' * (y + u)) / L_v, 's', (W_iter * Alpha) * (alpha / lambda / L_v));

                % Projection step
                phi_hat = Data_norm - L_norm * x_new_tran * Phi;
                Phi_i = (norms((U' * phi_hat * TBF_norm).').').^2;
                Tp_norm = phi_hat*T_norm';
                N_Phi = phi_hat - phi_hat * TBF_norm;
                N_f = norm(N_Phi, 'fro')^2;

                if norm(phi_hat, 'fro')^2 <= betta_tilda
                    x_new_er = zeros(size(x_new_tran));
                else
                    x_new_er = Proj_Flat_Hyper_Ellips(max(betta_tilda - N_f, 0), Phi_i, phi_hat, U, D, K_U, Tp_norm, numItProj);
                end
                x_new = x_new_tran + x_new_er;

                % FISTA acceleration
                t_new = 0.5 * (1 + sqrt(1 + 4 * t_old^2));
                x_tild = x_new + ((t_old - 1) / t_new) * (x_new - x_old);

                if norm(x_new - x_old, 'fro') / norm(x_new, 'fro') < epsilon
                    x_cond = false;
                end
                x_old = x_new;
                t_old = t_new;
            end

            % y-update
            y_new = wthresh(V * x_new - u, 's', W_d_iter / lambda);

            % Convergence check
            prim_res = norm(y_new - V * x_new, 'fro');
            dual_res = norm(lambda * (V' * (y - y_new)), 'fro');
            e_prim = sqrt(numel(y)) * e_abs + e_rel * max(sum(norms(V * x_new)), sum(norms(y_new)));
            e_dual = sqrt(numel(x_new)) * e_abs + e_rel * sum(norms(V' * y_new));

            if prim_res <= e_prim && dual_res <= e_dual
                y_cond = false;
            end

            y = y_new;
            u = u + y - V * x_new;

            % Dynamic lambda update
%             if prim_res > mu * dual_res
%                 lambda = lambda * tau;
%                 u = u / tau;
%             elseif dual_res > mu * prim_res
%                 lambda = lambda / tau;
%                 u = u * tau;
%             end
        end

        % Reweighting
        J_sol(:, :, weight_it, i_alpha) = x_new;
        Y_sol(:, :, weight_it, i_alpha) = y;

        J_n = reshape(x_new, [nd, nSources_loc, numTBFs]);
        Ab_J = squeeze(norms(J_n)) + 1e-20;
        Y_n = reshape(V * x_new, [nd, nEdges/nd, numTBFs]);
        Ab_Y = squeeze(norms(Y_n)) + 1e-20;

        W_old = W_iter;
        W_tr = 1 ./ (Ab_J ./ repmat(max(Ab_J), [nSources_loc, 1]) + eps_reweight);
%         W_iter = zeros(nSources_total, numTBFs);
%         for d = 1:nd
%             W_iter(d:nd:end, :) = W_tr;
%         end
        W_iter = repelem(W_tr,nd,1);
        
        W_d_old = W_d_iter;
        W_d_tr = 1 ./ (Ab_Y ./ repmat(max(Ab_Y), [nEdges/nd, 1]) + eps_reweight);
%         W_d_iter = zeros(nEdges, numTBFs);
%         for d = 1:nd
%             W_d_iter(d:nd:end, :) = W_d_tr;
%         end
        W_d_iter = repelem(W_d_tr,nd,1);

        % Cost function
        Residual = Data_norm - L_norm * x_new * Phi;
        cost = 0.5 * norm(Residual, 'fro')^2 + alpha * sum(sum(Ab_Y)) + sum(sum(Ab_J));
        costlist = [costlist, cost];

        % Stopping criterion
        if norm(W_iter - W_old, 'fro') / norm(W_old, 'fro') < epsilon && ...
           norm(W_d_iter - W_d_old, 'fro') / norm(W_d_old, 'fro') < epsilon
            stop_itr = true;
        end

        x_0_iter = x_new;
    end

    Parameters.LearnedAlpha(i_alpha) = alpha;
    Parameters.CostHistory{i_alpha} = costlist;
    x = wthresh(x_new, 's', (W_iter * Alpha) * (alpha / lambda / L_v));
    y = wthresh(V * x_new - u, 's', W_d_iter / lambda);
    Source_Estimate{i_alpha} = x * Phi;
    Parameters.InverseOperator{i_alpha} = x;
    Parameters.EdgeEstimate{i_alpha} = y;
    if num_alpha == 1
        Source_Estimate = x * Phi;
        Parameters.LearnedAlpha = alpha;
        Parameters.InverseOperator = x;
        Parameters.EdgeEstimate = y;
    end
    

    
end

%% --- Final Outputs ---

Parameters.LearnedBettaTilda = betta_tilda;
Parameters.LearnedLambda = lambda;
Parameters.LearnedNoiseCovariance = Cnoise;
Parameters.TemporalBasisFunctions = Phi;
Parameters.NumIterationsPerformed = weight_it;
Parameters.OptionsPassed = p.Results;



end

function [E] = Proj_Flat_Hyper_Ellips(betta_tilda, Phi_i, phi_hat, U, D, K_U, Tp_norm, num_it)
% Performs projection onto a flat hyper-ellipsoid for FAST-IRES.
% Inputs: betta_tilda (noise constraint), Phi_i (normed projection), phi_hat (residual),
%         U, D (SVD of L*L'), K_U (U'*L), Tp_norm (TBF projection), num_it (iterations).
% Output: E (correction term for J).
myfun = @(lambda,Phi,D_0,Betta) sum(Phi./(1+lambda*D_0).^2) - Betta;

Phi       = Phi_i;
D_0       = diag(D);
Betta     = betta_tilda;
Init_zero = 0; 

fun   = @(lambda) myfun(lambda,Phi,D_0,Betta);
lambda_ini = max(fzero(fun,Init_zero),0);

E = (K_U.'*diag(lambda_ini./(1+lambda_ini*diag(D)))*U.')*Tp_norm;
end

function n = norms(X)
% Compute vector norms for each dipole orientation
    n = sqrt(sum(X.^2, 1));
end

function out = Soft(th, input)
% Soft thresholding function
    out = sign(input) .* max(abs(input) - th, 0);
end