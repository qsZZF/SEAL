function [Source_eLORETA, Param_eLORETA] = seal_eLORETA(Data, L, varargin)
%SEAL_ELORETA Computes the exact LORETA (eLORETA) estimate using an iterative approach.
%   [Source_eLORETA, Param_eLORETA] = seal_eLORETA(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Implements the iterative eLORETA algorithm to find the optimal source
%   weighting matrix W, and then computes the eLORETA inverse solution.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured EEG/MEG data.
%       L (double matrix): Nchannels x Nsources original lead field matrix (K in paper).
%
%   Optional Name-Value Pair Inputs:
%       'RegularizationParameter' (double): Alpha (alpha) regularization parameter for sensor space.
%                                           Default: 1e-4 (typically small for eLORETA).
%       'NoiseCovariance' (double matrix): Noise covariance (C_noise). If provided and not identity,
%                                          data and lead field are whitened. The regularization
%                                          alpha is then applied to identity in whitened space.
%                                          Default: eye(Nchannels).
%       'NumOrientations' (integer): Number of orientations per source location (e.g., 1 or 3).
%                                    Default: 1.
%       'MaxIterations' (integer): Maximum number of iterations for W convergence. Default: 100.
%       'Tolerance' (double): Convergence tolerance for W. Default: 1e-6.
%   Outputs:
%       Source_eLORETA (double matrix): Nsources x Ntimepoints eLORETA source estimates.
%       Param_eLORETA (struct):
%           .InverseOperator (double matrix): The final eLORETA inverse operator (T_eLORETA)
%                                             that maps original Data to Source_eLORETA.
%           .WeightMatrix_W (sparse block-diag matrix): Converged eLORETA weight matrix W.
%                                                      (Diagonal blocks W_j)
%           .NumIterationsPerformed (integer)
%           .WhiteningMatrix (matrix): iW used.
%           .RegularizationAlpha (double): Alpha value used.
%           .OptionsPassed (struct)
%
%   See also: SEAL_MNE.

%   Author: FengZhao

%% Input Parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true; % Allow other params to be passed but not used by this core function

defaultRegParamAlpha = 1e-4;
if size(L,1) > 0, defaultNoiseCov = eye(size(L,1)); else, defaultNoiseCov = []; end
defaultNumOrientations = 1;
defaultMaxIter = 100;
defaultTol = 1e-6;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

addParameter(p, 'RegularizationParameter', defaultRegParamAlpha, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x) && ismember(x,[1,3]));
addParameter(p, 'MaxIterations', defaultMaxIter, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'Tolerance', defaultTol, @(x) isnumeric(x) && isscalar(x) && x>0);

try
    parse(p, Data, L, varargin{:});
catch ME
    disp('Error parsing inputs for seal_eLORETA:');
    rethrow(ME);
end

Param_eLORETA.OptionsPassed = p.Results;
Data_orig = p.Results.Data;
L_orig = p.Results.L;
alpha_reg = p.Results.RegularizationParameter;
C_noise = p.Results.NoiseCovariance;
nd = p.Results.NumOrientations;
max_iter = p.Results.MaxIterations;
tolerance = p.Results.Tolerance;

[Nchannels, ~] = size(Data_orig);
[Nchannels_L, Nsources_total] = size(L_orig);

if isempty(C_noise), C_noise = eye(Nchannels); Param_eLORETA.OptionsPassed.NoiseCovariance = C_noise; end
if Nchannels ~= Nchannels_L, error('seal_eLORETA:DimensionMismatch', 'Data and L must have the same number of channels.'); end
if mod(Nsources_total, nd) ~= 0, error('seal_eLORETA:NumOrientationsMismatch', 'Total sources in L not divisible by NumOrientations.'); end
Nsources_locations = Nsources_total / nd;

%% Noise Covariance Whitening
Data_proc = Data_orig;

iW = speye(Nchannels); % Whitening matrix
L_w = L_orig;          % Whitened Lead Field (K in paper's whitened equations)
Data_w = Data_proc;      % Whitened Data (Phi in paper's whitened equations)

if ~isequal(C_noise, eye(Nchannels))
    try
        [U_c, S_c_diag_vec] = eig((C_noise+C_noise')/2); S_c_diag_vec = diag(S_c_diag_vec);
        eig_tol_c = Nchannels*eps(max(S_c_diag_vec)); S_c_diag_vec(S_c_diag_vec < eig_tol_c) = eig_tol_c;
        iW_dense = diag(1./sqrt(S_c_diag_vec)) * U_c';
        iW = sparse(iW_dense); % Keep sparse if possible
        
        L_w = iW * L_orig;
        Data_w = iW * Data_proc;
    catch ME_w
        warning('seal_eLORETA:WhiteningFailed', 'Whitening failed: %s. Using unwhitened data/lead field.', ME_w.message);
        iW = speye(Nchannels); L_w = L_orig; Data_w = Data_proc;
    end
end
alpha_reg = max(alpha_reg * trace(L_w * L_w')/Nchannels, 1e-6);
Param_eLORETA.WhiteningMatrix = iW;
Param_eLORETA.RegularizationAlpha = alpha_reg;

%% Iterative eLORETA Algorithm
if nd == 1
    % Initialize W as Identity
    W_inv = ones(1, Nsources_locations);
    W_old = W_inv;
    for iter = 1:max_iter
        M_inv = bsxfun(@times, L_w, W_inv) * L_w' + alpha_reg * eye(Nchannels);
        if rank(full(M_inv)) < Nchannels % Check rank before inv
            M = pinv(M_inv);
        else
            M = inv(M_inv);
        end
        W = sqrt(sum(L_w .* (M * L_w), 1));
        delta_W = norm(W - W_old) / (norm(W_old) + eps);
        W_inv = 1 ./ W;
        W_old = W;
        if iter > 1 && delta_W < tolerance
            fprintf('eLORETA converged at iteration %d (delta_W = %g)\n', iter, delta_W);
            break;
        end
        if iter == max_iter
            warning('seal_eLORETA:MaxIterReached', 'eLORETA W convergence reached max iterations (%d). Delta_W = %g.', max_iter, delta_W);
        end
    end
    W_inv = W_inv';
    Param_eLORETA.WeightMatrix_W = sparse(diag(W));
elseif nd > 1
    % W is block-diagonal. W_j are nd x nd blocks.
    % Initialize W as Identity;
    W_inv = ones(Nsources_total, nd);
    W_old = W_inv;
    W = W_old;
    for iter = 1:max_iter
        KW_inv_K_T = zeros(Nchannels, Nchannels);
        % Construct KW_inv_iter_K_T = L_w * W_inv_iter * L_w'
        for i_loc = 1:Nsources_locations
            curIdx = (i_loc-1)*nd + (1:nd);
            KW_inv_K_T = KW_inv_K_T + L_w(:,curIdx) * W_inv(curIdx,:) * L_w(:,curIdx)';
        end
        M_inv = KW_inv_K_T + alpha_reg * eye(Nchannels);
        if rank(full(M_inv)) < Nchannels % Check rank before inv
            M = pinv(M_inv);
        else
            M = inv(M_inv);
        end
        % Update W_j blocks
        for j_loc = 1:Nsources_locations
            curIdx = (i_loc-1)*nd + (1:nd);
            Term_for_sqrtm = L_w(:,curIdx)' * M * L_w(:,curIdx);
            Term_for_sqrtm = (Term_for_sqrtm + Term_for_sqrtm')/2; % Ensure symmetry for sqrtm
            try
                W_j = sqrtm(Term_for_sqrtm);
            catch ME_sqrtm
                warning('seal_eLORETA:sqrtmFailed', 'sqrtm failed for block %d in iter %d: %s. Using regularized pseudo-inverse sqrt.', i_loc, iter, ME_sqrtm.message);
                % Fallback for sqrtm if matrix is not SPD enough
                [V_sqrt, D_sqrt] = eig(Term_for_sqrtm);
                d_sqrt = diag(D_sqrt);
                d_sqrt(d_sqrt < 0) = eps; % Ensure non-negative before sqrt
                W_j = V_sqrt * diag(sqrt(d_sqrt)) * V_sqrt';
            end
            W(curIdx,:) = W_j;
            if rank(W_j) < nd
                W_inv(curIdx,:) = pinv(W_j);
            else
                W_inv(curIdx,:) = inv(W_j);
            end
        end
        delta_W = norm((W - W_old), 'fro') / norm(W_old, 'fro') < tolerance;
        if iter > 1 && delta_W < tolerance
            fprintf('eLORETA converged at iteration %d (delta_W = %g)\n', iter, delta_W);
            break;
        end
        if iter == max_iter
            warning('seal_eLORETA:MaxIterReached', 'eLORETA W convergence reached max iterations (%d). Delta_W = %g.', max_iter, delta_W);
        end
        
    end
    % Construct block-diagonal W for output
    row_idx = []; col_idx = []; val_w = [];
    for i_loc = 1:Nsources_locations
        curIdx = (i_loc-1)*nd + (1:nd);
        [r_loc, c_loc, v_loc] = find(W(curIdx,:));
        row_idx = [row_idx; curIdx(r_loc)'];
        col_idx = [col_idx; curIdx(c_loc)'];
        val_w = [val_w; v_loc];
    end
    Param_eLORETA.WeightMatrix_W = sparse(row_idx, col_idx, val_w, Nsources_total, Nsources_total);
    
end

Param_eLORETA.NumIterationsPerformed = iter;




%% Construct Final eLORETA Inverse Operator (T_eLORETA from paper Eq. 37)

% T_w = W_inv * L_w' * inv(K_w * W_inv * K_w' + alpha_reg * I) (use the latest W_inv for simplicity)
% Pre-allocate T_w (Nsources_total x Nchannels)
T_w = sparse(Nsources_total, Nchannels);
for i_loc = 1:Nsources_locations
    curIdx = (i_loc-1)*nd + (1:nd);
    W_inv_j_block = W_inv(curIdx,:); % nd x nd    
    T_w_block = W_inv_j_block * L_w(:,curIdx)' * M; % nd x Nchannels
    T_w(curIdx, :) = T_w_block;
end

%% Source Estimation
Source_eLORETA = T_w * Data_w;

%% Final Inverse Operator (maps original Data_orig to Source_eLORETA)
Param_eLORETA.InverseOperator = T_w * iW;

end
