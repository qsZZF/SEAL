function [Source_MNE, Param_MNE] = seal_MNE(Data,L,varargin)
%SEAL_MNE Computes the L2 Minimum Norm Estimate (MNE).
%   [Source, Param] = seal_MNE(Data, L, 'ParameterName', ParameterValue, ...)
%
%   This function implements the basic L2 Minimum Norm Estimate for
%   electromagnetic source imaging.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured EEG/MEG data.
%       L (double matrix): Nchannels x Nsources matrix representing the
%                          lead field (gain matrix). Assumes sources are single orientation
%                          unless 'NumOrientations' is specified.
%
%   Optional Name-Value Pair Inputs:
%       'RegularizationParameter' (double or char/string): regularization value for Tikhonov regularization 
%                                           lambda (the value used in seal_MNE is lambda^2*trace(L*L')/trace(Noise_Cov)).
%                                           If 'auto', attempts to estimate from data.
%                                           Default: 0.1
%       'SourceCovariance' (double matrix or char/string):
%                                           - Nsources x Nsources prior source covariance matrix (R_s).
%                                           - 'none' or 'depth_weighted'. If 'depth_weighted', computes diagonal 
%                                             depth weighting, else use eye(Nsources). 
%                                           Default:'depth_weighted'.
%       'DepthWeightExp'(float): Coefficent for depth weights q, R_s_(i,i) = norm(L_(:,i),2)^q. 
%                                           Default: 0.5 
%       'NoiseCovariance' (double matrix): Nchannels x Nchannels noise covariance
%                                          matrix (C). Default: eye(Nchannels).
%       'NumOrientations' (integer): Number of orientations per source location (e.g., 1 for fixed, 3 for free).
%                                    If > 1, L should be Nchannels x (Nsources_locations * NumOrientations).
%                                    SourceCovariance 'auto' will produce a diagonal spatial covariance
%                                    which is then kron'd with eye(NumOrientations).
%                                    Default: 1.
%
%   Outputs:
%       Source (double matrix): Nsources x Ntimepoints matrix of estimated source activity.
%       Param (struct): Structure containing parameters and intermediate results:
%           .InverseOperator (double matrix): Nsources x Nchannels inverse operator (M_final)
%                                             to be applied to original (unwhitened) data.
%           .RegularizationParameter (double): The actual regularization value used.
%           .SourceCovariance (double matrix): The Nsources x Nsources source covariance (R_s) used.
%           .NoiseCovariance (double matrix): The Nchannels x Nchannels input noise covariance.
%           .WhiteningMatrix (double matrix): The whitening matrix (iW) applied to data and leadfield.
%                                             If NoiseCovariance is identity, this is identity.
%           .OptionsPassed (struct): Struct of all parsed input options.
%
%   Formulation:
%       If C is noise covariance, iW is whitening matrix (e.g., iW = chol(inv(C))').
%       Data_w = iW * Data
%       L_w = iW * L
%       Inverse Operator for whitened data: M_w = R_s * L_w' * inv(L_w * R_s * L_w' + lambda^2 * I)
%       Source = M_w * Data_w
%       Final Inverse Operator (for original data): M_final = M_w * iW
%
%   Example:
%       Nch = 32; Nsrc_loc = 100; Nt = 200; nd = 1; Nsrc = Nsrc_loc*nd;
%       Data = randn(Nch, Nt);
%       L = randn(Nch, Nsrc);
%       Cnoise = diag(rand(Nch,1)*0.5 + 0.1);
%
%       [SrcEst, P] = seal_MNE(Data, L, 'NoiseCovariance', Cnoise, ...
%                                       'RegularizationParameter', 0.05, ...
%                                       'SourceCovariance', 'auto', ...
%                                       'NumOrientations', nd);

%       Visualize
%
%   See also: 

%   Author: FengZhao


%% Check Input
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = false; % Set to true if you expect other params from a higher-level wrapper
% Define default values
defaultRegParam = 0.1;
defaultSourceCov = 'depth_weighted'; % Default to automatic depth weighting
defaultNoiseCov = eye(size(L, 1));
defaultNumOrientations = 1;
defaultDepthWeightExp = 1;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

addParameter(p, 'RegularizationParameter', defaultRegParam);
addParameter(p, 'DepthWeightExp', defaultDepthWeightExp, @(x) isnumeric(x));
addParameter(p, 'SourceCovariance', defaultSourceCov); % Can be matrix or string
addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x) && x>=1);

try
    parse(p, Data, L, varargin{:});
catch ME_solve
    disp('Error parsing inputs for seal_MNE:');
    rethrow(ME_solve);
end
%% Assign parsed inputs to local variables
Data_orig = p.Results.Data;
L_orig = p.Results.L;
lambda_param_input = p.Results.RegularizationParameter;
Rs_input = p.Results.SourceCovariance;
C_noise = p.Results.NoiseCovariance;
nd = p.Results.NumOrientations;
DepthExp = p.Results.DepthWeightExp;
Param_MNE.OptionsPassed = p.Results; % Store all parsed options

[Nchannels, Ntimepoints] = size(Data_orig);
[Nchannels_L, Nsources_total] = size(L_orig);

if Nchannels ~= Nchannels_L
    error('seal_MNE:DimensionMismatch', 'Data and L must have the same number of channels.');
end
if mod(Nsources_total, nd) ~= 0
    error('seal_MNE:NumOrientationsMismatch', 'Total number of sources in L is not divisible by NumOrientations.');
end
Nsources_locations = Nsources_total / nd;
%% Noise Covariance Whitening
Param_MNE.NoiseCovariance = C_noise;
iW = eye(Nchannels); % Default to identity (no whitening)
L_w = L_orig;
Data_w = Data_orig;

if ~isequal(C_noise, eye(Nchannels))
    try
        % Check if C_noise_input is positive definite
        [~, chol_flag] = chol(C_noise);
        if chol_flag ~= 0
            warning('seal_MNE:NoiseCovNotPosDef', 'NoiseCovariance is not positive definite. Using pseudo-inverse for whitening matrix.');
            % SVD based whitening for non-positive definite or ill-conditioned C
            [U_c, S_c, ~] = svd(C_noise);
            s_c = diag(S_c);
            tol_c = max(size(C_noise)) * eps(max(s_c));
            rank_c = sum(s_c > tol_c);
            s_inv_sqrt_c = zeros(size(s_c));
            s_inv_sqrt_c(1:rank_c) = 1./sqrt(s_c(1:rank_c));
            iW = diag(s_inv_sqrt_c) * U_c'; % More robust whitening matrix
        else
            [U_c, S_c_diag] = eig((C_noise_input + C_noise_input')/2); % Ensure symmetry
            S_c_diag = diag(S_c_diag);
            S_c_diag(S_c_diag < eps) = eps; % Floor small/negative eigenvalues
            iW = diag(1./sqrt(S_c_diag)) * U_c';
        end        
        L_w = iW * L_orig;
        Data_w = iW * Data_orig;
    catch ME_whitening
        warning('seal_MNE:WhiteningFailed', 'Failed to compute whitening matrix. Using identity. Error: %s', ME_whitening.message);
        iW = eye(Nchannels);
        L_w = L_orig;
        Data_w = Data_orig;
    end
end
Param_MNE.WhiteningMatrix = iW;
%% 2. Source Covariance (R_s)
if isnumeric(Rs_input) && ismatrix(Rs_input)
    % User provided Rs directly
    if (nd == 1 && isequal(size(Rs_input), [Nsources_total, Nsources_total])) || ...
            (nd > 1 && (isequal(size(Rs_input), [Nsources_locations, Nsources_locations]) || isequal(size(Rs_input), [Nsources_total, Nsources_total])))
        Rs_spatial_or_full = Rs_input;
    else
        error('seal_MNE:SourceCovDimMismatch', 'User-provided SourceCovariance has incorrect dimensions.');
    end
    
    if nd > 1 && isequal(size(Rs_spatial_or_full), [Nsources_locations, Nsources_locations])
        Rs = kron(Rs_spatial_or_full, eye(nd));
    else
        Rs = Rs_spatial_or_full; % Assumed to be full Nsources_total x Nsources_total
    end
    % Calculate weighted L_w
    Rc = chol(Rs,'lower');
    L_w = L_w*Rc;
    
elseif ischar(Rs_input) || isstring(Rs_input)
    switch lower(Rs_input)
        case 'none'
            Rs = speye(Nsources_total);
        case 'depth_weighted'         
            Lw_norm_sq = sum(L_orig.^2,1);   % based on original L         
            if nd == 1               
                col_norms_Lw = sqrt((Lw_norm_sq)).^DepthExp;
            else % nd > 1
                col_norms_Lw = sqrt(sum(reshape(Lw_norm_sq, nd, []),1)).^DepthExp;               
            end
            low_bound = max(col_norms_Lw)./100; % avoid too small values
            col_norms_Lw(col_norms_Lw < low_bound) = low_bound;
            Rs_spatial_diag_weights = max(col_norms_Lw) ./ col_norms_Lw;
            Rs_spatial = spdiags(Rs_spatial_diag_weights.^2);
            
            if nd > 1
                Rs = kron(Rs_spatial, eye(nd));
            else
                Rs = Rs_spatial;
            end
            % Calculate weighted L_w
            L_w = bsxfun(@times,L_w,Rs_spatial_diag_weights);
        otherwise
            error('seal_MNE:InvalidSourceCovString', 'Unknown string for SourceCovariance.');
    end
    Rc = speye(Nsources_total);
else
    error('seal_MNE:InvalidSourceCovType', 'SourceCovariance must be a matrix or a valid string.');
end
Param_MNE.SourceCovariance = Rs;

%% 3. Regularization Parameter (lambda^2 for the whitened problem)
% The parameter is often referred to as lambda^2 in MNE equations like
% inv(L_w*R_s*L_w' + lambda^2*I). Here, lambda_param_input is lambda.
if ischar(lambda_param_input) || isstring(lambda_param_input)
    switch lower(lambda_param_input)
        case 'auto'
        % learn lambda from data under Gaussian distribution    
            Max_iter = 100;
            evidence = zeros(Max_iter,1);
            gamma = 1;
            cost_old = 0;
            for iter = 1 : Max_iter
                Cov_B = eye(Nchannels) + gamma * (L_w * L_w)';
                gamma = norm(gamma * L_w' / Cov_B * Data_w,'fro')^2 / (Ntimepoints * trace(gamma * L_w' /Cov_B * L_w));
                cost = -(Ntimepoints * log(det(Cov_B)) + trace(Data_w * Data_w' / Cov_B) + Ntimepoints * Nchannels * log(2 * pi)) / 2;
                MSE = (cost - cost_old)/cost;
                cost_old = cost;
                evidence(iter) = cost;
                if abs(MSE) < 1e-5
                    break;
                end
            end
            lambda_sq_used = 1/gamma;
            
        otherwise
            error('seal_MNE:InvalidRegParamString', 'Unknown string for RegularizationParameter.');
    end
elseif isnumeric(lambda_param_input) && isscalar(lambda_param_input) && lambda_param_input >= 0
    lambda_sq_used = lambda_param_input^2 * trace(L_w * L_w') / Nchannels; % Square the input lambda, and scaled
else
    error('seal_MNE:InvalidRegParam', 'RegularizationParameter must be a non-negative scalar or a valid string.');
end
    Param_MNE.RegularizationParameter = sqrt(lambda_sq_used); % Store actually used lambda
%% 4. MNE Inverse Operator for Whitened Data (M_w)
% M_w = R_s * L_w' * inv(L_w * R_s * L_w' + lambda_sq_used * I_nchannels)
Gram = L_w * L_w' + lambda_sq_used * eye(Nchannels);
try
   if rank(Gram) < Nchannels
        warning('seal_MNE:RankDeficientSystemInv', 'System matrix for inversion (L_w*R_s*L_w'' + eff_lambda_sq*I) is rank deficient. Using pseudo-inverse.');
        M_w = Rc * L_w' * pinv(Gram);
   else
        M_w = Rc * L_w' / Gram; % Solves M_w * Term_inv = Term_mult for M_w
   end
catch ME_solve
    error('seal_MNE:SolverError', 'Error solving for inverse operator M_w: %s', ME_solve.message);
end

%% 5. Source Estimate
Source_MNE = M_w * Data_w;

%% 6. Final Inverse Operator (to be applied to original, unwhitened data)
Param_MNE.InverseOperator = M_w * iW;    

end