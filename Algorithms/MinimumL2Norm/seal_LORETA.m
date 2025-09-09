function [Source_LORETA, Param_LORETA] = seal_LORETA(Data, L, sourceSpaceInfo, varargin)
%SEAL_LORETA Computes the LORETA estimate.
%   [Source_LORETA, Param_LORETA] = seal_LORETA(Data, L, sourceSpaceInfo, 'ParameterName', ParameterValue, ...)
%
%   LORETA incorporates a spatial smoothness constraint (Laplacian operator)
%   into the source covariance matrix Rs.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints data.
%       L (double matrix): Nchannels x Nsources original lead field.
%       sourceSpaceInfo (struct): Struct with source space geometry. Required:
%           .VertConn (sparse matrix): Nloc x Nloc vertex connectivity/adjacency.
%           (Optional: .VertLoc for distance-based Laplacian variants, not used in this version)
%
%   Optional Name-Value Pair Inputs:
%       'RegularizationParameter' (double): For the main MNE solve. Default: 0.1.
%       'NoiseCovariance' (double matrix): Default: eye(Nchannels).
%       'NumOrientations' (integer): Default: 1.
%       'DepthWeightingExponent' (double): Exponent for depth weights combined with Laplacian.
%                                          Default: 1 (i.e., weights are (1/norm_L)^1).
%                                          Set to 0 for no additional depth weighting beyond Laplacian.
%       'LaplacianRegFactor' (double): Small regularization for (M'*M + factor*I)
%                                      Default: 1e-4.
%       'SourceCovariance' (double matrix or char/string):
%                                           - Nsources x Nsources prior source covariance matrix (R_s).
%                                           - 'none' or 'depth_weighted'. If 'depth_weighted', computes diagonal
%                                             depth weighting, else use eye(Nsources).
%                                           Default:'depth_weighted'.
%       % Other parameters accepted by seal_MNE can be passed.
%
%   Outputs:
%       Source_LORETA (double matrix): Nsources x Ntimepoints LORETA estimates.
%       Param_LORETA (struct): Parameters and results:
%           .InverseOperator (double matrix): Operator mapping original Data to LORETA Source.
%           .BaseMNE_Parameters (struct): The 'Param' struct from seal_MNE.
%           .Rs_LORETA (sparse matrix): The constructed LORETA source covariance.
%           .OptionsPassed (struct): Parsed input options for LORETA.
%
%   See also: SEAL_MNE, SEAL_LAURA, RUNESIANALYSIS.


%% Input Parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;

defaultRegParam = 0.1;
defaultNoiseCov = eye(size(L, 1));
defaultNumOrientations = 1;
defaultDepthWeightExp = 1; % (1/norm)^1
defaultLaplacianRegFactor = 1e-2;
defaultSourceCov = 'depth_weighted'; % Default to automatic depth weighting

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);

addParameter(p, 'RegularizationParameter', defaultRegParam);
addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'DepthWeightingExponent', defaultDepthWeightExp, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'LaplacianRegFactor', defaultLaplacianRegFactor, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'SourceCovariance', defaultSourceCov);
try
    parse(p, Data, L, sourceSpaceInfo, varargin{:});
catch ME
    disp('Error parsing inputs for seal_LORETA:');
    rethrow(ME);
end

Param_LORETA.OptionsPassed = p.Results;
srcInfo = p.Results.sourceSpaceInfo;
nd = p.Results.NumOrientations;
DepthExp = p.Results.DepthWeightingExponent;
lap_reg_factor = p.Results.LaplacianRegFactor;
Rs_input = p.Results.SourceCovariance;
if ~isfield(srcInfo, 'VertConn')
    error('seal_LORETA:MissingVertConn', 'sourceSpaceInfo must contain .VertConn for LORETA.');
end

Nsources_total = size(L,2);
Nsources_locations = Nsources_total / nd;
if size(srcInfo.VertConn,1) ~= Nsources_locations
    error('seal_LORETA:VertConnDimMismatch', 'VertConn dimensions do not match Nsources_locations.');
end

%% 1. Construct LORETA Source Covariance (Rs_LORETA)

VertConn = srcInfo.VertConn;
% graph Laplacian D - A
B_laplacian = spdiags(sum(VertConn, 2), 0, Nsources_locations, Nsources_locations) - VertConn;
BB_lap = B_laplacian' * B_laplacian;
% Regularize BB_lap to ensure invertibility if needed
BB_lap_reg = BB_lap + lap_reg_factor * trace(BB_lap)/Nsources_locations * speye(Nsources_locations);

%% 2. Determine Source Covariance
% W_depth elements are w_j_depth (e.g., 1/norm of L_orig columns)
if isnumeric(Rs_input) && ismatrix(Rs_input)
    % User provided Rs directly
    if (nd == 1 && isequal(size(Rs_input), [Nsources_total, Nsources_total])) || ...
            (nd > 1 && (isequal(size(Rs_input), [Nsources_locations, Nsources_locations]) || isequal(size(Rs_input), [Nsources_total, Nsources_total])))
        Rs_spatial_or_full = Rs_input;
    else
        error('seal_LORETA:SourceCovDimMismatch', 'User-provided SourceCovariance has incorrect dimensions.');
    end
    
    if nd > 1 && isequal(size(Rs_spatial_or_full), [Nsources_locations, Nsources_locations])
        Rs = kron(Rs_spatial_or_full, eye(nd));
    else
        Rs = Rs_spatial_or_full; % Assumed to be full Nsources_total x Nsources_total
    end
    
elseif ischar(Rs_input) || isstring(Rs_input)
    switch lower(Rs_input)
        case 'none'
            Rs = BB_lap_reg\speye(size(BB_lap_reg));
        case 'depth_weighted'
            Lw_norm_sq = sum(L.^2, 1);   % based on original L
            if nd == 1
                col_norms_Lw = sqrt((Lw_norm_sq)).^DepthExp;
            else % nd > 1
                col_norms_Lw = sqrt(sum(reshape(Lw_norm_sq, nd, []),1)).^DepthExp;
            end
            low_bound = max(col_norms_Lw)./100; % avoid too small values
            col_norms_Lw(col_norms_Lw < low_bound) = low_bound;
            Rs_spatial_diag_weights = max(col_norms_Lw) ./ col_norms_Lw;
            Rs_weight = speye(size(BB_lap_reg));
            Rs_weight = spdiags(Rs_spatial_diag_weights', 0 , Rs_weight);
            % Combine: Rs = Rs_depth_spatial * pinv(BB_lap_reg) * Rs_depth_spatial
            Rs = Rs_weight / BB_lap_reg * Rs_weight;        
        otherwise
            error('seal_LORETA:InvalidSourceCovString', 'Unknown string for SourceCovariance.');
    end
    if nd > 1
        Rs = kron(Rs, eye(nd));
    end
else
    error('seal_LORETA:InvalidSourceCovType', 'SourceCovariance must be a matrix or a valid string.');
end

Param_LORETA.Rs_LORETA = Rs;

%% 3. Call seal_MNE with the LORETA source covariance
mneOpts = {'SourceCovariance', Rs, ...
    'RegularizationParameter', p.Results.RegularizationParameter, ...
    'NoiseCovariance', p.Results.NoiseCovariance, ...
    'NumOrientations', nd};

unmatchedArgs = namedargs2cell(p.Unmatched);
mneOpts = [mneOpts, unmatchedArgs];

[Source_LORETA, Param_mne] = seal_MNE(Data, L, mneOpts{:});
Param_LORETA.BaseMNE_Parameters = Param_mne;
Param_LORETA.InverseOperator = Param_mne.InverseOperator;
end
