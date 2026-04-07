function [Source_LOGURA, Param_LOGURA] = seal_LOGURA(Data, L, sourceSpaceInfo, varargin)
%SEAL_LOGURA Computes a smooth-kernel-based  transformation to standard MNE.
%   [Source_LOGURA, Param_LOGURA] = seal_LOGURA(Data, L, sourceSpaceInfo, varargin)
%
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints data.
%       L (double matrix): Nchannels x Nsources original lead field.
%       sourceSpaceInfo (struct): Must contain .VertConn, .VertLoc (Nsources * 3).
%
%   Optional Name-Value Pair Inputs:
%       'RegularizationParameter' (double): Lambda. Default: 0.1.
%       'NoiseCovariance' (double matrix): Default: eye(Nchannels).
%       'NumOrientations' (integer): Default: 1.
%       'DepthWeightingExponent' (double): Exponent for depth weights combined with Kernel.
%                                          Default: 1 (i.e., weights are (1/norm_L)^1).
%                                          Set to 0 for no additional depth weighting beyond Laplacian.
%       'Exponent_e_LOGURA' (double): Exponent e_i for LOGURA operator A. Default: 2.
%       % Other parameters for seal_MNE can be passed.
%
%   Outputs:
%       Source_LOGURA (double matrix): Nsources x Ntimepoints LOGURA estimates.
%       Param_LOGURA (struct):
%           .InverseOperator (double matrix): Maps original Data to LOGURA Source.
%           .BaseMNE_Parameters (struct): The 'Param' struct from seal_MNE.
%           .Rs_LOGURA (sparse matrix): The constructed LOGURA source covariance.
%           .OptionsPassed (struct)
%


%% Input Parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;

defaultRegParam = 0.1; % Lambda for LOGURA objective
defaultNoiseCov = eye(size(L,1));
defaultNumOrientations = 1;
defaultDepthWeightExp = 1; % (1/norm)^1 for depth weighting
defaultExponent_e_LOGURA = 4; 
defaultSourceCov = 'depth_weighted'; % Default to automatic depth weighting

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);

addParameter(p, 'RegularizationParameter', defaultRegParam);
addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'Exponent_e_LOGURA', defaultExponent_e_LOGURA, @isscalar);
addParameter(p, 'SourceCovariance', defaultSourceCov);
addParameter(p, 'DepthWeightingExponent', defaultDepthWeightExp, @(x) isnumeric(x) && isscalar(x));

try
    parse(p, Data, L, sourceSpaceInfo, varargin{:});
catch ME
    disp('Error parsing inputs for seal_LAURA:');
    rethrow(ME);
end

Param_LOGURA.OptionsPassed = p.Results;
Nsources_total = size(L, 2);
nd = p.Results.NumOrientations;
Nsources_locations = Nsources_total / nd;
Rs_input = p.Results.SourceCovariance;
DepthExp = p.Results.DepthWeightingExponent;
srcInfo = p.Results.sourceSpaceInfo;

if size(srcInfo.VertConn,1) ~= Nsources_locations || size(srcInfo.Vertices,1) ~= Nsources_locations
    error('seal_LAURA_Rs:SourceInfoDimMismatch', 'VertConn/VertLoc dimensions do not match Nsources_locations.');
end
%% 1. Compute LOGURA Operator 
VertConn = p.Results.sourceSpaceInfo.VertConn;
VertLoc = p.Results.sourceSpaceInfo.Vertices;
exponent_e = p.Results.Exponent_e_LOGURA;


dist_mat = inf(Nsources_locations,Nsources_locations);
for iter = 1:Nsources_locations
    t_dist = sqrt(sum((repmat(VertLoc(iter,:),Nsources_locations,1) - VertLoc).^2,2));    
    dist_mat(VertConn(:,iter),iter) = t_dist(VertConn(:,iter));
end
D_min = min(min(dist_mat));
dist_mat = exp(-(dist_mat./(2*sqrt(D_min))).^exponent_e);
dist_mat(isinf(dist_mat)) = 0;
dist_mat = sparse(dist_mat+eye(Nsources_locations));

LAU = dist_mat*dist_mat';

% LAU = LAU + 1e-4 * trace(LAU)/Nsources_locations * speye(Nsources_locations);
%% 2. Determine Source Covariance
if isnumeric(Rs_input) && ismatrix(Rs_input)
    % User provided Rs directly
    if (nd == 1 && isequal(size(Rs_input), [Nsources_total, Nsources_total])) || ...
            (nd > 1 && (isequal(size(Rs_input), [Nsources_locations, Nsources_locations]) || isequal(size(Rs_input), [Nsources_total, Nsources_total])))
        Rs_spatial_or_full = Rs_input;
    else
        error('seal_LOGURA:SourceCovDimMismatch', 'User-provided SourceCovariance has incorrect dimensions.');
    end
    
    if nd > 1 && isequal(size(Rs_spatial_or_full), [Nsources_locations, Nsources_locations])
        Rs = kron(Rs_spatial_or_full, eye(nd));
    else
        Rs = Rs_spatial_or_full; % Assumed to be full Nsources_total x Nsources_total
    end
    
elseif ischar(Rs_input) || isstring(Rs_input)
    switch lower(Rs_input)
        case 'none'
            Rs = LAU;
        case 'depth_weighted'
            Lw_norm_sq = sum(L.^2,1);   % based on original L
            if nd == 1
                col_norms_Lw = sqrt((Lw_norm_sq)).^DepthExp;
            else % nd > 1
                col_norms_Lw = sqrt(sum(reshape(Lw_norm_sq, nd, []),1)).^DepthExp;
            end
            low_bound = max(col_norms_Lw)./100; % avoid too small values
            col_norms_Lw(col_norms_Lw < low_bound) = low_bound;
            Rs_spatial_diag_weights = max(col_norms_Lw) ./ col_norms_Lw;
            Rs_weight = speye(size(LAU));
            Rs_weight = spdiags(Rs_spatial_diag_weights', 0 , Rs_weight);
            Rs = Rs_weight * LAU * Rs_weight;
        otherwise
            error('seal_LOGURA:InvalidSourceCovString', 'Unknown string for SourceCovariance.');
    end
    if nd > 1
        Rs = kron(Rs, eye(nd));
    end
else
    error('seal_LOGURA:InvalidSourceCovType', 'SourceCovariance must be a matrix or a valid string.');
end

Param_LOGURA.Rs_LAURA = Rs;

%% 3. Call seal_MNE to solve
mneOpts = {'SourceCovariance', Rs, ...
    'RegularizationParameter', p.Results.RegularizationParameter, ...
    'NoiseCovariance', p.Results.NoiseCovariance, ...
    'NumOrientations', nd};

unmatchedArgs = namedargs2cell(p.Unmatched);
mneOpts = [mneOpts, unmatchedArgs];

[Source_LOGURA, Param_mne] = seal_MNE(Data, L, mneOpts{:});
Param_LOGURA.Param_MNE = Param_mne;
Param_LOGURA.InverseOperator = Param_mne.InverseOperator;

end
