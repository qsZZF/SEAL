function [Source_LAURA, Param_LAURA] = seal_LAURA(Data, L, sourceSpaceInfo, varargin)
%SEAL_LAURA Computes LAURA using a transformation to standard MNE.
%   [Source_LAURA, Param_LAURA] = seal_LAURA(Data, L, sourceSpaceInfo, 'param', value, ...)
%
%   Based on Grave de Peralta Menendez et al. (2001).
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints data.
%       L (double matrix): Nchannels x Nsources original lead field.
%       sourceSpaceInfo (struct): Must contain .VertConn, .VertLoc (Nsources * 3).
%
%   Optional Name-Value Pair Inputs:
%       'RegularizationParameter' (double): Lambda for LAURA. Default: 0.1.
%       'NoiseCovariance' (double matrix): Default: eye(Nchannels).
%       'NumOrientations' (integer): Default: 1.
%       'DepthWeightingExponent' (double): Exponent for depth weights combined with Laplacian.
%                                          Default: 2 (i.e., weights are (1/norm_L)^2).
%                                          Set to 0 for no additional depth weighting beyond Laplacian.
%       'Exponent_e_LAURA' (double): Exponent e_i for LAURA operator A. Default: 2.
%       'MaxNeighborsN_LAURA' (integer): Max neighbors N for LAURA operator A. Default: 9 for 2D cortex.
%       % Other parameters for seal_MNE can be passed.
%
%   Outputs:
%       Source_LAURA (double matrix): Nsources x Ntimepoints LAURA estimates.
%       Param_LAURA (struct):
%           .InverseOperator (double matrix): Maps original Data to LAURA Source.
%           .BaseMNE_Parameters (struct): The 'Param' struct from seal_MNE.
%           .Rs_LAURA (sparse matrix): The constructed LAURA source covariance.
%           .OptionsPassed (struct)
%
%   Reference: [1] Grave de Peralta Menendez et al. (2001), Brain Topography, 14(2).

%   Copyright 2025 Your Name/Institution

%% Input Parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;

defaultRegParam = 0.1; % Lambda for LAURA objective
defaultNoiseCov = eye(size(L,1));
defaultNumOrientations = 1;
defaultDepthWeightExp = 2; % (1/norm)^2 for depth weighting
defaultExponent_e_LAURA = 2; % From refer [1]
defaultMaxNeighborsN_LAURA = 9; % In refer [1] 26 is used for 3D volume. Here set to 9 for 2D cortex by default
defaultSourceCov = 'depth_weighted'; % Default to automatic depth weighting

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);

addParameter(p, 'RegularizationParameter', defaultRegParam);
addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
addParameter(p, 'NumOrientations', defaultNumOrientations, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'Exponent_e_LAURA', defaultExponent_e_LAURA, @isscalar);
addParameter(p, 'MaxNeighborsN_LAURA', defaultMaxNeighborsN_LAURA, @isscalar);
addParameter(p, 'SourceCovariance', defaultSourceCov);
addParameter(p, 'DepthWeightingExponent', defaultDepthWeightExp, @(x) isnumeric(x) && isscalar(x));

try
    parse(p, Data, L, sourceSpaceInfo, varargin{:});
catch ME
    disp('Error parsing inputs for seal_LAURA:');
    rethrow(ME);
end

Param_LAURA.OptionsPassed = p.Results;
Nsources_total = size(L, 2);
nd = p.Results.NumOrientations;
Nsources_locations = Nsources_total / nd;
Rs_input = p.Results.SourceCovariance;
DepthExp = p.Results.DepthWeightingExponent;

if size(srcInfo.VertConn,1) ~= Nsources_locations || size(srcInfo.VertLoc,1) ~= Nsources_locations
    error('seal_LAURA_Rs:SourceInfoDimMismatch', 'VertConn/VertLoc dimensions do not match Nsources_locations.');
end
%% 1. Compute LAURA Operator (M_LAURA_op_full = WA)
VertConn = p.Results.sourceSpaceInfo.VertConn;
VertLoc = p.Results.sourceSpaceInfo.VertLoc;
exponent_e = p.Results.Exponent_e_LAURA;
N_maxNeighbor = p.Results.MaxNeighborsN_LAURA;

A_spatial = spalloc(Nsources_locations, Nsources_locations, Nsources_locations * (N_maxNeighbor + 1)); % Preallocate sparse
% Calculate distance-based smooth operator A
for i_loc = 1:Nsources_locations
    neighbor_indices = find(VertConn(i_loc, :));
    N_i_actual_neighbors = length(neighbor_indices);
    
    if N_i_actual_neighbors == 0
        A_spatial(i_loc, i_loc) = 1; % Or some other sensible default for isolated points
        continue;
    end   
    dist_vec = sqrt(sum((VertLoc(neighbor_indices,:) - VertLoc(i_loc,:)).^2,2));
    dist_neg_pow = dist_vec.^(-exponent_e);
    A_spatial(i_loc, neighbor_indices) = -dist_neg_pow;
    A_spatial(i_loc, i_loc) = (N_maxNeighbor / N_i_actual_neighbors)*sum(dist_neg_pow);
end

%% 2. Determine Source Covariance
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
            Rs = pinv(full(A_spatial * A_spatial'));
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
            % Combine: Rs = Rs_depth_spatial * pinv(A_spatial*A_spatial') * Rs_depth_spatial
            Rs = spdiags(Rs_spatial_diag_weights) / full(A_spatial * A_spatial') * spdiags(Rs_spatial_diag_weights);
        otherwise
            error('seal_LORETA:InvalidSourceCovString', 'Unknown string for SourceCovariance.');
    end
    if nd > 1
        Rs = kron(Rs, eye(nd));
    end
else
    error('seal_LORETA:InvalidSourceCovType', 'SourceCovariance must be a matrix or a valid string.');
end

Param_LAURA.Rs_LAURA = Rs;

%% 3. Call seal_MNE to solve
mneOpts = {'SourceCovariance', Rs, ...
    'RegularizationParameter', p.Results.RegularizationParameter, ...
    'NoiseCovariance', p.Results.NoiseCovariance, ...
    'NumOrientations', nd};

unmatchedArgs = namedargs2cell(p.Unmatched);
mneOpts = [mneOpts, unmatchedArgs];

[Source_LAURA, Param_mne] = seal_MNE(Data, L, mneOpts{:});
Param_LAURA.Param_MNE = Param_mne;
Param_LAURA.InverseOperator = Param_mne.InverseOperator;

end
