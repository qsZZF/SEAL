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
%       'DepthWeightingExponent' (double): Exponent for depth weights combined with LAURA constraint.
%                                          Default: 1 (i.e., weights are (1/norm_L)^1).
%                                          Set to 0 for no additional depth weighting beyond Laplacian.
%       'Exponent_e_LAURA' (double): Exponent e_i for LAURA operator A. Default: 2.
%       'MaxNeighborsN_LAURA' (integer): Max neighbors N for LAURA operator A. Default: 9 for 2D cortex.
%       'LauReg' (double): Tikhonov regularisation factor applied to LAU = A'*A,
%                          which is otherwise singular (A has row sums of order
%                          zero when N_maxNeighbor matches actual neighbour counts).
%                          LAU <- LAU + LauReg * trace(LAU)/N * I. Default: 1e-4.
%       % Other parameters for seal_MNE can be passed.
%
%   Source covariance construction:
%       A_ij = -(d_ij / d_min)^{-e} for i a neighbour of j (d_min = min
%           neighbour distance; dividing by d_min makes the operator
%           invariant to the units of VertLoc).
%       A_ii = (N_max / n_neighbours_i) * sum_j |A_ij|.
%       LAU  = A'*A + LauReg * trace(A'*A)/N * I.
%       Rs   = inv(LAU)                              (SourceCovariance='none'), or
%       Rs   = diag(w) * inv(LAU) * diag(w)          ('depth_weighted'),
%       with w_i = max(||L_i||^p) / ||L_i||^p, p = DepthWeightingExponent.
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

%% Input Parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;

defaultRegParam = 0.1; % Lambda for LAURA objective
defaultNoiseCov = eye(size(L,1));
defaultNumOrientations = 1;
defaultDepthWeightExp = 1; % (1/norm)^1 for depth weighting
defaultExponent_e_LAURA = 2; % From refer [1]
defaultMaxNeighborsN_LAURA = 9; % In refer [1] 26 is used for 3D volume. Here set to 9 for 2D cortex by default
defaultSourceCov = 'depth_weighted'; % Default to automatic depth weighting
defaultLauReg = 1e-4;  % Tikhonov regularization for LAU = A'*A, which is
                       % singular/near-singular in theory (A has near-zero
                       % row sums). Without this, inv(LAU) is ill-posed.

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
addParameter(p, 'LauReg', defaultLauReg, @(x) isnumeric(x) && isscalar(x) && x >= 0);

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
srcInfo = p.Results.sourceSpaceInfo;

if size(srcInfo.VertConn,1) ~= Nsources_locations || size(srcInfo.Vertices,1) ~= Nsources_locations
    error('seal_LAURA_Rs:SourceInfoDimMismatch', 'VertConn/VertLoc dimensions do not match Nsources_locations.');
end
%% 1. Compute LAURA Operator (M_LAURA_op_full = WA)
VertConn = p.Results.sourceSpaceInfo.VertConn;
VertLoc = p.Results.sourceSpaceInfo.Vertices;
exponent_e = p.Results.Exponent_e_LAURA;
N_maxNeighbor = p.Results.MaxNeighborsN_LAURA;
LauReg = p.Results.LauReg;

% Build dist_mat(i, j) = ||VertLoc(i,:) - VertLoc(j,:)|| if i is a neighbour
% of j (and i != j); inf otherwise. We use find() on each column so that
% VertConn can be either a logical or a (sparse) numeric 0/1 adjacency
% matrix -- VertConn(:,iter) used directly as a subscript is unreliable
% across these two types.
dist_mat = inf(Nsources_locations, Nsources_locations);
for iter = 1:Nsources_locations
    nb = find(VertConn(:, iter));
    nb(nb == iter) = [];   % exclude self-loops so D_min > 0 is guaranteed
    if isempty(nb)
        continue;
    end
    t_dist = sqrt(sum((VertLoc(nb, :) - VertLoc(iter, :)).^2, 2));
    dist_mat(nb, iter) = t_dist;
end

% Non-dimensionalise the distances by the minimum neighbour distance so that
% the constructed operator is invariant to the units of VertLoc (mm vs m).
% Dividing by D_min (not sqrt(D_min)) is what makes d_ij/D_min dimensionless;
% the sqrt in the original code left d^{1/2} units in and caused the whole
% operator to rescale by a factor of scale^e under a units change.
D_min = min(dist_mat(isfinite(dist_mat)));
if isempty(D_min) || D_min <= 0
    error('seal_LAURA:InvalidMinDistance', ...
        'Could not determine a positive minimum neighbour distance. Check VertConn/Vertices.');
end
dist_mat = dist_mat ./ D_min;
dist_mat = -dist_mat.^(-exponent_e);   % A_ij = -(d_ij/D_min)^{-e} on neighbours
dist_mat(isinf(dist_mat)) = 0;         % non-neighbours (were inf) -> 0
% Diagonal weighting from Grave de Peralta et al. (2001): scale so that
% diag(A) = (N_max / n_neighbours_i) * sum_j |A_ij|. For locations whose
% realised neighbour count equals N_maxNeighbor this makes A have row sum 0
% (a proper difference operator); otherwise the scaling differs from 1.
neighbor_count = sum(dist_mat < 0, 2);
neighbor_count(neighbor_count == 0) = 1;   % avoid 0/0 for isolated vertices
dist_mat = dist_mat + diag((N_maxNeighbor ./ neighbor_count) .* sum(-dist_mat, 2));
dist_mat = sparse(dist_mat);

% A'*A is singular in exact arithmetic whenever A has a zero row-sum (the
% constant vector lies in its null space), and is badly conditioned more
% generally. Apply a small Tikhonov regulariser so that inv(LAU) is
% well-defined numerically. Matches the equivalent LaplacianRegFactor used
% by seal_LORETA.
LAU = dist_mat * dist_mat';
LAU = LAU + LauReg * trace(LAU)/Nsources_locations * speye(Nsources_locations);

%% 2. Determine Source Covariance
if isnumeric(Rs_input) && ismatrix(Rs_input)
    % User provided Rs directly
    if (nd == 1 && isequal(size(Rs_input), [Nsources_total, Nsources_total])) || ...
            (nd > 1 && (isequal(size(Rs_input), [Nsources_locations, Nsources_locations]) || isequal(size(Rs_input), [Nsources_total, Nsources_total])))
        Rs_spatial_or_full = Rs_input;
    else
        error('seal_LAURA:SourceCovDimMismatch', 'User-provided SourceCovariance has incorrect dimensions.');
    end
    
    if nd > 1 && isequal(size(Rs_spatial_or_full), [Nsources_locations, Nsources_locations])
        Rs = kron(Rs_spatial_or_full, eye(nd));
    else
        Rs = Rs_spatial_or_full; % Assumed to be full Nsources_total x Nsources_total
    end
    
elseif ischar(Rs_input) || isstring(Rs_input)
    switch lower(Rs_input)
        case 'none'
            Rs = LAU\speye(size(LAU));
%             Rs = LAU;
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
            % Combine: Rs = Rs_depth_spatial * pinv(BB_lap_reg) * Rs_depth_spatial
            Rs = Rs_weight / LAU * Rs_weight;
        otherwise
            error('seal_LAURA:InvalidSourceCovString', 'Unknown string for SourceCovariance.');
    end
    if nd > 1
        Rs = kron(Rs, eye(nd));
    end
else
    error('seal_LAURA:InvalidSourceCovType', 'SourceCovariance must be a matrix or a valid string.');
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