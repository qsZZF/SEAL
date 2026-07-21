function [Source_Estimate, Parameters] = seal_SmoothChampagne(Data, L, sourceSpaceInfo, varargin)
%SEAL_SMOOTHCHAMPAGNE Smooth Champagne for distributed E/MEG sources.
%   [S, PARAMETERS] = SEAL_SMOOTHCHAMPAGNE(DATA, L, SOURCESPACEINFO, ...)
%   implements the Smooth Champagne model of Cai et al. (2020).  The
%   physical source is written as S = H*Z, where H is a local cortical
%   smoothing kernel and Z has one learned variance per spatial tile.
%
%   DATA is Nchannels-by-Ntime and L is Nchannels-by-(Nlocations*nd).
%   SOURCESPACEINFO must contain VertConn and may contain Vertices.
%
%   Name-value options
%     NoiseCovariance          Fixed sensor noise covariance (default []).
%     UpdateNoiseCovariance    Learn diagonal residual covariance (default false).
%     KernelWidth              Paper parameter m (default 2).
%     KernelPower              Paper parameter p (default 4).
%     KernelThreshold          Discard H entries smaller than this (default 0.1).
%     TileSize                 Target number of locations per graph tile (default 6).
%     TileAssignments          Optional Nlocations-by-1 positive integer labels.
%     MaxIterations            Maximum fixed-point iterations (default 100).
%     Tolerance                Relative objective tolerance (default 1e-6).
%     NumOrientations          Orientations per location, 1 or 3 (default 1).
%     CovarianceRegularization Relative diagonal loading (default 1e-10).
%     Verbose                  Print iteration diagnostics (default true).
%
%   Reference
%     C. Cai et al., "Robust Empirical Bayesian Reconstruction of
%     Distributed Sources for Electromagnetic Brain Imaging," IEEE TMI,
%     39(3), 567-577, 2020. doi:10.1109/TMI.2019.2932290

p = inputParser;
p.CaseSensitive = false;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);
addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'UpdateNoiseCovariance', false, ...
    @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
addParameter(p, 'KernelWidth', 2, @(x) isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'KernelPower', 4, @(x) isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'KernelThreshold', 0.1, ...
    @(x) isscalar(x) && isfinite(x) && x > 0 && x < 1);
addParameter(p, 'TileSize', 6, @(x) isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'TileAssignments', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'MaxIterations', 100, @(x) isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'Tolerance', 1e-6, @(x) isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'NumOrientations', 1, @(x) isscalar(x) && ismember(x, [1, 3]));
addParameter(p, 'CovarianceRegularization', 1e-10, ...
    @(x) isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'Verbose', true, ...
    @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parse(p, Data, L, sourceSpaceInfo, varargin{:});

opts = p.Results;
Data = double(Data);
L = double(L);
verbose = logical(opts.Verbose);
updateNoise = logical(opts.UpdateNoiseCovariance);
maxIter = max(1, round(opts.MaxIterations));
nd = round(opts.NumOrientations);

if any(~isfinite(Data(:))) || any(~isfinite(L(:)))
    error('seal_SmoothChampagne:NonFiniteInput', ...
        'Data and L must contain only finite values.');
end

[nChannels, nTime] = size(Data);
if size(L, 1) ~= nChannels
    error('seal_SmoothChampagne:DimensionMismatch', ...
        'Data and L must have the same number of sensor rows.');
end

nSourceComponents = size(L, 2);
if mod(nSourceComponents, nd) ~= 0
    error('seal_SmoothChampagne:OrientationMismatch', ...
        'The number of lead-field columns must be divisible by NumOrientations.');
end
nLocations = nSourceComponents / nd;

if ~isfield(sourceSpaceInfo, 'VertConn') || isempty(sourceSpaceInfo.VertConn)
    error('seal_SmoothChampagne:MissingVertConn', ...
        'sourceSpaceInfo.VertConn is required to build the local kernel and tiles.');
end
adjacency = prepare_connectivity(sourceSpaceInfo.VertConn, nLocations);
vertices = [];
if isfield(sourceSpaceInfo, 'Vertices') && ~isempty(sourceSpaceInfo.Vertices)
    vertices = double(sourceSpaceInfo.Vertices);
    if ~isequal(size(vertices), [nLocations, 3]) || any(~isfinite(vertices(:)))
        error('seal_SmoothChampagne:InvalidVertices', ...
            'Vertices must be a finite Nlocations-by-3 matrix.');
    end
end

[Hlocation, kernelInfo] = build_local_kernel(adjacency, vertices, ...
    opts.KernelWidth, opts.KernelPower, opts.KernelThreshold);
H = kron(Hlocation, speye(nd));
G = L * H;

[tileAssignments, tileIndices] = build_tiles(adjacency, opts.TileSize, ...
    opts.TileAssignments, nLocations);
nTiles = numel(tileIndices);
componentTiles = repelem(tileAssignments, nd);

Cy = (Data * Data.') / nTime;
Cy = (Cy + Cy.') / 2;
Cnoise = initialize_noise_covariance(opts.NoiseCovariance, Cy, nChannels, ...
    opts.CovarianceRegularization);

initialVariance = trace(Cy) / max(trace(G * G.'), eps);
if ~isfinite(initialVariance) || initialVariance <= 0
    initialVariance = 1;
end
regionVariances = initialVariance * ones(nTiles, 1);

costHistory = nan(maxIter, 1);
relativeChange = inf;
converged = false;
noiseFloor = max(trace(Cnoise) / max(nChannels, 1) * 1e-8, eps);

if verbose
    fprintf('Running Smooth Champagne (nd=%d, tiles=%d, T=%d)...\n', ...
        nd, nTiles, nTime);
end

for iter = 1:maxIter
    componentVariances = regionVariances(componentTiles);
    weightedG = bsxfun(@times, G, componentVariances.');
    SigmaY = weightedG * G.' + Cnoise;
    SigmaY = stabilize_spd(SigmaY, opts.CovarianceRegularization);

    SigmaYInvG = solve_spd(SigmaY, G);
    latentFilter = bsxfun(@times, componentVariances, SigmaYInvG.');
    Z = latentFilter * Data;

    % Equation (12): regional resolution weights from the current evidence
    % covariance.  Grouping is done at the location level, then expanded to
    % the requested number of dipole orientations.
    componentResolution = sum(G .* SigmaYInvG, 1).';
    resolutionWeights = accumarray(componentTiles, componentResolution, ...
        [nTiles, 1], @sum, 0) / nd;
    resolutionWeights = max(real(resolutionWeights), eps);

    % Equation (11): fixed-point update for one variance per spatial tile.
    componentEnergy = sum(Z.^2, 2);
    tileEnergy = accumarray(componentTiles, componentEnergy, [nTiles, 1], @sum, 0);
    regionVariances = sqrt(max(tileEnergy ./ ...
        (nTime * nd * resolutionWeights), eps));
    regionVariances = min(max(regionVariances, 1e-16), 1e16);

    if updateNoise
        residual = Data - G * Z;
        newNoise = diag(max(mean(residual.^2, 2), noiseFloor));
        Cnoise = stabilize_spd(newNoise, opts.CovarianceRegularization);
    end

    % The data covariance objective is the negative log evidence up to an
    % additive constant.  It also gives a scale-independent stopping rule.
    componentVariances = regionVariances(componentTiles);
    weightedG = bsxfun(@times, G, componentVariances.');
    SigmaY = weightedG * G.' + Cnoise;
    SigmaY = stabilize_spd(SigmaY, opts.CovarianceRegularization);
    SigmaYInv = solve_spd(SigmaY, eye(nChannels));
    cost = real(trace(Cy * SigmaYInv) + spd_logdet(SigmaY));
    costHistory(iter) = cost;

    if iter > 1
        relativeChange = abs(costHistory(iter) - costHistory(iter - 1)) / ...
            max(1, abs(costHistory(iter - 1)));
    end

    if verbose && (iter == 1 || mod(iter, 10) == 0 || relativeChange < opts.Tolerance)
        activeTiles = sum(regionVariances > max(regionVariances) * 1e-3);
        fprintf('  Iter %d: cost %.6g, relative change %.3g, active tiles %d\n', ...
            iter, cost, relativeChange, activeTiles);
    end

    if ~isfinite(cost)
        warning('seal_SmoothChampagne:NonFiniteCost', ...
            'The objective became non-finite at iteration %d.', iter);
        break;
    end
    if iter > 1 && relativeChange < opts.Tolerance
        converged = true;
        break;
    end
end

costHistory = costHistory(1:iter);
componentVariances = regionVariances(componentTiles);
weightedG = bsxfun(@times, G, componentVariances.');
SigmaY = stabilize_spd(weightedG * G.' + Cnoise, opts.CovarianceRegularization);
SigmaYInvG = solve_spd(SigmaY, G);
latentFilter = bsxfun(@times, componentVariances, SigmaYInvG.');
Latent_Estimate = latentFilter * Data;
InverseOperator = H * latentFilter;
Source_Estimate = H * Latent_Estimate;

Parameters = struct();
Parameters.Algorithm = 'Smooth-Cham';
Parameters.OptionsPassed = opts;
Parameters.SpatialKernel = Hlocation;
Parameters.KernelInfo = kernelInfo;
Parameters.ModifiedLeadfield = G;
Parameters.TileAssignments = tileAssignments;
Parameters.TileIndices = tileIndices;
Parameters.LearnedRegionVariances = regionVariances;
Parameters.RegionResolutionWeights = resolutionWeights;
Parameters.LearnedNoiseCovariance = Cnoise;
Parameters.LatentSourceEstimate = Latent_Estimate;
Parameters.InverseOperator = InverseOperator;
Parameters.CostHistory = costHistory;
Parameters.NumIterationsPerformed = iter;
Parameters.Converged = converged;
Parameters.RelativeChange = relativeChange;
Parameters.NumOrientations = nd;
Parameters.Reference = 'doi:10.1109/TMI.2019.2932290';

end

function adjacency = prepare_connectivity(vertConn, nLocations)
if ~isnumeric(vertConn) && ~islogical(vertConn)
    error('seal_SmoothChampagne:InvalidVertConn', ...
        'VertConn must be numeric or logical.');
end
if ~isequal(size(vertConn), [nLocations, nLocations])
    error('seal_SmoothChampagne:VertConnDimensionMismatch', ...
        'VertConn must be Nlocations-by-Nlocations.');
end
if any(~isfinite(double(vertConn(:))))
    error('seal_SmoothChampagne:InvalidVertConn', ...
        'VertConn must contain finite values.');
end

adjacency = spones(sparse(vertConn));
adjacency = spones(adjacency | adjacency.');
adjacency = adjacency - spdiags(diag(adjacency), 0, nLocations, nLocations);
end

function [H, info] = build_local_kernel(adjacency, vertices, m, p, threshold)
nLocations = size(adjacency, 1);
hasVertices = ~isempty(vertices);

if hasVertices
    [row, col] = find(triu(adjacency, 1));
    edgeLengths = sqrt(sum((vertices(row, :) - vertices(col, :)).^2, 2));
    edgeLengths = edgeLengths(edgeLengths > 0 & isfinite(edgeLengths));
    if isempty(edgeLengths)
        baseDistance = 1;
    else
        % The paper's D is the regular voxel spacing.  Cortical meshes are
        % irregular, so the median adjacent-edge length is its robust analogue.
        baseDistance = median(edgeLengths);
    end
else
    baseDistance = 1;
end

supportDistance = sqrt(m) * baseDistance * (-log(threshold))^(1 / p);
maxHops = max(1, ceil(supportDistance / max(baseDistance, eps)) + 1);

rowCells = cell(nLocations, 1);
colCells = cell(nLocations, 1);
valueCells = cell(nLocations, 1);
for iLocation = 1:nLocations
    [candidates, hopDistance] = local_neighborhood(adjacency, iLocation, maxHops);
    if hasVertices
        distance = sqrt(sum((vertices(candidates, :) - vertices(iLocation, :)).^2, 2));
    else
        distance = hopDistance * baseDistance;
    end

    values = exp(-(distance ./ (sqrt(m) * baseDistance)).^p);
    keep = values >= threshold;
    candidates = candidates(keep);
    values = values(keep);

    if any(candidates == iLocation)
        kernelCols = candidates;
        kernelValues = values;
    else
        kernelCols = [iLocation; candidates];
        kernelValues = [1; values];
    end
    rowCells{iLocation} = repmat(iLocation, numel(kernelCols), 1);
    colCells{iLocation} = kernelCols;
    valueCells{iLocation} = kernelValues;
end

H = sparse(cell2mat(rowCells), cell2mat(colCells), cell2mat(valueCells), ...
    nLocations, nLocations);
H = max(H, H.');
H = H + spdiags(1 - diag(H), 0, nLocations, nLocations);

info = struct();
info.Width_m = m;
info.Power_p = p;
info.Threshold = threshold;
info.BaseDistance = baseDistance;
info.SupportDistance = supportDistance;
info.MaxHops = maxHops;
info.NumNonzeros = nnz(H);
info.UsesVertices = hasVertices;
end

function [visitedIndices, hopDistance] = local_neighborhood(adjacency, seed, maxHops)
nLocations = size(adjacency, 1);
visited = false(nLocations, 1);
distance = inf(nLocations, 1);
frontier = seed;
visited(seed) = true;
distance(seed) = 0;

for hop = 1:maxHops
    neighborMask = full(any(adjacency(frontier, :), 1)).';
    nextFrontier = find(neighborMask & ~visited);
    if isempty(nextFrontier)
        break;
    end
    visited(nextFrontier) = true;
    distance(nextFrontier) = hop;
    frontier = nextFrontier;
end

visitedIndices = find(visited);
hopDistance = distance(visitedIndices);
end

function [assignment, tileIndices] = build_tiles(adjacency, tileSize, requestedAssignment, nLocations)
if ~isempty(requestedAssignment)
    assignment = round(double(requestedAssignment(:)));
    if numel(assignment) ~= nLocations || any(~isfinite(assignment)) || any(assignment < 1)
        error('seal_SmoothChampagne:InvalidTileAssignments', ...
            'TileAssignments must contain one positive integer label per location.');
    end
    [~, ~, assignment] = unique(assignment, 'stable');
    nTiles = max(assignment);
    tileIndices = arrayfun(@(r) find(assignment == r), (1:nTiles).', ...
        'UniformOutput', false);
    return;
end

tileSize = max(1, round(tileSize));
assignment = zeros(nLocations, 1);
remaining = true(nLocations, 1);
tileIndices = {};
tileNumber = 0;
while any(remaining)
    seed = find(remaining, 1, 'first');
    tile = grow_tile(adjacency, seed, remaining, tileSize);
    tileNumber = tileNumber + 1;
    assignment(tile) = tileNumber;
    remaining(tile) = false;
    tileIndices{tileNumber, 1} = tile; %#ok<AGROW>
end
end

function tile = grow_tile(adjacency, seed, available, tileSize)
tile = seed;
selected = false(size(available));
selected(seed) = true;
frontier = seed;

while numel(tile) < tileSize && ~isempty(frontier)
    neighborMask = full(any(adjacency(frontier, :), 1)).';
    candidates = find(neighborMask & available & ~selected);
    if isempty(candidates)
        break;
    end
    nToTake = min(tileSize - numel(tile), numel(candidates));
    newNodes = candidates(1:nToTake);
    selected(newNodes) = true;
    tile = [tile; newNodes]; %#ok<AGROW>
    frontier = newNodes;
end
end

function Cnoise = initialize_noise_covariance(inputCovariance, Cy, nChannels, relativeLoading)
if isempty(inputCovariance)
    noiseScale = max(trace(Cy) / max(nChannels, 1), eps);
    Cnoise = 1e-2 * noiseScale * eye(nChannels);
else
    Cnoise = double(inputCovariance);
    if ~isequal(size(Cnoise), [nChannels, nChannels]) || any(~isfinite(Cnoise(:)))
        error('seal_SmoothChampagne:NoiseCovarianceDimensionMismatch', ...
            'NoiseCovariance must be a finite Nchannels-by-Nchannels matrix.');
    end
end
Cnoise = stabilize_spd(Cnoise, relativeLoading);
end

function Areg = stabilize_spd(A, relativeLoading)
A = full((A + A.') / 2);
n = size(A, 1);
scale = trace(A) / max(n, 1);
if ~isfinite(scale) || scale <= 0
    scale = norm(A, 'fro') / max(n, 1);
end
if ~isfinite(scale) || scale <= 0
    scale = 1;
end
jitter = max(relativeLoading * scale, eps(scale));
Areg = A + jitter * eye(n);
[~, flag] = chol(Areg);
attempt = 0;
while flag ~= 0 && attempt < 8
    jitter = 10 * jitter;
    Areg = A + jitter * eye(n);
    [~, flag] = chol(Areg);
    attempt = attempt + 1;
end
if flag ~= 0
    [V, D] = eig(A);
    values = max(real(diag(D)), jitter);
    Areg = real(V * diag(values) * V.');
    Areg = (Areg + Areg.') / 2;
end
end

function X = solve_spd(A, B)
[R, flag] = chol((A + A.') / 2);
if flag == 0
    X = R \ (R.' \ B);
else
    X = pinv(A) * B;
end
end

function value = spd_logdet(A)
[R, flag] = chol((A + A.') / 2);
if flag == 0
    value = 2 * sum(log(max(diag(R), realmin('double'))));
else
    values = eig((A + A.') / 2);
    value = sum(log(max(real(values), realmin('double'))));
end
end
