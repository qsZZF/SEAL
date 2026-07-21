function [Source_Estimate, Parameters] = seal_IRES(Data, L, sourceSpaceInfo, varargin)
%SEAL_IRES Iteratively reweighted edge sparsity minimization for E/MEG.
%   [S, PARAMETERS] = SEAL_IRES(DATA, L, SOURCESPACEINFO, ...) implements
%   the IRES formulation of Sohrabpour et al. (2016):
%
%       min ||Wd*V*S||_1 + Alpha*||Ws*S||_1
%       s.t. ||Data - L*S||_{Sigma^-1}^2 <= NoiseBound.
%
%   V is a cortical incidence matrix.  The outer loop updates source and
%   edge weights as 1 / (amplitude + epsilon).  For NumOrientations = 3,
%   the VIRES extension uses orientation-vector magnitudes in both terms.
%   Columns of DATA are solved independently but vectorized together; no
%   temporal regularizer is imposed.
%
%   Required inputs
%     Data              Nchannels-by-Ntime data matrix.
%     L                 Nchannels-by-(Nlocations*nd) lead field matrix.
%     sourceSpaceInfo   Struct with VertConn or EdgeOperator.
%
%   Name-value options
%     Alpha                 Source sparsity balance (default 0.05).
%     MaxIterations         Maximum reweighting iterations (default 4).
%     MaxADMMIterations     Maximum ADMM iterations per reweighting step.
%     Tolerance             Relative outer-source tolerance (default 1e-3).
%     ADMMTolerance         ADMM residual tolerance (default 1e-4).
%     ProgressInterval      ADMM iterations between progress messages (25).
%     Rho                   Shared ADMM penalty parameter (default 100).
%     OverRelaxation        ADMM relaxation in [1, 2) (default 1.5).
%     NoiseCovariance       Optional sensor noise covariance.
%     PrewhitenNoise        Whiten with noise covariance or robust scale.
%     NoiseBound            Squared residual bound per time point ([] uses
%                           the NoiseConfidence chi-square quantile).
%     NoiseConfidence       Family-wise confidence for the default bound
%                           across all time points (default .99).
%     WeightEpsilon         Relative reweighting epsilon (default 1e-3).
%     NormalizeLeadfield    Normalize each location lead-field block.
%     NumOrientations       1 (IRES) or 3 (VIRES extension).
%     Verbose               Print concise iteration diagnostics.
%     CancelCheck           Optional zero-input function that returns true
%                           when the caller requests cancellation.
%
%   Reference
%     A. Sohrabpour et al., NeuroImage 142, 27-42, 2016.
%     doi:10.1016/j.neuroimage.2016.05.064

p = inputParser;
p.CaseSensitive = false;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(p, 'sourceSpaceInfo', @isstruct);
addParameter(p, 'Alpha', 0.05, @(x) isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'MaxIterations', 4, @(x) isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'MaxADMMIterations', 200, @(x) isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'Tolerance', 1e-3, @(x) isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'ADMMTolerance', 1e-4, @(x) isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'ProgressInterval', 25, ...
    @(x) isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
addParameter(p, 'Rho', 100, @(x) isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'OverRelaxation', 1.5, ...
    @(x) isscalar(x) && isfinite(x) && x >= 1 && x < 2);
addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'PrewhitenNoise', true, @is_logical_like);
addParameter(p, 'NoiseBound', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'NoiseConfidence', 0.99, ...
    @(x) isscalar(x) && isfinite(x) && x > 0 && x < 1);
addParameter(p, 'WeightEpsilon', 1e-3, ...
    @(x) isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'NormalizeLeadfield', true, @is_logical_like);
addParameter(p, 'NumOrientations', 1, @(x) isscalar(x) && ismember(x, [1, 3]));
addParameter(p, 'Verbose', true, @is_logical_like);
addParameter(p, 'CancelCheck', [], @(x) isempty(x) || isa(x, 'function_handle'));
parse(p, Data, L, sourceSpaceInfo, varargin{:});

opts = p.Results;
Data = double(Data);
L = double(L);
nd = round(double(opts.NumOrientations));
maxOuter = max(1, round(double(opts.MaxIterations)));
maxADMM = max(1, round(double(opts.MaxADMMIterations)));
alpha = double(opts.Alpha);
rho = double(opts.Rho);
overRelaxation = double(opts.OverRelaxation);
outerTolerance = double(opts.Tolerance);
admmTolerance = double(opts.ADMMTolerance);
verbose = local_to_logical(opts.Verbose);

if any(~isfinite(Data(:))) || any(~isfinite(L(:)))
    error('seal_IRES:NonFiniteInput', 'Data and L must contain only finite values.');
end

[nChannels, nTime] = size(Data);
if size(L, 1) ~= nChannels
    error('seal_IRES:DimensionMismatch', ...
        'Data and L must have the same number of sensor rows.');
end

nComponents = size(L, 2);
if mod(nComponents, nd) ~= 0
    error('seal_IRES:OrientationMismatch', ...
        'Lead-field columns must be divisible by NumOrientations.');
end
nLocations = nComponents / nd;

if verbose
    initializationTimer = tic;
    fprintf('Preparing IRES: whitening data and normalizing lead field...\n');
    drawnow limitrate;
end
[DataW, LW, whiteningInfo] = whiten_inputs(Data, L, opts.NoiseCovariance, ...
    local_to_logical(opts.PrewhitenNoise));
local_check_cancel(opts.CancelCheck);
[LW, leadfieldScale] = normalize_leadfield_blocks(LW, nLocations, nd, ...
    local_to_logical(opts.NormalizeLeadfield));
if verbose
    fprintf('IRES: constructing the cortical edge operator...\n');
    drawnow limitrate;
end
[V, edgeInfo] = resolve_edge_operator(sourceSpaceInfo, nLocations, nd);
local_check_cancel(opts.CancelCheck);
nEdgeComponents = size(V, 1);
nEdges = nEdgeComponents / nd;

if mod(nEdgeComponents, nd) ~= 0
    error('seal_IRES:EdgeOperatorRows', ...
        'The expanded edge operator must have a multiple of NumOrientations rows.');
end

[noiseBoundSquared, perTimeNoiseConfidence] = resolve_noise_bound( ...
    opts.NoiseBound, opts.NoiseConfidence, whiteningInfo.ConstraintScale, ...
    nChannels, nTime);
noiseRadius = sqrt(noiseBoundSquared);

% The Woodbury system avoids forming the dense L''*L matrix.
if verbose
    fprintf('IRES: factorizing the %d-by-%d sparse source system...\n', nComponents, nComponents);
    drawnow limitrate;
end
A0 = rho * (V' * V + speye(nComponents));
[R0, factorFlag] = chol(A0, 'lower');
local_check_cancel(opts.CancelCheck);
if factorFlag ~= 0
    error('seal_IRES:FactorizationFailure', ...
        'Unable to factor the ADMM source-space system.');
end

solveA0 = @(B) R0' \ (R0 \ B);
if verbose
    fprintf('IRES: constructing and factorizing the %d-by-%d sensor Woodbury system...\n', ...
        nChannels, nChannels);
    drawnow limitrate;
end
A0InvLt = solveA0(LW.');
sensorSystem = (eye(nChannels) / rho) + LW * A0InvLt;
sensorSystem = (sensorSystem + sensorSystem.') / 2;
[Rs, sensorFlag] = chol(sensorSystem, 'lower');
if sensorFlag ~= 0
    sensorSystem = sensorSystem + max(eps(norm(sensorSystem, 'fro')), 1e-12) * eye(nChannels);
    [Rs, sensorFlag] = chol(sensorSystem, 'lower');
    if sensorFlag ~= 0
        error('seal_IRES:SensorFactorizationFailure', ...
            'Unable to factor the Woodbury sensor system.');
    end
end
local_check_cancel(opts.CancelCheck);

S = zeros(nComponents, nTime);
sourceWeights = ones(nLocations, nTime);
edgeWeights = ones(nEdges, nTime);
outerHistory = repmat(empty_outer_history(), maxOuter, 1);
admmHistory = cell(maxOuter, 1);
converged = false;
relativeChange = inf;

if verbose
    fprintf(['Running IRES (locations=%d, edges=%d, nd=%d, T=%d; ', ...
        'reweights=%d, ADMM/reweight=%d; initialization %.1fs)...\n'], ...
        nLocations, nEdges, nd, nTime, maxOuter, maxADMM, toc(initializationTimer));
    drawnow limitrate;
end

for outerIter = 1:maxOuter
    local_check_cancel(opts.CancelCheck);
    previousS = S;
    if verbose
        outerTimer = tic;
        fprintf('  Reweight %d/%d: starting weighted ADMM...\n', outerIter, maxOuter);
        drawnow limitrate;
    end
    [S, admmInfo] = solve_weighted_ires_step(S, DataW, LW, V, ...
        sourceWeights, edgeWeights, alpha, rho, noiseRadius, maxADMM, ...
        admmTolerance, overRelaxation, solveA0, A0InvLt, Rs, verbose, ...
        outerIter, maxOuter, opts.ProgressInterval, opts.CancelCheck);

    edgeEstimate = V * S;
    sourceMagnitude = group_magnitude(S, nd);
    edgeMagnitude = group_magnitude(edgeEstimate, nd);
    sourcePenalty = sourceWeights .* sourceMagnitude;
    edgePenalty = edgeWeights .* edgeMagnitude;
    objective = alpha * sum(sourcePenalty(:)) + sum(edgePenalty(:));
    residualSquared = sum((DataW - LW * S).^2, 1);

    if outerIter > 1
        relativeChange = norm(S - previousS, 'fro') / max(norm(S, 'fro'), eps);
    end

    [nextSourceWeights, nextEdgeWeights, epsilonAbsolute, weightReference] = ...
        update_weights(sourceMagnitude, edgeMagnitude, opts.WeightEpsilon);

    outerHistory(outerIter).Objective = objective;
    outerHistory(outerIter).RelativeSourceChange = relativeChange;
    outerHistory(outerIter).MaximumResidualSquared = max(residualSquared);
    outerHistory(outerIter).ConstraintSatisfied = all( ...
        residualSquared <= noiseBoundSquared .* (1 + 10 * admmTolerance));
    outerHistory(outerIter).WeightEpsilonAbsolute = epsilonAbsolute;
    outerHistory(outerIter).WeightReference = weightReference;
    outerHistory(outerIter).ADMMConverged = admmInfo.Converged;
    outerHistory(outerIter).ADMMIterations = admmInfo.Iterations;
    admmHistory{outerIter} = admmInfo;

    if verbose
        fprintf(['  Reweight %d: objective %.6g, relative change %.3g, ', ...
            'ADMM %d (%d), max residual^2 %.6g, elapsed %.1fs\n'], outerIter, objective, ...
            relativeChange, admmInfo.Iterations, admmInfo.Converged, ...
            outerHistory(outerIter).MaximumResidualSquared, toc(outerTimer));
        drawnow limitrate;
    end

    if outerIter > 1 && relativeChange <= outerTolerance && admmInfo.Converged
        converged = true;
        break;
    end

    sourceWeights = nextSourceWeights;
    edgeWeights = nextEdgeWeights;
end

outerHistory = outerHistory(1:outerIter);
admmHistory = admmHistory(1:outerIter);
Source_Estimate = restore_leadfield_scale(S, leadfieldScale, nd);
finalResidualSquared = sum((DataW - LW * S).^2, 1);

Parameters = struct();
Parameters.Algorithm = 'IRES';
Parameters.OptionsPassed = opts;
Parameters.OptionsPassed.PrewhitenNoise = local_to_logical(opts.PrewhitenNoise);
Parameters.OptionsPassed.NormalizeLeadfield = local_to_logical(opts.NormalizeLeadfield);
Parameters.EdgeOperator = V;
Parameters.EdgeOperatorInfo = edgeInfo;
Parameters.SourceWeights = sourceWeights;
Parameters.EdgeWeights = edgeWeights;
Parameters.WhitenedSourceEstimate = S;
Parameters.WhitenedEdgeEstimate = V * S;
Parameters.NoiseBoundSquared = noiseBoundSquared;
Parameters.PerTimeNoiseConfidence = perTimeNoiseConfidence;
Parameters.FinalResidualSquared = finalResidualSquared;
Parameters.DataConstraintSatisfied = all( ...
    finalResidualSquared <= noiseBoundSquared .* (1 + 10 * admmTolerance));
Parameters.OuterHistory = outerHistory;
Parameters.ADMMHistory = admmHistory;
Parameters.NumIterationsPerformed = outerIter;
Parameters.Converged = converged;
Parameters.RelativeChange = relativeChange;
Parameters.NoiseWhitening = whiteningInfo;
Parameters.LeadfieldNormalization = struct( ...
    'Applied', local_to_logical(opts.NormalizeLeadfield), ...
    'Scale', leadfieldScale, ...
    'ReturnedOriginalUnits', local_to_logical(opts.NormalizeLeadfield));
Parameters.NumOrientations = nd;
Parameters.Reference = 'doi:10.1016/j.neuroimage.2016.05.064';

end

function [S, info] = solve_weighted_ires_step(S, Data, L, V, sourceWeights, ...
    edgeWeights, alpha, rho, noiseRadius, maxIter, tolerance, relaxation, ...
    solveA0, A0InvLt, Rs, verbose, outerIter, maxOuter, progressInterval, cancelCheck)

nComponents = size(S, 1);
nTime = size(S, 2);
nEdgeComponents = size(V, 1);
nd = nComponents / size(sourceWeights, 1);

Y = V * S;
Z = S;
Q = project_residual_ball(L * S, Data, noiseRadius);
Ud = zeros(nEdgeComponents, nTime);
Us = zeros(nComponents, nTime);
Uq = zeros(size(Data));

info = struct('Iterations', 0, 'Converged', false, 'PrimalResidual', inf, ...
    'DualResidual', inf, 'PrimalTolerance', inf, 'DualTolerance', inf);
stepTimer = tic;

for iter = 1:maxIter
    local_check_cancel(cancelCheck);
    Yold = Y;
    Zold = Z;
    Qold = Q;

    rhs = rho * (V' * (Y - Ud) + (Z - Us) + L' * (Q - Uq));
    A0InvRhs = solveA0(rhs);
    correction = Rs' \ (Rs \ (L * A0InvRhs));
    S = A0InvRhs - A0InvLt * correction;

    edgeValue = V * S;
    sourceValue = S;
    dataValue = L * S;
    relaxedEdge = relaxation * edgeValue + (1 - relaxation) * Yold;
    relaxedSource = relaxation * sourceValue + (1 - relaxation) * Zold;
    relaxedData = relaxation * dataValue + (1 - relaxation) * Qold;

    Y = group_soft_threshold(relaxedEdge + Ud, edgeWeights / rho, nd);
    Z = group_soft_threshold(relaxedSource + Us, alpha * sourceWeights / rho, nd);
    Q = project_residual_ball(relaxedData + Uq, Data, noiseRadius);

    primalEdge = edgeValue - Y;
    primalSource = sourceValue - Z;
    primalData = dataValue - Q;
    Ud = Ud + relaxedEdge - Y;
    Us = Us + relaxedSource - Z;
    Uq = Uq + relaxedData - Q;

    primalResidual = sqrt(norm(primalEdge, 'fro')^2 + ...
        norm(primalSource, 'fro')^2 + norm(primalData, 'fro')^2);
    dualEdge = rho * V' * (Y - Yold);
    dualSource = rho * (Z - Zold);
    dualData = rho * L' * (Q - Qold);
    dualResidual = sqrt(norm(dualEdge, 'fro')^2 + ...
        norm(dualSource, 'fro')^2 + norm(dualData, 'fro')^2);

    primalTolerance = sqrt((nEdgeComponents + nComponents + size(Data, 1)) * nTime) * ...
        tolerance + tolerance * max( ...
        sqrt(norm(V * S, 'fro')^2 + norm(S, 'fro')^2 + norm(L * S, 'fro')^2), ...
        sqrt(norm(Y, 'fro')^2 + norm(Z, 'fro')^2 + norm(Q, 'fro')^2));
    dualTolerance = sqrt(nComponents * nTime) * tolerance + tolerance * rho * ...
        sqrt(norm(V' * Ud, 'fro')^2 + norm(Us, 'fro')^2 + norm(L' * Uq, 'fro')^2);

    info.Iterations = iter;
    info.PrimalResidual = primalResidual;
    info.DualResidual = dualResidual;
    info.PrimalTolerance = primalTolerance;
    info.DualTolerance = dualTolerance;

    if primalResidual <= primalTolerance && dualResidual <= dualTolerance
        info.Converged = true;
        if verbose
            fprintf(['    IRES ADMM %d/%d (reweight %d/%d): r %.3g/%.3g, ', ...
                's %.3g/%.3g, converged (%.1fs)\n'], iter, maxIter, ...
                outerIter, maxOuter, primalResidual, primalTolerance, ...
                dualResidual, dualTolerance, toc(stepTimer));
            drawnow limitrate;
        end
        break;
    end

    if verbose && (iter == 1 || mod(iter, progressInterval) == 0)
        fprintf(['    IRES ADMM %d/%d (reweight %d/%d): r %.3g/%.3g, ', ...
            's %.3g/%.3g (%.1fs)\n'], iter, maxIter, outerIter, maxOuter, ...
            primalResidual, primalTolerance, dualResidual, dualTolerance, toc(stepTimer));
        drawnow limitrate;
    end
end

if verbose && ~info.Converged
    fprintf('    ADMM stopped at %d iterations (r=%.3g, s=%.3g).\n', ...
        info.Iterations, info.PrimalResidual, info.DualResidual);
end
end

function local_check_cancel(cancelCheck)
drawnow limitrate;
if ~isempty(cancelCheck) && cancelCheck()
    error('SEAL:Cancelled', 'Source-imaging run was cancelled by the user.');
end
end

function Q = project_residual_ball(candidate, data, radii)
residual = candidate - data;
residualNorm = sqrt(sum(residual.^2, 1));
scale = min(1, radii ./ max(residualNorm, eps));
Q = data + residual .* scale;
end

function output = group_soft_threshold(input, thresholds, nd)
nRows = size(input, 1);
nTime = size(input, 2);
if mod(nRows, nd) ~= 0
    error('seal_IRES:GroupSizeMismatch', ...
        'Input rows must be divisible by NumOrientations.');
end
nGroups = nRows / nd;
if ~isequal(size(thresholds), [nGroups, nTime])
    error('seal_IRES:ThresholdSizeMismatch', ...
        'Group threshold matrix must be Ngroups-by-Ntime.');
end

inputTensor = reshape(input, nd, nGroups, nTime);
magnitude = reshape(sqrt(sum(inputTensor.^2, 1)), nGroups, nTime);
shrinkage = max(1 - thresholds ./ max(magnitude, eps), 0);
outputTensor = inputTensor .* reshape(shrinkage, 1, nGroups, nTime);
output = reshape(outputTensor, nRows, nTime);
end

function magnitude = group_magnitude(input, nd)
nRows = size(input, 1);
nTime = size(input, 2);
if mod(nRows, nd) ~= 0
    error('seal_IRES:GroupSizeMismatch', ...
        'Input rows must be divisible by NumOrientations.');
end
nGroups = nRows / nd;
inputTensor = reshape(input, nd, nGroups, nTime);
magnitude = reshape(sqrt(sum(inputTensor.^2, 1)), nGroups, nTime);
end

function [sourceWeights, edgeWeights, epsilonAbsolute, reference] = ...
    update_weights(sourceMagnitude, edgeMagnitude, relativeEpsilon)

reference = max([sourceMagnitude(:); edgeMagnitude(:)]);
if ~isfinite(reference) || reference <= 0
    reference = 1;
end
epsilonAbsolute = max(relativeEpsilon * reference, eps(reference));

% This is the paper's 1/(amplitude + epsilon) rule, scaled by a common
% positive factor so the constrained objective remains numerically stable.
sourceWeights = reference ./ (sourceMagnitude + epsilonAbsolute);
edgeWeights = reference ./ (edgeMagnitude + epsilonAbsolute);
end

function [bounds, perTimeConfidence] = resolve_noise_bound(option, confidence, constraintScale, ...
    nChannels, nTime)
if isempty(option)
    % The paper's constraint is defined for one time point.  The SEAL
    % matrix interface enforces it at every column, so use a Bonferroni
    % adjustment to retain the requested family-wise confidence.
    perTimeConfidence = 1 - (1 - confidence) / max(nTime, 1);
    defaultBound = local_chi2inv(perTimeConfidence, nChannels) * constraintScale^2;
    bounds = defaultBound * ones(1, nTime);
else
    perTimeConfidence = NaN;
    bounds = double(option(:)).';
    if isscalar(bounds)
        bounds = repmat(bounds, 1, nTime);
    elseif numel(bounds) ~= nTime
        error('seal_IRES:NoiseBoundSize', ...
            'NoiseBound must be scalar or contain one squared bound per time point.');
    end
    if any(~isfinite(bounds)) || any(bounds <= 0)
        error('seal_IRES:InvalidNoiseBound', ...
            'NoiseBound must contain finite positive squared residual bounds.');
    end
end
end

function [DataOut, LOut, info] = whiten_inputs(Data, L, inputCovariance, doWhiten)
DataOut = Data;
LOut = L;
info = struct('Applied', false, 'Method', 'none', 'NoiseScale', 1, ...
    'ConstraintScale', 1, 'ConditionNumber', NaN, 'Rank', []);

channelScale = robust_channel_scale(Data);
validScale = channelScale(isfinite(channelScale) & channelScale > 0);
if isempty(validScale)
    robustScale = 1;
else
    robustScale = median(validScale);
end

if isempty(inputCovariance)
    if doWhiten
        DataOut = Data / robustScale;
        LOut = L / robustScale;
        info.Applied = true;
        info.Method = 'robust_global_scale';
        info.NoiseScale = robustScale;
        info.ConstraintScale = 1;
    else
        info.Method = 'robust_scale_not_applied';
        info.NoiseScale = robustScale;
        info.ConstraintScale = robustScale;
    end
    return;
end

nChannels = size(Data, 1);
covariance = double(inputCovariance);
if ~isequal(size(covariance), [nChannels, nChannels]) || ...
        any(~isfinite(covariance(:)))
    error('seal_IRES:NoiseCovarianceMismatch', ...
        'NoiseCovariance must be a finite Nchannels-by-Nchannels matrix.');
end
covariance = (covariance + covariance.') / 2;
noiseScale = sqrt(max(trace(covariance) / max(nChannels, 1), eps));

if ~doWhiten
    info.Method = 'covariance_not_applied';
    info.NoiseScale = noiseScale;
    info.ConstraintScale = noiseScale;
    return;
end

[U, S, ~] = svd(covariance, 'econ');
values = real(diag(S));
if isempty(values) || values(1) <= 0
    error('seal_IRES:InvalidNoiseCovariance', ...
        'NoiseCovariance must have positive finite energy.');
end
floorValue = max(values(1) * 1e-12, eps(values(1)));
rankTolerance = max(size(covariance)) * eps(values(1));
rankNoise = sum(values > rankTolerance);
values(values < floorValue) = floorValue;
whitener = diag(1 ./ sqrt(values)) * U.';
DataOut = whitener * Data;
LOut = whitener * L;

info.Applied = true;
info.Method = 'noise_covariance';
info.NoiseScale = noiseScale;
info.ConstraintScale = 1;
info.ConditionNumber = values(1) / values(end);
info.Rank = rankNoise;
end

function scale = robust_channel_scale(Data)
centered = Data - median(Data, 2);
scale = median(abs(centered), 2) / 0.6745;
scale(~isfinite(scale) | scale <= 0) = NaN;
end

function [Lout, scale] = normalize_leadfield_blocks(Lin, nLocations, nd, doNormalize)
Lout = Lin;
scale = ones(nLocations, 1);
if ~doNormalize
    return;
end

for location = 1:nLocations
    index = (location - 1) * nd + (1:nd);
    blockNorm = norm(Lin(:, index), 'fro');
    if ~isfinite(blockNorm) || blockNorm <= 0
        error('seal_IRES:InvalidLeadfieldBlock', ...
            'Lead-field block %d has zero or invalid norm.', location);
    end
    scale(location) = blockNorm;
    Lout(:, index) = Lin(:, index) / blockNorm;
end
end

function Sout = restore_leadfield_scale(Sin, scale, nd)
Sout = Sin;
for location = 1:numel(scale)
    index = (location - 1) * nd + (1:nd);
    Sout(index, :) = Sin(index, :) / scale(location);
end
end

function [V, info] = resolve_edge_operator(sourceSpaceInfo, nLocations, nd)
nComponents = nLocations * nd;
info = struct('Mode', '', 'NumLocationEdges', [], 'NumComponentEdges', []);

if isfield(sourceSpaceInfo, 'EdgeOperator') && ~isempty(sourceSpaceInfo.EdgeOperator)
    Vinput = sparse(double(sourceSpaceInfo.EdgeOperator));
    if size(Vinput, 2) == nComponents
        V = Vinput;
        info.Mode = 'component_edge_operator';
    elseif size(Vinput, 2) == nLocations
        V = kron(Vinput, speye(nd));
        info.Mode = 'location_edge_operator_expanded';
    else
        error('seal_IRES:EdgeOperatorSize', ...
            'EdgeOperator must have %d location or %d component columns.', ...
            nLocations, nComponents);
    end
elseif isfield(sourceSpaceInfo, 'VertConn') && ~isempty(sourceSpaceInfo.VertConn)
    Vlocation = incidence_from_connectivity(sourceSpaceInfo.VertConn, nLocations);
    V = kron(Vlocation, speye(nd));
    info.Mode = 'vertconn_incidence_expanded';
else
    error('seal_IRES:MissingSourceGraph', ...
        'sourceSpaceInfo must contain VertConn or EdgeOperator for IRES.');
end

if mod(size(V, 1), nd) ~= 0
    error('seal_IRES:EdgeOperatorRows', ...
        'EdgeOperator rows must be divisible by NumOrientations.');
end
info.NumLocationEdges = size(V, 1) / nd;
info.NumComponentEdges = size(V, 1);
end

function V = incidence_from_connectivity(connectivity, nLocations)
if ~isequal(size(connectivity), [nLocations, nLocations])
    error('seal_IRES:VertConnSizeMismatch', ...
        'VertConn must be Nlocations-by-Nlocations.');
end
if any(~isfinite(double(connectivity(:))))
    error('seal_IRES:InvalidVertConn', 'VertConn must contain finite values.');
end

A = spones(sparse(connectivity));
A = spones(A | A.');
A = A - spdiags(diag(A), 0, nLocations, nLocations);
[row, col] = find(triu(A, 1));
nEdges = numel(row);
if nEdges == 0
    error('seal_IRES:EmptyGraph', ...
        'VertConn must contain at least one off-diagonal edge.');
end
edgeIndex = (1:nEdges).';
V = sparse([edgeIndex; edgeIndex], [row; col], ...
    [-ones(nEdges, 1); ones(nEdges, 1)], nEdges, nLocations);
end

function entry = empty_outer_history()
entry = struct('Objective', NaN, 'RelativeSourceChange', inf, ...
    'MaximumResidualSquared', NaN, 'ConstraintSatisfied', false, ...
    'WeightEpsilonAbsolute', NaN, 'WeightReference', NaN, ...
    'ADMMConverged', false, 'ADMMIterations', 0);
end

function value = local_chi2inv(probability, degreesOfFreedom)
if exist('chi2inv', 'file') == 2
    value = chi2inv(probability, degreesOfFreedom);
else
    value = 2 * gammaincinv(probability, degreesOfFreedom / 2);
end
end

function tf = is_logical_like(value)
tf = (islogical(value) && isscalar(value)) || ...
    (isnumeric(value) && isscalar(value)) || ...
    (ischar(value) || (isstring(value) && isscalar(value)));
end

function value = local_to_logical(value)
if islogical(value)
    value = value(1);
elseif isnumeric(value)
    value = value(1) ~= 0;
else
    value = strtrim(char(string(value)));
    value = any(strcmpi(value, {'true', '1', 'yes', 'on'}));
end
value = logical(value);
end
