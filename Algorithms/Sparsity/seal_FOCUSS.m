function [Source, Parameters] = seal_FOCUSS(Data, LeadField, varargin)
%SEAL_FOCUSS Recursive weighted minimum-norm electromagnetic source imaging.
%
%   [Source, Parameters] = seal_FOCUSS(Data, LeadField, ...) estimates
%   spatially sparse sources with FOCUSS (Gorodnitsky et al., 1995).  The
%   original algorithm recursively solves a weighted minimum-norm problem,
%   using the previous source magnitude to update the diagonal weights.
%
%   This implementation uses one spatial weight per source location.  It is
%   therefore rotation invariant for free-orientation lead fields and shares
%   the support across time samples, which extends the original single-map
%   formulation to SEAL's spatio-temporal matrix interface.
%
%   Name-value options:
%       'NoiseCovariance'          Sensor noise covariance matrix ([])
%       'PrewhitenNoise'           Prewhiten when covariance is supplied (true)
%       'RegularizationParameter'  Relative Tikhonov regularization (0.005)
%       'RegularizationMode'       'fixed_initial' or 'adaptive_legacy'
%       'MaxIterations'            Maximum recursive updates (40)
%       'Tolerance'                Relative source-change tolerance (1e-4)
%       'PruneThreshold'           Relative weight pruning threshold (1e-6)
%       'WeightFloor'              Minimum recursive weight (1e-12)
%       'DepthWeighting'           'none' or 'leadfield_norm' (none)
%       'WeightMode'               'compound' or 'simple'
%       'NumOrientations'          Components per location (1 or 3)
%       'Verbose'                  Print iteration progress (true)

    parser = inputParser;
    parser.FunctionName = mfilename;
    addParameter(parser, 'NoiseCovariance', [], @(x) isempty(x) || ...
        (isnumeric(x) && ismatrix(x) && size(x, 1) == size(x, 2)));
    addParameter(parser, 'PrewhitenNoise', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(parser, 'RegularizationParameter', 0.005, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(parser, 'RegularizationMode', 'fixed_initial', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'MaxIterations', 40, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'Tolerance', 1e-4, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'PruneThreshold', 1e-6, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 1);
    addParameter(parser, 'WeightFloor', 1e-12, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'DepthWeighting', 'none', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'WeightMode', 'compound', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'NumOrientations', 1, ...
        @(x) isnumeric(x) && isscalar(x) && any(x == [1, 3]));
    addParameter(parser, 'Verbose', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    parse(parser, varargin{:});
    opt = parser.Results;

    validateattributes(Data, {'numeric'}, {'2d', 'nonempty', 'finite'}, mfilename, 'Data', 1);
    validateattributes(LeadField, {'numeric'}, {'2d', 'nonempty', 'finite'}, mfilename, 'LeadField', 2);
    if size(Data, 1) ~= size(LeadField, 1)
        error('seal_FOCUSS:DimensionMismatch', ...
            'Data and LeadField must have the same number of sensor channels.');
    end

    nComponents = size(LeadField, 2);
    nOrientations = opt.NumOrientations;
    if mod(nComponents, nOrientations) ~= 0
        error('seal_FOCUSS:InvalidOrientationCount', ...
            'LeadField columns must be divisible by NumOrientations.');
    end
    nLocations = nComponents / nOrientations;

    depthMode = lower(char(opt.DepthWeighting));
    if ~ismember(depthMode, {'leadfield_norm', 'none'})
        error('seal_FOCUSS:InvalidDepthWeighting', ...
            'DepthWeighting must be ''none'' or ''leadfield_norm''.');
    end
    regularizationMode = lower(char(opt.RegularizationMode));
    if ~ismember(regularizationMode, {'fixed_initial', 'adaptive_legacy'})
        error('seal_FOCUSS:InvalidRegularizationMode', ...
            'RegularizationMode must be ''fixed_initial'' or ''adaptive_legacy''.');
    end
    weightMode = lower(char(opt.WeightMode));
    if ~ismember(weightMode, {'compound', 'simple'})
        error('seal_FOCUSS:InvalidWeightMode', ...
            'WeightMode must be ''compound'' or ''simple''.');
    end

    [DataWork, LeadFieldWork, whitening] = localPrewhiten(Data, LeadField, ...
        opt.NoiseCovariance, logical(opt.PrewhitenNoise));

    [depthWeights, depthInfo] = localDepthWeights(LeadFieldWork, nOrientations, ...
        depthMode, opt.WeightFloor);
    recursiveWeights = ones(nLocations, 1);
    locationWeights = localNormalizeWeights(depthWeights .* recursiveWeights, opt.WeightFloor);
    fixedLambda = localFixedRegularization(LeadFieldWork, locationWeights, ...
        nOrientations, opt.RegularizationParameter, opt.WeightFloor);

    [Source, inverseOperator, activeLocations, lambda] = localWeightedMinimumNorm( ...
        DataWork, LeadFieldWork, locationWeights, nOrientations, ...
        fixedLambda, regularizationMode, opt.RegularizationParameter, ...
        opt.PruneThreshold, opt.WeightFloor);

    relativeChangeHistory = nan(opt.MaxIterations, 1);
    costHistory = nan(opt.MaxIterations, 1);
    regularizationHistory = nan(opt.MaxIterations, 1);
    activeLocationHistory = nan(opt.MaxIterations, 1);
    converged = false;

    for iter = 1:opt.MaxIterations
        sourceMagnitude = localLocationMagnitude(Source, nOrientations);
        maxMagnitude = max(sourceMagnitude);
        if ~isfinite(maxMagnitude) || maxMagnitude <= 0
            warning('seal_FOCUSS:ZeroSolution', ...
                'FOCUSS reached a zero or non-finite source estimate; stopping iteration.');
            break;
        end

        normalizedMagnitude = max(sourceMagnitude ./ maxMagnitude, opt.WeightFloor);
        switch weightMode
            case 'compound'
                recursiveWeights = recursiveWeights .* normalizedMagnitude;
            case 'simple'
                recursiveWeights = normalizedMagnitude;
        end
        recursiveWeights = localNormalizeWeights(recursiveWeights, opt.WeightFloor);
        newLocationWeights = localNormalizeWeights(depthWeights .* recursiveWeights, opt.WeightFloor);

        [newSource, inverseOperator, activeLocations, lambda] = localWeightedMinimumNorm( ...
            DataWork, LeadFieldWork, newLocationWeights, nOrientations, ...
            fixedLambda, regularizationMode, opt.RegularizationParameter, ...
            opt.PruneThreshold, opt.WeightFloor);

        relativeChange = norm(newSource - Source, 'fro') / max(norm(Source, 'fro'), opt.WeightFloor);
        residual = DataWork - LeadFieldWork * newSource;
        costHistory(iter) = norm(residual, 'fro')^2 + ...
            lambda * localWeightedEnergy(newSource, newLocationWeights, nOrientations, opt.WeightFloor);
        relativeChangeHistory(iter) = relativeChange;
        regularizationHistory(iter) = lambda;
        activeLocationHistory(iter) = numel(activeLocations);

        Source = newSource;
        locationWeights = newLocationWeights;

        if logical(opt.Verbose)
            fprintf('FOCUSS iter %d: rel_change = %.3e, active_source = %d, lambda = %.3e\n', ...
                iter, relativeChange, numel(activeLocations), lambda);
        end
        if relativeChange <= opt.Tolerance
            converged = true;
            break;
        end
    end

    nIterations = find(~isnan(relativeChangeHistory), 1, 'last');
    if isempty(nIterations)
        nIterations = 0;
    end
    relativeChangeHistory = relativeChangeHistory(1:nIterations);
    costHistory = costHistory(1:nIterations);
    regularizationHistory = regularizationHistory(1:nIterations);
    activeLocationHistory = activeLocationHistory(1:nIterations);

    Parameters = struct();
    Parameters.Method = 'FOCUSS';
    Parameters.Reference = 'Gorodnitsky IF, George JS, Rao BD (1995), Electroencephalography and Clinical Neurophysiology 95:231-251.';
    Parameters.InverseOperator = inverseOperator;
    Parameters.NumOrientations = nOrientations;
    Parameters.DepthWeighting = depthMode;
    Parameters.DepthWeightingInfo = depthInfo;
    Parameters.WeightMode = weightMode;
    Parameters.StaticDepthWeights = depthWeights;
    Parameters.RecursiveWeights = recursiveWeights;
    Parameters.LocationWeights = locationWeights;
    Parameters.ActiveLocations = activeLocations;
    Parameters.RegularizationParameter = opt.RegularizationParameter;
    Parameters.RegularizationMode = regularizationMode;
    Parameters.FixedRegularization = fixedLambda;
    Parameters.FinalRegularization = lambda;
    Parameters.RegularizationHistory = regularizationHistory;
    Parameters.RelativeChangeHistory = relativeChangeHistory;
    Parameters.CostHistory = costHistory;
    Parameters.ActiveLocationHistory = activeLocationHistory;
    Parameters.NumIterationsPerformed = nIterations;
    Parameters.Converged = converged;
    Parameters.NoiseWhitening = whitening;
end

function [DataWork, LeadFieldWork, info] = localPrewhiten(Data, LeadField, noiseCovariance, usePrewhitening)
    nChannels = size(Data, 1);
    info = struct('Applied', false, 'Matrix', eye(nChannels));
    DataWork = Data;
    LeadFieldWork = LeadField;
    if isempty(noiseCovariance) || ~usePrewhitening
        return;
    end
    if ~isequal(size(noiseCovariance), [nChannels, nChannels])
        error('seal_FOCUSS:InvalidNoiseCovariance', ...
            'NoiseCovariance must be a square matrix with one row per sensor channel.');
    end

    noiseCovariance = (noiseCovariance + noiseCovariance.') / 2;
    [factor, flag] = chol(noiseCovariance, 'lower');
    if flag ~= 0
        jitter = max(trace(noiseCovariance) / max(nChannels, 1), 1) * 1e-10;
        [factor, flag] = chol(noiseCovariance + jitter * eye(nChannels), 'lower');
    end
    if flag ~= 0
        error('seal_FOCUSS:InvalidNoiseCovariance', ...
            'NoiseCovariance must be positive definite for prewhitening.');
    end

    whiteningMatrix = factor \ eye(nChannels);
    DataWork = whiteningMatrix * Data;
    LeadFieldWork = whiteningMatrix * LeadField;
    info = struct('Applied', true, 'Matrix', whiteningMatrix);
end

function [weights, info] = localDepthWeights(LeadField, nOrientations, mode, weightFloor)
    nLocations = size(LeadField, 2) / nOrientations;
    info = struct('Mode', mode, 'InitialWeights', ones(nLocations, 1));
    switch mode
        case 'none'
            weights = ones(nLocations, 1);
            info.InitialWeights = weights;
        case 'leadfield_norm'
            blockNorm = zeros(nLocations, 1);
            for location = 1:nLocations
                componentIndex = (location - 1) * nOrientations + (1:nOrientations);
                blockNorm(location) = norm(LeadField(:, componentIndex), 'fro');
            end
            weights = 1 ./ max(blockNorm, weightFloor);
            weights = localNormalizeWeights(weights, weightFloor);
            info.InitialWeights = weights;
    end
    weights = localNormalizeWeights(weights, weightFloor);
end

function weights = localNormalizeWeights(weights, weightFloor)
    weights = max(real(weights(:)), weightFloor);
    scale = max(weights);
    if ~isfinite(scale) || scale <= 0
        weights = ones(size(weights));
    else
        weights = weights / scale;
    end
end

function lambda = localFixedRegularization(LeadField, locationWeights, nOrientations, ...
        regularizationParameter, weightFloor)
    nChannels = size(LeadField, 1);
    componentWeights = repelem(locationWeights, nOrientations);
    weightedLeadField = LeadField .* componentWeights.';
    scale = trace(weightedLeadField * weightedLeadField.') / max(nChannels, 1);
    lambda = max(regularizationParameter * scale, weightFloor);
end

function [Source, inverseOperator, activeLocations, lambda] = localWeightedMinimumNorm( ...
        Data, LeadField, locationWeights, nOrientations, fixedLambda, ...
        regularizationMode, regularizationParameter, pruneThreshold, weightFloor)
    nChannels = size(LeadField, 1);
    nComponents = size(LeadField, 2);
    maxWeight = max(locationWeights);
    activeMask = locationWeights >= maxWeight * pruneThreshold;
    if ~any(activeMask)
        [~, strongestLocation] = max(locationWeights);
        activeMask(strongestLocation) = true;
    end
    activeLocations = find(activeMask);
    activeComponents = reshape((activeLocations(:).' - 1) * nOrientations + (1:nOrientations).', [], 1);
    componentWeights = repelem(locationWeights(activeLocations), nOrientations);

    activeLeadField = LeadField(:, activeComponents);
    weightedLeadField = activeLeadField .* componentWeights.';
    sensorCovariance = weightedLeadField * weightedLeadField.';
    if strcmp(regularizationMode, 'adaptive_legacy')
        signalScale = trace(sensorCovariance) / max(nChannels, 1);
        lambda = max(regularizationParameter * signalScale, weightFloor);
    else
        lambda = fixedLambda;
    end
    sensorCovariance = sensorCovariance + lambda * eye(nChannels);

    [factor, flag] = chol((sensorCovariance + sensorCovariance.') / 2, 'lower');
    if flag ~= 0
        jitter = max(trace(sensorCovariance) / max(nChannels, 1), 1) * 1e-10;
        factor = chol(sensorCovariance + jitter * eye(nChannels), 'lower');
    end
    weightedTranspose = factor' \ (factor \ weightedLeadField);
    activeInverse = weightedTranspose.' .* componentWeights;

    Source = zeros(nComponents, size(Data, 2));
    Source(activeComponents, :) = activeInverse * Data;
    inverseOperator = zeros(nComponents, nChannels);
    inverseOperator(activeComponents, :) = activeInverse;
end

function magnitude = localLocationMagnitude(Source, nOrientations)
    nLocations = size(Source, 1) / nOrientations;
    sourceBlocks = reshape(Source, nOrientations, nLocations, []);
    magnitude = sqrt(reshape(sum(sum(sourceBlocks .^ 2, 1), 3), nLocations, 1));
end

function energy = localWeightedEnergy(Source, locationWeights, nOrientations, weightFloor)
    nLocations = size(Source, 1) / nOrientations;
    sourceBlocks = reshape(Source, nOrientations, nLocations, []);
    locationEnergy = reshape(sum(sum(sourceBlocks .^ 2, 1), 3), nLocations, 1);
    energy = sum(locationEnergy ./ max(locationWeights .^ 2, weightFloor));
end
