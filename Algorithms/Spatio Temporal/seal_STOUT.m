function [Source, Parameters] = seal_STOUT(Data, L, varargin)
%SEAL_STOUT Spectrally constrained source imaging using biological profiles.
%
%   [Source, Parameters] = seal_STOUT(Data, L, ...) implements the
%   multimodal source-imaging framework of Luan et al. (2015).  A weighted
%   minimum-norm estimate is first obtained from MEG/EEG.  Each cortical
%   source's baseline-to-active spectral change is then compared to a
%   biological reference spectrum, producing the relationship score xi in
%   [0, 2].  The source amplitudes are constrained with one of:
%
%       CoCo: xi^2 * s
%       ExCo: s^xi
%       HyCo: alpha^(-1) * CoCo + alpha * ExCo, alpha = SNR / Gamma
%
%   For free orientations, the constraints operate on the orientation-vector
%   magnitude and preserve the estimated orientation at each time point.
%
%   Name-value options:
%       'ReferenceSpectrum'         ECoG spectral-change vector or 'motor'
%                                   for the paper's motor beta/gamma template
%       'SamplingRate'              Sampling rate in Hz (200)
%       'BaselineSamples'           Baseline sample indices ([] uses first third)
%       'ActiveSamples'             Active sample indices ([] uses remaining samples)
%       'ConstraintMode'            'hyco', 'coco', or 'exco' ('hyco')
%       'SNR'                       Linear signal-to-noise ratio for HyCo ([] estimates it)
%       'HybridGamma'               HyCo transition value (10)
%       'RegularizationParameter'   Relative wMNE regularization (0.05)
%       'NoiseCovariance'           Optional sensor noise covariance ([])
%       'PrewhitenNoise'            Prewhiten when covariance is supplied (true)
%       'NumOrientations'           Components per location: 1 or 3 (1)
%       'Verbose'                   Print concise diagnostics (true)

    parser = inputParser;
    parser.FunctionName = mfilename;
    parser.CaseSensitive = false;
    addRequired(parser, 'Data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addRequired(parser, 'L', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addParameter(parser, 'ReferenceSpectrum', 'motor', ...
        @(x) isnumeric(x) || ischar(x) || isstring(x));
    addParameter(parser, 'SamplingRate', 200, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'BaselineSamples', [], @(x) isempty(x) || ...
        (isnumeric(x) && isvector(x) && all(isfinite(x))));
    addParameter(parser, 'ActiveSamples', [], @(x) isempty(x) || ...
        (isnumeric(x) && isvector(x) && all(isfinite(x))));
    addParameter(parser, 'ConstraintMode', 'hyco', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'SNR', [], @(x) isempty(x) || ...
        (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0));
    addParameter(parser, 'HybridGamma', 10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'RegularizationParameter', 0.05, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(parser, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(parser, 'PrewhitenNoise', true, @localIsLogicalLike);
    addParameter(parser, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && any(x == [1, 3]));
    addParameter(parser, 'ExponentFloor', 1e-8, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
    addParameter(parser, 'Verbose', true, @localIsLogicalLike);
    parse(parser, Data, L, varargin{:});
    opt = parser.Results;

    Data = double(Data);
    L = double(L);
    if any(~isfinite(Data(:))) || any(~isfinite(L(:)))
        error('seal_STOUT:NonFiniteInput', 'Data and L must contain only finite values.');
    end
    if size(Data, 1) ~= size(L, 1)
        error('seal_STOUT:DimensionMismatch', ...
            'Data and L must have the same number of sensor rows.');
    end

    nTime = size(Data, 2);
    nComponents = size(L, 2);
    nOrientations = round(double(opt.NumOrientations));
    if mod(nComponents, nOrientations) ~= 0
        error('seal_STOUT:OrientationMismatch', ...
            'Lead-field columns must be divisible by NumOrientations.');
    end
    nLocations = nComponents / nOrientations;
    [baselineSamples, activeSamples] = localResolveWindows( ...
        opt.BaselineSamples, opt.ActiveSamples, nTime);
    mode = lower(char(opt.ConstraintMode));
    if ~ismember(mode, {'coco', 'exco', 'hyco'})
        error('seal_STOUT:InvalidConstraintMode', ...
            'ConstraintMode must be ''coco'', ''exco'', or ''hyco''.');
    end

    [DataWork, LWork, whiteningInfo] = localPrewhiten(Data, L, ...
        opt.NoiseCovariance, localToLogical(opt.PrewhitenNoise));
    [initialSource, inverseOperator, depthWeights, lambda] = localWMNE( ...
        DataWork, LWork, nOrientations, opt.RegularizationParameter);
    sourceMagnitude = localLocationMagnitude(initialSource, nOrientations);
    [spectralChange, frequencies] = localSpectralChange(sourceMagnitude, ...
        baselineSamples, activeSamples, opt.SamplingRate);
    referenceSpectrum = localReferenceSpectrum(opt.ReferenceSpectrum, frequencies);
    relationshipScore = localRelationshipScore(spectralChange, referenceSpectrum);

    if isempty(opt.SNR)
        fittedData = LWork * initialSource;
        residual = DataWork - fittedData;
        snr = norm(fittedData, 'fro') / max(norm(residual, 'fro'), eps);
    else
        snr = double(opt.SNR);
    end
    snrDb = 20 * log10(max(snr, eps));
    alpha = max(snr / double(opt.HybridGamma), eps);
    [Source, constraintGain] = localApplyConstraint(initialSource, sourceMagnitude, ...
        relationshipScore, mode, alpha, opt.ExponentFloor, nOrientations);

    if localToLogical(opt.Verbose)
        fprintf(['STOUT (%s): locations = %d, SNR = %.3f (%.2f dB), alpha = %.3f, ', ...
            'mean relationship = %.3f\n'], mode, nLocations, snr, snrDb, alpha, mean(relationshipScore));
    end

    Parameters = struct();
    Parameters.Method = 'STOUT';
    Parameters.Reference = ['Luan F, Zhang Z, Zhang T (2015). A multimodal framework for integrating ', ...
        'biological spectral characteristics and MEG/EEG in brain-source imaging. ICNC.'];
    Parameters.DOI = '10.1109/ICNC.2015.7377983';
    Parameters.InitialSource = initialSource;
    Parameters.InverseOperator = inverseOperator;
    Parameters.DepthWeights = depthWeights;
    Parameters.RegularizationParameter = opt.RegularizationParameter;
    Parameters.FinalRegularization = lambda;
    Parameters.ConstraintMode = mode;
    Parameters.RelationshipScore = relationshipScore;
    Parameters.ConstraintGain = constraintGain;
    Parameters.SpectralChange = spectralChange;
    Parameters.ReferenceSpectrum = referenceSpectrum;
    Parameters.Frequencies = frequencies;
    Parameters.BaselineSamples = baselineSamples;
    Parameters.ActiveSamples = activeSamples;
    Parameters.SNR = snr;
    Parameters.SNRdB = snrDb;
    Parameters.HybridAlpha = alpha;
    Parameters.NumOrientations = nOrientations;
    Parameters.NoiseWhitening = whiteningInfo;
end

function [baselineSamples, activeSamples] = localResolveWindows(baselineSamples, activeSamples, nTime)
    if isempty(baselineSamples) && isempty(activeSamples)
        split = max(2, floor(nTime / 3));
        baselineSamples = 1:split;
        activeSamples = (split + 1):nTime;
    elseif isempty(baselineSamples)
        activeSamples = localValidateSamples(activeSamples, nTime, 'ActiveSamples');
        baselineSamples = setdiff(1:nTime, activeSamples);
    elseif isempty(activeSamples)
        baselineSamples = localValidateSamples(baselineSamples, nTime, 'BaselineSamples');
        activeSamples = setdiff(1:nTime, baselineSamples);
    else
        baselineSamples = localValidateSamples(baselineSamples, nTime, 'BaselineSamples');
        activeSamples = localValidateSamples(activeSamples, nTime, 'ActiveSamples');
    end
    if numel(baselineSamples) < 2 || numel(activeSamples) < 2
        error('seal_STOUT:WindowTooShort', ...
            'BaselineSamples and ActiveSamples must each contain at least two samples.');
    end
end

function samples = localValidateSamples(samples, nTime, name)
    samples = unique(round(double(samples(:).')));
    if any(samples < 1) || any(samples > nTime)
        error('seal_STOUT:InvalidSamples', '%s must be valid sample indices.', name);
    end
end

function [DataWork, LWork, info] = localPrewhiten(Data, L, noiseCovariance, usePrewhitening)
    nChannels = size(Data, 1);
    DataWork = Data;
    LWork = L;
    info = struct('Applied', false, 'Matrix', eye(nChannels));
    if isempty(noiseCovariance) || ~usePrewhitening
        return;
    end
    if ~isequal(size(noiseCovariance), [nChannels, nChannels])
        error('seal_STOUT:InvalidNoiseCovariance', ...
            'NoiseCovariance must be square with one row per sensor channel.');
    end
    noiseCovariance = double((noiseCovariance + noiseCovariance.') / 2);
    [factor, flag] = chol(noiseCovariance, 'lower');
    if flag ~= 0
        jitter = max(trace(noiseCovariance) / max(nChannels, 1), 1) * 1e-10;
        [factor, flag] = chol(noiseCovariance + jitter * eye(nChannels), 'lower');
    end
    if flag ~= 0
        error('seal_STOUT:InvalidNoiseCovariance', ...
            'NoiseCovariance must be positive definite for prewhitening.');
    end
    whiteningMatrix = factor \ eye(nChannels);
    DataWork = whiteningMatrix * Data;
    LWork = whiteningMatrix * L;
    info = struct('Applied', true, 'Matrix', whiteningMatrix);
end

function [Source, inverseOperator, locationWeights, lambda] = localWMNE(Data, L, nOrientations, regularizationParameter)
    nChannels = size(L, 1);
    nLocations = size(L, 2) / nOrientations;
    locationWeights = zeros(nLocations, 1);
    for location = 1:nLocations
        componentIndex = (location - 1) * nOrientations + (1:nOrientations);
        locationWeights(location) = 1 / max(norm(L(:, componentIndex), 'fro'), eps);
    end
    locationWeights = locationWeights / max(locationWeights);
    componentWeights = repelem(locationWeights, nOrientations);
    weightedLeadField = L .* componentWeights.';
    sensorSystem = weightedLeadField * weightedLeadField.';
    lambda = max(double(regularizationParameter) * ...
        trace(sensorSystem) / max(nChannels, 1), eps);
    sensorSystem = sensorSystem + lambda * eye(nChannels);
    [factor, flag] = chol((sensorSystem + sensorSystem.') / 2, 'lower');
    if flag ~= 0
        error('seal_STOUT:WMNEFactorizationFailure', ...
            'Unable to factor the wMNE sensor-space system.');
    end
    weightedTranspose = factor' \ (factor \ weightedLeadField);
    inverseOperator = weightedTranspose.' .* componentWeights;
    Source = inverseOperator * Data;
end

function magnitude = localLocationMagnitude(Source, nOrientations)
    nLocations = size(Source, 1) / nOrientations;
    blocks = reshape(Source, nOrientations, nLocations, []);
    magnitude = sqrt(reshape(sum(blocks .^ 2, 1), nLocations, []));
end

function [change, frequencies] = localSpectralChange(sourceMagnitude, baselineSamples, activeSamples, samplingRate)
    nTime = size(sourceMagnitude, 2);
    nFFT = 2 ^ nextpow2(max(numel(baselineSamples), numel(activeSamples)));
    nFFT = min(max(nFFT, 8), 2 ^ nextpow2(nTime));
    baselinePSD = localLogPSD(sourceMagnitude(:, baselineSamples), nFFT);
    activePSD = localLogPSD(sourceMagnitude(:, activeSamples), nFFT);
    change = activePSD - baselinePSD;
    frequencies = (0:(nFFT / 2)) * (double(samplingRate) / nFFT);
end

function spectrum = localLogPSD(signal, nFFT)
    signal = signal - mean(signal, 2);
    transform = fft(signal, nFFT, 2);
    power = abs(transform(:, 1:(nFFT / 2 + 1))) .^ 2 / max(size(signal, 2), 1);
    averagePower = mean(power, 2);
    spectrum = log(max(power ./ max(averagePower, eps), eps));
end

function reference = localReferenceSpectrum(referenceOption, frequencies)
    nFrequency = numel(frequencies);
    if isnumeric(referenceOption)
        reference = double(referenceOption(:).');
        if isempty(reference) || any(~isfinite(reference))
            error('seal_STOUT:InvalidReferenceSpectrum', ...
                'ReferenceSpectrum must be a non-empty finite vector.');
        end
        if numel(reference) ~= nFrequency
            sourceAxis = linspace(0, 1, numel(reference));
            targetAxis = linspace(0, 1, nFrequency);
            reference = interp1(sourceAxis, reference, targetAxis, 'linear', 'extrap');
        end
    else
        if ~strcmpi(char(referenceOption), 'motor')
            error('seal_STOUT:InvalidReferenceSpectrum', ...
                'ReferenceSpectrum must be a numeric ECoG profile or ''motor''.');
        end
        normalizedFrequency = frequencies / max(frequencies(end), eps);
        reference = zeros(1, nFrequency);
        reference(normalizedFrequency >= 0.10 & normalizedFrequency <= 0.30) = -1;
        reference(normalizedFrequency >= 0.60 & normalizedFrequency <= 0.80) = 1;
        reference = localSmoothVector(reference, 3);
    end
    reference = reference - mean(reference);
end

function output = localSmoothVector(input, halfWidth)
    kernel = ones(1, 2 * halfWidth + 1) / (2 * halfWidth + 1);
    output = conv(input, kernel, 'same');
end

function relationship = localRelationshipScore(spectralChange, referenceSpectrum)
    referenceNorm = norm(referenceSpectrum, 2);
    spectralNorm = sqrt(sum(spectralChange .^ 2, 2));
    numerator = spectralChange * referenceSpectrum.';
    relationship = 1 + numerator ./ max(spectralNorm * referenceNorm, eps);
    relationship = min(max(real(relationship), 0), 2);
end

function [Source, gain] = localApplyConstraint(initialSource, sourceMagnitude, xi, mode, alpha, exponentFloor, nOrientations)
    maxMagnitude = max(sourceMagnitude(:));
    magnitudeScale = max(maxMagnitude, eps);
    normalizedMagnitude = max(sourceMagnitude / magnitudeScale, double(exponentFloor));
    cocoMagnitude = (xi .^ 2) .* sourceMagnitude;
    excoMagnitude = magnitudeScale * (normalizedMagnitude .^ xi);
    switch mode
        case 'coco'
            targetMagnitude = cocoMagnitude;
        case 'exco'
            targetMagnitude = excoMagnitude;
        case 'hyco'
            targetMagnitude = cocoMagnitude / alpha + alpha * excoMagnitude;
    end
    gain = targetMagnitude ./ max(sourceMagnitude, magnitudeScale * double(exponentFloor));
    nLocations = size(sourceMagnitude, 1);
    sourceBlocks = reshape(initialSource, nOrientations, nLocations, []);
    Source = reshape(sourceBlocks .* reshape(gain, 1, nLocations, []), size(initialSource));
end

function value = localToLogical(value)
    if islogical(value)
        value = value(1);
    elseif isnumeric(value)
        value = value(1) ~= 0;
    else
        value = strcmpi(char(string(value)), 'true');
    end
end

function valid = localIsLogicalLike(value)
    valid = islogical(value) || (isnumeric(value) && isscalar(value)) || ischar(value) || isstring(value);
end
