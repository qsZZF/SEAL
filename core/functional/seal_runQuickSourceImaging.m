function out = seal_runQuickSourceImaging(packagePath, varargin)
%SEAL_RUNQUICKSOURCEIMAGING Run one-step source imaging from a packaged MAT file.
%   OUT = seal_runQuickSourceImaging(PACKAGEPATH) loads a .mat file that
%   contains EEG data, Cortex, Leadfield, and optional metadata, runs source
%   imaging with default settings, saves the source matrix, and opens a
%   cortex visualization.
%
%   Recognized top-level fields include:
%     data/Data/F/EEG, Cortex/cortex/BrainModel/Surface,
%     Leadfield/LeadField/L/Gain/Forward, srate/fs/SamplingRate,
%     NoiseCovariance/NoiseCov/noise_cov, Algorithm, Parameters.
%
%   Optional name-value pairs:
%     Algorithm      - default 'MNE'
%     Parameters     - struct of algorithm parameters
%     NoiseStrategy  - 'brainstorm' (default), 'identity', 'baseline',
%                      'data_cov', 'precomputed', or 'pass_through'
%     SaveDir        - output folder; default is the package folder
%     OpenFigure     - true/false; default true
%     ResultPrefix   - output filename prefix; default package basename
%     DepthWeightExp - MNE depth weighting exponent; default 0.5

if nargin < 1 || isempty(packagePath)
    error('SEAL:QuickSource:MissingFile', 'A MAT package path is required.');
end

packagePath = char(packagePath);
if ~isfile(packagePath)
    error('SEAL:QuickSource:FileNotFound', 'File not found: %s', packagePath);
end

rawPackage = load(packagePath, '-mat');
package = unwrapSingletonStruct(rawPackage);

opts = defaultQuickOptions();
opts = mergeStructs(opts, extractOptionsFromPackage(package));
opts = mergeStructs(opts, parseNameValueOptions(varargin{:}));

quickData = parseQuickPackage(package, packagePath);
if ~isempty(quickData.Algorithm) && ~isUserOption(varargin, 'Algorithm')
    opts.Algorithm = quickData.Algorithm;
end
if ~isempty(fieldnames(quickData.Parameters))
    opts.Parameters = mergeStructs(quickData.Parameters, opts.Parameters);
end

algorithm = canonicalAlgorithmName(opts.Algorithm);
parameters = fillAlgorithmDefaults(algorithm, opts.Parameters);
if strcmpi(algorithm, 'STOUT') && ~isempty(quickData.Srate) && ...
        (~isfield(parameters, 'SamplingRate') || isempty(parameters.SamplingRate))
    parameters.SamplingRate = quickData.Srate;
end

data = double(quickData.Data);
leadfield = double(quickData.Leadfield);
cortex = quickData.Cortex;

opts.SaveDir = char(string(opts.SaveDir));
if isempty(opts.SaveDir)
    opts.SaveDir = fileparts(packagePath);
end
if ~isfolder(opts.SaveDir)
    mkdir(opts.SaveDir);
end

[data, leadfield, quickData.NoiseCovariance] = alignDataAndLeadfield( ...
    data, leadfield, quickData.NoiseCovariance);

[nChannels, nTimes, nTrials, allData] = reshapeDataForRun(data);
[nOrientations, nSources] = inferOrientationCount(leadfield, cortex);

if size(leadfield, 1) ~= nChannels
    error('SEAL:QuickSource:DimensionMismatch', ...
        ['Data channels (%d) and leadfield rows (%d) do not match. ', ...
         'If your leadfield is transposed, save it as channels x sources.'], ...
        nChannels, size(leadfield, 1));
end
if nSources ~= size(cortex.Vertices, 1)
    error('SEAL:QuickSource:SourceMismatch', ...
        'Leadfield implies %d source locations, but Cortex has %d vertices.', ...
        nSources, size(cortex.Vertices, 1));
end

[W, R, whiteningInfo] = computeQuickWhitener(allData, opts.NoiseStrategy, ...
    quickData.NoiseCovariance);
leadfieldW = W * leadfield;

sourceEstimate = zeros(nSources, nTimes, nTrials);
algorithmInfo = cell(1, nTrials);
for iTrial = 1:nTrials
    if nTrials == 1
        trialData = data;
    else
        trialData = data(:, :, iTrial);
    end

    trialDataW = W * trialData;
    [rawSource, algorithmInfo{iTrial}] = runQuickAlgorithm( ...
        algorithm, trialDataW, leadfieldW, cortex, R, nOrientations, ...
        parameters, opts, quickData.DepthLeadfield);
    sourceEstimate(:, :, iTrial) = collapseOrientations(rawSource, nOrientations, nSources);
end

if nTrials == 1
    sourceEstimate = sourceEstimate(:, :, 1);
    algorithmInfo = algorithmInfo{1};
end

displayVector = buildDisplayVector(sourceEstimate);
quickSourceInfo = buildQuickSourceInfo(packagePath, algorithm, parameters, ...
    quickData, whiteningInfo, nOrientations, size(sourceEstimate));
resultPath = saveQuickResult(opts.SaveDir, opts.ResultPrefix, packagePath, ...
    algorithm, sourceEstimate, displayVector, cortex, quickSourceInfo, ...
    algorithmInfo);

figureHandles = [];
if toLogical(opts.OpenFigure)
    figureHandles = PlotSource(displayVector, cortex, ...
        'ValueMode', 'brainstorm', 'Colormap', 'brainstorm');
    if isfield(figureHandles, 'h') && isvalid(figureHandles.h)
        set(figureHandles.h, 'Name', sprintf('SEAL Quick Source - %s', algorithm), ...
            'NumberTitle', 'off');
    end
end

out = struct();
out.SourceEstimate = sourceEstimate;
out.DisplayVector = displayVector;
out.ResultPath = resultPath;
out.Algorithm = algorithm;
out.Parameters = parameters;
out.WhiteningInfo = whiteningInfo;
out.Package = quickSourceInfo.Package;
out.FigureHandles = figureHandles;
end

function opts = defaultQuickOptions()
opts = struct();
opts.Algorithm = 'MNE';
opts.Parameters = struct();
opts.NoiseStrategy = 'brainstorm';
opts.SaveDir = '';
opts.OpenFigure = true;
opts.ResultPrefix = '';
opts.DepthWeightExp = 0.5;
end

function opts = parseNameValueOptions(varargin)
opts = struct();
if isempty(varargin)
    return;
end
if mod(numel(varargin), 2) ~= 0
    error('SEAL:QuickSource:InvalidNameValue', ...
        'Options must be provided as name-value pairs.');
end
for ii = 1:2:numel(varargin)
    name = canonicalOptionName(varargin{ii});
    opts.(name) = varargin{ii + 1};
end
end

function tf = isUserOption(args, wantedName)
tf = false;
if mod(numel(args), 2) ~= 0
    return;
end
wantedName = canonicalOptionName(wantedName);
for ii = 1:2:numel(args)
    if strcmp(canonicalOptionName(args{ii}), wantedName)
        tf = true;
        return;
    end
end
end

function name = canonicalOptionName(name)
name = char(string(name));
key = lower(strrep(name, '_', ''));
switch key
    case {'algorithm', 'algorithmname', 'alg'}
        name = 'Algorithm';
    case {'parameters', 'params', 'algorithmparameters', 'algparams'}
        name = 'Parameters';
    case {'noisestrategy', 'whitening', 'noisewhitening'}
        name = 'NoiseStrategy';
    case {'savedir', 'outputdir', 'resultdir', 'resultfolder'}
        name = 'SaveDir';
    case {'openfigure', 'showfigure', 'plot', 'visualize'}
        name = 'OpenFigure';
    case {'resultprefix', 'prefix'}
        name = 'ResultPrefix';
    case {'depthweightexp', 'depthweighting', 'depthweightingexponent'}
        name = 'DepthWeightExp';
    otherwise
        name = matlab.lang.makeValidName(name);
end
end

function algorithm = canonicalAlgorithmName(algorithm)
algorithm = char(string(algorithm));
try
    registry = seal_getAlgorithmRegistry();
    names = fieldnames(registry);
    idx = find(strcmpi(names, algorithm), 1);
    if ~isempty(idx)
        algorithm = names{idx};
    end
catch
end
end

function opts = extractOptionsFromPackage(package)
opts = struct();
[optionStruct, hasOptions] = getFieldByAliases(package, ...
    {'QuickSourceOptions', 'quickSourceOptions', 'QuickOptions', 'Options', 'options'});
if hasOptions && isstruct(optionStruct)
    opts = mergeStructs(opts, normalizeOptionStruct(optionStruct));
end

[algorithm, hasAlgorithm] = getFieldByAliases(package, ...
    {'Algorithm', 'algorithm', 'AlgorithmName', 'algorithmName'});
if hasAlgorithm && (ischar(algorithm) || isstring(algorithm))
    opts.Algorithm = char(string(algorithm));
end

[parameters, hasParameters] = getFieldByAliases(package, ...
    {'Parameters', 'parameters', 'Params', 'params', 'AlgorithmParameters'});
if hasParameters && isstruct(parameters)
    opts.Parameters = parameters;
end
end

function opts = normalizeOptionStruct(inputStruct)
opts = struct();
names = fieldnames(inputStruct);
for ii = 1:numel(names)
    opts.(canonicalOptionName(names{ii})) = inputStruct.(names{ii});
end
end

function merged = mergeStructs(base, override)
merged = base;
if isempty(override) || ~isstruct(override)
    return;
end
names = fieldnames(override);
for ii = 1:numel(names)
    name = names{ii};
    if isfield(merged, name) && isstruct(merged.(name)) && isstruct(override.(name))
        merged.(name) = mergeStructs(merged.(name), override.(name));
    else
        merged.(name) = override.(name);
    end
end
end

function package = unwrapSingletonStruct(package)
if ~isstruct(package)
    return;
end
fields = fieldnames(package);
if numel(fields) == 1 && isstruct(package.(fields{1}))
    package = package.(fields{1});
end
end

function quickData = parseQuickPackage(package, packagePath)
quickData = struct();
quickData.Data = [];
quickData.Cortex = [];
quickData.Leadfield = [];
quickData.Srate = [];
quickData.Time = [];
quickData.NoiseCovariance = [];
quickData.DepthLeadfield = [];
quickData.Algorithm = '';
quickData.Parameters = struct();
quickData.PackagePath = packagePath;
quickData.FieldMap = struct();

[dataObj, hasData, dataField] = getFieldByAliases(package, ...
    {'data', 'Data', 'DataFile', 'EEG', 'F', 'signal', 'Signal', ...
     'matrix', 'X', 'ERP', 'ERPset', 'EEGData'});
if ~hasData
    error('SEAL:QuickSource:MissingData', ...
        'No EEG data field found. Expected data/Data/F/EEG/signal/X.');
end
[quickData.Data, dataMeta] = extractDataMatrix(dataObj);
quickData.FieldMap.Data = dataField;

[cortexObj, hasCortex, cortexField] = getFieldByAliases(package, ...
    {'Cortex', 'cortex', 'CortexModel', 'BrainModel', 'Surface', ...
     'surface', 'HeadModel', 'headModel'});
if ~hasCortex
    error('SEAL:QuickSource:MissingCortex', ...
        'No Cortex field found. Expected Cortex/cortex/BrainModel/Surface.');
end
quickData.Cortex = normalizeCortex(cortexObj);
quickData.FieldMap.Cortex = cortexField;

[leadfieldObj, hasLeadfield, leadfieldField] = getFieldByAliases(package, ...
    {'Leadfield', 'leadfield', 'LeadField', 'leadField', 'L', ...
     'Gain', 'Gain_Matrix', 'GainMatrix', 'Forward', 'G', 'LF'});
if ~hasLeadfield
    error('SEAL:QuickSource:MissingLeadfield', ...
        'No Leadfield field found. Expected Leadfield/L/Gain/Forward.');
end
[quickData.Leadfield, leadfieldMeta] = extractLeadfieldMatrix(leadfieldObj);
quickData.DepthLeadfield = leadfieldMeta.DepthLeadfield;
quickData.FieldMap.Leadfield = leadfieldField;

[srate, hasSrate, srateField] = getFieldByAliases(package, ...
    {'srate', 'Srate', 'SamplingRate', 'SampleRate', 'fs', 'sfreq'});
if ~hasSrate && isfield(dataMeta, 'Srate')
    srate = dataMeta.Srate;
    hasSrate = true;
    srateField = sprintf('%s.%s', dataField, dataMeta.SrateField);
end
if hasSrate
    quickData.Srate = double(srate(1));
    quickData.FieldMap.Srate = srateField;
end

[timeValue, hasTime, timeField] = getFieldByAliases(package, ...
    {'Time', 'time', 'Times', 'times', 'TimeVector', 't'});
if hasTime
    quickData.Time = timeValue;
    quickData.FieldMap.Time = timeField;
    if isempty(quickData.Srate) && isnumeric(timeValue) && numel(timeValue) > 1
        dt = double(timeValue(2) - timeValue(1));
        if isfinite(dt) && dt ~= 0
            quickData.Srate = 1 / dt;
            quickData.FieldMap.Srate = timeField;
        end
    end
end

[noiseCov, hasNoiseCov, noiseField] = getFieldByAliases(package, ...
    {'NoiseCovariance', 'noiseCovariance', 'NoiseCov', 'noise_cov', ...
     'NoiseCovMat', 'Cov_n', 'CovNoise'});
if ~hasNoiseCov && isfield(dataMeta, 'NoiseCovariance')
    noiseCov = dataMeta.NoiseCovariance;
    hasNoiseCov = true;
    noiseField = sprintf('%s.NoiseCovariance', dataField);
end
if hasNoiseCov
    quickData.NoiseCovariance = extractNoiseCovariance(noiseCov);
    quickData.FieldMap.NoiseCovariance = noiseField;
end

[algorithm, hasAlgorithm] = getFieldByAliases(package, ...
    {'Algorithm', 'algorithm', 'AlgorithmName', 'algorithmName'});
if hasAlgorithm && (ischar(algorithm) || isstring(algorithm))
    quickData.Algorithm = char(string(algorithm));
end

[parameters, hasParameters] = getFieldByAliases(package, ...
    {'Parameters', 'parameters', 'Params', 'params', 'AlgorithmParameters'});
if hasParameters && isstruct(parameters)
    quickData.Parameters = parameters;
end
end

function [value, found, matchedName] = getFieldByAliases(s, aliases)
value = [];
found = false;
matchedName = '';
if ~isstruct(s)
    return;
end
names = fieldnames(s);
lowerNames = cellfun(@lower, names, 'UniformOutput', false);
for ii = 1:numel(aliases)
    idx = find(strcmp(lowerNames, lower(char(aliases{ii}))), 1);
    if ~isempty(idx)
        matchedName = names{idx};
        value = s.(matchedName);
        found = true;
        return;
    end
end

for ii = 1:numel(aliases)
    alias = lower(char(aliases{ii}));
    if numel(alias) < 3
        continue;
    end
    idx = find(cellfun(@(name) ~isempty(strfind(name, alias)), lowerNames), 1);
    if ~isempty(idx)
        matchedName = names{idx};
        value = s.(matchedName);
        found = true;
        return;
    end
end
end

function [data, meta] = extractDataMatrix(dataObj)
meta = struct();
dataObj = unwrapSingletonStruct(dataObj);
if isnumeric(dataObj) && (ismatrix(dataObj) || ndims(dataObj) == 3)
    data = dataObj;
    return;
end
if ~isstruct(dataObj)
    error('SEAL:QuickSource:InvalidData', ...
        'The Data field must be a numeric matrix/3-D array or a struct containing one.');
end

[srate, hasSrate, srateField] = getFieldByAliases(dataObj, ...
    {'srate', 'Srate', 'SamplingRate', 'SampleRate', 'fs', 'sfreq'});
if hasSrate
    meta.Srate = srate;
    meta.SrateField = srateField;
end

[noiseCov, hasNoiseCov] = getFieldByAliases(dataObj, ...
    {'NoiseCovariance', 'noiseCovariance', 'NoiseCov', 'NoiseCovMat'});
if hasNoiseCov
    meta.NoiseCovariance = extractNoiseCovariance(noiseCov);
end

[data, hasData] = getFieldByAliases(dataObj, ...
    {'data', 'Data', 'F', 'EEG', 'signal', 'Signal', 'matrix', 'X', 'ERP'});
if hasData && isnumeric(data) && (ismatrix(data) || ndims(data) == 3)
    return;
end

fields = fieldnames(dataObj);
for ii = 1:numel(fields)
    candidate = dataObj.(fields{ii});
    if isnumeric(candidate) && (ismatrix(candidate) || ndims(candidate) == 3)
        data = candidate;
        return;
    end
end

error('SEAL:QuickSource:InvalidData', ...
    'Could not identify a numeric EEG matrix inside the Data field.');
end

function cortex = normalizeCortex(cortexObj)
cortexObj = unwrapSingletonStruct(cortexObj);
if ~isstruct(cortexObj)
    error('SEAL:QuickSource:InvalidCortex', ...
        'Cortex must be a struct containing Vertices and Faces.');
end

[vertices, hasVertices] = getFieldByAliases(cortexObj, ...
    {'Vertices', 'vertices', 'Vertex', 'vert', 'Vert', 'pos', 'Position', 'nodes'});
[faces, hasFaces] = getFieldByAliases(cortexObj, ...
    {'Faces', 'faces', 'Face', 'face', 'tri', 'Triangles', 'polygons'});
if ~hasVertices || ~hasFaces
    error('SEAL:QuickSource:InvalidCortex', ...
        'Cortex must contain Vertices and Faces (or recognized aliases).');
end

cortex = cortexObj;
vertices = double(vertices);
faces = double(faces);
if size(vertices, 2) ~= 3 && size(vertices, 1) == 3
    vertices = vertices';
end
if size(faces, 2) ~= 3 && size(faces, 1) == 3
    faces = faces';
end
if min(faces(:)) == 0
    faces = faces + 1;
end
cortex.Vertices = vertices;
cortex.Faces = faces;

[vertConn, hasVertConn] = getFieldByAliases(cortexObj, ...
    {'VertConn', 'vertConn', 'VertexConnectivity', 'Connectivity'});
if hasVertConn && ~isempty(vertConn)
    cortex.VertConn = vertConn;
else
    cortex.VertConn = buildVertexConnectivity(faces, size(vertices, 1));
end
end

function conn = buildVertexConnectivity(faces, nVertices)
faces = double(faces);
i = [faces(:, 1); faces(:, 2); faces(:, 3); faces(:, 2); faces(:, 3); faces(:, 1)];
j = [faces(:, 2); faces(:, 3); faces(:, 1); faces(:, 1); faces(:, 2); faces(:, 3)];
conn = sparse(i, j, true(size(i)), nVertices, nVertices);
conn = conn | conn';
end

function [leadfield, meta] = extractLeadfieldMatrix(leadfieldObj)
meta = struct('DepthLeadfield', []);
leadfieldObj = unwrapSingletonStruct(leadfieldObj);
if isnumeric(leadfieldObj) && ismatrix(leadfieldObj)
    leadfield = leadfieldObj;
    return;
end
if ~isstruct(leadfieldObj)
    error('SEAL:QuickSource:InvalidLeadfield', ...
        'Leadfield must be a numeric matrix or a struct containing one.');
end

[depthLeadfield, hasDepthLeadfield] = getFieldByAliases(leadfieldObj, ...
    {'DepthLeadfield', 'L_free_avgref', 'Gain_free_avgref', ...
     'Leadfield_free_avgref', 'L_free', 'Gain_raw'});
if hasDepthLeadfield && isnumeric(depthLeadfield) && ismatrix(depthLeadfield)
    meta.DepthLeadfield = double(depthLeadfield);
end

[leadfield, hasLeadfield] = getFieldByAliases(leadfieldObj, ...
    {'Leadfield', 'leadfield', 'LeadField', 'L', 'Gain', 'Gain_Matrix', ...
     'GainMatrix', 'Forward', 'G', 'LF', 'matrix'});
if hasLeadfield && isnumeric(leadfield) && ismatrix(leadfield)
    return;
end

fields = fieldnames(leadfieldObj);
for ii = 1:numel(fields)
    candidate = leadfieldObj.(fields{ii});
    if isnumeric(candidate) && ismatrix(candidate)
        leadfield = candidate;
        return;
    end
end

error('SEAL:QuickSource:InvalidLeadfield', ...
    'Could not identify a numeric leadfield matrix.');
end

function cov = extractNoiseCovariance(obj)
cov = [];
obj = unwrapSingletonStruct(obj);
if isnumeric(obj) && ismatrix(obj)
    cov = double(obj);
    return;
end
if isstruct(obj)
    [value, found] = getFieldByAliases(obj, ...
        {'NoiseCovariance', 'NoiseCov', 'noiseCovariance', 'NoiseCovMat', 'Cov'});
    if found
        cov = extractNoiseCovariance(value);
        return;
    end
    fields = fieldnames(obj);
    for ii = 1:numel(fields)
        value = obj.(fields{ii});
        if isnumeric(value) && ismatrix(value) && size(value, 1) == size(value, 2)
            cov = double(value);
            return;
        end
    end
end
end

function [data, leadfield, noiseCov] = alignDataAndLeadfield(data, leadfield, noiseCov)
nChannels = size(data, 1);
if size(leadfield, 1) ~= nChannels && size(leadfield, 2) == nChannels
    leadfield = leadfield';
end
if ~isempty(noiseCov)
    if size(noiseCov, 1) ~= size(noiseCov, 2)
        warning('SEAL:QuickSource:InvalidNoiseCovariance', ...
            'Noise covariance is not square and will be ignored.');
        noiseCov = [];
    elseif size(noiseCov, 1) ~= nChannels
        warning('SEAL:QuickSource:NoiseCovarianceMismatch', ...
            'Noise covariance size does not match data channels and will be ignored.');
        noiseCov = [];
    end
end
end

function [nChannels, nTimes, nTrials, allData] = reshapeDataForRun(data)
if ndims(data) == 3
    [nChannels, nTimes, nTrials] = size(data);
    allData = reshape(data, nChannels, nTimes * nTrials);
elseif ismatrix(data)
    [nChannels, nTimes] = size(data);
    nTrials = 1;
    allData = data;
else
    error('SEAL:QuickSource:InvalidDataDimension', ...
        'Data must be channels x time or channels x time x trials.');
end
end

function [nOrientations, nSources] = inferOrientationCount(leadfield, cortex)
nVertices = size(cortex.Vertices, 1);
nColumns = size(leadfield, 2);
if nColumns == nVertices
    nOrientations = 1;
    nSources = nVertices;
elseif nColumns == 3 * nVertices
    nOrientations = 3;
    nSources = nVertices;
else
    error('SEAL:QuickSource:LeadfieldCortexMismatch', ...
        'Leadfield columns (%d) must equal Cortex vertices (%d) or 3x vertices.', ...
        nColumns, nVertices);
end
end

function params = fillAlgorithmDefaults(algorithm, params)
if nargin < 2 || isempty(params) || ~isstruct(params)
    params = struct();
end
try
    registry = seal_getAlgorithmRegistry();
catch
    registry = struct();
end
if ~isfield(registry, algorithm) || ~isfield(registry.(algorithm), 'params')
    return;
end
specs = registry.(algorithm).params;
for ii = 1:numel(specs)
    fieldName = specs{ii}{1};
    defaultValue = specs{ii}{4};
    if ~isfield(params, fieldName) || isempty(params.(fieldName))
        params.(fieldName) = defaultValue;
    end
end
end

function [W, RWhite, info] = computeQuickWhitener(X, strategy, noiseCov)
nChannels = size(X, 1);
info = struct('strategy', strategy, 'baselineSamples', 0, ...
    'condNumber', NaN, 'rankDeficient', false, 'source', '');
strategy = lower(char(string(strategy)));

switch strategy
    case {'brainstorm', 'brainstorm_reg'}
        if ~isempty(noiseCov)
            [W, ~, rankNoise, condNumber] = brainstormRegWhitener(noiseCov, 0.1);
            RWhite = eye(nChannels);
            info.source = 'package_noise_cov_brainstorm_reg';
            info.noiseReg = 0.1;
            info.rankNoise = rankNoise;
            info.rankDeficient = rankNoise < nChannels;
            info.condNumber = condNumber;
            return;
        end
        [W, RWhite, info] = baselineWhitener(X, info, 'brainstorm_fallback_baseline');

    case 'precomputed'
        if isempty(noiseCov)
            W = eye(nChannels);
            RWhite = eye(nChannels);
            info.source = 'fallback_identity_no_precomputed_cov';
            return;
        end
        [W, RWhite, info] = covarianceWhitener(noiseCov, info);
        info.source = 'package_noise_cov';

    case 'baseline'
        [W, RWhite, info] = baselineWhitener(X, info, 'baseline');

    case 'data_cov'
        Xc = X - mean(X, 2);
        covData = (Xc * Xc') / max(size(X, 2) - 1, 1);
        [W, RWhite, info] = covarianceWhitener(covData, info);
        info.source = 'data_cov';

    case {'identity', 'none'}
        W = eye(nChannels);
        RWhite = eye(nChannels);
        info.source = 'identity';

    case 'pass_through'
        W = eye(nChannels);
        RWhite = [];
        info.source = 'pass_through';

    otherwise
        warning('SEAL:QuickSource:UnknownNoiseStrategy', ...
            'Unknown NoiseStrategy "%s"; using identity.', strategy);
        W = eye(nChannels);
        RWhite = eye(nChannels);
        info.strategy = 'identity';
        info.source = 'fallback_identity_unknown_strategy';
end
end

function [W, RWhite, info] = baselineWhitener(X, info, sourceName)
nChannels = size(X, 1);
baselineLen = min(90, size(X, 2));
info.baselineSamples = baselineLen;
info.source = sourceName;
if baselineLen < 2
    W = eye(nChannels);
    RWhite = eye(nChannels);
    info.source = 'fallback_identity_short_baseline';
    return;
end
if baselineLen < nChannels
    info.rankDeficient = true;
end
baseline = X(:, 1:baselineLen);
baseline = baseline - mean(baseline, 2);
covNoise = (baseline * baseline') / max(baselineLen - 1, 1);
[W, RWhite, info] = covarianceWhitener(covNoise, info);
end

function [W, RWhite, info] = covarianceWhitener(covNoise, info)
nChannels = size(covNoise, 1);
covNoise = (double(covNoise) + double(covNoise)') / 2;
[U, S] = svd(covNoise, 'econ');
singularValues = diag(S);
if isempty(singularValues) || ~isfinite(singularValues(1)) || singularValues(1) <= 0
    W = eye(nChannels);
    RWhite = eye(nChannels);
    info.source = 'fallback_identity_invalid_cov';
    return;
end
sMax = singularValues(1);
rankTol = max(size(covNoise)) * eps(sMax);
rankNoise = sum(singularValues > rankTol);
lowerBound = sMax * 1e-6;
singularValues(singularValues < lowerBound) = lowerBound;
W = diag(1 ./ sqrt(singularValues)) * U';
RWhite = eye(nChannels);
info.rankNoise = rankNoise;
info.rankDeficient = info.rankDeficient || rankNoise < nChannels;
info.condNumber = sMax / singularValues(end);
end

function [W, covReg, rankNoise, condNumber] = brainstormRegWhitener(covNoise, noiseReg)
covNoise = (double(covNoise) + double(covNoise)') / 2;
[U, S] = svd(covNoise, 'econ');
eigVals = diag(S);
if isempty(eigVals) || eigVals(1) <= 0
    error('SEAL:QuickSource:InvalidNoiseCovariance', ...
        'Noise covariance has no positive eigenvalue.');
end
noiseStd = sqrt(eigVals);
tol = length(noiseStd) * eps(single(noiseStd(1)));
rankNoise = sum(noiseStd > tol);
UKeep = U(:, 1:rankNoise);
noiseStd = noiseStd(1:rankNoise);
ridgeFactor = mean(eigVals) * noiseReg;
covReg = UKeep * diag(noiseStd .^ 2) * UKeep' + ridgeFactor * eye(size(covNoise, 1));
W = UKeep * diag(1 ./ sqrt(noiseStd .^ 2 + ridgeFactor)) * UKeep';
regEig = eig((covReg + covReg') / 2);
posEig = regEig(regEig > 0);
condNumber = max(posEig) / min(posEig);
end

function [source, algorithmInfo] = runQuickAlgorithm(algorithm, X, L, cortex, R, ...
    nOrientations, params, opts, depthLeadfield)
algorithmInfo = struct();
algorithm = char(string(algorithm));
RForAlgorithm = R;
if isempty(RForAlgorithm)
    RForAlgorithm = eye(size(L, 1));
end

switch algorithm
    case 'MNE'
        if ~isfield(params, 'RegularizationParameter')
            params.RegularizationParameter = 1 / 3;
        end
        if ~isfield(params, 'sourceCovariance') || isempty(params.sourceCovariance)
            params.sourceCovariance = 'brainstorm_depth';
        end
        [source, algorithmInfo] = seal_MNE(X, L, ...
            'RegularizationParameter', params.RegularizationParameter, ...
            'SourceCovariance', params.sourceCovariance, ...
            'DepthWeightExp', opts.DepthWeightExp, ...
            'DepthWeightLimit', 10, ...
            'DepthLeadfield', depthLeadfield, ...
            'NoiseCovariance', RForAlgorithm, ...
            'NumOrientations', nOrientations);

    case 'sLORETA'
        [source, algorithmInfo] = seal_sLORETA(X, L, ...
            'RegularizationParameter', params.RegularizationParameter, ...
            'NoiseCovariance', RForAlgorithm, ...
            'NumOrientations', nOrientations);

    case 'dSPM'
        [source, algorithmInfo] = seal_dSPM(X, L, ...
            'RegularizationParameter', params.RegularizationParameter, ...
            'NoiseCovariance', RForAlgorithm, ...
            'NumOrientations', nOrientations);

    case 'LCMV'
        [source, algorithmInfo] = seal_LCMV(X, L, ...
            'Type', params.Type, ...
            'RegularizationParameter', params.RegularizationParameter, ...
            'ScalarBeamformer', toLogical(params.ScalarBeamformer), ...
            'NoiseCovariance', RForAlgorithm, ...
            'NumOrientations', nOrientations);

    otherwise
        [source, algorithmInfo] = runAlgorithmByConvention( ...
            algorithm, X, L, cortex, RForAlgorithm, nOrientations, params);
end
end

function [source, algorithmInfo] = runAlgorithmByConvention(algorithm, X, L, cortex, ...
    R, nOrientations, params)
functionName = algorithm;
if exist(functionName, 'file') ~= 2
    functionName = ['seal_' algorithm];
end
if exist(functionName, 'file') ~= 2
    error('SEAL:QuickSource:UnknownAlgorithm', ...
        'Algorithm "%s" was not found on the MATLAB path.', algorithm);
end

params.NoiseCovariance = R;
params.NumOrientations = nOrientations;
paramCell = structToNameValue(params);
nInputs = abs(nargin(functionName));
needsVertConn = any(strcmpi(functionName, {'seal_dMAP', 'seal_dSPN', 'seal_TS_Champagne'}));
needsCortex = any(strcmpi(functionName, ...
    {'seal_BESTIES', 'seal_SmoothChampagne', 'seal_IRES', 'seal_VSSI_Lp', 'seal_VSSI_GGD'}));
if needsVertConn
    if ~isfield(cortex, 'VertConn') || isempty(cortex.VertConn)
        error('SEAL:QuickSource:MissingVertConn', ...
            'Algorithm "%s" requires Cortex.VertConn.', algorithm);
    end
    thirdArg = cortex.VertConn;
elseif needsCortex || nInputs >= 4
    thirdArg = cortex;
end

if needsVertConn || needsCortex || nInputs >= 4
    [source, algorithmInfo] = feval(functionName, X, L, thirdArg, paramCell{:});
else
    [source, algorithmInfo] = feval(functionName, X, L, paramCell{:});
end
end

function cells = structToNameValue(s)
skip = {'Algorithm', 'AlgorithmName', 'name', 'ResultFolderPath', 'srate'};
names = fieldnames(s);
names = setdiff(names, skip, 'stable');
cells = cell(1, 2 * numel(names));
for ii = 1:numel(names)
    cells{2 * ii - 1} = names{ii};
    cells{2 * ii} = s.(names{ii});
end
end

function value = toLogical(value)
if islogical(value)
    return;
end
if isnumeric(value)
    value = value ~= 0;
    return;
end
value = strcmpi(char(string(value)), 'true');
end

function result = collapseOrientations(source, nOrientations, nSources)
if size(source, 1) == nSources
    result = source;
    return;
end
if nOrientations == 1
    error('SEAL:QuickSource:UnexpectedSourceRows', ...
        'Algorithm returned %d source rows; expected %d.', size(source, 1), nSources);
end
if size(source, 1) ~= nOrientations * nSources
    error('SEAL:QuickSource:UnexpectedSourceRows', ...
        'Algorithm returned %d source rows; expected %d or %d.', ...
        size(source, 1), nSources, nOrientations * nSources);
end
result = reshape(source, nOrientations, [], size(source, 2));
result = squeeze(sqrt(sum(result .^ 2, 1)));
end

function displayVector = buildDisplayVector(sourceEstimate)
if ndims(sourceEstimate) == 2
    flattened = sourceEstimate;
else
    flattened = reshape(sourceEstimate, size(sourceEstimate, 1), []);
end
displayVector = sqrt(mean(abs(flattened) .^ 2, 2));
end

function info = buildQuickSourceInfo(packagePath, algorithm, parameters, quickData, ...
    whiteningInfo, nOrientations, sourceSize)
info = struct();
info.Package = packagePath;
info.CreatedDate = datetime('now');
info.Algorithm = algorithm;
info.Parameters = parameters;
info.NoiseWhitening = whiteningInfo;
info.NumOrientations = nOrientations;
info.SourceSize = sourceSize;
info.Srate = quickData.Srate;
info.FieldMap = quickData.FieldMap;
info.DataSize = size(quickData.Data);
info.LeadfieldSize = size(quickData.Leadfield);
info.CortexVertexCount = size(quickData.Cortex.Vertices, 1);
end

function resultPath = saveQuickResult(saveDir, resultPrefix, packagePath, algorithm, ...
    sourceEstimate, displayVector, cortex, quickSourceInfo, algorithmInfo)
resultPrefix = char(string(resultPrefix));
if isempty(resultPrefix)
    [~, baseName] = fileparts(packagePath);
else
    baseName = resultPrefix;
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS_FFF');
safeAlgorithm = regexprep(algorithm, '[^\w\-]', '_');
fileName = sprintf('%s_%s_quick_source_%s.mat', baseName, safeAlgorithm, timestamp);
resultPath = fullfile(saveDir, fileName);

Cortex = cortex; %#ok<NASGU>
QuickSourceInfo = quickSourceInfo; %#ok<NASGU>
AlgorithmInfo = algorithmInfo; %#ok<NASGU>
try
    save(resultPath, 'sourceEstimate', 'displayVector', 'Cortex', ...
        'QuickSourceInfo', 'AlgorithmInfo', '-mat');
catch ME
    if contains(ME.message, '2GB') || contains(lower(ME.message), 'large')
        save(resultPath, 'sourceEstimate', 'displayVector', 'Cortex', ...
            'QuickSourceInfo', 'AlgorithmInfo', '-v7.3');
    else
        rethrow(ME);
    end
end
end
