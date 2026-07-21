function [Source, Parameters] = seal_VSSI_Lp(Data, L, sourceSpaceInfo, varargin)
%SEAL_VSSI_LP Variation sparse source imaging with an Lp edge penalty.
%
%   [Source, Parameters] = seal_VSSI_Lp(Data, L, sourceSpaceInfo, ...)
%   implements VSSI-Lp (Peng et al., 2024):
%
%       min_S ||Data - L*S||_F^2 + lambda * ||V*S||_p^p,  0 < p <= 1,
%
%   where V is the source-space incidence matrix formed from cortical
%   connectivity.  The non-convex edge-domain Lp proximal is solved with
%   generalized soft thresholding (GST) inside an ADMM iteration.
%
%   sourceSpaceInfo may be a cortex struct containing .VertConn or
%   .EdgeOperator, or a location-by-location VertConn matrix directly.
%   Fixed orientation (NumOrientations=1) follows the paper formulation.
%   For free orientations, SEAL expands V per component and applies the
%   same Lp variation penalty componentwise.
%
%   Name-value options:
%       'RegularizationParameter'  lambda, the direct edge penalty weight (0.05)
%       'P'                        Lp exponent in (0, 1] (0.6)
%       'Rho'                      ADMM penalty parameter (1)
%       'MaxIterations'            Maximum ADMM iterations (500)
%       'Tolerance'                Relative ADMM tolerance (1e-4)
%       'AbsoluteTolerance'        Absolute ADMM tolerance (1e-6)
%       'ProgressInterval'         ADMM iterations between progress messages (25)
%       'AdaptiveRho'              Residual-balanced rho updates (true)
%       'GSTIterations'            GST fixed-point iterations (3)
%       'SourceRidge'              Relative numerical ridge for V''*V (1e-8)
%       'NoiseCovariance'          Optional sensor noise covariance ([])
%       'PrewhitenNoise'           Prewhiten if covariance is supplied (true)
%       'NormalizeProblem'         Normalize data/leadfield solver coordinates (true)
%       'NumOrientations'          Components per location: 1 or 3 (1)
%       'Verbose'                  Print concise ADMM diagnostics (true)
%       'CancelCheck'              Optional zero-input cancellation function ([])

    parser = inputParser;
    parser.FunctionName = mfilename;
    parser.CaseSensitive = false;
    addRequired(parser, 'Data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addRequired(parser, 'L', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addRequired(parser, 'sourceSpaceInfo', @(x) isstruct(x) || isnumeric(x) || islogical(x));
    addParameter(parser, 'RegularizationParameter', 0.05, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(parser, 'P', 0.6, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1);
    addParameter(parser, 'Rho', 1, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'MaxIterations', 500, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'Tolerance', 1e-4, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'AbsoluteTolerance', 1e-6, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'ProgressInterval', 25, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'AdaptiveRho', true, @localIsLogicalLike);
    addParameter(parser, 'RhoUpdateInterval', 10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'RhoBalance', 10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
    addParameter(parser, 'RhoScale', 2, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
    addParameter(parser, 'GSTIterations', 3, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'SourceRidge', 1e-8, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(parser, 'PrewhitenNoise', true, @localIsLogicalLike);
    addParameter(parser, 'NormalizeProblem', true, @localIsLogicalLike);
    addParameter(parser, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && any(x == [1, 3]));
    addParameter(parser, 'Verbose', true, @localIsLogicalLike);
    addParameter(parser, 'CancelCheck', [], @(x) isempty(x) || isa(x, 'function_handle'));
    parse(parser, Data, L, sourceSpaceInfo, varargin{:});
    opt = parser.Results;

    Data = double(Data);
    L = double(L);
    if any(~isfinite(Data(:))) || any(~isfinite(L(:)))
        error('seal_VSSI_Lp:NonFiniteInput', 'Data and L must contain only finite values.');
    end
    if size(Data, 1) ~= size(L, 1)
        error('seal_VSSI_Lp:DimensionMismatch', ...
            'Data and L must have the same number of sensor rows.');
    end

    nTime = size(Data, 2);
    nComponents = size(L, 2);
    nOrientations = round(double(opt.NumOrientations));
    if mod(nComponents, nOrientations) ~= 0
        error('seal_VSSI_Lp:OrientationMismatch', ...
            'Lead-field columns must be divisible by NumOrientations.');
    end
    nLocations = nComponents / nOrientations;

    verbose = localToLogical(opt.Verbose);
    if verbose
        initializationTimer = tic;
        fprintf('Preparing VSSI-Lp: whitening data and lead field...\n');
        drawnow limitrate;
    end
    [DataWork, LWork, whiteningInfo] = localPrewhiten(Data, L, ...
        opt.NoiseCovariance, localToLogical(opt.PrewhitenNoise));
    localCheckCancel(opt.CancelCheck);
    [DataSolve, LSolve, coordinateInfo] = localNormalizeSolveCoordinates( ...
        DataWork, LWork, nLocations, nOrientations, localToLogical(opt.NormalizeProblem));
    if verbose && coordinateInfo.Applied
        fprintf(['VSSI-Lp: normalized solver coordinates (data RMS %.3e, ', ...
            'lead-field block scale %.3e, source scale %.3e).\n'], ...
            coordinateInfo.DataScale, coordinateInfo.LeadfieldScale, coordinateInfo.SourceScale);
        drawnow limitrate;
    end
    if verbose
        fprintf('VSSI-Lp: constructing the cortical edge operator...\n');
        drawnow limitrate;
    end
    [V, edgeInfo] = localResolveEdgeOperator(sourceSpaceInfo, nLocations, nOrientations);
    localCheckCancel(opt.CancelCheck);
    nEdgeComponents = size(V, 1);
    if nEdgeComponents == 0
        error('seal_VSSI_Lp:EmptyEdgeOperator', ...
            'The source-space graph must contain at least one edge.');
    end

    lambda = double(opt.RegularizationParameter);
    pNorm = double(opt.P);
    rho = double(opt.Rho);
    maxIterations = round(double(opt.MaxIterations));
    relTolerance = double(opt.Tolerance);
    absTolerance = double(opt.AbsoluteTolerance);
    adaptiveRho = localToLogical(opt.AdaptiveRho);
    progressInterval = round(double(opt.ProgressInterval));
    if verbose
        fprintf('VSSI-Lp: assembling the weighted source system...\n');
        drawnow limitrate;
    end
    VtV = V' * V;
    dataGradient = 2 * (LSolve' * DataSolve);
    localCheckCancel(opt.CancelCheck);

    if verbose
        fprintf('VSSI-Lp: factorizing the %d-by-%d sparse source system...\n', ...
            nComponents, nComponents);
        drawnow limitrate;
    end
    [sourceSolver, ridgeValue] = localBuildSourceSolver(LSolve, VtV, rho, opt.SourceRidge);
    localCheckCancel(opt.CancelCheck);
    Source = zeros(nComponents, nTime);
    variation = zeros(nEdgeComponents, nTime);
    scaledDual = zeros(nEdgeComponents, nTime);

    objectiveHistory = nan(maxIterations, 1);
    dataFitHistory = nan(maxIterations, 1);
    variationPenaltyHistory = nan(maxIterations, 1);
    primalResidualHistory = nan(maxIterations, 1);
    dualResidualHistory = nan(maxIterations, 1);
    relativeChangeHistory = nan(maxIterations, 1);
    rhoHistory = nan(maxIterations, 1);
    converged = false;

    if verbose
        fprintf(['Running VSSI-Lp (locations=%d, edges=%d, nd=%d, T=%d, p=%.2f; ', ...
            'ADMM=%d, GST/ADMM=%d; initialization %.1fs)...\n'], ...
            nLocations, edgeInfo.NumEdges, nOrientations, nTime, pNorm, ...
            maxIterations, opt.GSTIterations, toc(initializationTimer));
        drawnow limitrate;
    end

    admmTimer = tic;
    for iter = 1:maxIterations
        localCheckCancel(opt.CancelCheck);
        previousSource = Source;
        previousVariation = variation;

        sourceRhs = dataGradient + rho * (V' * (variation - scaledDual));
        Source = localSolveSourceSystem(sourceSolver, LSolve, sourceRhs);
        sourceVariation = V * Source;

        variation = localGST(sourceVariation + scaledDual, pNorm, lambda / rho, opt.GSTIterations);
        scaledDual = scaledDual + sourceVariation - variation;

        primalResidual = norm(sourceVariation - variation, 'fro');
        dualResidual = rho * norm(V' * (variation - previousVariation), 'fro');
        primalTolerance = sqrt(nEdgeComponents * nTime) * absTolerance + ...
            relTolerance * max(norm(sourceVariation, 'fro'), norm(variation, 'fro'));
        dualTolerance = sqrt(nComponents * nTime) * absTolerance + ...
            relTolerance * rho * norm(V' * scaledDual, 'fro');
        relativeChange = norm(Source - previousSource, 'fro') / max(norm(previousSource, 'fro'), eps);

        dataFit = norm(DataSolve - LSolve * Source, 'fro')^2;
        variationPenalty = sum(abs(sourceVariation(:)).^pNorm);
        objectiveHistory(iter) = dataFit + lambda * variationPenalty;
        dataFitHistory(iter) = dataFit;
        variationPenaltyHistory(iter) = variationPenalty;
        primalResidualHistory(iter) = primalResidual;
        dualResidualHistory(iter) = dualResidual;
        relativeChangeHistory(iter) = relativeChange;
        rhoHistory(iter) = rho;

        if verbose && (iter == 1 || mod(iter, progressInterval) == 0 || ...
                (primalResidual <= primalTolerance && dualResidual <= dualTolerance))
            fprintf(['VSSI-Lp iter %d: objective = %.6g, primal = %.3e/%.3e, ', ...
                'dual = %.3e/%.3e, rho = %.3e (%.1fs)\n'], ...
                iter, objectiveHistory(iter), primalResidual, primalTolerance, ...
                dualResidual, dualTolerance, rho, toc(admmTimer));
            drawnow limitrate;
        end

        if primalResidual <= primalTolerance && dualResidual <= dualTolerance
            converged = true;
            break;
        end

        if adaptiveRho && mod(iter, opt.RhoUpdateInterval) == 0
            oldRho = rho;
            if primalResidual > opt.RhoBalance * max(dualResidual, eps)
                rho = rho * opt.RhoScale;
            elseif dualResidual > opt.RhoBalance * max(primalResidual, eps)
                rho = rho / opt.RhoScale;
            end
            if rho ~= oldRho
                scaledDual = scaledDual * (oldRho / rho);
                if verbose
                    fprintf('  VSSI-Lp: rho changed to %.3e; rebuilding the source system...\n', rho);
                    drawnow limitrate;
                end
                [sourceSolver, ridgeValue] = localBuildSourceSolver(LSolve, VtV, rho, opt.SourceRidge);
                localCheckCancel(opt.CancelCheck);
            end
        end
    end

    nIterations = find(~isnan(objectiveHistory), 1, 'last');
    if isempty(nIterations)
        nIterations = 0;
    end

    normalizedSource = Source;
    Source = coordinateInfo.SourceScale * Source;

    Parameters = struct();
    Parameters.Method = 'VSSI-Lp';
    Parameters.Reference = ['Peng S, Qi F, Yu H, Liu K (2024). EEG Extended Source Imaging ', ...
        'with Variation Sparsity and Lp-Norm Constraint. CICAI 2023, pp. 500-511.'];
    Parameters.DOI = '10.1007/978-981-99-9119-8_45';
    Parameters.EdgeOperator = V;
    Parameters.EdgeInfo = edgeInfo;
    Parameters.NumOrientations = nOrientations;
    Parameters.OrientationPenalty = 'componentwise_edge_lp';
    Parameters.RegularizationParameter = lambda;
    Parameters.P = pNorm;
    Parameters.SourceRidge = ridgeValue;
    Parameters.NormalizedSourceEstimate = normalizedSource;
    Parameters.SolverCoordinateNormalization = coordinateInfo;
    Parameters.ObjectiveHistory = objectiveHistory(1:nIterations);
    Parameters.DataFitHistory = dataFitHistory(1:nIterations);
    Parameters.VariationPenaltyHistory = variationPenaltyHistory(1:nIterations);
    Parameters.PrimalResidualHistory = primalResidualHistory(1:nIterations);
    Parameters.DualResidualHistory = dualResidualHistory(1:nIterations);
    Parameters.RelativeChangeHistory = relativeChangeHistory(1:nIterations);
    Parameters.RhoHistory = rhoHistory(1:nIterations);
    Parameters.FinalRho = rho;
    Parameters.NumIterationsPerformed = nIterations;
    Parameters.Converged = converged;
    Parameters.FinalVariation = V * Source;
    Parameters.FinalScaledDual = scaledDual;
    Parameters.NoiseWhitening = whiteningInfo;
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
        error('seal_VSSI_Lp:InvalidNoiseCovariance', ...
            'NoiseCovariance must be square with one row per sensor channel.');
    end
    noiseCovariance = double((noiseCovariance + noiseCovariance.') / 2);
    [factor, flag] = chol(noiseCovariance, 'lower');
    if flag ~= 0
        jitter = max(trace(noiseCovariance) / max(nChannels, 1), 1) * 1e-10;
        [factor, flag] = chol(noiseCovariance + jitter * eye(nChannels), 'lower');
    end
    if flag ~= 0
        error('seal_VSSI_Lp:InvalidNoiseCovariance', ...
            'NoiseCovariance must be positive definite for prewhitening.');
    end
    whiteningMatrix = factor \ eye(nChannels);
    DataWork = whiteningMatrix * Data;
    LWork = whiteningMatrix * L;
    info = struct('Applied', true, 'Matrix', whiteningMatrix);
end

function [V, info] = localResolveEdgeOperator(sourceSpaceInfo, nLocations, nOrientations)
    nComponents = nLocations * nOrientations;
    info = struct('Mode', '', 'NumEdges', 0);
    if isstruct(sourceSpaceInfo)
        if isfield(sourceSpaceInfo, 'EdgeOperator') && ~isempty(sourceSpaceInfo.EdgeOperator)
            candidate = sourceSpaceInfo.EdgeOperator;
            if size(candidate, 2) == nLocations
                V = kron(sparse(candidate), speye(nOrientations));
                info.Mode = 'edge_operator_expanded';
                info.NumEdges = size(candidate, 1);
            elseif size(candidate, 2) == nComponents
                V = sparse(candidate);
                if mod(size(V, 1), nOrientations) ~= 0
                    error('seal_VSSI_Lp:EdgeOperatorRows', ...
                        'Expanded EdgeOperator rows must be divisible by NumOrientations.');
                end
                info.Mode = 'edge_operator';
                info.NumEdges = size(V, 1) / nOrientations;
            else
                error('seal_VSSI_Lp:EdgeOperatorSize', ...
                    'EdgeOperator must have %d or %d columns.', nLocations, nComponents);
            end
        elseif isfield(sourceSpaceInfo, 'VertConn') && ~isempty(sourceSpaceInfo.VertConn)
            [Vlocation, edgeCount] = localIncidenceFromVertConn(sourceSpaceInfo.VertConn, nLocations);
            V = kron(Vlocation, speye(nOrientations));
            info.Mode = 'vertconn_incidence_expanded';
            info.NumEdges = edgeCount;
        else
            error('seal_VSSI_Lp:MissingGraph', ...
                'sourceSpaceInfo must contain EdgeOperator or VertConn.');
        end
    else
        [Vlocation, edgeCount] = localIncidenceFromVertConn(sourceSpaceInfo, nLocations);
        V = kron(Vlocation, speye(nOrientations));
        info.Mode = 'vertconn_incidence_expanded';
        info.NumEdges = edgeCount;
    end
    if any(~isfinite(nonzeros(V)))
        error('seal_VSSI_Lp:InvalidEdgeOperator', 'The edge operator must contain finite values.');
    end
    V = sparse(V);
end

function [V, nEdges] = localIncidenceFromVertConn(VertConn, nLocations)
    if ~isequal(size(VertConn), [nLocations, nLocations])
        error('seal_VSSI_Lp:VertConnSize', ...
            'VertConn must be %d-by-%d.', nLocations, nLocations);
    end
    if any(~isfinite(double(nonzeros(VertConn))))
        error('seal_VSSI_Lp:InvalidVertConn', 'VertConn must contain finite values.');
    end
    adjacency = spones(sparse(VertConn));
    adjacency = adjacency - spdiags(diag(adjacency), 0, nLocations, nLocations);
    adjacency = spones(adjacency + adjacency.');
    [row, col] = find(tril(adjacency, -1));
    nEdges = numel(row);
    edgeIndex = (1:nEdges).';
    V = sparse([edgeIndex; edgeIndex], [row; col], ...
        [ones(nEdges, 1); -ones(nEdges, 1)], nEdges, nLocations);
end

function [solver, ridgeValue] = localBuildSourceSolver(L, VtV, rho, ridgeFactor)
    nComponents = size(L, 2);
    nChannels = size(L, 1);
    diagonalScale = max(full(diag(VtV)));
    ridgeValue = max(ridgeFactor * max(rho * diagonalScale, 1), eps(max(rho * diagonalScale, 1)));
    baseSystem = rho * VtV + ridgeValue * speye(nComponents);
    [baseFactor, flag] = chol(baseSystem, 'lower');
    if flag ~= 0
        error('seal_VSSI_Lp:SourceFactorizationFailure', ...
            'Unable to factor the source-space ADMM system.');
    end
    baseInverseLeadField = baseFactor' \ (baseFactor \ L.');
    sensorSystem = 0.5 * eye(nChannels) + L * baseInverseLeadField;
    sensorSystem = (sensorSystem + sensorSystem.') / 2;
    [sensorFactor, flag] = chol(sensorSystem, 'lower');
    if flag ~= 0
        jitter = max(eps(norm(sensorSystem, 'fro')), 1e-12);
        [sensorFactor, flag] = chol(sensorSystem + jitter * eye(nChannels), 'lower');
    end
    if flag ~= 0
        error('seal_VSSI_Lp:SensorFactorizationFailure', ...
            'Unable to factor the Woodbury sensor-space system.');
    end
    solver = struct('BaseFactor', baseFactor, ...
        'BaseInverseLeadField', baseInverseLeadField, ...
        'SensorFactor', sensorFactor);
end

function Source = localSolveSourceSystem(solver, L, rhs)
    baseSolution = solver.BaseFactor' \ (solver.BaseFactor \ rhs);
    sensorCorrection = solver.SensorFactor' \ (solver.SensorFactor \ (L * baseSolution));
    Source = baseSolution - solver.BaseInverseLeadField * sensorCorrection;
end

function output = localGST(input, pNorm, tau, nIterations)
    if tau == 0
        output = input;
        return;
    end
    magnitude = abs(input);
    if pNorm == 1
        output = sign(input) .* max(magnitude - tau, 0);
        return;
    end

    fixedPoint = (2 * tau * (1 - pNorm))^(1 / (2 - pNorm));
    threshold = fixedPoint + tau * pNorm * fixedPoint^(pNorm - 1);
    active = magnitude > threshold;
    output = zeros(size(input));
    if ~any(active(:))
        return;
    end
    values = magnitude(active);
    for iter = 1:nIterations
        values = max(values, eps);
        values = magnitude(active) - tau * pNorm * values.^(pNorm - 1);
        values = max(values, 0);
    end
    output(active) = sign(input(active)) .* values;
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
valid = (islogical(value) || (isnumeric(value) && isscalar(value)) || ischar(value) || isstring(value));
end

function [DataOut, LOut, info] = localNormalizeSolveCoordinates(Data, L, nLocations, nOrientations, doNormalize)
    info = struct('Applied', false, 'DataScale', 1, 'LeadfieldScale', 1, ...
        'SourceScale', 1);
    DataOut = Data;
    LOut = L;
    if ~doNormalize
        return;
    end

    dataScale = norm(Data, 'fro') / sqrt(max(numel(Data), 1));
    if ~isfinite(dataScale) || dataScale <= realmin('double')
        error('seal_VSSI_Lp:ZeroData', 'Data must have non-zero finite energy.');
    end
    blockScales = zeros(nLocations, 1);
    for location = 1:nLocations
        index = (location - 1) * nOrientations + (1:nOrientations);
        blockScales(location) = norm(L(:, index), 'fro');
    end
    blockScales = blockScales(isfinite(blockScales) & blockScales > realmin('double'));
    if isempty(blockScales)
        error('seal_VSSI_Lp:ZeroLeadfield', 'Lead field must contain a non-zero finite block.');
    end
    leadfieldScale = median(blockScales);
    sourceScale = dataScale / leadfieldScale;
    DataOut = Data / dataScale;
    LOut = L * (sourceScale / dataScale);
    info = struct('Applied', true, 'DataScale', dataScale, ...
        'LeadfieldScale', leadfieldScale, 'SourceScale', sourceScale);
end

function localCheckCancel(cancelCheck)
drawnow limitrate;
if ~isempty(cancelCheck) && cancelCheck()
    error('SEAL:Cancelled', 'Source-imaging run was cancelled by the user.');
end
end
