function [Source, Parameters] = seal_VSSI_GGD(Data, L, sourceSpaceInfo, varargin)
%SEAL_VSSI_GGD Variation sparse source imaging with a GGD prior.
%
%   [Source, Parameters] = seal_VSSI_GGD(Data, L, sourceSpaceInfo, ...)
%   estimates extended sources using a generalized Gaussian distribution
%   (GGD) prior on edge-domain source differences.  The solver alternates
%   between an ADMM/IRLS source update and a conditional-mean posterior
%   covariance update:
%
%       gamma_e = p * lambda^(-p) * (|| (V*S)_e ||_2^2 + xi_e)^(p/2 - 1)
%
%   V is the cortical incidence operator and xi approximates the diagonal
%   of V*Sigma*V'.  This follows the VSSI conditional-mean framework and
%   the later VSSI-GGD generalized Gaussian prior.
%
%   sourceSpaceInfo may be a cortex struct containing .VertConn or
%   .EdgeOperator, or a location-by-location VertConn matrix directly.
%   Fixed orientation (NumOrientations=1) matches the published model.
%   For free orientations, SEAL expands V per component and applies the
%   same GGD edge prior to each orientation component.
%
%   Name-value options:
%       'SparseWeight'             Paper-compatible lambda scale (0.01)
%       'Lambda'                   Direct GGD scale override ([] derives it)
%       'P'                        GGD shape parameter in (0, 2] (0.7)
%       'MaxOuterIterations'       Maximum conditional-mean updates (10)
%       'MaxADMMIterations'        Maximum ADMM iterations per outer loop (400)
%       'MaxIRLSIterations'        Edge-domain IRLS updates per ADMM step (10)
%       'Tolerance'                Relative outer surrogate-cost tolerance (1e-3)
%       'ADMMTolerance'            Relative ADMM residual tolerance (1e-5)
%       'ProgressInterval'         ADMM iterations between progress messages (25)
%       'Rho'                      Initial ADMM penalty parameter (1)
%       'AdaptiveRho'              Residual-balanced rho updates (true)
%       'HutchinsonProbes'         Probes for diag(V*Sigma*V'') estimate (8)
%       'NoiseCovariance'          Optional sensor noise covariance ([])
%       'PrewhitenNoise'           Prewhiten if covariance is supplied (true)
%       'NormalizeProblem'         Normalize data/leadfield solver coordinates (true)
%       'NumOrientations'          Components per location: 1 or 3 (1)
%       'Verbose'                  Print concise outer-loop diagnostics (true)
%       'CancelCheck'              Optional zero-input cancellation function ([])

    parser = inputParser;
    parser.FunctionName = mfilename;
    parser.CaseSensitive = false;
    addRequired(parser, 'Data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addRequired(parser, 'L', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addRequired(parser, 'sourceSpaceInfo', @(x) isstruct(x) || isnumeric(x) || islogical(x));
    addParameter(parser, 'SparseWeight', 0.01, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'Lambda', [], @(x) isempty(x) || ...
        (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'P', 0.7, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 2);
    addParameter(parser, 'MaxOuterIterations', 10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'MaxADMMIterations', 400, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'MaxIRLSIterations', 10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'Tolerance', 1e-3, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'ADMMTolerance', 1e-5, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'ProgressInterval', 25, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'AbsoluteADMMTolerance', 1e-7, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'Rho', 1, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'AdaptiveRho', true, @localIsLogicalLike);
    addParameter(parser, 'RhoUpdateInterval', 10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'RhoBalance', 10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
    addParameter(parser, 'RhoScale', 2, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
    addParameter(parser, 'Relaxation', 0.6, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 2);
    addParameter(parser, 'HutchinsonProbes', 8, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    addParameter(parser, 'ProbeSeed', 1291, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x));
    addParameter(parser, 'SourceRidge', 1e-8, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'InitialMNERidge', 1e-6, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'VarianceFloor', 1e-12, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'WeightFloorRatio', 1e-6, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
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
        error('seal_VSSI_GGD:NonFiniteInput', 'Data and L must contain only finite values.');
    end
    if size(Data, 1) ~= size(L, 1)
        error('seal_VSSI_GGD:DimensionMismatch', ...
            'Data and L must have the same number of sensor rows.');
    end

    nTime = size(Data, 2);
    nComponents = size(L, 2);
    nOrientations = round(double(opt.NumOrientations));
    if mod(nComponents, nOrientations) ~= 0
        error('seal_VSSI_GGD:OrientationMismatch', ...
            'Lead-field columns must be divisible by NumOrientations.');
    end
    nLocations = nComponents / nOrientations;

    verbose = localToLogical(opt.Verbose);
    if verbose
        initializationTimer = tic;
        fprintf('Preparing VSSI-GGD: whitening data and lead field...\n');
        drawnow limitrate;
    end
    [DataWork, LWork, whiteningInfo] = localPrewhiten(Data, L, ...
        opt.NoiseCovariance, localToLogical(opt.PrewhitenNoise));
    localCheckCancel(opt.CancelCheck);
    [DataSolve, LSolve, coordinateInfo] = localNormalizeSolveCoordinates( ...
        DataWork, LWork, nLocations, nOrientations, localToLogical(opt.NormalizeProblem));
    if verbose && coordinateInfo.Applied
        fprintf(['VSSI-GGD: normalized solver coordinates (data RMS %.3e, ', ...
            'lead-field block scale %.3e, source scale %.3e).\n'], ...
            coordinateInfo.DataScale, coordinateInfo.LeadfieldScale, coordinateInfo.SourceScale);
        drawnow limitrate;
    end
    if verbose
        fprintf('VSSI-GGD: constructing the cortical edge operator...\n');
        drawnow limitrate;
    end
    [V, edgeInfo] = localResolveEdgeOperator(sourceSpaceInfo, nLocations, nOrientations);
    localCheckCancel(opt.CancelCheck);
    nEdgeComponents = size(V, 1);
    if nEdgeComponents == 0
        error('seal_VSSI_GGD:EmptyEdgeOperator', ...
            'The source-space graph must contain at least one edge.');
    end

    pNorm = double(opt.P);
    if verbose
        fprintf('VSSI-GGD: calibrating the GGD penalty scale...\n');
        drawnow limitrate;
    end
    lambda = localResolveLambda(opt.Lambda, opt.SparseWeight, V, LSolve, DataSolve);
    localCheckCancel(opt.CancelCheck);
    if verbose
        fprintf('VSSI-GGD: computing the initial wMNE estimate...\n');
        drawnow limitrate;
    end
    Source = localInitialMNE(DataSolve, LSolve, opt.InitialMNERidge);
    localCheckCancel(opt.CancelCheck);
    variation = V * Source;
    scaledDual = zeros(nEdgeComponents, nTime);
    rowEnergy = sum(variation .^ 2, 2);
    xi = max(mean(rowEnergy), opt.VarianceFloor) * ones(nEdgeComponents, 1);
    gamma = localGGDPrecision(variation, xi, lambda, pNorm, ...
        opt.WeightFloorRatio, opt.VarianceFloor);
    rho = double(opt.Rho);

    maxOuter = round(double(opt.MaxOuterIterations));
    outerHistory = repmat(localEmptyOuterHistory(), maxOuter, 1);
    converged = false;
    previousCost = inf;

    if verbose
        fprintf(['Running VSSI-GGD (locations=%d, edges=%d, nd=%d, T=%d, p=%.2f; ', ...
            'outer=%d, ADMM/outer=%d, IRLS/ADMM=%d; initialization %.1fs)...\n'], ...
            nLocations, edgeInfo.NumEdges, nOrientations, nTime, pNorm, ...
            maxOuter, opt.MaxADMMIterations, opt.MaxIRLSIterations, toc(initializationTimer));
        drawnow limitrate;
    end

    for outerIter = 1:maxOuter
        localCheckCancel(opt.CancelCheck);
        if verbose
            outerTimer = tic;
            fprintf('  VSSI-GGD outer %d/%d: starting ADMM...\n', outerIter, maxOuter);
            drawnow limitrate;
        end
        [Source, variation, scaledDual, rho, admmInfo] = localRunADMM( ...
            DataSolve, LSolve, V, Source, variation, scaledDual, xi, lambda, ...
            pNorm, rho, opt, outerIter, maxOuter, opt.CancelCheck);

        edgeVariation = V * Source;
        gamma = localGGDPrecision(edgeVariation, xi, lambda, pNorm, ...
            opt.WeightFloorRatio, opt.VarianceFloor);
        [posteriorSolver, ~] = localBuildWeightedSolver( ...
            LSolve, V, gamma, opt.SourceRidge);
        localCheckCancel(opt.CancelCheck);
        if verbose
            fprintf('    estimating posterior edge variance with %d Hutchinson probes...\n', ...
                opt.HutchinsonProbes);
            drawnow limitrate;
        end
        xi = localEstimateEdgeVariance(V, posteriorSolver, LSolve, nTime, ...
            opt.HutchinsonProbes, opt.ProbeSeed + outerIter - 1, opt.VarianceFloor);
        localCheckCancel(opt.CancelCheck);

        dataFit = 0.5 * norm(DataSolve - LSolve * Source, 'fro')^2;
        ggdPenalty = sum((sqrt(sum(edgeVariation .^ 2, 2) + xi) / lambda) .^ pNorm);
        surrogateCost = dataFit + ggdPenalty;
        relativeCostChange = abs(surrogateCost - previousCost) / max(abs(surrogateCost), eps);
        previousCost = surrogateCost;

        outerHistory(outerIter).SurrogateCost = surrogateCost;
        outerHistory(outerIter).DataFit = dataFit;
        outerHistory(outerIter).GGDPenalty = ggdPenalty;
        outerHistory(outerIter).RelativeCostChange = relativeCostChange;
        outerHistory(outerIter).Rho = rho;
        outerHistory(outerIter).ADMMIterations = admmInfo.NumIterations;
        outerHistory(outerIter).PrimalResidual = admmInfo.PrimalResidual;
        outerHistory(outerIter).DualResidual = admmInfo.DualResidual;
        outerHistory(outerIter).MeanGamma = mean(gamma);
        outerHistory(outerIter).MeanXi = mean(xi);

        if verbose
            fprintf(['VSSI-GGD outer %d: surrogate = %.6g, rel_change = %.3e, ', ...
                'ADMM = %d, rho = %.3e, elapsed %.1fs\n'], outerIter, surrogateCost, ...
                relativeCostChange, admmInfo.NumIterations, rho, toc(outerTimer));
            drawnow limitrate;
        end
        if outerIter > 1 && relativeCostChange <= opt.Tolerance
            converged = true;
            break;
        end
    end

    nOuter = outerIter;
    % The conditional mean of the final Gaussian posterior is the output.
    [posteriorSolver, ridgeValue] = localBuildWeightedSolver(LSolve, V, gamma, opt.SourceRidge);
    Source = localSolveWeightedSystem(posteriorSolver, LSolve, LSolve.' * DataSolve);
    normalizedSource = Source;
    Source = coordinateInfo.SourceScale * Source;

    Parameters = struct();
    Parameters.Method = 'VSSI-GGD';
    Parameters.Reference = ['Liu K, Peng S, Liang C, Yu ZL, Xiao B, Wang G, Wu W (2024). ', ...
        'VSSI-GGD: A Variation Sparse EEG Source Imaging Approach Based on Generalized Gaussian Distribution.'];
    Parameters.DOI = '10.1109/TNSRE.2024.3383452';
    Parameters.EdgeOperator = V;
    Parameters.EdgeInfo = edgeInfo;
    Parameters.NumOrientations = nOrientations;
    Parameters.OrientationPenalty = 'componentwise_edge_ggd';
    Parameters.Lambda = lambda;
    Parameters.SparseWeight = opt.SparseWeight;
    Parameters.P = pNorm;
    Parameters.Gamma = gamma;
    Parameters.Xi = xi;
    Parameters.SourceRidge = ridgeValue;
    Parameters.NormalizedSourceEstimate = normalizedSource;
    Parameters.SolverCoordinateNormalization = coordinateInfo;
    Parameters.OuterHistory = outerHistory(1:nOuter);
    Parameters.NumOuterIterations = nOuter;
    Parameters.Converged = converged;
    Parameters.FinalRho = rho;
    Parameters.FinalVariation = V * Source;
    Parameters.NoiseWhitening = whiteningInfo;
end

function lambda = localResolveLambda(lambdaOption, sparseWeight, V, L, Data)
    if ~isempty(lambdaOption)
        lambda = double(lambdaOption);
        return;
    end
    edgeCorrelation = sqrt(sum((V * (L.' * Data)) .^ 2, 2));
    scale = max(edgeCorrelation);
    lambda = 1 / (double(sparseWeight) * max(scale, eps));
end

function Source = localInitialMNE(Data, L, ridgeFactor)
    nChannels = size(L, 1);
    sensorCovariance = L * L.';
    ridge = max(double(ridgeFactor) * max(trace(sensorCovariance) / max(nChannels, 1), 1), eps);
    sensorCovariance = sensorCovariance + ridge * eye(nChannels);
    [factor, flag] = chol((sensorCovariance + sensorCovariance.') / 2, 'lower');
    if flag ~= 0
        error('seal_VSSI_GGD:InitialMNEFailure', ...
            'Unable to factor the initial minimum-norm sensor system.');
    end
    Source = L.' * (factor' \ (factor \ Data));
end

function [Source, U, scaledDual, rho, info] = localRunADMM( ...
        Data, L, V, Source, U, scaledDual, xi, lambda, pNorm, rho, opt, outerIter, maxOuter, cancelCheck)
    nComponents = size(L, 2);
    nEdges = size(V, 1);
    nTime = size(Data, 2);
    dataGradient = L.' * Data;
    maxADMM = round(double(opt.MaxADMMIterations));
    maxIRLS = round(double(opt.MaxIRLSIterations));
    relativeTolerance = double(opt.ADMMTolerance);
    absoluteTolerance = double(opt.AbsoluteADMMTolerance);
    adaptiveRho = localToLogical(opt.AdaptiveRho);
    verbose = localToLogical(opt.Verbose);
    progressInterval = round(double(opt.ProgressInterval));
    [sourceSolver, ~] = localBuildWeightedSolver(L, V, rho * ones(nEdges, 1), opt.SourceRidge);
    localCheckCancel(cancelCheck);

    converged = false;
    primalResidual = inf;
    dualResidual = inf;
    primalTolerance = inf;
    dualTolerance = inf;
    admmTimer = tic;
    for iter = 1:maxADMM
        localCheckCancel(cancelCheck);
        previousU = U;
        rhs = dataGradient + rho * (V.' * (U - scaledDual));
        Source = localSolveWeightedSystem(sourceSolver, L, rhs);
        sourceVariation = V * Source;
        relaxedVariation = opt.Relaxation * sourceVariation + (1 - opt.Relaxation) * previousU;

        U = localGGDIRLSProx(relaxedVariation + scaledDual, xi, lambda, pNorm, ...
            rho, maxIRLS, opt.WeightFloorRatio, opt.VarianceFloor);
        scaledDual = scaledDual + relaxedVariation - U;

        primalResidual = norm(sourceVariation - U, 'fro');
        dualResidual = rho * norm(V.' * (U - previousU), 'fro');
        primalTolerance = sqrt(nEdges * nTime) * absoluteTolerance + ...
            relativeTolerance * max(norm(sourceVariation, 'fro'), norm(U, 'fro'));
        dualTolerance = sqrt(nComponents * nTime) * absoluteTolerance + ...
            relativeTolerance * rho * norm(V.' * scaledDual, 'fro');
        if primalResidual <= primalTolerance && dualResidual <= dualTolerance
            converged = true;
            if verbose
                fprintf(['    VSSI-GGD ADMM %d/%d (outer %d/%d): r %.3g/%.3g, ', ...
                    's %.3g/%.3g, converged (%.1fs)\n'], iter, maxADMM, ...
                    outerIter, maxOuter, primalResidual, primalTolerance, ...
                    dualResidual, dualTolerance, toc(admmTimer));
                drawnow limitrate;
            end
            break;
        end

        if verbose && (iter == 1 || mod(iter, progressInterval) == 0)
            fprintf(['    VSSI-GGD ADMM %d/%d (outer %d/%d): r %.3g/%.3g, ', ...
                's %.3g/%.3g, rho %.3e (%.1fs)\n'], iter, maxADMM, ...
                outerIter, maxOuter, primalResidual, primalTolerance, ...
                dualResidual, dualTolerance, rho, toc(admmTimer));
            drawnow limitrate;
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
                [sourceSolver, ~] = localBuildWeightedSolver( ...
                    L, V, rho * ones(nEdges, 1), opt.SourceRidge);
                localCheckCancel(cancelCheck);
            end
        end
    end

    info = struct('NumIterations', iter, 'Converged', converged, ...
        'PrimalResidual', primalResidual, 'DualResidual', dualResidual, ...
        'PrimalTolerance', primalTolerance, 'DualTolerance', dualTolerance);
end

function U = localGGDIRLSProx(target, xi, lambda, pNorm, rho, maxIterations, weightFloorRatio, varianceFloor)
    nEdges = size(target, 1);
    nTime = size(target, 2);
    U = zeros(nEdges, nTime);
    weights = ones(nEdges, 1);
    for iter = 1:maxIterations
        previousU = U;
        U = (rho ./ (rho + 4 * weights)) .* target;
        weights = localGGDPrecision(U, xi, lambda, pNorm, weightFloorRatio, varianceFloor);
        relativeChange = norm(U - previousU, 'fro') / max(norm(U, 'fro'), eps);
        if relativeChange <= 1e-6
            break;
        end
    end
end

function gamma = localGGDPrecision(variation, xi, lambda, pNorm, weightFloorRatio, varianceFloor)
    energy = sum(variation .^ 2, 2) + max(xi, varianceFloor);
    gamma = pNorm * lambda^(-pNorm) * max(energy, varianceFloor) .^ ((pNorm - 2) / 2);
    floorValue = max(max(gamma) * weightFloorRatio, eps);
    gamma = max(gamma, floorValue);
end

function [solver, ridgeValue] = localBuildWeightedSolver(L, V, edgePrecision, ridgeFactor)
    nComponents = size(L, 2);
    nChannels = size(L, 1);
    edgePrecision = max(double(edgePrecision(:)), eps);
    baseSystem = V.' * (spdiags(edgePrecision, 0, numel(edgePrecision), numel(edgePrecision)) * V);
    diagonalScale = max(full(diag(baseSystem)));
    ridgeValue = max(double(ridgeFactor) * max(diagonalScale, 1), eps(max(diagonalScale, 1)));
    baseSystem = baseSystem + ridgeValue * speye(nComponents);
    [baseFactor, flag] = chol(baseSystem, 'lower');
    if flag ~= 0
        error('seal_VSSI_GGD:SourceFactorizationFailure', ...
            'Unable to factor the source-space precision matrix.');
    end
    baseInverseLeadField = baseFactor' \ (baseFactor \ L.');
    sensorSystem = eye(nChannels) + L * baseInverseLeadField;
    sensorSystem = (sensorSystem + sensorSystem.') / 2;
    [sensorFactor, flag] = chol(sensorSystem, 'lower');
    if flag ~= 0
        jitter = max(eps(norm(sensorSystem, 'fro')), 1e-12);
        [sensorFactor, flag] = chol(sensorSystem + jitter * eye(nChannels), 'lower');
    end
    if flag ~= 0
        error('seal_VSSI_GGD:SensorFactorizationFailure', ...
            'Unable to factor the sensor-space Woodbury system.');
    end
    solver = struct('BaseFactor', baseFactor, ...
        'BaseInverseLeadField', baseInverseLeadField, ...
        'SensorFactor', sensorFactor);
end

function Source = localSolveWeightedSystem(solver, L, rhs)
    baseSolution = solver.BaseFactor' \ (solver.BaseFactor \ rhs);
    sensorCorrection = solver.SensorFactor' \ (solver.SensorFactor \ (L * baseSolution));
    Source = baseSolution - solver.BaseInverseLeadField * sensorCorrection;
end

function xi = localEstimateEdgeVariance(V, solver, L, nTime, nProbes, seed, varianceFloor)
    nEdges = size(V, 1);
    stream = RandStream('mt19937ar', 'Seed', double(seed));
    probes = sign(randn(stream, nEdges, round(double(nProbes))));
    probes(probes == 0) = 1;
    posteriorSamples = localSolveWeightedSystem(solver, L, V.' * probes);
    diagonalEstimate = sum(probes .* (V * posteriorSamples), 2) / size(probes, 2);
    xi = nTime * max(real(diagonalEstimate), double(varianceFloor));
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
        error('seal_VSSI_GGD:InvalidNoiseCovariance', ...
            'NoiseCovariance must be square with one row per sensor channel.');
    end
    noiseCovariance = double((noiseCovariance + noiseCovariance.') / 2);
    [factor, flag] = chol(noiseCovariance, 'lower');
    if flag ~= 0
        jitter = max(trace(noiseCovariance) / max(nChannels, 1), 1) * 1e-10;
        [factor, flag] = chol(noiseCovariance + jitter * eye(nChannels), 'lower');
    end
    if flag ~= 0
        error('seal_VSSI_GGD:InvalidNoiseCovariance', ...
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
                    error('seal_VSSI_GGD:EdgeOperatorRows', ...
                        'Expanded EdgeOperator rows must be divisible by NumOrientations.');
                end
                info.Mode = 'edge_operator';
                info.NumEdges = size(V, 1) / nOrientations;
            else
                error('seal_VSSI_GGD:EdgeOperatorSize', ...
                    'EdgeOperator must have %d or %d columns.', nLocations, nComponents);
            end
        elseif isfield(sourceSpaceInfo, 'VertConn') && ~isempty(sourceSpaceInfo.VertConn)
            [Vlocation, edgeCount] = localIncidenceFromVertConn(sourceSpaceInfo.VertConn, nLocations);
            V = kron(Vlocation, speye(nOrientations));
            info.Mode = 'vertconn_incidence_expanded';
            info.NumEdges = edgeCount;
        else
            error('seal_VSSI_GGD:MissingGraph', ...
                'sourceSpaceInfo must contain EdgeOperator or VertConn.');
        end
    else
        [Vlocation, edgeCount] = localIncidenceFromVertConn(sourceSpaceInfo, nLocations);
        V = kron(Vlocation, speye(nOrientations));
        info.Mode = 'vertconn_incidence_expanded';
        info.NumEdges = edgeCount;
    end
    if any(~isfinite(nonzeros(V)))
        error('seal_VSSI_GGD:InvalidEdgeOperator', 'The edge operator must contain finite values.');
    end
    V = sparse(V);
end

function [V, nEdges] = localIncidenceFromVertConn(VertConn, nLocations)
    if ~isequal(size(VertConn), [nLocations, nLocations])
        error('seal_VSSI_GGD:VertConnSize', ...
            'VertConn must be %d-by-%d.', nLocations, nLocations);
    end
    if any(~isfinite(double(nonzeros(VertConn))))
        error('seal_VSSI_GGD:InvalidVertConn', 'VertConn must contain finite values.');
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

function history = localEmptyOuterHistory()
    history = struct('SurrogateCost', 0, 'DataFit', 0, 'GGDPenalty', 0, ...
        'RelativeCostChange', inf, 'Rho', 0, 'ADMMIterations', 0, ...
        'PrimalResidual', inf, 'DualResidual', inf, 'MeanGamma', 0, 'MeanXi', 0);
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
        error('seal_VSSI_GGD:ZeroData', 'Data must have non-zero finite energy.');
    end
    blockScales = zeros(nLocations, 1);
    for location = 1:nLocations
        index = (location - 1) * nOrientations + (1:nOrientations);
        blockScales(location) = norm(L(:, index), 'fro');
    end
    blockScales = blockScales(isfinite(blockScales) & blockScales > realmin('double'));
    if isempty(blockScales)
        error('seal_VSSI_GGD:ZeroLeadfield', 'Lead field must contain a non-zero finite block.');
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
