function [Source_Estimate, Parameters] = seal_STBFSI(Data, L, varargin)
%SEAL_STBFSI Source imaging based on spatio-temporal basis functions.
%   [Source, Param] = seal_STBFSI(Data, L, sourceSpaceInfo, ...)
%   [Source, Param] = seal_STBFSI(Data, L, VertConn, ...)
%   [Source, Param] = seal_STBFSI(Data, L, 'SpatialBasis', A, ...)
%
%   This is a MATLAB implementation of the SI-STBF model in:
%   Liu et al. (2019), Bayesian Electromagnetic Spatio-Temporal Imaging of
%   Extended Sources Based on Matrix Factorization, IEEE TBME.
%
%   The model uses S = A*Theta*Phi and B = L*A*Theta*Phi + noise, where Phi
%   contains temporal basis functions and A contains spatial basis functions.
%   If a spatial basis is not supplied, this function builds local graph
%   diffusion bases from SourceSpaceInfo.VertConn. The paper's DDP bases can
%   be reproduced by passing them directly through 'SpatialBasis'.

    [posSourceInfo, posVertConn, nvArgs, inputInfo] = preparseStbfsiInputs(varargin{:});

    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'SourceSpaceInfo', posSourceInfo, @isstruct);
    addParameter(p, 'VertConn', posVertConn, @(x) isempty(x) || isnumeric(x) || islogical(x) || issparse(x));
    addParameter(p, 'SpatialBasis', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'InitialTBFs', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NumTBFs', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'MaxSpatialBases', 120, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NeighborOrder', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'DiffusionSigma', 0.6, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'DiffusionOrder', 8, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'SeedSelection', 'mne', @(x) ischar(x) || isstring(x));
    addParameter(p, 'MaxIterations', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'PruneAlphas', 1e-4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'PruneGammas', 1e-4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'PruneSpatialBases', true, @(x) isLogicalLike(x));
    addParameter(p, 'GammaUpdate', 'em', @(x) ischar(x) || isstring(x));
    addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'PrewhitenNoise', true, @(x) isLogicalLike(x));
    addParameter(p, 'NoiseObservationMatrix', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x)); %#ok<NVREPL>
    addParameter(p, 'Regularization', 1e-8, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'CovarianceMode', 'auto', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FullCovarianceLimit', 450, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Verbose', true, @(x) isLogicalLike(x));

    try
        parse(p, Data, L, nvArgs{:});
    catch ME
        disp('Error parsing inputs for seal_STBFSI:');
        rethrow(ME);
    end
    opts = p.Results;

    opts.NumOrientations = round(opts.NumOrientations);
    opts.NumTBFs = round(opts.NumTBFs);
    opts.MaxSpatialBases = round(opts.MaxSpatialBases);
    opts.NeighborOrder = round(opts.NeighborOrder);
    opts.DiffusionOrder = round(opts.DiffusionOrder);
    opts.MaxIterations = round(opts.MaxIterations);
    opts.FullCovarianceLimit = round(opts.FullCovarianceLimit);
    opts.PruneSpatialBases = parseLogical(opts.PruneSpatialBases, 'PruneSpatialBases');
    opts.PrewhitenNoise = parseLogical(opts.PrewhitenNoise, 'PrewhitenNoise');
    opts.Verbose = parseLogical(opts.Verbose, 'Verbose');
    opts.GammaUpdate = lower(char(opts.GammaUpdate));
    opts.CovarianceMode = lower(char(opts.CovarianceMode));

    if ~strcmp(opts.GammaUpdate, 'em')
        warning('seal_STBFSI:UnsupportedGammaUpdate', ...
            'Only the SI-STBF EM gamma update is implemented. Falling back to ''em''.');
        opts.GammaUpdate = 'em';
    end
    if ~any(strcmp(opts.CovarianceMode, {'auto', 'full', 'diagonal'}))
        error('seal_STBFSI:InvalidCovarianceMode', ...
            'CovarianceMode must be ''auto'', ''full'', or ''diagonal''.');
    end

    if isempty(opts.VertConn) && isfield(opts.SourceSpaceInfo, 'VertConn')
        opts.VertConn = opts.SourceSpaceInfo.VertConn;
        inputInfo.ExplicitVertConn = true;
    end

    Data0 = double(opts.Data);
    L0 = double(opts.L);
    [nChan, nTime] = size(Data0);
    if size(L0, 1) ~= nChan
        error('seal_STBFSI:DataLeadfieldMismatch', ...
            'Data rows (%d) must match L rows (%d).', nChan, size(L0, 1));
    end
    if any(~isfinite(Data0(:))) || any(~isfinite(L0(:)))
        error('seal_STBFSI:InvalidInput', 'Data and L must contain finite values.');
    end

    nd = opts.NumOrientations;
    nSource = size(L0, 2);
    if mod(nSource, nd) ~= 0
        error('seal_STBFSI:NumOrientationsMismatch', ...
            'size(L,2) (%d) must be divisible by NumOrientations (%d).', nSource, nd);
    end
    nLoc = nSource / nd;

    [DataW, LW, Cnoise, noiseInfo] = prepareNoiseModel( ...
        Data0, L0, opts.NoiseCovariance, opts.PrewhitenNoise);
    opts.NoiseCovariance = Cnoise;

    [A, spatialInfo] = resolveSpatialBasis(DataW, LW, opts, nLoc, nSource, nd, inputInfo);
    if isempty(A) || size(A, 2) == 0
        error('seal_STBFSI:EmptySpatialBasis', 'No spatial basis functions were generated.');
    end

    F = LW * A;
    keepGood = sqrt(sum(F.^2, 1)).' > eps;
    if ~all(keepGood)
        A = A(:, keepGood);
        F = F(:, keepGood);
        spatialInfo.RemovedZeroForwardBases = find(~keepGood);
    else
        spatialInfo.RemovedZeroForwardBases = [];
    end
    Rbasis = size(F, 2);
    spatialKeepList = (1:Rbasis).';

    [PhiMean, temporalInfo] = resolveTemporalBasis(opts.InitialTBFs, DataW, opts.NumTBFs, nTime);
    K = size(PhiMean, 1);
    temporalKeepList = (1:K).';

    ThetaMean = zeros(Rbasis, K);
    thetaCovDiag = zeros(Rbasis, K);
    thetaCovTrace = zeros(K, 1);
    SigmaPhi = eye(K) * 1e-6;
    alpha = ones(K, 1);
    gamma = ones(Rbasis, 1);
    costHistory = nan(opts.MaxIterations, 1);
    converged = false;
    usedDiagonalCovariance = false;

    if opts.Verbose
        fprintf('Running SI-STBF (nd=%d, K=%d, spatial bases=%d)...\n', nd, K, Rbasis);
    end

    previousCost = inf;
    for iter = 1:opts.MaxIterations
        FtF = full(F' * F);
        FtF = (FtF + FtF') / 2;
        FtB = F' * DataW;

        Ephi = PhiMean * PhiMean' + nTime * SigmaPhi;
        Ephi = full((Ephi + Ephi') / 2);

        for k = 1:K
            otherTheta = ThetaMean * Ephi(:, k) - ThetaMean(:, k) * Ephi(k, k);
            rhsTheta = FtB * PhiMean(k, :)' - FtF * otherTheta;
            PrecTheta = Ephi(k, k) * FtF;
            diagIdxTheta = 1:Rbasis + 1:numel(PrecTheta);
            PrecTheta(diagIdxTheta) = PrecTheta(diagIdxTheta) + gamma(:).';

            [thetaK, covDiagK, covTraceK, modeUsed] = posteriorMeanAndMoments( ...
                PrecTheta, rhsTheta, FtF, opts.CovarianceMode, opts.FullCovarianceLimit, opts.Regularization);
            ThetaMean(:, k) = thetaK;
            thetaCovDiag(:, k) = covDiagK;
            thetaCovTrace(k) = covTraceK;
            usedDiagonalCovariance = usedDiagonalCovariance || strcmp(modeUsed, 'diagonal');
        end

        Etheta = ThetaMean' * FtF * ThetaMean + diag(thetaCovTrace);
        Etheta = full((Etheta + Etheta') / 2);
        PrecPhi = Etheta;
        PrecPhi(1:K + 1:end) = PrecPhi(1:K + 1:end) + alpha(:).';
        SigmaPhi = invertSpd(PrecPhi, opts.Regularization);
        PhiMean = SigmaPhi * (ThetaMean' * FtB);

        phiEnergy = sum(PhiMean.^2, 2) + nTime * diag(SigmaPhi);
        alpha = nTime ./ max(phiEnergy, realmin);

        thetaEnergy = sum(ThetaMean.^2, 2) + sum(thetaCovDiag, 2);
        gamma = K ./ max(thetaEnergy, realmin);

        coeff = ThetaMean * PhiMean;
        residual = DataW - F * coeff;
        cost = 0.5 * sum(residual(:).^2) ...
            + 0.5 * sum(gamma .* thetaEnergy) ...
            + 0.5 * sum(alpha .* phiEnergy);
        costHistory(iter) = cost;

        if isfinite(previousCost)
            relChange = abs(previousCost - cost) / max(1, abs(previousCost));
        else
            relChange = inf;
        end
        if opts.Verbose
            fprintf('SI-STBF iteration %d: cost %.6g, rel change %.3g, K=%d, R=%d\n', ...
                iter, cost, relChange, K, Rbasis);
        end
        if iter > 1 && relChange < opts.Tolerance
            converged = true;
            break;
        end
        previousCost = cost;

        [PhiMean, SigmaPhi, alpha, ThetaMean, thetaCovDiag, temporalKeepList] = pruneTemporalBases( ...
            PhiMean, SigmaPhi, alpha, ThetaMean, thetaCovDiag, temporalKeepList, opts.PruneAlphas);
        K = size(PhiMean, 1);
        if numel(thetaCovTrace) ~= K
            thetaCovTrace = thetaCovTrace(1:K);
        end

        if opts.PruneSpatialBases && iter > 1
            [A, F, gamma, ThetaMean, thetaCovDiag, spatialKeepList] = pruneSpatialBases( ...
                A, F, gamma, ThetaMean, thetaCovDiag, spatialKeepList, opts.PruneGammas);
            Rbasis = size(F, 2);
        end
    end

    nIter = find(~isnan(costHistory), 1, 'last');
    if isempty(nIter)
        nIter = 0;
        costHistory = [];
    else
        costHistory = costHistory(1:nIter);
    end

    Source_Estimate = full(A * (ThetaMean * PhiMean));

    Parameters = struct();
    Parameters.Algorithm = 'SI-STBF';
    Parameters.Reference = 'Liu et al. 2019, IEEE TBME, DOI:10.1109/TBME.2018.2890291';
    Parameters.SourceModel = 'S = A*Theta*Phi';
    Parameters.SpatialBasis = A;
    Parameters.EffectiveForward = F;
    Parameters.ThetaMean = ThetaMean;
    Parameters.ThetaCovarianceDiag = thetaCovDiag;
    Parameters.TemporalBasisFunctions = PhiMean;
    Parameters.SigmaPhi = SigmaPhi;
    Parameters.Alpha = alpha;
    Parameters.Gamma = gamma;
    Parameters.ActiveSpatialBasis = spatialKeepList;
    Parameters.TemporalBasisKeepList = temporalKeepList;
    Parameters.CostHistory = costHistory;
    Parameters.NumIterationsPerformed = nIter;
    Parameters.Converged = converged;
    Parameters.NoiseWhitening = noiseInfo;
    Parameters.SpatialBasisInfo = spatialInfo;
    Parameters.TemporalBasisInfo = temporalInfo;
    Parameters.InputParsing = inputInfo;
    Parameters.OptionsPassed = opts;
    Parameters.PosteriorCovarianceMode = opts.CovarianceMode;
    Parameters.UsedDiagonalCovarianceApproximation = usedDiagonalCovariance;
end

function [sourceInfo, vertConn, nvArgs, inputInfo] = preparseStbfsiInputs(varargin)
    sourceInfo = struct();
    vertConn = [];
    nvArgs = varargin;
    inputInfo = struct( ...
        'PositionalSourceSpaceInfo', false, ...
        'ExplicitVertConn', false, ...
        'ExplicitSpatialBasis', false);

    if isempty(nvArgs)
        return;
    end

    first = nvArgs{1};
    if ischar(first) || isstring(first)
        return;
    end

    nvArgs(1) = [];
    if isstruct(first)
        sourceInfo = first;
        inputInfo.PositionalSourceSpaceInfo = true;
        if isfield(first, 'VertConn')
            vertConn = first.VertConn;
            inputInfo.ExplicitVertConn = true;
        end
    else
        vertConn = first;
        inputInfo.ExplicitVertConn = true;
    end
    inputInfo.ExplicitSpatialBasis = hasNameValue(nvArgs, 'SpatialBasis');
end

function tf = hasNameValue(args, name)
    tf = false;
    for ii = 1:2:numel(args)
        if ischar(args{ii}) || isstring(args{ii})
            if strcmpi(char(args{ii}), name)
                tf = true;
                return;
            end
        end
    end
end

function tf = isLogicalLike(x)
    tf = islogical(x) || (isnumeric(x) && isscalar(x)) || ischar(x) || isstring(x);
end

function tf = parseLogical(x, name)
    if islogical(x)
        tf = x;
    elseif isnumeric(x) && isscalar(x)
        tf = x ~= 0;
    elseif ischar(x) || isstring(x)
        x = lower(strtrim(char(x)));
        if any(strcmp(x, {'true', '1', 'yes', 'on'}))
            tf = true;
        elseif any(strcmp(x, {'false', '0', 'no', 'off'}))
            tf = false;
        else
            error('seal_STBFSI:InvalidLogical', '%s must be logical-compatible.', name);
        end
    else
        error('seal_STBFSI:InvalidLogical', '%s must be logical-compatible.', name);
    end
end

function [DataW, LW, Cnoise, noiseInfo] = prepareNoiseModel(Data, L, Cnoise, prewhiten)
    nChan = size(Data, 1);
    noiseInfo = struct('Applied', false, 'Method', 'none', 'Regularization', 0);
    if isempty(Cnoise)
        Cnoise = eye(nChan);
        DataW = Data;
        LW = L;
        return;
    end
    if size(Cnoise, 1) ~= nChan || size(Cnoise, 2) ~= nChan
        error('seal_STBFSI:NoiseCovarianceSize', ...
            'NoiseCovariance must be %d x %d.', nChan, nChan);
    end
    Cnoise = double((Cnoise + Cnoise') / 2);
    if ~prewhiten
        DataW = Data;
        LW = L;
        return;
    end

    [V, D] = eig(Cnoise);
    d = real(diag(D));
    floorVal = max(max(d) * 1e-12, eps);
    d(d < floorVal) = floorVal;
    W = diag(1 ./ sqrt(d)) * V';
    DataW = W * Data;
    LW = W * L;
    noiseInfo = struct('Applied', true, 'Method', 'eig', 'Regularization', floorVal);
end

function [A, spatialInfo] = resolveSpatialBasis(Data, L, opts, nLoc, nSource, nd, inputInfo)
    spatialInfo = struct();
    spatialInfo.Mode = '';
    spatialInfo.NumLocations = nLoc;
    spatialInfo.NumOrientations = nd;
    spatialInfo.CustomSpatialBasis = inputInfo.ExplicitSpatialBasis;

    if ~isempty(opts.SpatialBasis)
        [A, customInfo] = sanitizeSpatialBasis(opts.SpatialBasis, nLoc, nSource, nd);
        spatialInfo.Mode = 'custom';
        spatialInfo.CustomInfo = customInfo;
        return;
    end

    if ~isempty(opts.VertConn)
        Adj = sanitizeAdjacency(opts.VertConn, nLoc);
        [A_loc, graphInfo] = buildGraphSpatialBasis(Data, L, Adj, opts, nLoc, nd);
        A = kron(A_loc, speye(nd));
        spatialInfo.Mode = 'graph-diffusion';
        spatialInfo.GraphInfo = graphInfo;
        spatialInfo.Note = ['Default graph-diffusion bases approximate the SI-STBF spatial prior. ', ...
            'Pass paper DDP bases through SpatialBasis for exact DDP components.'];
    else
        [A_loc, seedInfo] = buildIdentitySpatialBasis(Data, L, opts, nLoc, nd);
        A = kron(A_loc, speye(nd));
        spatialInfo.Mode = 'identity-seeds';
        spatialInfo.SeedInfo = seedInfo;
        spatialInfo.Note = 'No VertConn supplied; using selected source-location identity bases.';
    end

    A = normalizeColumns(A);
end

function [A, info] = sanitizeSpatialBasis(SpatialBasis, nLoc, nSource, nd)
    B = sparse(double(SpatialBasis));
    info = struct('InputRows', size(B, 1), 'InputColumns', size(B, 2));
    if size(B, 1) == nSource
        A = B;
        info.ExpandedOrientations = false;
    elseif size(B, 1) == nLoc
        A = kron(B, speye(nd));
        info.ExpandedOrientations = true;
    else
        error('seal_STBFSI:SpatialBasisSize', ...
            'SpatialBasis rows must match nLocations (%d) or nSources (%d).', nLoc, nSource);
    end
    A = normalizeColumns(A);
    colNorms = sqrt(sum(A.^2, 1)).';
    keep = colNorms > eps;
    A = A(:, keep);
    info.RemovedZeroColumns = find(~keep);
end

function Adj = sanitizeAdjacency(VertConn, nLoc)
    if size(VertConn, 1) ~= nLoc || size(VertConn, 2) ~= nLoc
        error('seal_STBFSI:VertConnSize', ...
            'VertConn must be %d x %d source-location adjacency.', nLoc, nLoc);
    end
    Adj = spones(sparse(VertConn));
    Adj = spones(Adj + Adj');
    Adj = Adj - spdiags(diag(Adj), 0, nLoc, nLoc);
    Adj = spones(Adj);
end

function [A_loc, info] = buildGraphSpatialBasis(Data, L, Adj, opts, nLoc, nd)
    maxBases = min(opts.MaxSpatialBases, nLoc);
    Reach = graphReach(Adj, opts.NeighborOrder);
    score = computeSeedScores(Data, L, nLoc, nd, opts.SeedSelection, opts.Regularization);
    seeds = selectSeeds(score, Reach, maxBases);

    deg = full(sum(Adj, 2));
    G = Adj - spdiags(deg, 0, nLoc, nLoc);
    A_loc = spalloc(nLoc, numel(seeds), max(1, numel(seeds)) * 16);
    for jj = 1:numel(seeds)
        support = find(Reach(seeds(jj), :));
        col = diffusionColumn(G, seeds(jj), support(:), opts.DiffusionSigma, opts.DiffusionOrder);
        A_loc(:, jj) = sparse(col);
    end
    A_loc = normalizeColumns(A_loc);

    info = struct();
    info.Seeds = seeds;
    info.SeedScores = score(seeds);
    info.NeighborOrder = opts.NeighborOrder;
    info.DiffusionSigma = opts.DiffusionSigma;
    info.DiffusionOrder = opts.DiffusionOrder;
    info.NumBasisLocations = size(A_loc, 2);
    info.SeedSelection = char(opts.SeedSelection);
end

function [A_loc, info] = buildIdentitySpatialBasis(Data, L, opts, nLoc, nd)
    maxBases = min(opts.MaxSpatialBases, nLoc);
    score = computeSeedScores(Data, L, nLoc, nd, opts.SeedSelection, opts.Regularization);
    [~, ord] = sort(score, 'descend');
    seeds = ord(1:maxBases);
    A_loc = sparse(seeds, 1:numel(seeds), 1, nLoc, numel(seeds));

    info = struct();
    info.Seeds = seeds(:);
    info.SeedScores = score(seeds);
    info.NumBasisLocations = size(A_loc, 2);
    info.SeedSelection = char(opts.SeedSelection);
end

function score = computeSeedScores(Data, L, nLoc, nd, seedMode, reg)
    seedMode = lower(char(seedMode));
    nSource = size(L, 2);
    switch seedMode
        case 'mne'
            C = L * L';
            lambda = reg * max(trace(C) / max(1, size(C, 1)), eps);
            if lambda == 0
                lambda = eps;
            end
            Kinv = L' / (C + lambda * eye(size(C)));
            S0 = Kinv * Data;
            sourceEnergy = sum(S0.^2, 2);
        case 'leadfield'
            sourceEnergy = sum(L.^2, 1).';
        case 'uniform'
            sourceEnergy = ones(nSource, 1);
        otherwise
            error('seal_STBFSI:UnknownSeedSelection', ...
                'SeedSelection must be ''mne'', ''leadfield'', or ''uniform''.');
    end

    score = zeros(nLoc, 1);
    for ii = 1:nLoc
        idx = (ii - 1) * nd + (1:nd);
        score(ii) = sum(sourceEnergy(idx));
    end
    score(~isfinite(score)) = 0;
    if all(score == 0)
        score(:) = 1;
    end
end

function Reach = graphReach(Adj, neighborOrder)
    nLoc = size(Adj, 1);
    Reach = speye(nLoc);
    if neighborOrder <= 0
        return;
    end
    Step = spones(Adj + speye(nLoc));
    for kk = 1:neighborOrder
        Reach = spones(Reach * Step);
    end
end

function seeds = selectSeeds(score, Reach, maxSeeds)
    [~, ord] = sort(score, 'descend');
    blocked = false(numel(score), 1);
    selected = zeros(maxSeeds, 1);
    count = 0;
    for ii = 1:numel(ord)
        cand = ord(ii);
        if ~blocked(cand)
            count = count + 1;
            selected(count) = cand;
            blocked(find(Reach(cand, :))) = true; %#ok<FNDSB>
            if count >= maxSeeds
                break;
            end
        end
    end
    if count < maxSeeds
        already = false(numel(score), 1);
        already(selected(1:count)) = true;
        for ii = 1:numel(ord)
            cand = ord(ii);
            if ~already(cand)
                count = count + 1;
                selected(count) = cand;
                if count >= maxSeeds
                    break;
                end
            end
        end
    end
    seeds = selected(1:count);
end

function col = diffusionColumn(G, seed, support, sigma, diffusionOrder)
    nLoc = size(G, 1);
    v = sparse(seed, 1, 1, nLoc, 1);
    col = v;
    term = v;
    for kk = 1:diffusionOrder
        term = (sigma / kk) * (G * term);
        col = col + term;
    end
    col = full(col);
    mask = false(nLoc, 1);
    mask(support) = true;
    col(~mask) = 0;
    col(col < 0) = 0;
    if norm(col) <= eps
        col = zeros(nLoc, 1);
        col(seed) = 1;
    else
        col = col / norm(col);
    end
end

function A = normalizeColumns(A)
    norms = sqrt(sum(A.^2, 1)).';
    keep = norms > eps;
    if any(keep)
        scale = ones(numel(norms), 1);
        scale(keep) = 1 ./ norms(keep);
        A = A * spdiags(scale, 0, numel(scale), numel(scale));
    end
end

function [Phi, info] = resolveTemporalBasis(initialPhi, Data, numTBFs, nTime)
    info = struct();
    if isempty(initialPhi)
        K = min([numTBFs, nTime, size(Data, 1)]);
        [~, ~, V] = svd(Data, 'econ');
        Phi = V(:, 1:K)';
        info.Mode = 'svd';
        info.NumRequested = numTBFs;
    else
        Phi = double(initialPhi);
        if size(Phi, 2) ~= nTime && size(Phi, 1) == nTime
            Phi = Phi';
            info.TransposedInput = true;
        else
            info.TransposedInput = false;
        end
        if size(Phi, 2) ~= nTime
            error('seal_STBFSI:InitialTBFsTimeMismatch', ...
                'InitialTBFs must have %d columns.', nTime);
        end
        info.Mode = 'user';
        info.NumRequested = size(Phi, 1);
    end
    if isempty(Phi) || size(Phi, 1) < 1
        error('seal_STBFSI:NoTemporalBasis', 'At least one temporal basis function is required.');
    end
    rowNorm = sqrt(sum(Phi.^2, 2));
    if any(rowNorm <= eps)
        error('seal_STBFSI:InvalidTemporalBasis', ...
            'Temporal basis rows must have non-zero energy.');
    end
    Phi = bsxfun(@rdivide, Phi, rowNorm);
    info.NumBasisFunctions = size(Phi, 1);
end

function [meanVec, covDiag, covTrace, modeUsed] = posteriorMeanAndMoments(Prec, rhs, FtF, mode, fullLimit, reg)
    n = size(Prec, 1);
    if strcmp(mode, 'auto')
        if n <= fullLimit
            modeUsed = 'full';
        else
            modeUsed = 'diagonal';
        end
    else
        modeUsed = mode;
    end

    Prec = full((Prec + Prec') / 2);
    if strcmp(modeUsed, 'diagonal')
        d = max(diag(Prec), realmin);
        covDiag = 1 ./ d;
        meanVec = covDiag .* rhs;
        covTrace = sum(diag(FtF) .* covDiag);
        return;
    end

    [R, jitter] = cholSafe(Prec, reg);
    meanVec = R \ (R' \ rhs);
    IR = R \ eye(n);
    covDiag = sum(IR.^2, 2);
    covTrace = sum(sum((FtF * IR) .* IR));
    covTrace = max(real(covTrace), 0);
    if jitter > 0
        covDiag = max(real(covDiag), realmin);
    end
end

function [Ainv, info] = invertSpd(A, reg)
    A = full((A + A') / 2);
    [R, jitter] = cholSafe(A, reg);
    IR = R \ eye(size(A, 1));
    Ainv = IR * IR';
    Ainv = (Ainv + Ainv') / 2;
    info = struct('Jitter', jitter);
end

function [R, jitter] = cholSafe(A, reg)
    n = size(A, 1);
    scale = max(mean(abs(diag(A))), 1);
    jitter = max(reg * scale, eps);
    [R, p] = chol(A);
    if p == 0
        jitter = 0;
        return;
    end
    for ii = 1:8
        [R, p] = chol(A + jitter * eye(n));
        if p == 0
            return;
        end
        jitter = jitter * 10;
    end
    warning('seal_STBFSI:CholeskyFailed', ...
        'Precision matrix was not SPD after jitter. Using eigenvalue repair.');
    A = (A + A') / 2;
    [V, D] = eig(A);
    d = real(diag(D));
    floorVal = max(max(abs(d)) * 1e-10, eps);
    d(d < floorVal) = floorVal;
    Afix = V * diag(d) * V';
    Afix = (Afix + Afix') / 2;
    R = chol(Afix + floorVal * eye(n));
    jitter = floorVal;
end

function [Phi, SigmaPhi, alpha, Theta, thetaCovDiag, keepList] = pruneTemporalBases( ...
        Phi, SigmaPhi, alpha, Theta, thetaCovDiag, keepList, threshold)
    if threshold <= 0 || numel(alpha) <= 1
        return;
    end
    relevance = 1 ./ max(alpha, realmin);
    maxRel = max(relevance);
    keep = relevance >= maxRel * threshold;
    if ~any(keep)
        [~, idx] = max(relevance);
        keep(idx) = true;
    end
    if all(keep)
        return;
    end
    Phi = Phi(keep, :);
    SigmaPhi = SigmaPhi(keep, keep);
    alpha = alpha(keep);
    Theta = Theta(:, keep);
    thetaCovDiag = thetaCovDiag(:, keep);
    keepList = keepList(keep);
end

function [A, F, gamma, Theta, thetaCovDiag, keepList] = pruneSpatialBases( ...
        A, F, gamma, Theta, thetaCovDiag, keepList, threshold)
    if threshold <= 0 || numel(gamma) <= 1
        return;
    end
    relevance = 1 ./ max(gamma, realmin);
    maxRel = max(relevance);
    keep = relevance >= maxRel * threshold;
    if ~any(keep)
        [~, idx] = max(relevance);
        keep(idx) = true;
    end
    if all(keep)
        return;
    end
    A = A(:, keep);
    F = F(:, keep);
    gamma = gamma(keep);
    Theta = Theta(keep, :);
    thetaCovDiag = thetaCovDiag(keep, :);
    keepList = keepList(keep);
end
