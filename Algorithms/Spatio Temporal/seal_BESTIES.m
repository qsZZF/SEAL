function [Source_BESTIES, Param_BESTIES] = seal_BESTIES(Data, L, varargin)
%SEAL_BESTIES Bayesian electromagnetic spatio-temporal imaging.
%   [Source, Param] = seal_BESTIES(Data, L, Phi, VertConn, ...)
%   [Source, Param] = seal_BESTIES(Data, L, sourceSpaceInfo, ...)
%   [Source, Param] = seal_BESTIES(Data, L, 'Phi', Phi, 'VertConn', VertConn, ...)
%
%   BESTIES models source activity as S = W*Phi + E, with spatial MRF
%   priors on both the temporal-basis coefficients W and the residual source
%   activity E. If Phi is not supplied, temporal bases are generated from
%   the data by SVD. Pass 'NumTBFs', 0 or an explicit empty Phi to run the
%   spatial-only MRF variant.
%
%   Inputs:
%       Data     - Nchannels x Ntimepoints measured data.
%       L        - Nchannels x Nsources lead field.
%       Phi      - K x Ntimepoints temporal basis functions.
%       VertConn - Nlocations x Nlocations source-space adjacency matrix.
%
%   Optional name-value inputs:
%       'NumOrientations'         - Orientations per location. Default: 1.
%       'Phi'                     - Temporal basis functions.
%       'VertConn'                - Source-space adjacency matrix.
%       'SourceSpaceInfo'         - Struct containing .VertConn.
%       'NumTBFs'                 - SVD bases to generate when Phi is absent.
%                                   Default: 4.
%       'TemporalBasisMode'       - 'svd' or 'none'. Default: 'svd'.
%       'NoiseCovariance'         - Sensor noise covariance. Default: eye.
%       'PrewhitenNoise'          - Whiten Data/L by NoiseCovariance before
%                                   inference. Default: true.
%       'MaxIterations'           - Maximum VB iterations. Default: 100.
%       'Tolerance'               - Convergence tolerance. Default: 1e-3.
%       'LearnTau'                - Learn MRF spatial parameter. Default: true.
%       'InitialTau'              - Initial tau in [0, 1). Default: 0.98.
%       'TauBounds'               - Clamp learned tau. Default: [0 0.99].
%       'ComputeFreeEnergy'       - Use free energy as cost. Default: true.

    %% 1. Input parsing and validation
    [posPhi, posVertConn, nvArgs, inputInfo] = preparseBestiesInputs(varargin{:});
    inputInfo.ExplicitPhi = inputInfo.ExplicitPhi || hasNameValue(nvArgs, 'Phi');
    inputInfo.ExplicitVertConn = inputInfo.ExplicitVertConn || hasNameValue(nvArgs, 'VertConn');

    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'Phi', posPhi, @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'VertConn', posVertConn, @(x) isempty(x) || isnumeric(x) || islogical(x) || issparse(x));
    addParameter(p, 'SourceSpaceInfo', struct(), @isstruct);

    addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NumTBFs', 4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'TemporalBasisMode', 'svd', @(x) ischar(x) || isstring(x));
    addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'PrewhitenNoise', true, @(x) isLogicalLike(x));
    addParameter(p, 'NoiseObservationMatrix', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x)); % accepted for GUI compatibility
    addParameter(p, 'MaxIterations', 100, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Tolerance', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'LearnTau', true, @(x) isLogicalLike(x));
    addParameter(p, 'InitialTau', 0.98, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
    addParameter(p, 'TauBounds', [0 0.99], @(x) isnumeric(x) && numel(x) == 2 && x(1) >= 0 && x(2) < 1 && x(1) <= x(2));
    addParameter(p, 'ComputeFreeEnergy', true, @(x) isLogicalLike(x));

    parse(p, Data, L, nvArgs{:});
    opts = p.Results;

    opts.LearnTau = parseLogical(opts.LearnTau, 'LearnTau');
    opts.ComputeFreeEnergy = parseLogical(opts.ComputeFreeEnergy, 'ComputeFreeEnergy');
    opts.PrewhitenNoise = parseLogical(opts.PrewhitenNoise, 'PrewhitenNoise');
    opts.NumOrientations = round(opts.NumOrientations);
    opts.NumTBFs = round(opts.NumTBFs);
    opts.MaxIterations = round(opts.MaxIterations);
    opts.TauBounds = double(opts.TauBounds(:)).';

    if isempty(opts.VertConn) && isfield(opts.SourceSpaceInfo, 'VertConn')
        opts.VertConn = opts.SourceSpaceInfo.VertConn;
        inputInfo.ExplicitVertConn = true;
    end
    if isempty(opts.VertConn)
        error('seal_BESTIES:MissingVertConn', ...
            'BESTIES requires VertConn, either positionally, by name-value, or in SourceSpaceInfo.VertConn.');
    end

    Data_orig = double(opts.Data);
    L_orig = double(opts.L);
    [nChan, nSnap] = size(Data_orig);
    if size(L_orig, 1) ~= nChan
        error('seal_BESTIES:DataLeadfieldMismatch', ...
            'Data rows (%d) must match L rows (%d).', nChan, size(L_orig, 1));
    end

    nd = opts.NumOrientations;
    nSource = size(L_orig, 2);
    if mod(nSource, nd) ~= 0
        error('seal_BESTIES:NumOrientationsMismatch', ...
            'size(L,2) (%d) must be divisible by NumOrientations (%d).', nSource, nd);
    end
    nLoc = nSource / nd;

    Adj_loc = sanitizeAdjacency(opts.VertConn, nLoc);
    [Data_orig, L_orig, C_noise, noiseInfo] = prepareNoiseModel( ...
        Data_orig, L_orig, opts.NoiseCovariance, opts.PrewhitenNoise);
    opts.NoiseCovariance = C_noise;

    [Phi_basis, basisInfo] = resolveTemporalBasis(opts.Phi, Data_orig, opts, inputInfo.ExplicitPhi);
    opts.Phi = Phi_basis;
    if isempty(Phi_basis)
        K = 0;
    else
        K = size(Phi_basis, 1);
    end

    %% 2. Initialization
    Cost = zeros(opts.MaxIterations, 1);
    Cost_old = 0;
    S_old = zeros(nSource, nSnap);

    W = zeros(nSource, K); W1 = W;
    V = zeros(nSource, nSnap); V1 = V;
    keeplist = (1:K)';
    c = ones(K, 1);

    % Row-normalized first-order MRF neighborhood matrix.
    deg = full(sum(Adj_loc, 2));
    deg(deg == 0) = 1;
    N_loc = spdiags(1 ./ deg, 0, nLoc, nLoc) * Adj_loc;

    % Expand source-space graph to orientation space. This assumes columns
    % in L are ordered as [loc1_ori1, loc1_ori2, ..., loc2_ori1, ...].
    if nd > 1
        N = kron(N_loc, speye(nd));
    else
        N = N_loc;
    end

    tau = min(max(opts.InitialTau, opts.TauBounds(1)), opts.TauBounds(2));
    M = speye(nSource) - tau * N;
    LL = L_orig / M;

    tr_phi = 0;
    if K > 0
        tr_phi = trace(Phi_basis' * Phi_basis);
    end
    init_prec = trace(LL * LL') * (tr_phi + nSnap) / max(nSnap * trace(C_noise), eps);
    if ~isfinite(init_prec) || init_prec <= 0
        init_prec = 1;
    end
    gamma_prec = init_prec * ones(nSource, 1);
    alpha_prec = init_prec * ones(nSource, 1);

    alpha_norm = norm(alpha_prec);
    c = c .* alpha_norm;
    alpha_prec = alpha_prec ./ alpha_norm;

    %% 3. Variational Bayes loop
    fprintf('Running BESTIES (nd=%d, K=%d)...\n', nd, K);
    for iter = 1:opts.MaxIterations

        if ~isempty(Phi_basis)
            index = find(abs(1 ./ c) > max(abs(1 ./ c)) * 1e-3);
            c = c(index);
            keeplist = keeplist(index);
            Phi_basis = Phi_basis(index, :);
            W = W(:, index);
            W1 = W1(:, index);
            K = numel(index);
            C2 = cell(numel(c), 1);
        end

        if opts.LearnTau
           M = speye(nSource) - tau * N;
           LL = L_orig / M;
        end

        % W step
        if ~isempty(Phi_basis)
            TrW = zeros(nSource, 1);
            for k = numel(c):-1:1
                x_phi = Phi_basis(k, :) * Phi_basis(k, :)';
                C2{k} = LL .* repmat(1 ./ (c(k) * alpha_prec)', nChan, 1) * LL' + C_noise / x_phi;
                otherIdx = setdiff(1:K, k);
                residual = (Data_orig - L_orig * V) * Phi_basis(k, :)' - ...
                    L_orig * sum((W(:, otherIdx) .* repmat(Phi_basis(k, :) * Phi_basis(otherIdx, :)', nSource, 1)), 2);
                Kernel = repmat(1 ./ (c(k) * alpha_prec), 1, nChan) .* LL' / C2{k};
                W1(:, k) = Kernel * residual / x_phi;
                if opts.LearnTau
                    G = M \ Kernel * LL;
                    TrW = TrW + sum(N .* G', 2);
                end
            end
            W = M \ W1;
        end

        % V/E step
        C1 = LL .* repmat(1 ./ gamma_prec', nChan, 1) * LL' + C_noise;
        if ~isempty(Phi_basis)
            V1 = repmat(1 ./ gamma_prec, 1, nChan) .* LL' / C1 * (Data_orig - L_orig * W * Phi_basis);
        else
            V1 = repmat(1 ./ gamma_prec, 1, nChan) .* LL' / C1 * Data_orig;
        end
        E = M \ (repmat(1 ./ gamma_prec, 1, nChan) .* LL' / C1 * LL);
        V = M \ V1;

        % Gamma step, grouped by location when multiple orientations exist.
        rou = diag(LL' / C1 * LL);
        if nd > 1
            rou_loc = sum(reshape(rou, nd, []), 1)';
            V1_sq_loc = sum(reshape(sum(V1.^2, 2), nd, []), 1)';
            gamma_loc = sqrt(nSnap * rou_loc ./ max(V1_sq_loc, 1e-12));
            gamma_prec = reshape(repmat(gamma_loc', nd, 1), [], 1);
        else
            gamma_prec = sqrt(nSnap * rou ./ max(sum(V1.^2, 2), 1e-12));
        end
        gamma_prec = sanitizePrecision(gamma_prec);

        % Alpha step
        if ~isempty(Phi_basis)
            kesi = zeros(nSource, 1);
            for k = 1:numel(c)
                kesi = kesi + diag(LL' / C2{k} * LL) / c(k);
            end
            if nd > 1
                kesi_loc = sum(reshape(kesi, nd, []), 1)';
                W1_sq_loc = sum(reshape(sum(W1 .* repmat(c', nSource, 1) .* W1, 2), nd, []), 1)';
                alpha_loc = sqrt(kesi_loc ./ max(W1_sq_loc, 1e-12));
                alpha_prec = reshape(repmat(alpha_loc', nd, 1), [], 1);
            else
                alpha_prec = sqrt(kesi ./ max(sum(W1 .* repmat(c', nSource, 1) .* W1, 2), 1e-12));
            end
            alpha_prec = sanitizePrecision(alpha_prec);

            % Resolve the c/alpha scaling ambiguity.
            alpha_norm = norm(alpha_prec);
            if alpha_norm > 0 && isfinite(alpha_norm)
                c = c .* alpha_norm;
                alpha_prec = alpha_prec ./ alpha_norm;
            end
            for k = 1:numel(c)
                omega = diag(LL' / C2{k} * LL)' * (1 ./ alpha_prec);
                c(k) = sqrt(omega / max(W1(:, k)' .* alpha_prec' * W1(:, k), 1e-12));
            end
            c = sanitizePrecision(c);
        end

        % Tau step
        if opts.LearnTau
             if ~isempty(Phi_basis)
                 W_mean = N * W;
             end
             V_mean = N * V;
             TrV = nSnap * sum(E' .* N, 2);

             if ~isempty(Phi_basis)
                 Num = sum(W .* repmat(c', nSource, 1) .* W_mean .* repmat(alpha_prec, 1, numel(c)), 2) - TrW + ...
                       sum(V .* V_mean .* repmat(gamma_prec, 1, nSnap), 2) - TrV;
                 Den = sum(repmat(c', nSource, 1) .* (W_mean .^ 2) .* repmat(alpha_prec, 1, numel(c)), 2) + ...
                       sum((V_mean .^ 2) .* repmat(gamma_prec, 1, nSnap), 2);
             else
                 Num = sum(V .* V_mean .* repmat(gamma_prec, 1, nSnap), 2) - TrV;
                 Den = sum((V_mean .^ 2) .* repmat(gamma_prec, 1, nSnap), 2);
             end
             tau = sum(Num) / max(sum(Den), 1e-12);
             tau = max(min(tau, opts.TauBounds(2)), opts.TauBounds(1));
        end

        if ~isempty(Phi_basis)
            S_est = W * Phi_basis + V;
        else
            S_est = V;
        end

        % Convergence check
        if ~opts.ComputeFreeEnergy
            MSE = norm(S_est - S_old, 'fro')^2 / max(norm(S_est, 'fro')^2, eps);
            S_old = S_est;
            Cost(iter) = MSE;
        else
            TOL_FE = 1e-14;
            alpha_fe = max(min(alpha_prec, 1 / TOL_FE), TOL_FE);
            gamma_fe = max(min(gamma_prec, 1 / TOL_FE), TOL_FE);

            temp = 0;
            logCWprior = 0;
            logCW = 0;
            if ~isempty(Phi_basis)
                for k = 1:K
                    temp = temp + c(k) * W1(:, k)' .* alpha_fe' * W1(:, k);
                    logCWprior = logCWprior + sum(log(c(k) * alpha_fe));

                    x_phi = Phi_basis(k, :) * Phi_basis(k, :)';
                    LALT = LL .* repmat(1 ./ (c(k) * alpha_fe)', nChan, 1) * LL';
                    logCW = logCW + robust_logdet(eye(nChan) / x_phi + LALT) + ...
                        nChan * log(x_phi) + sum(log(c(k) * alpha_fe));
                end
            end

            logCV = robust_logdet(eye(nChan) + LL .* repmat(1 ./ gamma_fe', nChan, 1) * LL') + sum(log(gamma_fe));
            logCVprior = sum(log(gamma_fe));

            FreeEnergy = -norm(Data_orig - L_orig * S_est, 'fro')^2 ...
                         - trace(V1' .* repmat(gamma_fe', nSnap, 1) * V1) ...
                         - temp ...
                         + (logCWprior - logCW) ...
                         + nSnap * (logCVprior - logCV);

            Cost(iter) = FreeEnergy;

            if iter > 1
                MSE = abs(FreeEnergy - Cost_old) / max(abs(FreeEnergy), eps);
            else
                MSE = inf;
            end
            Cost_old = FreeEnergy;
            if mod(iter, 10) == 0
                fprintf('BESTIES iteration: %d, tau: %.4f, Diff: %g\n', iter, tau, MSE);
            end
        end

        if abs(MSE) < opts.Tolerance
            fprintf('BESTIES Converged at iteration: %d, tau: %.4f, Diff: %g\n', iter, tau, MSE);
            Cost = Cost(1:iter);
            break;
        end
    end

    Source_BESTIES = S_est;
    Param_BESTIES.BasisCoefficients_W = W;
    Param_BESTIES.ResidualSource_E = V;
    Param_BESTIES.FinalPhi = Phi_basis;
    Param_BESTIES.TemporalBasisKeepList = keeplist;
    Param_BESTIES.TemporalBasisInfo = basisInfo;
    Param_BESTIES.SpatialPrecision_alpha = alpha_prec;
    Param_BESTIES.SourceErrorPrecision_gamma = gamma_prec;
    Param_BESTIES.TemporalPrecision_c = c;
    Param_BESTIES.FinalTau = tau;
    Param_BESTIES.CostFunction = Cost;
    Param_BESTIES.NoiseWhitening = noiseInfo;
    Param_BESTIES.InputParsing = inputInfo;
    Param_BESTIES.OptionsPassed = opts;
end

function [posPhi, posVertConn, nvArgs, inputInfo] = preparseBestiesInputs(varargin)
    posPhi = [];
    posVertConn = [];
    used = 0;
    inputInfo = struct( ...
        'ExplicitPhi', false, ...
        'ExplicitVertConn', false, ...
        'PositionalSourceSpaceInfo', false, ...
        'CallingStyle', 'name-value');

    if isempty(varargin)
        nvArgs = {};
        return;
    end

    first = varargin{1};
    if ~(ischar(first) || isstring(first))
        if isstruct(first)
            inputInfo.CallingStyle = 'sourceSpaceInfo';
            inputInfo.PositionalSourceSpaceInfo = true;
            if isfield(first, 'VertConn')
                posVertConn = first.VertConn;
                inputInfo.ExplicitVertConn = true;
            end
            used = 1;
        else
            inputInfo.CallingStyle = 'phi-vertconn';
            posPhi = first;
            inputInfo.ExplicitPhi = true;
            used = 1;
            if numel(varargin) >= 2 && ~(ischar(varargin{2}) || isstring(varargin{2}))
                posVertConn = varargin{2};
                inputInfo.ExplicitVertConn = true;
                used = 2;
            end
        end
    end
    nvArgs = varargin(used + 1:end);
end

function tf = hasNameValue(args, name)
    tf = false;
    for i = 1:2:numel(args)
        if ischar(args{i}) || isstring(args{i})
            if strcmpi(char(string(args{i})), name)
                tf = true;
                return;
            end
        end
    end
end

function tf = isLogicalLike(x)
    tf = (islogical(x) && isscalar(x)) || ...
         (isnumeric(x) && isscalar(x)) || ...
         ((ischar(x) || isstring(x)) && isscalar(string(x)));
end

function tf = parseLogical(x, name)
    if islogical(x) && isscalar(x)
        tf = x;
        return;
    end
    if isnumeric(x) && isscalar(x)
        tf = x ~= 0;
        return;
    end
    if ischar(x) || isstring(x)
        value = lower(strtrim(char(string(x))));
        switch value
            case {'true', 't', 'yes', 'y', 'on', '1'}
                tf = true;
                return;
            case {'false', 'f', 'no', 'n', 'off', '0'}
                tf = false;
                return;
        end
    end
    error('seal_BESTIES:InvalidLogical', '%s must be logical-compatible.', name);
end

function [Phi, info] = resolveTemporalBasis(Phi, Data, opts, explicitPhi)
    info = struct( ...
        'Generated', false, ...
        'Mode', char(string(opts.TemporalBasisMode)), ...
        'NumTBFsRequested', opts.NumTBFs, ...
        'NumTBFsUsed', 0);

    if isempty(Phi)
        if explicitPhi || opts.NumTBFs == 0 || strcmpi(info.Mode, 'none')
            Phi = [];
            return;
        end

        switch lower(info.Mode)
            case 'svd'
                [~, ~, V] = svd(Data, 'econ');
                k = min(opts.NumTBFs, size(V, 2));
                Phi = V(:, 1:k)';
                info.Generated = true;
                info.NumTBFsUsed = k;
            otherwise
                error('seal_BESTIES:UnsupportedTemporalBasisMode', ...
                    'TemporalBasisMode "%s" is not supported. Use "svd" or "none".', info.Mode);
        end
    else
        Phi = double(Phi);
        info.NumTBFsUsed = size(Phi, 1);
    end

    if ~isempty(Phi)
        if size(Phi, 2) ~= size(Data, 2)
            error('seal_BESTIES:PhiTimeMismatch', ...
                'Phi must have one column per time point: size(Phi,2)=%d, size(Data,2)=%d.', ...
                size(Phi, 2), size(Data, 2));
        end
        if any(~isfinite(Phi(:)))
            error('seal_BESTIES:InvalidPhi', 'Phi contains NaN or Inf.');
        end
        rowEnergy = sum(Phi .^ 2, 2);
        if any(rowEnergy <= eps(max(1, max(rowEnergy))))
            error('seal_BESTIES:InvalidPhi', 'Each row of Phi must have non-zero energy.');
        end
        info.NumTBFsUsed = size(Phi, 1);
    end
end

function Adj = sanitizeAdjacency(VertConn, nLoc)
    if ~isequal(size(VertConn), [nLoc, nLoc])
        error('seal_BESTIES:VertConnSize', ...
            'VertConn must be Nlocations x Nlocations, where Nlocations=size(L,2)/NumOrientations (%d).', nLoc);
    end
    Adj = double(sparse(VertConn ~= 0));
    Adj = Adj - spdiags(diag(Adj), 0, nLoc, nLoc);
    Adj = spones(Adj + Adj');
end

function [Data, L, C_noise, info] = prepareNoiseModel(Data, L, C_noise_in, prewhiten)
    nChan = size(Data, 1);
    info = struct( ...
        'Applied', false, ...
        'InputWasEmpty', isempty(C_noise_in), ...
        'Rank', nChan, ...
        'ConditionNumber', 1, ...
        'RegularizationFloor', 0);

    if isempty(C_noise_in)
        C_noise = eye(nChan);
        return;
    end
    if ~isequal(size(C_noise_in), [nChan, nChan])
        error('seal_BESTIES:NoiseCovarianceSize', ...
            'NoiseCovariance must be %d x %d.', nChan, nChan);
    end
    C_noise = (double(C_noise_in) + double(C_noise_in)') / 2;
    if any(~isfinite(C_noise(:)))
        error('seal_BESTIES:InvalidNoiseCovariance', 'NoiseCovariance contains NaN or Inf.');
    end

    identityScale = norm(C_noise - eye(nChan), 'fro') / max(norm(C_noise, 'fro'), eps);
    if ~prewhiten || identityScale < 1e-10
        return;
    end

    [U, S] = svd(C_noise, 'econ');
    s = diag(S);
    if isempty(s) || max(s) <= 0
        error('seal_BESTIES:InvalidNoiseCovariance', ...
            'NoiseCovariance must have at least one positive eigenvalue.');
    end
    sMax = max(s);
    floorVal = sMax * 1e-8;
    rankNoise = sum(s > max(numel(s) * eps(sMax), floorVal));
    s = max(s, floorVal);
    Wn = diag(1 ./ sqrt(s)) * U';
    Data = Wn * Data;
    L = Wn * L;
    C_noise = eye(nChan);

    info.Applied = true;
    info.Rank = rankNoise;
    info.ConditionNumber = sMax / min(s);
    info.RegularizationFloor = floorVal;
end

function x = sanitizePrecision(x)
    x = real(x);
    x(~isfinite(x) | x <= 0) = 1e12;
    minVal = min(x);
    if isfinite(minVal) && minVal > 0
        x(x > minVal * 1e12) = minVal * 1e12;
    end
end

function ld = robust_logdet(C)
% Computes log(det(C)) robustly to safely replace spm_logdet.
    C = (C + C') / 2;
    try
        R = chol(C);
        ld = 2 * sum(log(diag(R)));
    catch
        [~, S, ~] = svd(C);
        s_diag = diag(S);
        ld = sum(log(max(s_diag, 1e-12)));
    end
end
