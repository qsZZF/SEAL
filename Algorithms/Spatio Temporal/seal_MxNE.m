function [Source_MxNE, Param_MxNE] = seal_MxNE(Data, L, varargin)
%SEAL_MXNE Mixed-Norm Estimate for M/EEG source imaging.
%   [Source_MxNE, Param_MxNE] = seal_MxNE(Data, L, 'ParameterName', Value, ...)
%
%   This implementation follows the L21 MxNE formulation in
%   Gramfort, Kowalski, and Hamalainen (2012):
%
%       min_X 0.5 * ||M - G*X||_F^2 + lambda * sum_g ||X_g||_F
%
%   where each group g is one source location across all time samples. For
%   free-orientation leadfields, a group contains all orientations at that
%   location and all time samples.
%
%   The public interface keeps the previous SEAL/MNE-Python style options.
%   In particular, Alpha remains a percentage in [0, 100] of lambda_max:
%
%       lambda = (Alpha / 100) * max_g ||G_g' * M||_F
%
%   Optional inputs:
%       Alpha             Regularization percentage of lambda_max. Default: 15.
%       NumOrientations   Orientations per source location. Default: 1.
%       NoiseCovariance   Channel noise covariance. Default: eye(Nchannels).
%       DepthBiasComp     Struct with fields exponent and limit. Default: 0.5, 10.
%       Weights           Optional external source weights/vector. Default: [].
%       WeightsMin        Threshold for external weights. Default: 0.
%       Solver            'as_fista' or 'fista'. Default: 'as_fista'.
%       ActiveSetSize     Locations added per AS-FISTA pass. Default: 10.
%       SolverMaxIter     Maximum FISTA iterations. Default: 3000.
%       SolverTolerance   Duality-gap tolerance. Default: 1e-4.
%       Debias            Refit amplitudes on selected support. Default: false.
%       DebiasMaxIter     Max iterations for debiasing. Default: 1000.
%       DebiasTolerance   Debiasing tolerance. Default: 1e-7.
%       PCAWhitening      true/false or numeric PCA rank. Default: true.
%       PickOri           'none' returns magnitude per location for nd>1;
%                         'vector' returns all orientation components. Default: 'none'.
%       Verbose           Print progress information. Default: true.
%       DualGapFrequency  Duality-gap check interval. Default: 10.

%% Input parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = false;

defaultAlpha = 15.0;
if size(L, 1) > 0
    defaultNoiseCov = eye(size(L, 1));
else
    defaultNoiseCov = [];
end
defaultDepthBiasComp = struct('exponent', 0.5, 'limit', 10.0);

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

addParameter(p, 'Alpha', defaultAlpha, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 100);
addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'NoiseCovariance', defaultNoiseCov, @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'DepthBiasComp', defaultDepthBiasComp, @(x) isempty(x) || isstruct(x));
addParameter(p, 'Weights', [], @(x) isempty(x) || isnumeric(x) || isstruct(x));
addParameter(p, 'WeightsMin', 0.0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'Solver', 'as_fista', @(x) ischar(x) || isstring(x));
addParameter(p, 'ActiveSetSize', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'SolverMaxIter', 3000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'SolverTolerance', 1e-4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Debias', false, @isLogicalLike);
addParameter(p, 'DebiasMaxIter', 1000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'DebiasTolerance', 1e-7, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PCAWhitening', true, @(x) isLogicalLike(x) || (isnumeric(x) && isscalar(x) && x >= 1));
addParameter(p, 'PickOri', 'none', @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose', true, @isLogicalLike);
addParameter(p, 'DualGapFrequency', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);

parse(p, Data, L, varargin{:});

Data_orig = double(p.Results.Data);
L_orig = double(p.Results.L);
alpha_percent = double(p.Results.Alpha);
nd = round(double(p.Results.NumOrientations));
C_noise = p.Results.NoiseCovariance;
depth_params = p.Results.DepthBiasComp;
user_weights_input = p.Results.Weights;
user_weights_min = double(p.Results.WeightsMin);
solver_name = normalizeSolverName(p.Results.Solver);
active_set_size = max(1, round(double(p.Results.ActiveSetSize)));
solver_max_iter = max(1, round(double(p.Results.SolverMaxIter)));
solver_tol = double(p.Results.SolverTolerance);
perform_debias = localToLogical(p.Results.Debias);
debias_max_iter = max(1, round(double(p.Results.DebiasMaxIter)));
debias_tol = double(p.Results.DebiasTolerance);
pca_whitening_opt = p.Results.PCAWhitening;
pick_ori_opt = lower(char(string(p.Results.PickOri)));
verbose = localToLogical(p.Results.Verbose);
dual_gap_frequency = max(1, round(double(p.Results.DualGapFrequency)));
local_check_option(pick_ori_opt, {'none', 'vector'});

[Nchannels, Ntimepoints] = size(Data_orig);
[Nchan_L, Nsources_total] = size(L_orig);
if Nchannels ~= Nchan_L
    error('seal_MxNE:ChannelMismatch', ...
        'Data rows (%d) and leadfield rows (%d) must match.', Nchannels, Nchan_L);
end
if nd < 1 || mod(Nsources_total, nd) ~= 0
    error('seal_MxNE:NDMismatch', ...
        'Leadfield columns (%d) must be divisible by NumOrientations (%d).', ...
        Nsources_total, nd);
end
if ~all(isfinite(Data_orig(:))) || ~all(isfinite(L_orig(:)))
    error('seal_MxNE:InvalidInput', 'Data and L must contain only finite values.');
end
nSources_loc = Nsources_total / nd;

Param_MxNE = struct();
Param_MxNE.OptionsPassed = p.Results;
Param_MxNE.OptionsPassed.Debias = perform_debias;
Param_MxNE.OptionsPassed.Verbose = verbose;
Param_MxNE.OptionsPassed.Solver = solver_name;

%% 1. Prepare data and gain matrix
[Data_w, L_w, whiteningInfo] = local_whiten_data(Data_orig, L_orig, C_noise, verbose);
Param_MxNE.WhiteningMatrix = whiteningInfo.WhiteningMatrix;
Param_MxNE.NoiseWhitening = whiteningInfo;

[M_processed, Vh_pca, pcaInfo] = local_temporal_pca(Data_w, pca_whitening_opt, ...
    whiteningInfo.Rank, verbose);
Param_MxNE.PCA = pcaInfo;

[source_scaling, scalingInfo] = local_source_scaling( ...
    L_w, nSources_loc, nd, depth_params, user_weights_input);
Param_MxNE.SourceWeightingUsedForGain = source_scaling;
Param_MxNE.SourceScalingInfo = scalingInfo;

G_eff = bsxfun(@times, L_w, source_scaling');

[G_eff, mask_active_sources, mask_active_locations] = local_apply_weight_min( ...
    G_eff, source_scaling, nSources_loc, nd, user_weights_min, user_weights_input, verbose);
Param_MxNE.ActiveSourceMask_Initial = mask_active_sources;
Param_MxNE.ActiveLocationMask_Initial = mask_active_locations;

if mod(size(G_eff, 2), nd) ~= 0
    error('seal_MxNE:MaskBlockMismatch', ...
        'Active source mask does not respect orientation-block structure.');
end
n_locs_eff = size(G_eff, 2) / nd;

%% 2. Convert Alpha to lambda
GtM = G_eff' * M_processed;
[lambda_max, lambda_max_by_group] = local_l21_dual_norm(GtM, nd);
if ~isfinite(lambda_max) || lambda_max < eps
    lambda_max = 1.0;
end
lambda = (alpha_percent / 100.0) * lambda_max;

Param_MxNE.AlphaSolver = alpha_percent;
Param_MxNE.AlphaMaxScalingFactor = lambda_max / 100.0;
Param_MxNE.LambdaMax = lambda_max;
Param_MxNE.Lambda = lambda;
Param_MxNE.LambdaMaxByGroup = lambda_max_by_group;
Param_MxNE.Solver = solver_name;

%% 3. Solve L21 MxNE with FISTA / active-set FISTA
if verbose
    fprintf('MxNE: Solver=%s, Alpha=%.4g%%, Lambda=%.4g, nd=%d\n', ...
        solver_name, alpha_percent, lambda, nd);
end

switch solver_name
    case 'as_fista'
        [Z_reduced_proc, solverInfo] = local_as_fista_l21( ...
            M_processed, G_eff, lambda, nd, solver_max_iter, solver_tol, ...
            active_set_size, dual_gap_frequency, verbose);
    case 'fista'
        [Z_reduced_proc, solverInfo] = local_fista_l21( ...
            M_processed, G_eff, lambda, nd, solver_max_iter, solver_tol, ...
            dual_gap_frequency, verbose, []);
    otherwise
        error('seal_MxNE:InvalidSolver', 'Unsupported Solver: %s.', solver_name);
end

Param_MxNE.Converged = solverInfo.Converged;
Param_MxNE.ConvergedAt = solverInfo.Iterations;
Param_MxNE.E_path = solverInfo.Objective;
Param_MxNE.DualGapPath = solverInfo.DualGap;
Param_MxNE.SolverInfo = solverInfo;

active_set_reduced = local_compute_active_mask(Z_reduced_proc, nd);
if ~any(active_set_reduced)
    warning('seal_MxNE:NoActiveSourcesFound', ...
        'No active sources found. Alpha may be too high.');
    Source_MxNE = zerosForOutput(Nsources_total, nSources_loc, Ntimepoints, nd, pick_ori_opt);
    Param_MxNE.ActiveSetIndices_Solver = [];
    Param_MxNE.ActiveSetIndices_Locations = [];
    Param_MxNE.ActiveSetIndices_Components = [];
    return;
end

%% 4. Optional debiasing on selected support
Z_active_proc = Z_reduced_proc(active_set_reduced, :);
if perform_debias
    if verbose
        fprintf('MxNE: Refitting amplitudes on %d active component(s)...\n', ...
            size(Z_active_proc, 1));
    end
    G_active = G_eff(:, active_set_reduced);
    [Z_active_proc, debiasInfo] = local_debias_support_ls( ...
        M_processed, G_active, Z_active_proc, verbose);
    Param_MxNE.DebiasingWeights = [];
    Param_MxNE.DebiasingInfo = debiasInfo;
else
    Param_MxNE.DebiasingWeights = [];
    Param_MxNE.DebiasingInfo = struct('Applied', false, ...
        'Method', 'none', 'MaxIter', debias_max_iter, 'Tolerance', debias_tol);
end

Z_reduced_proc(:) = 0;
Z_reduced_proc(active_set_reduced, :) = Z_active_proc;

%% 5. Reconstruct full time series if temporal PCA was used
Z_reduced_full_time = Z_reduced_proc;
if ~isempty(Vh_pca)
    Z_reduced_full_time = Z_reduced_proc * Vh_pca;
end

%% 6. Map back to original source space and original source units
source_scaling_reduced = source_scaling(mask_active_sources);
X_reduced_full_time = bsxfun(@times, Z_reduced_full_time, source_scaling_reduced);

X_full_components = zeros(Nsources_total, Ntimepoints);
X_full_components(mask_active_sources, :) = X_reduced_full_time;

Param_MxNE.ActiveSetIndices_Solver = find(active_set_reduced);
Param_MxNE.ActiveSetIndices_Components = find(any(X_full_components, 2));
if nd > 1
    Param_MxNE.ActiveSetIndices_Locations = unique(ceil(Param_MxNE.ActiveSetIndices_Components / nd));
else
    Param_MxNE.ActiveSetIndices_Locations = Param_MxNE.ActiveSetIndices_Components;
end
Param_MxNE.SourceComponentSize = size(X_full_components);

%% 7. Format output
if nd > 1 && strcmpi(pick_ori_opt, 'none')
    Source_MxNE = local_collapse_orientations(X_full_components, nSources_loc, nd);
else
    Source_MxNE = X_full_components;
end

end

%% Solver helpers

function [X, info] = local_as_fista_l21(M, G, lambda, nd, max_iter, tol, active_set_size, gap_freq, verbose)
    [~, Ns] = size(G);
    nLoc = Ns / nd;
    X = zeros(Ns, size(M, 2));

    [lambda_max, score0] = local_l21_dual_norm(G' * M, nd);
    if lambda >= lambda_max || lambda_max < eps
        [gap0, obj0, ~] = local_dgap_l21(M, G, X, lambda, nd);
        info = emptySolverInfo(true, 0, obj0, gap0, 'as_fista');
        info.ActiveSetSizes = 0;
        return;
    end

    [~, order] = sort(score0, 'descend');
    nInitial = min(active_set_size, max(1, sum(score0 > lambda)));
    active_locs = false(nLoc, 1);
    active_locs(order(1:nInitial)) = true;

    total_iter = 0;
    outer_iter = 0;
    obj_path = [];
    gap_path = [];
    active_size_path = [];
    converged = false;
    last_gap = inf;
    last_obj = inf;

    if verbose
        fprintf('MxNE AS-FISTA: initial active locations=%d/%d\n', ...
            sum(active_locs), nLoc);
    end

    while total_iter < max_iter
        outer_iter = outer_iter + 1;
        active_comp = locsToComponents(active_locs, nd);
        remaining_iter = max_iter - total_iter;
        if all(active_locs)
            inner_iter_budget = remaining_iter;
        else
            inner_iter_budget = min(remaining_iter, max(50, ceil(max_iter / 10)));
        end

        [X_active, innerInfo] = local_fista_l21( ...
            M, G(:, active_comp), lambda, nd, inner_iter_budget, tol, ...
            gap_freq, verbose, X(active_comp, :));

        X(:) = 0;
        X(active_comp, :) = X_active;
        total_iter = total_iter + innerInfo.Iterations;

        [full_gap, full_obj, ~] = local_dgap_l21(M, G, X, lambda, nd);
        last_gap = full_gap;
        last_obj = full_obj;
        obj_path = [obj_path; innerInfo.Objective(:); full_obj]; %#ok<AGROW>
        gap_path = [gap_path; innerInfo.DualGap(:); full_gap]; %#ok<AGROW>
        active_size_path = [active_size_path; sum(active_locs)]; %#ok<AGROW>

        if verbose
            fprintf('  AS-FISTA pass %d | total iter=%d | active locs=%d | full gap=%.3e\n', ...
                outer_iter, total_iter, sum(active_locs), full_gap);
        end

        if full_gap <= tol
            converged = true;
            break;
        end

        R = M - G * X;
        [~, scores] = local_l21_dual_norm(G' * R, nd);
        scores(active_locs) = -inf;
        violations = scores - lambda;
        candidate_locs = find(violations > max(tol, 1e-10) * max(1, lambda));

        if isempty(candidate_locs)
            converged = full_gap <= tol || ...
                (innerInfo.Converged && full_gap <= max(tol, 1e-8) * max(1, abs(full_obj)));
            if converged || total_iter >= max_iter
                break;
            end
            % No new KKT violators were found; keep refining the current
            % active set with the remaining iteration budget.
            continue;
        end

        [~, add_order_local] = sort(violations(candidate_locs), 'descend');
        nAdd = min(active_set_size, numel(candidate_locs));
        add_locs = candidate_locs(add_order_local(1:nAdd));
        active_locs(add_locs) = true;

        if all(active_locs)
            % One final full pass is enough; the next loop solves the full problem.
            if verbose
                fprintf('  AS-FISTA active set reached the full source space.\n');
            end
        end
    end

    if ~converged
        warning('seal_MxNE:SolverMaxIter', ...
            'AS-FISTA reached the iteration limit (%d) with full duality gap %.3e.', ...
            max_iter, last_gap);
    end

    info = emptySolverInfo(converged, total_iter, last_obj, last_gap, 'as_fista');
    info.Objective = obj_path;
    info.DualGap = gap_path;
    info.ActiveSetSizes = active_size_path;
    info.OuterIterations = outer_iter;
    info.ActiveLocations = find(active_locs);
end

function [X, info] = local_fista_l21(M, G, lambda, nd, max_iter, tol, gap_freq, verbose, X0)
    [~, Ns] = size(G);
    if isempty(X0)
        X = zeros(Ns, size(M, 2));
    else
        X = X0;
    end
    Y = X;
    t = 1.0;
    L_lip = local_lipschitz_squared(G);
    step = 1.0 / L_lip;

    obj_path = [];
    gap_path = [];
    converged = false;
    last_obj = inf;
    last_gap = inf;
    last_iter = 0;

    if verbose
        fprintf('  FISTA: MaxIter=%d, Tol=%.1e, L=%.3e, Ncomp=%d\n', ...
            max_iter, tol, L_lip, Ns);
    end

    for iter = 1:max_iter
        last_iter = iter;
        X_old = X;

        Grad = G' * (G * Y - M);
        X = local_prox_l21(Y - step * Grad, lambda * step, nd);

        t_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * t^2));
        Y = X + ((t - 1.0) / t_new) * (X - X_old);
        t = t_new;

        check_gap = mod(iter, gap_freq) == 0 || iter == 1 || iter == max_iter;
        if check_gap
            [last_gap, last_obj, ~] = local_dgap_l21(M, G, X, lambda, nd);
            obj_path = [obj_path; last_obj]; %#ok<AGROW>
            gap_path = [gap_path; last_gap]; %#ok<AGROW>

            if verbose
                fprintf('    FISTA iter %4d | pobj=%.4e | dgap=%.3e\n', ...
                    iter, last_obj, last_gap);
            end

            if last_gap <= tol
                converged = true;
                break;
            end
        end

        rel_change = norm(X(:) - X_old(:)) / max(norm(X_old(:)), eps);
        if iter > 1 && rel_change < tol * 1e-2
            [last_gap, last_obj, ~] = local_dgap_l21(M, G, X, lambda, nd);
            obj_path = [obj_path; last_obj]; %#ok<AGROW>
            gap_path = [gap_path; last_gap]; %#ok<AGROW>
            converged = last_gap <= max(tol, tol * max(1, abs(last_obj)));
            if converged
                break;
            end
        end
    end

    info = emptySolverInfo(converged, last_iter, last_obj, last_gap, 'fista');
    info.Objective = obj_path;
    info.DualGap = gap_path;
    info.Lipschitz = L_lip;
end

function X = local_prox_l21(Y, beta, nd)
    if beta <= 0
        X = Y;
        return;
    end
    [Ntotal, T] = size(Y);
    if mod(Ntotal, nd) ~= 0
        error('seal_MxNE:GroupMismatch', ...
            'Rows (%d) must be divisible by NumOrientations (%d).', Ntotal, nd);
    end
    nLoc = Ntotal / nd;
    Yt = reshape(Y, nd, nLoc, T);
    group_norms = sqrt(sum(sum(Yt .^ 2, 1), 3));
    shrink = max(1 - beta ./ max(group_norms, eps), 0);
    Xt = bsxfun(@times, Yt, shrink);
    X = reshape(Xt, Ntotal, T);
end

function [gap, p_obj, d_obj] = local_dgap_l21(M, G, X, lambda, nd)
    GX = G * X;
    R = M - GX;
    nR2 = sum(R(:) .^ 2);
    p_obj = 0.5 * nR2 + lambda * local_l21_norm(X, nd);

    [dual_norm_l2inf, ~] = local_l21_dual_norm(G' * R, nd);
    if lambda <= 0
        scaling = 0;
    else
        scaling = min(lambda / max(dual_norm_l2inf, eps), 1.0);
    end
    d_obj = (scaling - 0.5 * scaling^2) * nR2 + scaling * sum(sum(R .* GX));
    gap = p_obj - d_obj;
    if ~isfinite(gap)
        gap = p_obj;
    end
    gap = max(gap, 0);
end

function val = local_l21_norm(X, nd)
    if isempty(X)
        val = 0;
        return;
    end
    if nd > 1
        Xt = reshape(X, nd, [], size(X, 2));
        val = sum(reshape(sqrt(sum(sum(Xt .^ 2, 1), 3)), [], 1));
    else
        val = sum(sqrt(sum(X .^ 2, 2)));
    end
end

function [max_norm, group_norms] = local_l21_dual_norm(GtR, nd)
    if isempty(GtR)
        max_norm = 0;
        group_norms = [];
        return;
    end
    if mod(size(GtR, 1), nd) ~= 0
        error('seal_MxNE:GroupMismatch', ...
            'Rows (%d) must be divisible by NumOrientations (%d).', size(GtR, 1), nd);
    end
    if nd > 1
        Zt = reshape(GtR, nd, [], size(GtR, 2));
        group_norms = reshape(sqrt(sum(sum(Zt .^ 2, 1), 3)), [], 1);
    else
        group_norms = sqrt(sum(GtR .^ 2, 2));
    end
    max_norm = max(group_norms);
end

function mask = local_compute_active_mask(X, nd)
    if isempty(X)
        mask = false(0, 1);
        return;
    end
    if nd == 1
        mask = any(abs(X) > 0, 2);
        return;
    end
    nLoc = size(X, 1) / nd;
    active_locs = false(nLoc, 1);
    for i = 1:nLoc
        idx = (i - 1) * nd + (1:nd);
        active_locs(i) = any(any(abs(X(idx, :)) > 0));
    end
    mask = repelem(active_locs, nd);
end

function comp_mask = locsToComponents(loc_mask, nd)
    comp_mask = repelem(loc_mask(:), nd);
end

function info = emptySolverInfo(converged, iter, obj, gap, solver)
    info = struct();
    info.Solver = solver;
    info.Converged = converged;
    info.Iterations = iter;
    info.Objective = obj;
    info.DualGap = gap;
end

%% Preprocessing helpers

function [Data_w, L_w, info] = local_whiten_data(Data, L, C_noise, verbose)
    nChannels = size(Data, 1);
    Data_w = Data;
    L_w = L;
    W = eye(nChannels);
    info = struct( ...
        'Applied', false, ...
        'Reason', 'identity_or_empty', ...
        'Rank', nChannels, ...
        'ConditionNumber', 1, ...
        'WhiteningMatrix', W);

    if isempty(C_noise)
        return;
    end
    if size(C_noise, 1) ~= size(C_noise, 2)
        error('seal_MxNE:InvalidNoiseCovariance', 'NoiseCovariance must be square.');
    end
    if size(C_noise, 1) ~= nChannels
        error('seal_MxNE:NoiseCovarianceMismatch', ...
            'NoiseCovariance size (%d) must match data channels (%d).', ...
            size(C_noise, 1), nChannels);
    end

    C_noise = full((double(C_noise) + double(C_noise)') * 0.5);
    if isIdentityLike(C_noise)
        return;
    end

    [U, S] = svd(C_noise, 'econ');
    eig_vals = diag(S);
    s_max = max(eig_vals);
    if isempty(eig_vals) || ~isfinite(s_max) || s_max <= 0
        error('seal_MxNE:InvalidNoiseCovariance', ...
            'NoiseCovariance must have a positive finite eigenvalue.');
    end
    rank_tol = max(size(C_noise)) * eps(s_max);
    rank_noise = sum(eig_vals > rank_tol);
    floor_val = s_max * 1e-12;
    eig_vals(eig_vals < floor_val) = floor_val;

    W = diag(1 ./ sqrt(eig_vals)) * U';
    Data_w = W * Data;
    L_w = W * L;

    info.Applied = true;
    info.Reason = 'noise_covariance';
    info.Rank = rank_noise;
    info.ConditionNumber = s_max / min(eig_vals);
    info.WhiteningMatrix = W;
    if verbose
        fprintf('MxNE: Applied noise whitening, rank=%d, cond=%.3g\n', ...
            rank_noise, info.ConditionNumber);
    end
end

function [M_proc, Vh_pca, info] = local_temporal_pca(M, pca_opt, rank_whitener, verbose)
    M_proc = M;
    Vh_pca = [];
    info = struct('Applied', false, 'Rank', size(M, 2), 'OriginalTimePoints', size(M, 2));

    if isLogicalLike(pca_opt)
        do_pca = localToLogical(pca_opt);
        requested_rank = rank_whitener;
    else
        do_pca = true;
        requested_rank = round(double(pca_opt));
    end

    if ~do_pca
        return;
    end

    [U, S, V] = svd(M, 'econ');
    max_rank = size(S, 1);
    requested_rank = min([requested_rank, max_rank, size(M, 2)]);
    if requested_rank >= size(M, 2)
        return;
    end

    M_proc = U(:, 1:requested_rank) * S(1:requested_rank, 1:requested_rank);
    Vh_pca = V(:, 1:requested_rank)';
    info.Applied = true;
    info.Rank = requested_rank;
    if verbose
        fprintf('MxNE: Temporal PCA reduced time samples from %d to %d.\n', ...
            size(M, 2), requested_rank);
    end
end

function [scaling, info] = local_source_scaling(L_w, nSources_loc, nd, depth_params, user_weights_input)
    loc_scaling = ones(nSources_loc, 1);
    info = struct();
    info.DepthApplied = false;
    info.UserWeightsApplied = false;
    info.LocationScaling = loc_scaling;

    if ~isempty(depth_params) && isfield(depth_params, 'exponent')
        exponent = double(depth_params.exponent);
        source_norms = zeros(nSources_loc, 1);
        for i = 1:nSources_loc
            idx = (i - 1) * nd + (1:nd);
            source_norms(i) = norm(L_w(:, idx), 'fro');
        end
        source_norms(source_norms < eps) = eps;
        depth_scaling = source_norms .^ (-exponent);
        if isfield(depth_params, 'limit') && ~isempty(depth_params.limit) && depth_params.limit > 0
            limit_val = min(depth_scaling) * (double(depth_params.limit) ^ 2);
            depth_scaling(depth_scaling > limit_val) = limit_val;
        end
        loc_scaling = loc_scaling .* depth_scaling;
        info.DepthApplied = exponent ~= 0;
        info.SourceNorms = source_norms;
        info.DepthScaling = depth_scaling;
    end

    user_scaling = [];
    if ~isempty(user_weights_input)
        if isstruct(user_weights_input) && isfield(user_weights_input, 'data')
            user_scaling = abs(max(user_weights_input.data, [], 2));
        elseif isnumeric(user_weights_input) && isvector(user_weights_input)
            user_scaling = abs(user_weights_input(:));
        end
        if ~isempty(user_scaling) && numel(user_scaling) == nSources_loc
            if max(user_scaling) > eps
                user_scaling = user_scaling / max(user_scaling);
            end
            loc_scaling = loc_scaling .* user_scaling;
            info.UserWeightsApplied = true;
            info.UserScaling = user_scaling;
        else
            warning('seal_MxNE:UserWeightsDim', ...
                'User weights ignored because their length does not match source locations.');
        end
    end

    info.LocationScaling = loc_scaling;
    scaling = repelem(loc_scaling, nd);
end

function [G_eff, mask_comp, mask_loc] = local_apply_weight_min( ...
    G_eff, source_scaling, nSources_loc, nd, weights_min, user_weights_input, verbose)
    mask_loc = true(nSources_loc, 1);
    if weights_min > 0 && ~isempty(user_weights_input)
        loc_scaling = source_scaling(1:nd:end);
        threshold = weights_min * (max(loc_scaling) + eps);
        mask_loc = loc_scaling > threshold;
        if ~any(mask_loc)
            error('seal_MxNE:NoSourcesAfterWeightMin', ...
                'No sources remain after applying WeightsMin.');
        end
        if verbose
            fprintf('MxNE: Reduced source space to %d/%d location(s) via WeightsMin.\n', ...
                sum(mask_loc), nSources_loc);
        end
    end
    mask_comp = locsToComponents(mask_loc, nd);
    G_eff = G_eff(:, mask_comp);
end

function L2 = local_lipschitz_squared(G)
    if isempty(G)
        L2 = 1.0;
        return;
    end
    try
        if numel(G) <= 2e6 || size(G, 2) <= 1000
            sigma = norm(G, 2);
        else
            sigma = normest(G);
        end
    catch
        sigma = norm(G, 2);
    end
    L2 = max(sigma ^ 2, eps);
end

%% Output and debias helpers

function [X_refit, info] = local_debias_support_ls(M_data, G_active_comp, X_active_comp, verbose)
% Refit amplitudes on the selected support by unpenalized least squares.
    info = struct('Applied', false, 'Method', 'support_least_squares', ...
        'Rank', 0, 'Ridge', 0, 'RelativeResidualBefore', NaN, ...
        'RelativeResidualAfter', NaN);
    X_refit = X_active_comp;
    if isempty(X_active_comp)
        return;
    end

    data_norm = max(norm(M_data, 'fro'), eps);
    info.RelativeResidualBefore = norm(M_data - G_active_comp * X_active_comp, 'fro') / data_norm;
    info.Rank = rank(G_active_comp);

    try
        if info.Rank == size(G_active_comp, 2)
            X_refit = G_active_comp \ M_data;
        else
            gram = G_active_comp' * G_active_comp;
            rhs = G_active_comp' * M_data;
            ridge_scale = max(trace(gram) / max(size(gram, 1), 1), 1);
            info.Ridge = ridge_scale * 1e-12;
            X_refit = (gram + info.Ridge * eye(size(gram))) \ rhs;
        end
    catch
        gram = G_active_comp' * G_active_comp;
        rhs = G_active_comp' * M_data;
        ridge_scale = max(trace(gram) / max(size(gram, 1), 1), 1);
        info.Ridge = ridge_scale * 1e-8;
        X_refit = (gram + info.Ridge * eye(size(gram))) \ rhs;
    end

    info.RelativeResidualAfter = norm(M_data - G_active_comp * X_refit, 'fro') / data_norm;
    info.Applied = true;
    if verbose
        fprintf('  Debiasing LS: rank=%d/%d, relres %.3g -> %.3g\n', ...
            info.Rank, size(G_active_comp, 2), ...
            info.RelativeResidualBefore, info.RelativeResidualAfter);
    end
end

function Source = local_collapse_orientations(X_full_components, nSources_loc, nd)
    Source = zeros(nSources_loc, size(X_full_components, 2));
    for i = 1:nSources_loc
        idx = (i - 1) * nd + (1:nd);
        Source(i, :) = sqrt(sum(X_full_components(idx, :) .^ 2, 1));
    end
end

function Source = zerosForOutput(Nsources_total, nSources_loc, nTime, nd, pick_ori_opt)
    if nd > 1 && strcmpi(pick_ori_opt, 'none')
        Source = zeros(nSources_loc, nTime);
    else
        Source = zeros(Nsources_total, nTime);
    end
end

%% Small utilities

function solver_name = normalizeSolverName(value)
    solver_name = lower(strrep(char(string(value)), '-', '_'));
    local_check_option(solver_name, {'as_fista', 'fista'});
end

function local_check_option(value, allowed)
    if ~ismember(value, allowed)
        error('seal_MxNE:InvalidOption', ...
            'Invalid value: %s. Allowed: %s', value, strjoin(allowed, ', '));
    end
end

function tf = isLogicalLike(value)
    tf = (islogical(value) && isscalar(value)) || ...
        (isnumeric(value) && isscalar(value)) || ...
        ischar(value) || ...
        (isstring(value) && isscalar(value));
end

function tf = localToLogical(value)
    if islogical(value)
        tf = value;
    elseif isnumeric(value)
        tf = value ~= 0;
    else
        value = lower(strtrim(char(string(value))));
        tf = any(strcmp(value, {'true', '1', 'yes', 'on'}));
    end
    tf = logical(tf(1));
end

function tf = isIdentityLike(A)
    n = size(A, 1);
    delta = A - speye(n);
    tf = norm(delta, 'fro') <= 1e-12 * max(1, norm(A, 'fro'));
end
