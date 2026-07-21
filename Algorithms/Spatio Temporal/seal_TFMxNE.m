function [Source_TFMxNE, Param_TFMxNE] = seal_TFMxNE(Data, L, varargin)
%SEAL_TFMXNE Time-frequency mixed-norm estimate for M/EEG source imaging.
%   [Source, Param] = seal_TFMxNE(Data, L, 'ParameterName', Value, ...)
%
%   This implementation follows the TF-MxNE model in Gramfort et al. (2013):
%
%       min_Z 0.5 * ||M - G * Z * Phi^H||_F^2
%             + lambda_space * ||Z||_21 + lambda_time * ||Z||_1
%
%   where Phi is represented by a short-time Fourier / Gabor transform.
%   For free-orientation source spaces, both penalties are grouped by
%   source location across the orientation components.

%% Input parsing
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = false;

addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

addParameter(p, 'TimeStep', 4, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'WindowSize', 64, @(x) isnumeric(x) && isscalar(x) && x >= 4);
addParameter(p, 'SpatialReg', 50, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'TemporalReg', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'MaxIterations', 200, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Tolerance', 1e-4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Lipschitz', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'LipschitzMaxIter', 30, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NumOrientations', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'DepthBiasComp', struct('exponent', 0.8, 'limit', 10.0), @(x) isempty(x) || isstruct(x));
addParameter(p, 'DepthBiasExponent', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
addParameter(p, 'DepthBiasLimit', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'LooseOrientationWeight', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'PickOri', 'none', @(x) ischar(x) || isstring(x));
addParameter(p, 'Debias', false, @isLogicalLike);
addParameter(p, 'Verbose', true, @isLogicalLike);

parse(p, Data, L, varargin{:});
opts = p.Results;

M_orig = double(opts.Data);
G_orig = double(opts.L);
nd = round(double(opts.NumOrientations));
time_step = round(double(opts.TimeStep));
window_size = round(double(opts.WindowSize));
max_iter = round(double(opts.MaxIterations));
tol = double(opts.Tolerance);
spatial_reg_percent = double(opts.SpatialReg);
temporal_reg_percent = double(opts.TemporalReg);
pick_ori = lower(char(string(opts.PickOri)));
perform_debias = localToLogical(opts.Debias);
verbose = localToLogical(opts.Verbose);

if ~ismember(pick_ori, {'none', 'vector'})
    error('seal_TFMxNE:InvalidPickOri', 'PickOri must be ''none'' or ''vector''.');
end
if mod(window_size, 2) ~= 0
    error('seal_TFMxNE:InvalidWindowSize', 'WindowSize must be even.');
end
if time_step > window_size
    error('seal_TFMxNE:InvalidTimeStep', 'TimeStep must be <= WindowSize.');
end
if mod(window_size - time_step, 2) ~= 0
    error('seal_TFMxNE:InvalidFrameGeometry', ...
        'WindowSize - TimeStep must be even for centered STFT padding.');
end

[n_channels, n_times] = size(M_orig);
[n_lf_channels, n_sources_total] = size(G_orig);
if n_channels ~= n_lf_channels
    error('seal_TFMxNE:ChannelMismatch', ...
        'Data rows (%d) and leadfield rows (%d) must match.', n_channels, n_lf_channels);
end
if nd < 1 || mod(n_sources_total, nd) ~= 0
    error('seal_TFMxNE:NDMismatch', ...
        'Leadfield columns (%d) must be divisible by NumOrientations (%d).', ...
        n_sources_total, nd);
end
if ~all(isfinite(M_orig(:))) || ~all(isfinite(G_orig(:)))
    error('seal_TFMxNE:InvalidInput', 'Data and L must contain only finite values.');
end
n_locations = n_sources_total / nd;

Param_TFMxNE = struct();
Param_TFMxNE.OptionsPassed = opts;
Param_TFMxNE.OptionsPassed.Debias = perform_debias;
Param_TFMxNE.OptionsPassed.Verbose = verbose;

%% 1. Noise whitening and source scaling
[M_w, G_w, whitening_info] = local_whiten_data(M_orig, G_orig, opts.NoiseCovariance, verbose);
Param_TFMxNE.NoiseWhitening = whitening_info;
Param_TFMxNE.WhiteningMatrix = whitening_info.WhiteningMatrix;

[source_scaling, scaling_info] = local_source_scaling( ...
    G_w, n_locations, nd, opts.DepthBiasComp, opts.DepthBiasExponent, ...
    opts.DepthBiasLimit, double(opts.LooseOrientationWeight));
G_eff = bsxfun(@times, G_w, source_scaling(:)');
Param_TFMxNE.SourceScaling = source_scaling;
Param_TFMxNE.SourceScalingInfo = scaling_info;

%% 2. TF dictionary dimensions and regularization scale
tf_probe = local_stft(M_w(1:min(1, n_channels), :), window_size, time_step);
n_freq = size(tf_probe, 2);
n_frames = size(tf_probe, 3);
n_coeff = n_freq * n_frames;

GtM = G_eff' * M_w;
[lambda_max, lambda_max_by_group] = local_l21_dual_norm_time(GtM, nd);
if ~isfinite(lambda_max) || lambda_max < eps
    lambda_max = 1.0;
end
lambda_space = (spatial_reg_percent / 100.0) * lambda_max;
lambda_time = (temporal_reg_percent / 100.0) * lambda_max;

Param_TFMxNE.LambdaMax = lambda_max;
Param_TFMxNE.LambdaMaxByGroup = lambda_max_by_group;
Param_TFMxNE.LambdaSpace = lambda_space;
Param_TFMxNE.LambdaTime = lambda_time;
Param_TFMxNE.TFShape = [n_freq, n_frames];

%% 3. Lipschitz constant
if isempty(opts.Lipschitz)
    if verbose
        fprintf('TFMxNE: Estimating Lipschitz constant...\n');
    end
    lipschitz_k = 1.1 * local_estimate_lipschitz( ...
        G_eff, n_times, n_freq, n_frames, time_step, window_size, ...
        round(double(opts.LipschitzMaxIter)), 1e-3);
else
    lipschitz_k = double(opts.Lipschitz);
end
lipschitz_k = max(lipschitz_k, eps);
Param_TFMxNE.Lipschitz = lipschitz_k;

mu_l21 = lambda_space / lipschitz_k;
mu_l1 = lambda_time / lipschitz_k;

%% 4. FISTA in the time-frequency coefficient domain
if verbose
    fprintf('TFMxNE: SpatialReg=%.4g%%, TemporalReg=%.4g%%, lambda_space=%.4g, lambda_time=%.4g\n', ...
        spatial_reg_percent, temporal_reg_percent, lambda_space, lambda_time);
    fprintf('TFMxNE: FISTA MaxIter=%d, Tol=%.1e, Ncoef=%d, nd=%d\n', ...
        max_iter, tol, n_coeff, nd);
end

Z = complex(zeros(n_sources_total, n_coeff));
Y = Z;
t_fista = 1.0;
objective_path = nan(max_iter, 1);
rel_change_path = nan(max_iter, 1);
converged = false;
last_error = inf;
last_iter = max_iter;

for iter = 1:max_iter
    Z_old = Z;

    X_y = local_istft_matrix(Y, n_freq, n_frames, time_step, n_times);
    residual = M_w - G_eff * real(X_y);
    grad_time = G_eff' * residual;
    grad_tf = reshape(local_stft(grad_time, window_size, time_step), n_sources_total, n_coeff);

    Z = Y + grad_tf / lipschitz_k;
    Z = local_prox_l1_tf(Z, mu_l1, nd);
    Z = local_prox_l21_tf(Z, mu_l21, nd);

    t_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * t_fista ^ 2));
    Y = Z + ((t_fista - 1.0) / t_new) * (Z - Z_old);
    t_fista = t_new;

    denom = max(norm(Z_old(:)), eps);
    last_error = norm(Z(:) - Z_old(:)) / denom;
    rel_change_path(iter) = last_error;

    if iter == 1 || mod(iter, 10) == 0 || iter == max_iter || last_error < tol
        X_z = real(local_istft_matrix(Z, n_freq, n_frames, time_step, n_times));
        residual_z = M_w - G_eff * X_z;
        objective_path(iter) = 0.5 * sum(residual_z(:) .^ 2) + ...
            lambda_space * local_l21_norm_tf(Z, nd) + ...
            lambda_time * local_l1_norm_tf(Z, nd);
        if verbose
            active_locs = sum(local_active_locations(Z, nd));
            fprintf('  FISTA iter %4d | pobj=%.4e | relchg=%.3e | active locs=%d\n', ...
                iter, objective_path(iter), last_error, active_locs);
        end
    end

    last_iter = iter;
    if iter > 1 && last_error < tol
        converged = true;
        break;
    end
end

objective_path = objective_path(~isnan(objective_path));
rel_change_path = rel_change_path(1:last_iter);

%% 5. Back-project to source time series and undo source scaling
X_eff = real(local_istft_matrix(Z, n_freq, n_frames, time_step, n_times));
X_components = bsxfun(@times, X_eff, source_scaling(:));

active_locations = local_active_locations(Z, nd);
if perform_debias && any(active_locations)
    [X_components, debias_info] = local_debias_location_scaling( ...
        M_w, G_w, X_components, active_locations, nd, verbose);
else
    debias_info = struct('Applied', false, 'Method', 'none', 'Scales', []);
end

if nd > 1 && strcmpi(pick_ori, 'none')
    Source_TFMxNE = local_collapse_orientations(X_components, n_locations, nd);
else
    Source_TFMxNE = X_components;
end

active_components = find(any(abs(X_components) > 0, 2));
if nd > 1
    active_location_indices = unique(ceil(active_components / nd));
else
    active_location_indices = active_components;
end

Param_TFMxNE.ActiveSet = active_locations;
Param_TFMxNE.ActiveSetIndices_Locations = active_location_indices;
Param_TFMxNE.ActiveSetIndices_Components = active_components;
Param_TFMxNE.SourceComponentSize = size(X_components);
Param_TFMxNE.Objective = objective_path;
Param_TFMxNE.RelChange = rel_change_path;
Param_TFMxNE.FinalError = last_error;
Param_TFMxNE.Converged = converged;
Param_TFMxNE.ConvergedAt = last_iter;
Param_TFMxNE.DebiasingInfo = debias_info;

end

%% Optimization helpers

function Z = local_prox_l1_tf(Z, beta, nd)
    if beta <= 0 || isempty(Z)
        return;
    end
    if nd == 1
        magnitudes = abs(Z);
        shrink = max(1 - beta ./ max(magnitudes, beta), 0);
        Z = Z .* shrink;
        return;
    end

    [n_total, n_coeff] = size(Z);
    n_locations = n_total / nd;
    Z3 = reshape(Z, nd, n_locations, n_coeff);
    group_norms = sqrt(sum(abs(Z3) .^ 2, 1));
    shrink = max(1 - beta ./ max(group_norms, beta), 0);
    Z3 = bsxfun(@times, Z3, shrink);
    Z = reshape(Z3, n_total, n_coeff);
end

function Z = local_prox_l21_tf(Z, beta, nd)
    if beta <= 0 || isempty(Z)
        return;
    end
    [n_total, n_coeff] = size(Z);
    if nd == 1
        group_norms = sqrt(sum(abs(Z) .^ 2, 2));
        shrink = max(1 - beta ./ max(group_norms, beta), 0);
        Z = bsxfun(@times, Z, shrink);
        return;
    end

    n_locations = n_total / nd;
    Z3 = reshape(Z, nd, n_locations, n_coeff);
    group_norms = reshape(sqrt(sum(sum(abs(Z3) .^ 2, 1), 3)), [], 1);
    shrink = max(1 - beta ./ max(group_norms, beta), 0);
    Z3 = bsxfun(@times, Z3, reshape(shrink, 1, n_locations, 1));
    Z = reshape(Z3, n_total, n_coeff);
end

function val = local_l1_norm_tf(Z, nd)
    if isempty(Z)
        val = 0;
        return;
    end
    if nd == 1
        val = sum(abs(Z(:)));
        return;
    end
    [n_total, n_coeff] = size(Z);
    n_locations = n_total / nd;
    Z3 = reshape(Z, nd, n_locations, n_coeff);
    tmp = sqrt(sum(abs(Z3) .^ 2, 1));
    val = sum(tmp(:));
end

function val = local_l21_norm_tf(Z, nd)
    if isempty(Z)
        val = 0;
        return;
    end
    if nd == 1
        val = sum(sqrt(sum(abs(Z) .^ 2, 2)));
        return;
    end
    [n_total, n_coeff] = size(Z);
    n_locations = n_total / nd;
    Z3 = reshape(Z, nd, n_locations, n_coeff);
    val = sum(reshape(sqrt(sum(sum(abs(Z3) .^ 2, 1), 3)), [], 1));
end

function active_locations = local_active_locations(Z, nd)
    if isempty(Z)
        active_locations = false(0, 1);
        return;
    end
    if nd == 1
        active_locations = any(abs(Z) > 0, 2);
        return;
    end
    [n_total, n_coeff] = size(Z);
    n_locations = n_total / nd;
    Z3 = reshape(Z, nd, n_locations, n_coeff);
    active_locations = reshape(any(any(abs(Z3) > 0, 1), 3), [], 1);
end

function [lambda_max, group_norms] = local_l21_dual_norm_time(GtM, nd)
    if isempty(GtM)
        lambda_max = 0;
        group_norms = [];
        return;
    end
    if nd == 1
        group_norms = sqrt(sum(GtM .^ 2, 2));
    else
        n_locations = size(GtM, 1) / nd;
        GtM3 = reshape(GtM, nd, n_locations, size(GtM, 2));
        group_norms = reshape(sqrt(sum(sum(GtM3 .^ 2, 1), 3)), [], 1);
    end
    lambda_max = max(group_norms);
end

function L_est = local_estimate_lipschitz(G, n_times, n_freq, n_frames, tstep, wsize, max_iter, tol)
    n_sources = size(G, 2);
    rng_state = rng;
    cleanup_obj = onCleanup(@() rng(rng_state)); %#ok<NASGU>
    rng(0, 'twister');

    Z = complex(randn(n_sources, n_freq * n_frames), randn(n_sources, n_freq * n_frames));
    Z = Z / max(norm(Z(:)), eps);
    L_old = 0;
    L_est = 1;
    for iter = 1:max_iter
        X = real(local_istft_matrix(Z, n_freq, n_frames, tstep, n_times));
        GX = G * X;
        AtA_time = G' * GX;
        W = reshape(local_stft(AtA_time, wsize, tstep), n_sources, n_freq * n_frames);
        L_est = norm(W(:));
        if ~isfinite(L_est) || L_est < eps
            L_est = 1;
            break;
        end
        W = W / L_est;
        if iter > 1 && abs(L_est - L_old) / max(L_old, eps) < tol
            break;
        end
        Z = W;
        L_old = L_est;
    end
end

%% Preprocessing and output helpers

function [Data_w, L_w, info] = local_whiten_data(Data, L, C_noise, verbose)
    n_channels = size(Data, 1);
    Data_w = Data;
    L_w = L;
    W = eye(n_channels);
    info = struct('Applied', false, 'Rank', n_channels, ...
        'ConditionNumber', 1, 'WhiteningMatrix', W);

    if isempty(C_noise)
        return;
    end
    if size(C_noise, 1) ~= size(C_noise, 2)
        error('seal_TFMxNE:InvalidNoiseCovariance', 'NoiseCovariance must be square.');
    end
    if size(C_noise, 1) ~= n_channels
        error('seal_TFMxNE:NoiseCovarianceMismatch', ...
            'NoiseCovariance size must match the number of data channels.');
    end
    C_noise = full((double(C_noise) + double(C_noise)') * 0.5);
    if norm(C_noise - eye(n_channels), 'fro') <= 1e-12 * max(1, norm(C_noise, 'fro'))
        return;
    end

    [U, S] = svd(C_noise, 'econ');
    eig_vals = diag(S);
    s_max = max(eig_vals);
    if isempty(eig_vals) || ~isfinite(s_max) || s_max <= 0
        error('seal_TFMxNE:InvalidNoiseCovariance', ...
            'NoiseCovariance must have a positive finite eigenvalue.');
    end
    rank_tol = max(size(C_noise)) * eps(s_max);
    info.Rank = sum(eig_vals > rank_tol);
    floor_val = s_max * 1e-12;
    eig_vals(eig_vals < floor_val) = floor_val;
    W = diag(1 ./ sqrt(eig_vals)) * U';

    Data_w = W * Data;
    L_w = W * L;
    info.Applied = true;
    info.ConditionNumber = s_max / min(eig_vals);
    info.WhiteningMatrix = W;
    if verbose
        fprintf('TFMxNE: Applied noise whitening, rank=%d, cond=%.3g\n', ...
            info.Rank, info.ConditionNumber);
    end
end

function [scaling, info] = local_source_scaling( ...
    L_w, n_locations, nd, depth_params, depth_exp_override, depth_limit_override, rho)
    if ~ismember(nd, [1, 3])
        error('seal_TFMxNE:UnsupportedOrientationCount', ...
            'TFMxNE currently supports NumOrientations equal to 1 or 3.');
    end

    exponent = 0.8;
    limit = 10.0;
    if ~isempty(depth_params)
        if isfield(depth_params, 'exponent') && ~isempty(depth_params.exponent)
            exponent = double(depth_params.exponent);
        end
        if isfield(depth_params, 'limit') && ~isempty(depth_params.limit)
            limit = double(depth_params.limit);
        end
    end
    if ~isempty(depth_exp_override)
        exponent = double(depth_exp_override);
    end
    if ~isempty(depth_limit_override)
        limit = double(depth_limit_override);
    end

    loc_norms = zeros(n_locations, 1);
    for i_loc = 1:n_locations
        idx = (i_loc - 1) * nd + (1:nd);
        loc_norms(i_loc) = norm(L_w(:, idx), 'fro');
    end
    loc_norms(loc_norms < eps) = eps;
    loc_scaling = loc_norms .^ (-exponent);
    if limit > 0
        max_scaling = min(loc_scaling) * (limit ^ 2);
        loc_scaling(loc_scaling > max_scaling) = max_scaling;
    end

    scaling = repelem(loc_scaling, nd);
    orientation_scaling = ones(n_locations * nd, 1);
    if nd == 3
        orientation_scaling = repmat([1; rho; rho], n_locations, 1);
        scaling = scaling .* orientation_scaling;
    end

    info = struct('DepthExponent', exponent, 'DepthLimit', limit, ...
        'LocationNorms', loc_norms, 'LocationScaling', loc_scaling, ...
        'LooseOrientationWeight', rho, 'OrientationScaling', orientation_scaling);
end

function [X_components, info] = local_debias_location_scaling(M, G, X_components, active_locations, nd, verbose)
    n_locations = numel(active_locations);
    active_idx = find(active_locations);
    n_active = numel(active_idx);
    info = struct('Applied', false, 'Method', 'location_scaling_ls', ...
        'Scales', ones(n_locations, 1), 'ActiveLocations', active_idx);
    if n_active == 0
        return;
    end

    B = zeros(numel(M), n_active);
    for ii = 1:n_active
        loc = active_idx(ii);
        comp = (loc - 1) * nd + (1:nd);
        contribution = G(:, comp) * X_components(comp, :);
        B(:, ii) = contribution(:);
    end
    rhs = M(:);
    gram = B' * B;
    scales_active = (gram + 1e-12 * max(trace(gram), 1) * eye(n_active)) \ (B' * rhs);
    scales_active = max(real(scales_active), 1.0);

    for ii = 1:n_active
        loc = active_idx(ii);
        comp = (loc - 1) * nd + (1:nd);
        X_components(comp, :) = X_components(comp, :) * scales_active(ii);
    end
    info.Scales(active_idx) = scales_active;
    info.Applied = true;
    if verbose
        fprintf('TFMxNE: Debiasing location scales range %.3g to %.3g\n', ...
            min(scales_active), max(scales_active));
    end
end

function Source = local_collapse_orientations(X_components, n_locations, nd)
    Source = zeros(n_locations, size(X_components, 2));
    for i_loc = 1:n_locations
        idx = (i_loc - 1) * nd + (1:nd);
        Source(i_loc, :) = sqrt(sum(X_components(idx, :) .^ 2, 1));
    end
end

%% STFT helpers

function X = local_stft(x, wsize, tstep)
    if isempty(x)
        X = [];
        return;
    end
    [n_rows, n_times] = size(x);
    n_step = ceil(n_times / tstep);
    n_freq = wsize / 2 + 1;
    X = complex(zeros(n_rows, n_freq, n_step));

    win = sin((0.5:wsize-0.5) / wsize * pi);
    win2 = win .^ 2;
    frame_len = (n_step - 1) * tstep + wsize;
    swin = zeros(1, frame_len);
    for t = 1:n_step
        idx = (t - 1) * tstep + (1:wsize);
        swin(idx) = swin(idx) + win2;
    end
    swin = sqrt(wsize * max(swin, eps));

    pad = (wsize - tstep) / 2;
    xp = zeros(n_rows, frame_len);
    xp(:, pad + (1:n_times)) = x;

    for t = 1:n_step
        idx = (t - 1) * tstep + (1:wsize);
        wwin = repmat(win ./ swin(idx), n_rows, 1);
        frame = xp(:, idx) .* wwin;
        fframe = fft(frame, [], 2);
        X(:, :, t) = fframe(:, 1:n_freq);
    end
end

function x = local_istft(X, tstep, n_times)
    if isempty(X)
        x = [];
        return;
    end
    [n_rows, n_freq, n_step] = size(X);
    wsize = 2 * (n_freq - 1);
    frame_len = (n_step - 1) * tstep + wsize;
    x = zeros(n_rows, frame_len);

    win = sin((0.5:wsize-0.5) / wsize * pi);
    win2 = win .^ 2;
    swin = zeros(1, frame_len);
    for t = 1:n_step
        idx = (t - 1) * tstep + (1:wsize);
        swin(idx) = swin(idx) + win2;
    end
    swin = sqrt(max(swin, eps) / wsize);

    fframe = complex(zeros(n_rows, wsize));
    for t = 1:n_step
        fframe(:, 1:n_freq) = X(:, :, t);
        if n_freq > 2
            fframe(:, n_freq+1:end) = conj(X(:, n_freq-1:-1:2, t));
        end
        frame = real(ifft(fframe, [], 2));
        idx = (t - 1) * tstep + (1:wsize);
        wwin = repmat(win ./ swin(idx), n_rows, 1);
        x(:, idx) = x(:, idx) + frame .* wwin;
    end

    pad = (wsize - tstep) / 2;
    x = x(:, pad + (1:n_times));
end

function X_time = local_istft_matrix(Z, n_freq, n_frames, tstep, n_times)
    if isempty(Z)
        X_time = [];
        return;
    end
    n_rows = size(Z, 1);
    X_tf = reshape(Z, n_rows, n_freq, n_frames);
    X_time = local_istft(X_tf, tstep, n_times);
end

%% Small utilities

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
