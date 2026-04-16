function [Source_TFMxNE, Param_TFMxNE] = seal_TFMxNE(Data, L, varargin)
%SEAL_TFMXNE Time-Frequency Mixed-Norm Estimate (STOUT / TF-MxNE).
%   [Source, Param] = seal_TFMxNE(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Solves the spatio-temporal inverse problem using the Fast Iterative 
%   Shrinkage-Thresholding Algorithm (FISTA). It projects the data into the 
%   time-frequency domain (STFT) and applies an L2,1 norm (spatial sparsity) 
%   and an L1 norm (time-frequency sparsity).
%
%   Based on Gramfort et al. (2013) / Castano-Candamil (2015) STOUT port.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints measured data.
%       L (double matrix): Nchannels x Nsources lead field.
%
%   Optional Name-Value Pair Inputs:
%       'SpatialReg' (double): % spatial regularization (sreg). Default: 90.
%       'TemporalReg' (double): % temporal regularization (treg). Default: 30.
%       'TimeStep' (integer): STFT time shift (a). Default: 4.
%       'WindowSize' (integer): STFT frequency bins (m). Default: 80.
%       'MaxIterations' (integer): FISTA maximum iterations. Default: 50.
%       'Tolerance' (double): Stopping criterion. Default: 1e-3.
%       'Lipschitz' (double): Lipschitz constant. If empty, it is estimated.
%       'NumOrientations' (integer): Orientations per dipole. Default: 1.
%
%   Outputs:
%       Source_TFMxNE (double matrix): Reconstructed source activity.
%       Param_TFMxNE (struct): Contains solver metadata (Active Set, Lipschitz, etc.)

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    
    addParameter(p, 'TimeStep', 4, @isscalar);
    addParameter(p, 'WindowSize', 80, @isscalar);
    addParameter(p, 'SpatialReg', 90, @isscalar);
    addParameter(p, 'TemporalReg', 30, @isscalar);
    addParameter(p, 'MaxIterations', 50, @isscalar);
    addParameter(p, 'Tolerance', 1e-3, @isscalar);
    addParameter(p, 'Lipschitz', [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'NumOrientations', 1, @isscalar);

    parse(p, Data, L, varargin{:});
    opts = p.Results;

    y = opts.Data;
    L = opts.L;
    nd = opts.NumOrientations;
    [Nchannels, Ntimepoints] = size(y);    
    Nsources_total = size(L, 2);

    %% 2. Initialization and STFT Setup
    normL = sum(L.^2, 1);
    if nd > 1
        normL = sum(reshape(normL, nd, []));
        weights = 1 ./ sqrt(normL);
        weights = repmat(weights, nd, 1);
        weights = weights(:);
    else
        weights = 1 ./ sqrt(normL);
    end
    thr = 0.01 * max(weights);
    weights(weights < thr) = thr;
    L = bsxfun(@times, L, weights');
    c = stft(y, opts.WindowSize, opts.TimeStep);
    T_stft = size(c, 3);
    K_stft = size(c, 2);
    L_ori = L;
    % Leadfield normalization
    tempGY = L' * y;
    if nd > 1
        aux = sum(reshape(tempGY.^2', [], Nsources_total/nd)', 2);
    else
        aux = tempGY.^2;
    end
    basepar = 0.01 * sqrt(max(aux(:)));
    L = L / basepar;

    %% 3. Lipschitz Constant Estimation
    if isempty(opts.Lipschitz)
        fprintf('Estimating Lipschitz constant... \n');
        lipschitz_k = 1.1 * estimate_lipschitz(y, L, 1e-3, opts.TimeStep, opts.WindowSize);
    else
        lipschitz_k = opts.Lipschitz;
    end
    
    mu_lc = opts.SpatialReg / lipschitz_k;
    lambda_lc = opts.TemporalReg / lipschitz_k;

    %% 4. FISTA Optimization Loop
    fprintf('Running FISTA algorithm...\n');
    R = y;
    active_set = logical(sparse(1, Nsources_total));
    Y_time_as = [];
    Y_as = [];
    
    Z = sparse(0, K_stft * T_stft);
    Y = sparse(Nsources_total, K_stft * T_stft);
    tau = 1;
    error = inf;
    
    for i = 1:opts.MaxIterations
        Z_0 = Z;
        active_set0 = active_set;
        
        % Gradient Step & Proximal Operators
        if ~isempty(Y_time_as) && sum(active_set) < Nchannels
            GTR = L' * R / lipschitz_k;
            A = GTR;
            A(Y_as,:) = A(Y_as,:) + Y_time_as;
            
            [~, active_set_l21] = prox_l21(A, mu_lc, nd);
            
            temp = stft(GTR(active_set_l21, :), opts.WindowSize, opts.TimeStep);
            temp = reshape(temp, sum(active_set_l21), []);
            Baux = Y(active_set_l21, :) + temp;
            
            [Z, active_set_l1] = prox_l1(Baux, lambda_lc, nd);
            
            active_set_l21(active_set_l21) = active_set_l1;
            active_set_l1 = active_set_l21;            
        else
            temp = stft(R, opts.WindowSize, opts.TimeStep);
            temp = reshape(temp, Nchannels, []);
            Y = Y + L' * temp / lipschitz_k;
            [Z, active_set_l1] = prox_l1(Y, lambda_lc, nd);
        end
        
        [Z, active_set_l21] = prox_l21(Z, mu_lc, nd);
        active_set = active_set_l1;
        active_set(active_set_l1) = active_set_l21;
        
        % Check Convergence
        if ~any(active_set ~= active_set0) && i > 1
            error = norm(abs(Z - Z_0)) / (norm(abs(Z_0)) + eps);
        else
            error = inf;
        end
        
        % FISTA Momentum Update (tau)
        if i < opts.MaxIterations
            tau_0 = tau;
            tau = 0.5 + sqrt(1 + 4 * tau_0^2) / 2;
            Y = sparse(Nsources_total, K_stft * T_stft);
            dt = (tau_0 - 1) / tau;
            Y(active_set, :) = (1 + dt) * Z;
            Y(active_set0, :) = Y(active_set0, :) - dt * Z_0;
            
            Y_as = active_set0(:) | active_set(:);
            act_Y_as = find(Y_as);
            
            temp = reshape(full(Y(act_Y_as, :)), length(act_Y_as), K_stft, T_stft);
            Y_time_as = istft(temp, opts.TimeStep, Ntimepoints);
            
            if isempty(Y_time_as)
                R = y;
            else
                R = y - L(:, act_Y_as) * Y_time_as(:, 1:Ntimepoints);
            end
        end
        
        if error < opts.Tolerance || (sum(full(active_set)) == 0 && i > 1)
            fprintf('FISTA Converged at iteration %d.\n', i);
            break
        end
    end

    %% 5. Transform Back to Time Domain
    temp = reshape(full(Z), [], K_stft, T_stft);
    temp = istft(temp, opts.TimeStep, Ntimepoints);
    
    if ~isempty(temp)
        J_recf = zeros(Nsources_total, size(temp, 2));
        J_recf(active_set, :) = temp;
        J_recf = J_recf(:, 1:Ntimepoints);
        
        % Rescale
        Jf = zeros(Nsources_total, Ntimepoints);
        Jf(1:Nsources_total, :) = J_recf;
        Jf = Jf * norm(y, 2) / (norm(L_ori * Jf, 2) + eps);
    else
        Jf = zeros(Nsources_total, Ntimepoints);
    end

    Source_TFMxNE = Jf;
    Param_TFMxNE.ActiveSet = active_set;
    Param_TFMxNE.Lipschitz = lipschitz_k;
    Param_TFMxNE.FinalError = error;
end


function k = estimate_lipschitz(y, L, tol, a, M)
    Nt = size(y, 2);
    Nd = size(L, 2);
    iv = ones(Nd, Nt);
    v = stft(iv, M, a);
    l = 1e100;
    l_old = 0;
    
    for i = 1:100
        l_old = l;
        aux = istft(v, a, Nt);
        iv = real(aux);
        Lv = L * iv;
        LtLv = L' * Lv;
        w = stft(LtLv, M, a);
        l = max(abs(w(:)));
        v = w / l;
        if abs(l - l_old) / l_old < tol
            break
        end
    end
    k = l;
end

function [Y, active_set] = prox_l21(Y, mu, n_orient)
    n_pos = size(Y, 1) / n_orient;
    rows_norm = sqrt(sum(reshape(abs(Y).^2', [], n_pos)', 2));
    shrink = max(1 - mu ./ max(rows_norm, mu), 0);
    active_set = (shrink > 0);
    shrink = shrink(active_set);
    
    if n_orient > 1
        active_set = repmat(active_set(:)', n_orient, 1);
        active_set = active_set(:)';
    end
    
    temp = repmat(shrink(:)', n_orient, 1);
    temp = temp(:);
    Y = Y(active_set, :) .* repmat(temp, 1, size(Y, 2));
end

function [Y, active_set] = prox_l1(Y, lambda, n_orient)
    n_pos = size(Y, 1) / n_orient;
    norms = sqrt(sum(reshape((abs(Y).^2), n_orient, []), 1));
    shrink = max(1 - lambda ./ max(norms, lambda), 0);
    shrink = reshape(shrink', n_pos, []);
    active_set = sum(shrink, 2) > 0;
    shrink = shrink(active_set, :);
    
    active_set = repmat(active_set(:)', n_orient, 1);
    active_set = active_set(:)';
    Y = Y(active_set, :);
    
    if ~isempty(Y)
        for i = 1:n_orient
            Y(i:n_orient:end, :) = Y(i:n_orient:end, :) .* shrink;
        end
    end
end

function X = stft(x, wsize, tstep)
    if isempty(x), X = []; return; end
    [Nc, Nt] = size(x);
    n_step = ceil(Nt / tstep);
    n_freq = wsize / 2 + 1;
    X = zeros(Nc, n_freq, n_step);
    
    win = sin((0.5:wsize-0.5) / wsize * pi);
    win2 = win.^2;
    
    swin = zeros(1, (n_step - 1) * tstep + wsize);
    for t = 1:n_step
        idx = (t-1)*tstep+1 : (t-1)*tstep+wsize;
        swin(idx) = swin(idx) + win2;
    end
    swin = sqrt(wsize * swin);
    
    xp = zeros(Nc, wsize + (n_step-1)*tstep);
    xp(:, (wsize-tstep)/2+1 : (wsize-tstep)/2+Nt) = x;
    
    for t = 1:n_step
        idx = (t-1)*tstep+1 : (t-1)*tstep+wsize;
        wwin = repmat(win ./ swin(idx), Nc, 1);
        frame = xp(:, idx) .* wwin;
        fframe = fft(frame, [], 2);
        X(:,:,t) = fframe(:, 1:n_freq);
    end
end

function x = istft(X, tstep, Tx)
    if isempty(X), x = []; return; end
    [Nc, n_win, n_step] = size(X);
    wsize = 2 * (n_win - 1);
    Nt = n_step * tstep;
    x = zeros(Nc, Nt + wsize - tstep);
    
    win = sin((0.5:wsize-0.5) / wsize * pi);
    win2 = win.^2;
    
    swin = zeros(1, (n_step - 1) * tstep + wsize);
    for t = 1:n_step
        idx = (t-1)*tstep+1 : (t-1)*tstep+wsize;
        swin(idx) = swin(idx) + win2;
    end
    swin = sqrt(swin / wsize);
    
    fframe = zeros(Nc, n_win + wsize/2 - 1);
    for t = 1:n_step
        fframe(:, 1:n_win) = X(:,:,t);
        fframe(:, n_win+1:end) = conj(X(:, wsize/2:-1:2, t));
        frame = ifft(fframe, [], 2, 'symmetric');
        
        idx = (t-1)*tstep+1 : (t-1)*tstep+wsize;
        wwin = repmat(win ./ swin(idx), Nc, 1);
        x(:, idx) = x(:, idx) + real(conj(frame) .* wwin);
    end
    x = x(:, (wsize-tstep)/2 : (wsize-tstep)/2+Tx+1);
    x = x(:,1:Tx);
end