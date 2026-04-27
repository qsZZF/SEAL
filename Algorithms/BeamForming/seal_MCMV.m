function [Source_MCMV, Param_MCMV] = seal_MCMV(Data, L, varargin)
%SEAL_MCMV Multi-Core Minimum Variance (MCMV) Beamformer (Moiseev 2011).
%   [Source, Param] = seal_MCMV(Data, L, 'Name', Value, ...)
%
%   Iteratively builds a joint MCMV beamformer that places null-constraints
%   between selected source locations to mitigate correlated-source bias.
%   Implements Moiseev et al. (NeuroImage 2011).
%
%   This rewrite supports two operating modes:
%     - Network mode (default for n <= MaxIterations): applies the iterative
%       nulling to all sources in L. Suitable when L is a pre-selected
%       network of interest (~tens of sources).
%     - Whole-brain mode (n > MaxIterations): selects the top
%       MaxIterations sources by single-source SNR, runs the iterative
%       MCMV on that subset, and assembles a sparse Nsources x Nchannels
%       inverse operator (zero rows for non-selected sources).
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints (can be [] if
%                             DataCovariance is provided).
%       L    (double matrix): Nchannels x Nsources lead-field matrix.
%
%   Optional Name-Value pairs:
%       'DataCovariance'           - Nch x Nch.   Default: computed from Data.
%       'NoiseCovariance'          - Nch x Nch.   Default: eye(Nch).
%       'EvokedCovariance'         - Nch x Nch.   If supplied, optimizes
%                                    SNR for evoked responses (Moiseev 2011).
%       'RegularizationParameter'  - Tikhonov shrinkage on data cov.
%                                    Default: 1e-3.
%       'NumOrientations'          - 1 (fixed) or 3 (free). Default: 1.
%       'MaxIterations'            - Hard cap on iterative selection. Use
%                                    Inf to process all sources (slow!).
%                                    Default: min(Nlocations, 30).
%       'AutoScaleNoise'           - If true and NoiseCov is identity-like,
%                                    rescale to match data covariance trace
%                                    (Moiseev's recommendation). Default: true.
%       'Verbose'                  - 0/1. Default: 0.
%
%   Outputs:
%       Source_MCMV (double matrix): Nsources x Ntimepoints.
%       Param_MCMV (struct):
%           .InverseOperator       - Nsources x Nchannels.
%           .OptimalOrientations   - nd x Nlocations.
%           .SourceOrder           - Order in which sources were selected.
%           .SelectedSources       - Indices of sources actually solved
%                                    via MCMV (others have zero rows).
%           .BeamSNR               - SNR values from the iterative loop.

    %% 1. Input parsing
    p = inputParser;
    p.CaseSensitive = false;
    addRequired(p, 'Data', @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addRequired(p, 'L',    @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'DataCovariance',          [],   @(x) isempty(x) || ismatrix(x));
    addParameter(p, 'NoiseCovariance',         [],   @(x) isempty(x) || ismatrix(x));
    addParameter(p, 'EvokedCovariance',        [],   @(x) isempty(x) || ismatrix(x));
    addParameter(p, 'RegularizationParameter', 1e-3, @isscalar);
    addParameter(p, 'NumOrientations',         1,    @isscalar);
    addParameter(p, 'MaxIterations',           [],   @(x) isempty(x) || isscalar(x));
    addParameter(p, 'AutoScaleNoise',          true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'Verbose',                 0,    @isscalar);
    parse(p, Data, L, varargin{:});

    Data_orig    = p.Results.Data;
    L_orig       = p.Results.L;
    C_y_user     = p.Results.DataCovariance;
    C_noise      = p.Results.NoiseCovariance;
    C_avg        = p.Results.EvokedCovariance;
    lambda       = p.Results.RegularizationParameter;
    nd           = p.Results.NumOrientations;
    auto_scale_n = logical(p.Results.AutoScaleNoise);
    verbose      = p.Results.Verbose;

    [Nchannels_L, Nsources_total] = size(L_orig);
    if mod(Nsources_total, nd) ~= 0
        error('seal_MCMV:NumOrientationsMismatch', ...
              'Total sources in L is not divisible by NumOrientations.');
    end
    Nlocs = Nsources_total / nd;

    if isempty(Data_orig)
        Nchannels = Nchannels_L;
    else
        Nchannels = size(Data_orig, 1);
        if Nchannels ~= Nchannels_L
            error('seal_MCMV:DimensionMismatch', ...
                  'Data and L must have same number of channels.');
        end
    end

    if isempty(C_noise),  C_noise = eye(Nchannels); end

    max_iter = p.Results.MaxIterations;
    if isempty(max_iter), max_iter = min(Nlocs, 30); end
    max_iter = min(max_iter, Nlocs);

    Is_Evoked = ~isempty(C_avg);

    if verbose
        fprintf('seal_MCMV: %d locations, %d orientations, max_iter=%d\n', ...
                Nlocs, nd, max_iter);
    end

    %% 2. Data covariance & inverse
    if isempty(C_y_user)
        if isempty(Data_orig)
            error('seal_MCMV:NoCovariance', ...
                  'Must provide Data or DataCovariance.');
        end
        Ntp = size(Data_orig, 2);
        Y0  = Data_orig - mean(Data_orig, 2);
        C_y = (Y0 * Y0') / max(Ntp - 1, 1);
    else
        C_y = C_y_user;
    end
    C_y   = (C_y + C_y') / 2;
    C_reg = C_y + lambda * (trace(C_y) / Nchannels) * eye(Nchannels);
    iR    = invSPD(C_reg);

    % Auto-scale noise covariance if user passed identity-like (Moiseev rec.)
    if auto_scale_n && isequal(C_noise, eye(Nchannels))
        C_noise = (Nchannels / max(trace(iR), eps)) * eye(Nchannels);
    end

    if Is_Evoked
        iR_working = iR * C_avg * iR;
    else
        iR_working = iR;
    end
    iR_N_iR = iR * C_noise * iR;

    %% 3. Single-source SNR initialization
    lstU = zeros(Nlocs, nd);    % optimal orientations [Nlocs x nd]
    SNR0 = zeros(Nlocs, 1);

    for i = 1:Nlocs
        idx = (i-1)*nd + (1:nd);
        h   = L_orig(:, idx);

        T = h' * iR_N_iR    * h;
        S = h' * iR_working * h;

        if nd > 1
            [V, Ev] = eig(S, T);
            [~, k]  = max(real(diag(Ev)));
            u       = real(V(:, k));
            u       = u / max(norm(u), eps);
        else
            u = 1;
        end
        lstU(i, :) = u(:).';

        hu      = h * u;
        Su_scal = hu' * iR_working * hu;
        Tu_scal = hu' * iR_N_iR    * hu;
        SNR0(i) = Su_scal / max(Tu_scal, eps);
    end

    %% 4. Iterative MCMV: pick top sources one by one
    selected   = false(Nlocs, 1);  % which sources are in the reference set
    order      = zeros(max_iter, 1);
    BeamSNR    = zeros(max_iter, 1);
    SNR_curr   = SNR0;             % live SNR estimates (updated each iter)

    % --- iteration 1: pick the strongest single source
    [snr1, i1] = max(SNR_curr);
    selected(i1) = true;
    order(1)     = i1;
    BeamSNR(1)   = snr1;

    if verbose
        fprintf('  iter 1/%d: picked src %d, SNR=%.3f\n', max_iter, i1, snr1);
    end

    % --- iterations 2..max_iter: re-evaluate non-selected sources under
    %     the constraint of zero leakage from selected ones
    for iIter = 2:max_iter
        ref_locs = find(selected);     % current reference set
        nRef     = numel(ref_locs);

        % Build Hr = [h_{r1}*u_{r1}, h_{r2}*u_{r2}, ...]  (Nch x nRef)
        Hr = zeros(Nchannels, nRef);
        for o = 1:nRef
            r   = ref_locs(o);
            idx = (r-1)*nd + (1:nd);
            Hr(:, o) = L_orig(:, idx) * lstU(r, :).';
        end

        iRNiR_Hr = iR_N_iR    * Hr;     % Nch x nRef
        iR_Hr    = iR_working * Hr;     % Nch x nRef
        Tr  = Hr' * iRNiR_Hr;           % nRef x nRef
        Sr  = Hr' * iR_Hr;              % nRef x nRef
        iTr = invSPD(Tr);
        iTrSrTr = iTr * Sr * iTr;       % cached for inner loop

        % Re-evaluate every non-selected location
        cand = find(~selected);
        cand_snr = -inf(numel(cand), 1);
        cand_u   = zeros(numel(cand), nd);

        for c = 1:numel(cand)
            i   = cand(c);
            idx = (i-1)*nd + (1:nd);
            h   = L_orig(:, idx);

            T   = h' * iR_N_iR    * h;
            S   = h' * iR_working * h;
            Tsr = h' * iRNiR_Hr;        % nd x nRef
            Ssr = h' * iR_Hr;           % nd x nRef

            % Effective S, T after constraining away the reference subspace
            D = Tsr * iTrSrTr * Tsr' - Tsr * iTr * Ssr' ...
                                    - Ssr * iTr * Tsr' + S;
            F = T - Tsr * iTr * Tsr';

            if nd > 1
                [V, Ev] = eig(D, F);
                [~, k]  = max(real(diag(Ev)));
                u       = real(V(:, k));
                u       = u / max(norm(u), eps);
            else
                u = 1;
            end
            cand_u(c, :) = u(:).';

            hu = h * u;
            Hu = [Hr, hu];                 %#ok<AGROW>  small, OK
            Su = Hu' * iR_working * Hu;
            Tu = Hu' * iR_N_iR    * Hu;
            iTu = invSPD(Tu);
            cand_snr(c) = trace(iTu * Su);
        end

        % Pick the candidate that maximizes joint SNR
        [snr_best, c_best] = max(cand_snr);
        i_best = cand(c_best);

        selected(i_best)    = true;
        lstU(i_best, :)     = cand_u(c_best, :);
        order(iIter)        = i_best;
        BeamSNR(iIter)      = snr_best;
        SNR_curr(i_best)    = snr_best;   % bookkeeping (not strictly needed)

        if verbose
            fprintf('  iter %d/%d: picked src %d, joint SNR=%.3f\n', ...
                    iIter, max_iter, i_best, snr_best);
        end
    end

    %% 5. Joint MCMV weights for the selected reference set
    sel_idx = order(1:max_iter);                        % size: max_iter
    Ht = zeros(max_iter, Nchannels);
    for k = 1:max_iter
        loc = sel_idx(k);
        idx = (loc-1)*nd + (1:nd);
        Ht(k, :) = lstU(loc, :) * L_orig(:, idx).';
    end

    %   W_sel = (Ht * iR * Ht')^{-1} * Ht * iR
    %   size: max_iter x Nchannels
    W_sel = invSPD(Ht * iR * Ht.') * Ht * iR;

    % Map back into a full Nlocs x Nchannels operator (zeros elsewhere)
    W_full = zeros(Nlocs, Nchannels);
    W_full(sel_idx, :) = W_sel;

    %% 6. Outputs
    Param_MCMV.InverseOperator     = W_full;
    Param_MCMV.OptimalOrientations = lstU.';        % nd x Nlocs
    Param_MCMV.SourceOrder         = sel_idx(:).';
    Param_MCMV.SelectedSources     = sel_idx(:).';
    Param_MCMV.BeamSNR             = BeamSNR(:).';
    Param_MCMV.OptionsPassed       = p.Results;

    if ~isempty(Data_orig)
        Source_MCMV = W_full * Data_orig;
    else
        Source_MCMV = [];
    end
end


%% ---------- Helper: robust SPD inverse via Cholesky ----------
function Am1 = invSPD(A)
    A = (A + A') * 0.5;
    [R, p] = chol(A);
    if p == 0
        Am1 = R \ (R.' \ eye(size(A)));
    else
        % Add small jitter and try again; fall back to pinv on failure
        jitter = 1e-12 * trace(A) / size(A, 1) * eye(size(A));
        [R, p] = chol(A + jitter);
        if p == 0
            Am1 = R \ (R.' \ eye(size(A)));
        else
            Am1 = pinv(A);
        end
    end
end