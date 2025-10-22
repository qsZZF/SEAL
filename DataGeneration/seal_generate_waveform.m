function waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, waveformParams)
%SEAL_GENERATE_WAVEFORM Generates various types of temporal waveforms with default parameters.
%   waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, waveformParams)
%
%   Inputs:
%       timeVector (double vector): 1xNtimepoints vector of time points (seconds).
%       samplingRate (double): Sampling rate in Hz.
%       waveformType (char/string): Type of waveform, e.g., 'erp', 'ersp',
%                                   'gaussian_pulse', 'spike', 'ar',
%                                   'resting_alpha', 'resting_pinknoise'.
%       waveformParams (struct): Parameters specific to the waveformType.
%                                If sub-fields are missing, defaults are used.
%                                An overall 'amplitude' field is generally expected for scaling.
%
%   Outputs:
%       waveform (double vector): 1xNtimepoints generated waveform.
%
%   Waveform-Specific Parameters in waveformParams (with defaults if not provided):
%       'erp':
%           .components (struct array): [{latency_s, amplitude_rel, width_s}, ...]
%                                       Default: single component at epoch_center, amp_rel=1, width=0.05s.
%           .amplitude (overall amplitude scale for the summed ERP). Default: 1.0.
%       'ersp':
%           .base_frequency_band_Hz. Default: 10 Hz (alpha).
%           .base_amplitude (amplitude of carrier before modulation). Default: 1.0.
%           .modulation (struct): Default: type='burst', latency=center, width=20% epoch, taper=0.25, min_rel_amp=0.1.
%               .type ('burst', 'invburst', 'ampmod')
%               .latency_s, .width_s, .taper (0-1)
%               .min_rel_amplitude (0-1)
%               .mod_frequency_Hz (for 'ampmod'). Default: 1 Hz.
%               .mod_phase_rad (for 'ampmod'). Default: 0.
%           .amplitude (final overall scaling factor for the ERSP waveform). Default: 1.0.
%       'gaussian_pulse':
%           .pulse_center_time_s. Default: center of timeVector.
%           .pulse_width_time_s (std dev). Default: 10% of epoch duration.
%           .amplitude. Default: 1.0.
%       'spike':
%           .latency_s (peak time). Default: center of timeVector.
%           .duration_s (approx total duration). Default: 0.07s (70ms).
%           .shape ('triangular', 'exponential_decay'). Default: 'exponential_decay'.
%           .decay_tau_s (for exponential_decay). Default: duration_s/4.
%           .amplitude. Default: 1.0.
%       'ar': (Autoregressive Model - placeholder, requires specific parameters)
%           .ar_coeffs (matrix): [nsignal x nsignal x n_lags] AR coefficients.
%                                (For single channel: [1 x 1 x n_lags] or [n_lags x 1]).
%           .innovation_variance (scalar): Variance of the innovation process. Default: 1.0.
%           .burn_in_samples (integer): Number of samples to discard for burn-in. Default: 1000.
%           .amplitude. Default: 1.0 (applied after generation and normalization).
%       'resting_alpha', 'resting_beta', etc.:
%           .center_frequency_Hz. Default: 10 (alpha), 20 (beta), 6 (theta), 2 (delta), 40 (gamma).
%           .bandwidth_Hz. Default: 4 (alpha), 8 (beta), 3 (theta), 2 (delta), 10 (gamma).
%           .burst_properties (optional struct): Default: duty_cycle=0.5, min_dur=0.5s, max_dur=1.5s.
%           .amplitude. Default: 1.0.
%       'resting_pinknoise':
%           .amplitude. Default: 1.0.
%
%   Example:
%       Fs = 100; t = 0:1/Fs:2-1/Fs; % 2 seconds
%       erp_params.components = struct('latency_s',{0.3,0.5},'amplitude_rel',{1,-0.5},'width_s',{0.05,0.08});
%       erp_params.amplitude = 5;
%       erp_wave = seal_generate_waveform(t, Fs, 'erp', erp_params);
%       figure; plot(t, erp_wave); title('Simulated ERP');
%
%       pulse_params.amplitude = 2; % other params will use defaults
%       pulse_wave = seal_generate_waveform(t, Fs, 'gaussian_pulse', pulse_params);

%   Author: FengZhao

    Ntimepoints = length(timeVector);
    waveform = zeros(1, Ntimepoints); % Initialize output
    epoch_duration_s = timeVector(end) - timeVector(1) + (1/samplingRate);
    epoch_center_s = timeVector(1) + epoch_duration_s / 2;

    % Ensure waveformParams is a struct
    if ~isstruct(waveformParams), waveformParams = struct(); end

    % Get overall amplitude scale, defaulting if not present
    if isfield(waveformParams, 'amplitude')
        amplitude_scale = waveformParams.amplitude;
    else
        amplitude_scale = 1.0;
        % For ERSP, base_amplitude is for carrier, main amplitude is for overall scaling
        if ~strcmpi(waveformType, 'ersp') && ~strcmpi(waveformType, 'ar') % AR handles its own scaling then applies amplitude
            warning('seal_generate_waveform:DefaultAmplitude', ...
                'waveformParams.amplitude not provided for type "%s". Using default 1.0.', waveformType);
        end
    end

    switch lower(waveformType)
        case 'erp'
            defaults.components = struct('type', 'gaussian', 'latency_s', epoch_center_s, 'amplitude_rel', 1.0, 'width_s', 0.05);
            params = merge_structs(waveformParams, defaults);
            
            if ~isstruct(params.components) || isempty(params.components)
                error('ERP type requires .components struct array in waveformParams.');
            end

            for i_comp = 1:length(params.components)
                comp = params.components(i_comp);
                comp_defaults.type = 'gaussian';
                comp_defaults.latency_s = epoch_center_s;
                comp_defaults.amplitude_rel = 1.0;
                % Defaults specific to component type
                if strcmpi(comp_defaults.type, 'gaussian')
                    comp_defaults.width_s = 0.05; % For Gaussian type
                elseif strcmpi(comp_defaults.type, 'gabor')
                    comp_defaults.frequency_Hz = 10; % For Gabor type
                    comp_defaults.std_dev_s = 0.05;  % Envelope width for Gabor
                    comp_defaults.phase_rad = 0;     % Phase of sine for Gabor
                end
                comp = merge_structs(comp, comp_defaults); % Apply general defaults first
                % Then apply type-specific defaults if type was user-specified but params missing
                if isfield(comp,'type') && strcmpi(comp.type, 'gaussian') && ~isfield(comp,'width_s')
                    comp.width_s = 0.05;
                end
                if isfield(comp,'type') && strcmpi(comp.type, 'gabor')
                    if ~isfield(comp,'frequency_Hz'), comp.frequency_Hz = 10; end
                    if ~isfield(comp,'std_dev_s'), comp.std_dev_s = 0.05; end
                    if ~isfield(comp,'phase_rad'), comp.phase_rad = 0; end
                end


                if strcmpi(comp.type, 'gaussian')
                    if ~isfield(comp,'width_s'), error('Gaussian ERP component needs .width_s'); end
                    waveform = waveform + (comp.amplitude_rel * exp(-((timeVector - comp.latency_s).^2) / (2 * comp.width_s^2)));
                elseif strcmpi(comp.type, 'gabor')
                    if ~isfield(comp,{'frequency_Hz','std_dev_s'}), error('Gabor ERP component needs .frequency_Hz and .std_dev_s'); end
                    gaussian_envelope = exp(-((timeVector - comp.latency_s).^2) / (2 * comp.std_dev_s^2));
                    sinusoid = sin(2 * pi * comp.frequency_Hz * (timeVector - comp.latency_s) + comp.phase_rad);
                    waveform = waveform + (comp.amplitude_rel * gaussian_envelope .* sinusoid);
                else
                    error('Unknown ERP component type: %s. Use "gaussian" or "gabor".', comp.type);
                end
            end
            waveform = waveform * amplitude_scale;

        case 'ersp'
            defaults.base_frequency_band_Hz = 10; % Alpha
            defaults.base_amplitude = 1.0;
            defaults.modulation = struct('type','burst', ...
                                         'latency_s', epoch_center_s, ...
                                         'width_s', epoch_duration_s * 0.2, ...
                                         'taper', 0.25, ...
                                         'min_rel_amplitude', 0.1, ...
                                         'mod_frequency_Hz', 1, ...
                                         'mod_phase_rad', 0);
            defaults.amplitude = 1.0; % Overall scaling for the final ERSP
            params = merge_structs(waveformParams, defaults);
            % Ensure modulation sub-struct has all its defaults if partially specified
            if isfield(waveformParams, 'modulation') && isstruct(waveformParams.modulation)
                params.modulation = merge_structs(waveformParams.modulation, defaults.modulation);
            end


            % Base signal generation
            if numel(params.base_frequency_band_Hz) == 1 
                base_freq = params.base_frequency_band_Hz;
                phase_val = 0; if isfield(params, 'base_phase_rad'), phase_val = params.base_phase_rad; end
                signal_base = sin(phase_val + 2*pi*base_freq*timeVector);
            elseif numel(params.base_frequency_band_Hz) == 4 
                freq_edges = params.base_frequency_band_Hz; 
                try
                    nyquist = samplingRate / 2;
                    freq_edges = sort(max(0, min(freq_edges, nyquist-eps))); 
                    if freq_edges(2) <= freq_edges(1), freq_edges(2) = freq_edges(1) + max(eps,nyquist/100); end
                    if freq_edges(4) <= freq_edges(3), freq_edges(3) = freq_edges(4) - max(eps,nyquist/100); end
                    if freq_edges(3) <= freq_edges(2), freq_edges(3) = freq_edges(2) + max(eps,nyquist/100); end
                    
                    order = 4; 
                    passband_norm = [freq_edges(2) freq_edges(3)] / nyquist;
                    if passband_norm(1) >= passband_norm(2) || any(passband_norm <=0) || any(passband_norm >=1)
                        warning('seal_generate_waveform:ERSPFilterBand', 'Problem with ERSP filter band. Using full band noise.');
                        signal_base_raw = randn(1, Ntimepoints + 2*round(samplingRate)) - 0.5; 
                        signal_base = signal_base_raw(round(samplingRate)+1:round(samplingRate)+Ntimepoints);
                    else
                        [b_filt, a_filt] = butter(order, passband_norm, 'bandpass');
                        signal_base_raw = randn(1, Ntimepoints + 2*round(samplingRate)) - 0.5; 
                        filtered_noise = filter(b_filt, a_filt, signal_base_raw);
                        signal_base = filtered_noise(round(samplingRate)+1:round(samplingRate)+Ntimepoints); 
                    end
                catch ME_filt
                     warning('seal_generate_waveform:ERSPFilterFail', 'ERSP filter design failed (%s). Using white noise.', ME_filt.message);
                     signal_base_raw = randn(1, Ntimepoints + 2*round(samplingRate)) - 0.5;
                     signal_base = signal_base_raw(round(samplingRate)+1:round(samplingRate)+Ntimepoints);
                end
                signal_base = signal_base / (max(abs(signal_base)) + eps); 
            else
                error('ERSP base_frequency_band_Hz must be scalar or 4-element vector.');
            end
            signal_base = params.base_amplitude * signal_base;

            % Modulation
            mod_params = params.modulation;
            mod_win = ones(1, Ntimepoints); 

            if ismember(lower(mod_params.type), {'burst', 'invburst'})
                latency_samples = floor((mod_params.latency_s)*samplingRate) + 1;
                width_samples = floor((mod_params.width_s)*samplingRate);
                if width_samples < 1, width_samples = 0; end
                
                temp_win = zeros(1, Ntimepoints);
                if width_samples > 0
                    tukey_w = tukeywin(width_samples, mod_params.taper)';
                    start_idx = latency_samples - floor(width_samples/2);
                    end_idx = start_idx + width_samples - 1;
                    
                    valid_tukey_start = 1; if start_idx < 1, valid_tukey_start = 1 - (start_idx-1); start_idx = 1; end
                    valid_tukey_end = width_samples; 
                    if end_idx > Ntimepoints, valid_tukey_end = width_samples - (end_idx-Ntimepoints); end_idx = Ntimepoints;end
                    
                    if start_idx <= end_idx && valid_tukey_start <= valid_tukey_end && valid_tukey_end > 0 && valid_tukey_start > 0
                         temp_win(start_idx:end_idx) = tukey_w(valid_tukey_start:valid_tukey_end);
                    end
                end

                if strcmpi(mod_params.type, 'invburst'), temp_win = 1 - temp_win; end
                
                min_w = min(temp_win); max_w = max(temp_win);
                if (max_w - min_w) > eps
                    mod_win = mod_params.min_rel_amplitude + (1-mod_params.min_rel_amplitude) * (temp_win - min_w) / (max_w - min_w);
                else 
                    if max_w > 0.5, mod_win = ones(1, Ntimepoints);
                    else, mod_win = ones(1, Ntimepoints) * mod_params.min_rel_amplitude;
                    end
                end
            elseif strcmpi(mod_params.type, 'ampmod')
                mod_sig = sin(mod_params.mod_phase_rad + 2*pi*mod_params.mod_frequency_Hz*timeVector);
                min_ms = min(mod_sig); max_ms = max(mod_sig);
                if (max_ms - min_ms) > eps
                    mod_win = mod_params.min_rel_amplitude + (1-mod_params.min_rel_amplitude) * (mod_sig - min_ms) / (max_ms - min_ms);
                else 
                    mod_win = ones(1, Ntimepoints) * (mod_params.min_rel_amplitude + 1)/2;
                end
            end
            waveform = signal_base .* mod_win;
            waveform = waveform * params.amplitude; % Final overall scaling from .amplitude

        case 'gaussian_pulse'
            defaults.pulse_center_time_s = epoch_center_s;
            defaults.pulse_width_time_s = epoch_duration_s / 10;
            params = merge_structs(waveformParams, defaults);
            waveform = amplitude_scale * exp(-((timeVector - params.pulse_center_time_s).^2) / (2 * params.pulse_width_time_s^2));
            
        case 'spike'
            defaults.latency_s = epoch_center_s;
            defaults.duration_s = 0.07; % 70ms
            defaults.shape = 'exponential_decay';
            defaults.decay_tau_s = defaults.duration_s / 4;
            params = merge_structs(waveformParams, defaults);
            if strcmpi(params.shape, 'exponential_decay') && ~isfield(waveformParams,'decay_tau_s') % if user specified shape but not tau
                params.decay_tau_s = params.duration_s / 4;
            end
            
            idx_peak = round(params.latency_s * samplingRate);
            if idx_peak < 1, idx_peak = 1; end
            if idx_peak > Ntimepoints, idx_peak = Ntimepoints; end
            half_duration_samples = round(params.duration_s * samplingRate / 2);
            
            switch params.shape
                case 'triangular'
                    start_idx = max(1, idx_peak - half_duration_samples);
                    end_idx = min(Ntimepoints, idx_peak + half_duration_samples);
                    temp_waveform = zeros(1,Ntimepoints);
                    if start_idx < idx_peak
                        temp_waveform(start_idx:idx_peak) = linspace(0, 1, idx_peak - start_idx + 1);
                    elseif start_idx == idx_peak 
                         temp_waveform(idx_peak) = 1;
                    end
                    if idx_peak < end_idx
                        temp_waveform(idx_peak+1:end_idx) = linspace(1 * (1 - 1/(end_idx-idx_peak + eps)), 0, end_idx - idx_peak);
                    elseif idx_peak == end_idx && start_idx == idx_peak 
                        % do nothing
                    end
                    waveform = amplitude_scale * temp_waveform;

                case 'exponential_decay'
                    time_rel_peak = timeVector - params.latency_s;
                    spike_shape_raw = exp(-abs(time_rel_peak) / params.decay_tau_s);
                    spike_shape_raw(abs(time_rel_peak) > params.duration_s/2 * 1.5) = 0; 
                    waveform = amplitude_scale * spike_shape_raw;
                otherwise
                    error('Unknown spike shape: %s', params.shape);
            end

        case 'ar' % Autoregressive model
            defaults.ar_coeffs = [0.5, -0.2]; % Example: AR(2) coeffs for single channel
            defaults.innovation_variance = 1.0;
            defaults.burn_in_samples = 1000;
            params = merge_structs(waveformParams, defaults);

            if ~isfield(params, 'ar_coeffs') || isempty(params.ar_coeffs)
                error('AR type requires .ar_coeffs in waveformParams.');
            end
            
            ar_coeffs = params.ar_coeffs;
            if size(ar_coeffs,1) > 1 && size(ar_coeffs,2) > 1 % Multivariate AR not supported for single channel output
                error('For single channel AR waveform, ar_coeffs should be [1x1xP] or [Px1] or [1xP].');
            end
            % Ensure ar_coeffs is a row vector for filter command [1, -a1, -a2, ...]
            if size(ar_coeffs,1) > 1, ar_coeffs = ar_coeffs'; end % make row
            filter_coeffs = [1, -ar_coeffs]; 
            order_p = length(ar_coeffs);

            total_samples_needed = Ntimepoints + params.burn_in_samples;
            innovation_noise = sqrt(params.innovation_variance) * randn(1, total_samples_needed);
            
            ar_signal_full = filter(1, filter_coeffs, innovation_noise);
            
            if length(ar_signal_full) > params.burn_in_samples
                ar_signal = ar_signal_full(params.burn_in_samples + 1 : params.burn_in_samples + Ntimepoints);
            else
                ar_signal = ar_signal_full(1:Ntimepoints); % Should not happen if burn_in is reasonable
                warning('seal_generate_waveform:ARBurnInTooLong', 'Burn-in samples might exceed generated AR signal length.');
            end
            
            % Normalize and scale
            if max(abs(ar_signal)) > eps
                waveform = amplitude_scale * (ar_signal / max(abs(ar_signal)));
            else
                waveform = zeros(1, Ntimepoints); % If signal is zero
            end


        case {'resting_alpha', 'resting_beta', 'resting_theta', 'resting_delta', 'resting_gamma'}
            type_str = lower(waveformType);
            defaults.center_frequency_Hz = 10; 
            defaults.bandwidth_Hz = 4;       
            if contains(type_str, 'beta'), defaults.center_frequency_Hz = 20; defaults.bandwidth_Hz = 8;
            elseif contains(type_str, 'theta'), defaults.center_frequency_Hz = 6; defaults.bandwidth_Hz = 3;
            elseif contains(type_str, 'delta'), defaults.center_frequency_Hz = 2; defaults.bandwidth_Hz = 2;
            elseif contains(type_str, 'gamma'), defaults.center_frequency_Hz = 40; defaults.bandwidth_Hz = 10; 
            end
            defaults.burst_properties = struct('duty_cycle', 0.5, 'min_duration_s', 0.5, 'max_duration_s', 1.5);
            params = merge_structs(waveformParams, defaults);
            if isfield(waveformParams, 'burst_properties') && isstruct(waveformParams.burst_properties)
                params.burst_properties = merge_structs(waveformParams.burst_properties, defaults.burst_properties);
            end

            pass_low = params.center_frequency_Hz - params.bandwidth_Hz/2;
            pass_high = params.center_frequency_Hz + params.bandwidth_Hz/2;
            if pass_low <= 0, pass_low = eps; end 
            if pass_high >= samplingRate/2, pass_high = samplingRate/2 - eps; end
            if pass_low >= pass_high, pass_low = pass_high - max(eps, pass_high/100); end % ensure pass_low < pass_high

            base_noise = seal_generate_pink_noise_channel(Ntimepoints);
            
            order = 4; nyquist = samplingRate / 2;
            passband_norm = [pass_low pass_high] / nyquist;
            
            if any(passband_norm <=0) || any(passband_norm >=1) || passband_norm(1) >= passband_norm(2)
                warning('seal_generate_waveform:RestingStateFilterBand', 'Problem with resting state filter band for %s. Using broadband pink noise.', waveformType);
                filtered_wave = base_noise;
            else
                [b_filt, a_filt] = butter(order, passband_norm, 'bandpass');
                filtered_wave = filter(b_filt, a_filt, base_noise);
            end

            if isfield(params, 'burst_properties') && isstruct(params.burst_properties) && ...
               isfield(params.burst_properties, 'duty_cycle') 
                bp = params.burst_properties;
                mod_envelope = zeros(1, Ntimepoints); 
                t_idx = 1;
                while t_idx <= Ntimepoints
                    if rand() < bp.duty_cycle 
                        burst_dur_s = bp.min_duration_s + rand() * (bp.max_duration_s - bp.min_duration_s);
                        burst_dur_samples = round(burst_dur_s * samplingRate);
                        end_burst_idx = min(Ntimepoints, t_idx + burst_dur_samples - 1);
                        mod_envelope(t_idx:end_burst_idx) = 1; 
                        t_idx = end_burst_idx + 1;
                    else 
                        gap_dur_s = bp.min_duration_s + rand() * (bp.max_duration_s - bp.min_duration_s); 
                        gap_dur_samples = round(gap_dur_s * samplingRate);
                        end_gap_idx = min(Ntimepoints, t_idx + gap_dur_samples -1);
                        t_idx = end_gap_idx + 1;
                    end
                end
                filtered_wave = filtered_wave .* mod_envelope;
            end
            waveform = amplitude_scale * (filtered_wave / (max(abs(filtered_wave)) + eps));

        case 'resting_pinknoise'
            pink_wave = seal_generate_pink_noise_channel(Ntimepoints);
            waveform = amplitude_scale * (pink_wave / (max(abs(pink_wave)) + eps)); 

        otherwise
            error('seal_generate_waveform:UnknownType', 'Unknown waveformType: %s', waveformType);
    end
end

% pink noise
function noise_channel = seal_generate_pink_noise_channel(n_samples)
    if n_samples == 0, noise_channel = []; return; end
    n_fft = 2^nextpow2(n_samples); 
    
    f = linspace(0, 1, n_fft/2 + 1); 
    inv_f_spectrum = zeros(size(f));
    inv_f_spectrum(f > 0) = 1./sqrt(f(f > 0)); 
    inv_f_spectrum(1) = 0; 

    random_phases = exp(1i * 2 * pi * rand(size(f)));
    half_spectrum = inv_f_spectrum .* random_phases;
    
    full_spectrum = zeros(1, n_fft);
    full_spectrum(1:n_fft/2+1) = half_spectrum;
    if n_fft/2+2 <= n_fft % Ensure indices are valid for fliplr
        full_spectrum(n_fft/2+2:end) = conj(fliplr(half_spectrum(2:end-1)));
    end

    noise_channel_temp = real(ifft(full_spectrum, n_fft));
    noise_channel = noise_channel_temp(1:n_samples); 
    noise_channel = noise_channel(:)'; 
end

% merge default and user-provided structs
function S_out = merge_structs(S_user, S_default)
    S_out = S_default;
    if ~isstruct(S_user) || isempty(S_user) 
        return;
    end
    user_fields = fieldnames(S_user);
    for i = 1:length(user_fields)
        field = user_fields{i};
        if isfield(S_default, field) && isstruct(S_default.(field)) && isstruct(S_user.(field))
            S_out.(field) = merge_structs(S_user.(field), S_default.(field));
            S_out.(field) = S_user.(field);
        else
            S_out.(field) = S_user.(field);
        end
    end
end
