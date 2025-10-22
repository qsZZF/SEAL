% DEMO_generate_all_waveforms.m
%
% This script demonstrates how to generate various types of temporal waveforms
% using the seal_generate_waveform function.

% clear; close all; clc;

% --- Common Parameters ---
samplingRate = 200; % Hz
duration_s = 2;     % seconds
timeVector = (0:1/samplingRate:duration_s - 1/samplingRate);
Ntimepoints = length(timeVector);

fprintf('Generating waveforms for %d timepoints at %d Hz...\n\n', Ntimepoints, samplingRate);

% --- Store waveforms and titles for plotting ---
generated_waveforms = {};
waveform_titles = {};

% --- 1. Event-Related Potential (ERP) ---
waveformType = 'erp';
erpParams.amplitude = 5; % Overall amplitude
erpParams.components = [
    struct('latency_s', 0.3, 'amplitude_rel', 1.0, 'width_s', 0.05), ... % P1-like
    struct('latency_s', 0.5, 'amplitude_rel', -0.7, 'width_s', 0.07), ...% N1-like
    struct('latency_s', 0.8, 'amplitude_rel', 1.5, 'width_s', 0.1)    % P2-like
];
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, erpParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (Amp: %.1f)', upper(waveformType), erpParams.amplitude);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end

% Generate Gabor-like waveform
erpParams.components = [
    struct('type', 'gabor', 'latency_s', 0.8, 'amplitude_rel', -0.8, 'frequency_Hz', 10, 'std_dev_s', 0.1, 'phase_rad', pi/2)
];
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, erpParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (Gabor, Amp: %.1f)', upper(waveformType), erpParams.amplitude);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end

% --- 2. Event-Related Spectral Perturbation (ERSP) - Alpha Burst ---
waveformType = 'ersp';
erspParams.base_frequency_band_Hz = 10; % Alpha
erspParams.base_amplitude = 1.5;      % Amplitude of the alpha carrier
erspParams.modulation = struct(...
    'type', 'burst', ...
    'latency_s', timeVector(round(Ntimepoints/2)), ...
    'width_s', 0.5, ...            % 500ms burst duration
    'taper', 0.25, ...
    'min_rel_amplitude', 0.05 ...   % Nearly off outside burst
);
erspParams.amplitude = 2.0; % Overall scaling of the final ERSP waveform
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, erspParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (Alpha Burst, Overall Amp: %.1f)', upper(waveformType), erspParams.amplitude);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end

% --- 3. Gaussian Pulse ---
waveformType = 'gaussian_pulse';
gpParams.amplitude = 3;
gpParams.pulse_center_time_s = 0.75;
gpParams.pulse_width_time_s = 0.1; % std dev in seconds
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, gpParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (Amp: %.1f, Center: %.2fs, Width: %.2fs)', upper(waveformType), gpParams.amplitude, gpParams.pulse_center_time_s, gpParams.pulse_width_time_s);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end

% --- 4. Spike ---
waveformType = 'spike';
spikeParams.amplitude = 10;
spikeParams.latency_s = 1.2;
spikeParams.duration_s = 0.08; % 80ms
spikeParams.shape = 'exponential_decay'; % or 'triangular'
% decay_tau_s will use default (duration_s/4)
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, spikeParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (Amp: %.1f, Lat: %.1fs, Shape: %s)', upper(waveformType), spikeParams.amplitude, spikeParams.latency_s, spikeParams.shape);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end

% --- 5. Autoregressive (AR) Model ---
waveformType = 'ar';
arParams.amplitude = 2.5;
arParams.ar_coeffs = [0.6, -0.3, 0.1]; % AR(3) model example coefficients (for a_1, a_2, a_3)
arParams.innovation_variance = 0.5;
arParams.burn_in_samples = 500;
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, arParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (AR(%d), Amp: %.1f)', upper(waveformType), length(arParams.ar_coeffs), arParams.amplitude);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end

% --- 6. Resting Alpha ---
waveformType = 'resting_alpha';
alphaParams.amplitude = 1.5;
alphaParams.center_frequency_Hz = 10; % Overriding default if needed
alphaParams.bandwidth_Hz = 3;       % Overriding default if needed
alphaParams.burst_properties = struct(...
    'duty_cycle', 0.6, ...
    'min_duration_s', 0.4, ...
    'max_duration_s', 1.2 ...
);
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, alphaParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (Amp: %.1f, Bursty)', upper(waveformType), alphaParams.amplitude);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end

% --- 7. Resting Pink Noise ---
waveformType = 'resting_pinknoise';
pinkParams.amplitude = 2;
try
    waveform = seal_generate_waveform(timeVector, samplingRate, waveformType, pinkParams);
    generated_waveforms{end+1} = waveform;
    waveform_titles{end+1} = sprintf('%s (Amp: %.1f)', upper(waveformType), pinkParams.amplitude);
catch ME
    fprintf('Error generating %s: %s\n', waveformType, ME.message);
end


% --- Plotting the Waveforms ---
num_waveforms = length(generated_waveforms);
if num_waveforms == 0
    disp('No waveforms were generated.');
    return;
end

figure('Name', 'Generated Waveforms Demo', 'NumberTitle', 'off', 'Position', [100, 100, 800, 200*ceil(num_waveforms/2)]);
cols = min(2, num_waveforms);
rows = ceil(num_waveforms/cols);

for i = 1:num_waveforms
    subplot(rows, cols, i);
    plot(timeVector, generated_waveforms{i});
    title(waveform_titles{i},'interpret','none');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    xlim([timeVector(1) timeVector(end)]);
end

if mod(num_waveforms, 2) == 1 && num_waveforms > 1 % Adjust layout if odd number and more than 1 plot
    fig_pos = get(gcf, 'Position');
    if rows > 1 % only adjust if it makes sense
      subplot(rows, cols, num_waveforms); % Get current axes
      current_ax_pos = get(gca, 'Position');
      % Make the last plot wider if it's alone in the last row
      % current_ax_pos(3) = current_ax_pos(3) * 1.8; % Increase width
      % set(gca, 'Position', current_ax_pos);
    end
end

sgtitle('Demonstration of Various Waveform Generations', 'FontSize', 14, 'FontWeight', 'bold');

disp('Demo script finished.');

% --- Helper function for this demo script (if not already on path) ---
% You would typically have seal_generate_waveform.m and its helpers in your MATLAB path.
% If seal_generate_waveform.m uses merge_structs as a separate helper, ensure it's available.
% For this demo, I'll assume seal_generate_waveform handles its own dependencies.