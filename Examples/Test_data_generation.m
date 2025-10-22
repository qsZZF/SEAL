%% Test Script for seal_generate_source_activity



% ---  Define Simulation Parameters ---
fprintf('Defining simulation parameters...\n');
cortex = Cortex;
% Time vector
Fs = 500; % Sampling frequency in Hz
t_start = -0.2; % seconds
t_end = 0.8; % seconds
time_vector = t_start : 1/Fs : t_end;

% Define Regions of Interest (ROIs)
ROIs = [];

% ROI 1: An ERP-like response in the "frontal" area
% This ROI uses the user-defined 'decay' parameter for spatial falloff.
roi1.seedvox = 150; % An arbitrary vertex index
roi1.extent = 10*1e-4; % Target patch area in cm^2
roi1.waveform_params.waveformtype = 'erp';
roi1.waveform_params.amplitude = 15e-9; % 15 nAm
roi1.waveform_params.decay = 10; % Amp drops to 50% at 10mm from seed
roi1.waveform_params.samplingRate = Fs;
erp_components(1) = struct('type', 'gaussian', 'latency_s', 0.100, 'amplitude_rel', 1.0, 'width_s', 0.02);
erp_components(2) = struct('type', 'gaussian', 'latency_s', 0.180, 'amplitude_rel', -0.5, 'width_s', 0.04);
roi1.waveform_params.component = erp_components;
ROIs = [ROIs, roi1];

% ROI 2: An oscillatory burst in the "parietal" area
% This ROI uses the fallback method for spatial falloff (no .decay field).
roi2.seedvox = 850; % A different vertex index
roi2.extent = 20*1e-4; % A larger patch
roi2.waveform_params.waveformtype = 'ersp';
roi2.waveform_params.samplingRate = Fs;
roi2.waveform_params.amplitude = 20e-9; % 20 nAm
roi2.waveform_params.base_frequency_band_Hz = 15; % Beta band
mod_params.type = 'burst';
mod_params.latency_s = 0.4;
mod_params.width_s = 0.3;
roi2.waveform_params.component = mod_params;
ROIs = [ROIs, roi2];

% --- 4. Run Source Generation ---
fprintf('Running seal_generate_source_activity...\n');
[S, active_indices] = seal_generate_source_activity(cortex, ROIs, time_vector);
fprintf('Simulation complete.\n');

% --- 5. Validation and Visualization ---
fprintf('Validating outputs and generating plots...\n');



% Figure 1: Visualize the generated patches on the cortex
figure('Name', 'Generated Simulation Patches', 'NumberTitle', 'off', 'Color', 'w');
h_patch = trisurf(cortex.Faces, cortex.Vertices(:,1), cortex.Vertices(:,2), cortex.Vertices(:,3), ...
    'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
axis equal;
lighting gouraud;
material dull;
camlight;
view(-90, 15);

% Plot ROI 1 patch in red
patch_verts1 = cortex.Vertices(active_indices{1}, :);
plot3(patch_verts1(:,1), patch_verts1(:,2), patch_verts1(:,3), '.r', 'MarkerSize', 20);

% Plot ROI 2 patch in blue
patch_verts2 = cortex.Vertices(active_indices{2}, :);
plot3(patch_verts2(:,1), patch_verts2(:,2), patch_verts2(:,3), '.b', 'MarkerSize', 20);

title('Spatial Location of Generated Patches');
legend({'Cortex', 'ROI 1 (ERP)', 'ROI 2 (Beta Burst)'});


% Figure 2: Visualize the source activity map at peak times
figure('Name', 'Source Activity Maps', 'NumberTitle', 'off', 'Color', 'w');

% Subplot 1: Peak of ROI 1 (ERP)
subplot(1, 2, 1);
peak_time_roi1 = 0.100; % Time of N100 peak
[~, time_idx1] = min(abs(time_vector - peak_time_roi1));
activity1 = S(:, time_idx1);

PlotSource(activity1,Cortex);

trisurf(cortex.Faces, cortex.Vertices(:,1), cortex.Vertices(:,2), cortex.Vertices(:,3), activity1, ...
    'EdgeColor', 'none');
axis equal; view(-90, 15); colormap('hot'); colorbar;
title(sprintf('Activity at t = %.3f s (ROI 1 Peak)', time_vector(time_idx1)));
material dull; camlight;

% Subplot 2: Peak of ROI 2 (Beta Burst)
subplot(1, 2, 2);
peak_time_roi2 = 0.400; % Time of burst center
[~, time_idx2] = min(abs(time_vector - peak_time_roi2));
activity2 = S(:, time_idx2);

trisurf(cortex.Faces, cortex.Vertices(:,1), cortex.Vertices(:,2), cortex.Vertices(:,3), activity2, ...
    'EdgeColor', 'none');
axis equal; view(-90, 15); colormap('hot'); colorbar;
title(sprintf('Activity at t = %.3f s (ROI 2 Peak)', time_vector(time_idx2)));
material dull; camlight;


% Figure 3: Plot the time courses of the seed vertices
figure('Name', 'Generated Time Courses at Seed Vertices', 'NumberTitle', 'off', 'Color', 'w');

% Time course for ROI 1 seed
subplot(2, 1, 1);
plot(time_vector, S(roi1.seedvox, :) * 1e9, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Time Course for ROI 1 (Seed Vertex %d)', roi1.seedvox));
xlabel('Time (s)');
ylabel('Amplitude (nAm)');
xlim([t_start t_end]);

% Time course for ROI 2 seed
subplot(2, 1, 2);
plot(time_vector, S(roi2.seedvox, :) * 1e9, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Time Course for ROI 2 (Seed Vertex %d)', roi2.seedvox));
xlabel('Time (s)');
ylabel('Amplitude (nAm)');
xlim([t_start t_end]);

fprintf('Test script finished successfully.\n');
