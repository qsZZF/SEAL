function [S, active_indices, metadata] = seal_generate_source_activity(cortex, ROIs, time_vector, correlated_setting, options)
%SEAL_GENERATE_SOURCE_ACTIVITY Generates simulated source activity on a cortical surface.
%   [S, active_indices] = seal_generate_source_activity(cortex, ROIs, time_vector)
%   creates a source activity matrix 'S' by simulating one or more active
%   regions (patches) on a cortical surface model.
%
%   Inputs:
%       cortex (struct):
%           A cortical surface model structure. It must contain the field:
%           .Vertices - [n_vertices x 3] matrix of vertex coordinates.
%           .Faces - [n_vertices x 3] matrix of connected vertex
%           .VertConn (optional) - [n_vertices x n_vertices] adjacent matrix
%
%       ROIs (struct array):
%           An array of structs, where each struct defines a Region of Interest (ROI).
%           Each struct must contain the following fields:
%           .seedvox (integer) - Vertex index for the center of the ROI.
%           .extent (double) - Active patch area in mm^2 by default.
%                              Set .extent_type='radius_mm' to interpret
%                              extent as a radius and convert it to area.
%           .patch_area_mm2 (double, optional) - Explicit active patch area.
%           .waveform_params (struct) - Parameters for the temporal waveform,
%                                       which are passed to the local
%                                       'generate_waveform' function.
%                                       e.g., params.waveformtype = 'erp';
%                                             params.samplingRate = 100;
%                                             params.component.amplitude = 10e-9;
%                                             params.component.frequency = 10;
%                                             params.component.phase = 0;
%                                             params.component.decay = 50 (mm);
%
%       time_vector (vector):
%           A 1xN vector defining the time points for the simulation in seconds.
%
%       correlated_setting (struct):
%           A simulation parameter structure. It must contain the field:
%           .type (char/string) - 'correlated'/'causal', for simulating
%           activities with specific inter-ROI correlation / causality
%           .Coef (matrix) - a matrix specifying correlation/causality
%           coefficients.
%   Outputs:
%       S (matrix):
%           The generated [n_vertices x n_time_points] source activity matrix.
%           Units depend on the amplitude set in waveform_params (e.g., Am).
%
%       active_indices (cell):
%           A cell array where each cell contains the vector of vertex indices
%           corresponding to an active patch defined in ROIs.
%
%   See also: seal_generate_waveform, seal_load_source_activity.

if nargin < 4
    correlated_setting = [];
end
if nargin < 5 || isempty(options)
    options = struct();
end

restore_rng = false;
if isfield(options, 'random_seed') && ~isempty(options.random_seed)
    old_rng = rng;
    rng(options.random_seed);
    restore_rng = true;
end

% --- Input Validation and Initialization ---
if ~isfield(cortex, 'Vertices')
    error('The cortex structure must contain a .Vertices field.');
end
if ~isstruct(ROIs) || ~all(isfield(ROIs, {'seedvox', 'extent', 'waveform_params'}))
    error('ROIs must be a struct array with fields: center, extent, and waveform_params.');
end
if nargin>3 && ~isempty(correlated_setting)
    simu_corr = 1;
else
    simu_corr = 0;
end

n_sources = size(cortex.Vertices, 1);
n_time_points = length(time_vector);

% Initialize the main source activity matrix with zeros
S = zeros(n_sources, n_time_points);
active_indices = cell(1, numel(ROIs));
spatial_weights = cell(1, numel(ROIs));
[~, VertArea] = tess_area1(cortex);
coord_to_mm = seal_coordinate_scale_to_mm(cortex.Vertices);
VertArea_mm2 = VertArea * coord_to_mm^2;

% --- Main Loop: Iterate Through Each ROI ---

simulated_series = zeros(numel(ROIs),n_time_points);

for i = 1:numel(ROIs)
    % Get parameters for the current ROI
    roi = ROIs(i);
    center_vertex_idx = roi.seedvox;
    patch_area_mm2 = seal_resolve_patch_area_mm2(roi);

    % --- 1. Spatial Patch Generation ---

    % Store the indices of the active patch
    patch_indices = PatchGenerate(center_vertex_idx, cortex.VertConn, VertArea_mm2, patch_area_mm2);
    active_indices{i} = patch_indices;
    
    
    % --- 2. Temporal Waveform Generation ---
    % Generate the time course for this ROI
    waveform_params = roi.waveform_params;
    waveform = seal_generate_waveform(time_vector, waveform_params.samplingRate, waveform_params.waveformtype, waveform_params);
    simulated_series(i,:) = waveform;

    % --- 3. Spatial Weight Calculation ---
    % Weight is always 1 at the seed vertex. The selected decay model controls
    % how the waveform attenuates across the generated patch.
    spatial_weights{i} = ones(length(patch_indices), 1);
    decay_model = 'gaussian';
    if isfield(waveform_params, 'decay_model') && ~isempty(waveform_params.decay_model)
        decay_model = lower(waveform_params.decay_model);
    end

    if isfield(waveform_params, 'decay') && ~isempty(waveform_params.decay) && ...
            ~strcmpi(decay_model, 'flat')
        center_coords = cortex.Vertices(center_vertex_idx, :);
        patch_coords = cortex.Vertices(patch_indices, :);
        patch_distances = sqrt(sum((patch_coords - center_coords).^2, 2)) * coord_to_mm;
        decay_distance = max(waveform_params.decay, eps);

        switch decay_model
            case 'gaussian'
                sigma2 = decay_distance.^2 / log(2);
                spatial_weights{i} = exp(-patch_distances.^2 / sigma2);
            case 'linear'
                spatial_weights{i} = max(0, 1 - patch_distances / decay_distance);
            case 'exponential'
                tau = decay_distance / log(2);
                spatial_weights{i} = exp(-patch_distances / tau);
            otherwise
                warning('seal_generate_source_activity:UnknownDecayModel', ...
                    'Unknown decay_model "%s". Falling back to Gaussian.', decay_model);
                sigma2 = decay_distance.^2 / log(2);
                spatial_weights{i} = exp(-patch_distances.^2 / sigma2);
        end
    end
    % Add the activity of the current patch to the main source matrix 'S'.
    % Using addition allows for the linear superposition of overlapping patches.
%     S(patch_indices, :) = S(patch_indices, :) + waveform;
end

if simu_corr == 1
    switch correlated_setting.type
        case 'correlated'
            corr_mat = seal_make_positive_definite(correlated_setting.Coef, numel(ROIs));
            simulated_series = seal_apply_correlation(simulated_series, corr_mat);
            
        case 'causal'
            simulated_series = seal_apply_causal_mixing(simulated_series, correlated_setting.Coef);
        otherwise
            warning('undefined simulation type (correlated/causal), use default settings');
            
    end
    
end
for i = 1:numel(ROIs)
    patch_indices = active_indices{i};
    waveform = simulated_series(i,:);
    weight = spatial_weights{i};
    waveform = repmat(waveform,length(patch_indices),1);
    S(patch_indices, :) = S(patch_indices, :) + weight.*waveform;
end

metadata = struct();
metadata.n_sources = n_sources;
metadata.n_time_points = n_time_points;
metadata.coordinate_unit_scale_to_mm = coord_to_mm;
metadata.random_seed = [];
if isfield(options, 'random_seed')
    metadata.random_seed = options.random_seed;
end
metadata.seed_vertices = arrayfun(@(x) x.seedvox, ROIs);
metadata.active_indices = active_indices;
metadata.spatial_weights = spatial_weights;
metadata.correlated_setting = correlated_setting;

if restore_rng
    rng(old_rng);
end

end

function coord_to_mm = seal_coordinate_scale_to_mm(vertices)
% Treat small-coordinate surfaces as meters and convert distances to mm.
if max(abs(vertices(:))) < 10
    coord_to_mm = 1000;
else
    coord_to_mm = 1;
end
end

function patch_area_mm2 = seal_resolve_patch_area_mm2(roi)
if isfield(roi, 'patch_area_mm2') && ~isempty(roi.patch_area_mm2)
    patch_area_mm2 = roi.patch_area_mm2;
    return;
end
if isfield(roi, 'extent_type') && any(strcmpi(roi.extent_type, {'radius', 'radius_mm'}))
    patch_area_mm2 = pi * roi.extent^2;
else
    patch_area_mm2 = roi.extent;
end
end

function corr_mat = seal_make_positive_definite(coef, n_roi)
corr_mat = double(coef);
if ~isequal(size(corr_mat), [n_roi, n_roi])
    error('correlated_setting.Coef must be nROI x nROI.');
end
corr_mat = (corr_mat + corr_mat') / 2;
corr_mat(1:n_roi+1:end) = 1;
[V, D] = eig(corr_mat);
eig_vals = max(diag(D), 1e-6);
corr_mat = V * diag(eig_vals) * V';
d = sqrt(diag(corr_mat));
corr_mat = corr_mat ./ (d * d');
end

function mixed = seal_apply_correlation(series, corr_mat)
mu = mean(series, 2);
sigma = std(series, 0, 2);
sigma(sigma < eps) = 1;
z = (series - mu) ./ sigma;
L = chol(corr_mat, 'lower');
mixed = L * z;
mixed = mixed .* sigma + mu;
end

function mixed = seal_apply_causal_mixing(series, coef)
n_roi = size(series, 1);
causal_mat = double(coef);
if ~isequal(size(causal_mat), [n_roi, n_roi])
    error('causal correlated_setting.Coef must be nROI x nROI.');
end
spec = max(abs(eig(causal_mat)));
if spec >= 0.95
    causal_mat = causal_mat * (0.95 / spec);
end
mixed = zeros(size(series));
mixed(:, 1) = series(:, 1);
for t = 2:size(series, 2)
    mixed(:, t) = series(:, t) + causal_mat * mixed(:, t - 1);
end
mu = mean(series, 2);
sigma = std(series, 0, 2);
out_sigma = std(mixed, 0, 2);
out_sigma(out_sigma < eps) = 1;
mixed = (mixed - mean(mixed, 2)) ./ out_sigma .* sigma + mu;
end


function [FaceArea, VertArea] = tess_area1(cortex)
% Compute the surface area associated with each face and each vertex.
% Brainstorm source: tess_area.m

% Compute the area of all the faces
r12 = cortex.Vertices(cortex.Faces(:,1),:);        % temporary holding
r13 = cortex.Vertices(cortex.Faces(:,3),:) - r12;  % negative of r31
r12 = cortex.Vertices(cortex.Faces(:,2),:) - r12;  % from 1 to 2
FaceArea = sqrt(sum(cross(r12,r13,2).^2, 2)) / 2;

% Compute the triangle area only if needed
if (nargout >= 2)
    % Build vertex-face connectivity matrix, with the area information
    nFaces = size(cortex.Faces,1);
    rowno = double([cortex.Faces(:,1); cortex.Faces(:,2); cortex.Faces(:,3)]);
    colno = [1:nFaces, 1:nFaces, 1:nFaces]';
    data  = [FaceArea; FaceArea; FaceArea];
    VertFacesArea = sparse(rowno,colno,data);

    % Compute the vertex area: 1/3 of each triangle involving this vertex
    VertArea = 1/3 * full(sum(VertFacesArea,2));
end
    
end

function Patch = PatchGenerate(Seed, VertConn, VertArea, AreaDef)
% Grows a patch of vertices from a seed until a defined area is reached.

Patch = Seed;
Area = sum(VertArea(Patch));
while Area <= AreaDef
    % Find the next layer of connected vertices
    newverts = tess_scout_swell(Patch, VertConn);
    % If there are no more neighbors, break the loop
    if isempty(newverts)
        break;
    end
    
    Nouter = union(Patch, newverts);
    Area = sum(VertArea(Nouter));

    if Area > AreaDef
        % If adding the whole layer overshoots the target area,
        % add vertices from the new layer one-by-one until the
        % area is just exceeded.
        Ndiff = setdiff(Nouter, Patch);
        for i = 1:numel(Ndiff)
             Patch = union(Patch, Ndiff(i));
             Area = sum(VertArea(Patch));
             if Area > AreaDef
                break;
             end
        end
    else
        % If it doesn't overshoot, add the whole layer and continue
        Patch = Nouter;
    end
end
end

function newverts = tess_scout_swell(iverts, vconn)
% TESS_SCOUT_SWELL: Enlarge a patch by appending the next set of adjacent vertices.
% Brainstorm source: tess_scout_swell.m

if (size(iverts,1) ~= 1)
  iverts = iverts(:)'; % ensure row vector
end

% Concatenate all vertex connections for all verts in the patch
newverts = find(max(vconn(iverts, :), [], 1));

% Extract unique set of indices, remove existing vertices
newverts = setdiff(newverts, iverts);
end
