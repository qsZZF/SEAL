function [S, active_indices] = seal_generate_source_activity(cortex, ROIs, time_vector, correlated_setting)
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
%           .extent (double) - Radius of the active patch in millimeters (mm).
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

% --- Main Loop: Iterate Through Each ROI ---

simulated_series = zeros(numel(ROIs),n_time_points);

for i = 1:numel(ROIs)
    % Get parameters for the current ROI
    roi = ROIs(i);
    center_vertex_idx = roi.seedvox;
    patch_extent = roi.extent;

    % --- 1. Spatial Patch Generation ---

    % Store the indices of the active patch
    patch_indices = PatchGenerate(center_vertex_idx, cortex.VertConn, VertArea, patch_extent);
    active_indices{i} = patch_indices;
    
    
    % --- 2. Temporal Waveform Generation ---
    % Generate the time course for this ROI
    waveform_params = roi.waveform_params;
    waveform = seal_generate_waveform(time_vector, waveform_params.samplingRate, waveform_params.waveformtype, waveform_params);
    simulated_series(i,:) = waveform;
    spatial_weights{i} = ones(size(waveform,1),1);
%     waveform = repmat(waveform,length(patch_indices),1);
    % --- 3. Activity Assignment with Gaussian Falloff ---
    % Calculate a standard deviation for the Gaussian spatial falloff.
    if isfield(waveform_params, 'decay')
     % Get coordinates of the center (seed) vertex
    center_coords = cortex.Vertices(center_vertex_idx, :);

    % Get coordinates of all vertices in the generated patch
    patch_coords = cortex.Vertices(patch_indices, :);

    % Calculate Euclidean distance from the center to all vertices in the patch
    patch_distances = sqrt(sum((patch_coords - center_coords).^2, 2));
    
    half_amplitude_dist = waveform_params.decay;
    sigma2 = half_amplitude_dist.^2/log(2);
    
    % Calculate the Gaussian weights for a smooth spatial falloff
    % The weight is 1 at the center and decreases towards the edge.
    gaussian_weights = exp(-patch_distances.^2 / sigma2);
    
    % Combine the spatial weights with the temporal waveform
    % The result is a matrix where each row is the scaled waveform for one vertex
    spatial_weights{i} = gaussian_weights;
%     waveform = gaussian_weights .* waveform;
    end
    % Add the activity of the current patch to the main source matrix 'S'.
    % Using addition allows for the linear superposition of overlapping patches.
%     S(patch_indices, :) = S(patch_indices, :) + waveform;
end

if simu_corr == 1
    switch correlated_setting.type
        case 'correlated'
            corr_mat = chol(correlated_setting.Coef)';
            simulated_series = orth(simulated_series')';
            simulated_series = corr_mat*simulated_series;
            
        case 'causal'
            causal_mat = correlated_setting.Coef;
            spec = svd(causal_mat); spec = max(abs(spec));
            if spec >= 1
                causal_mat = causal_mat*0.95 / spec;
            end
            simulated_series = causal_mat*simulated_series;
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