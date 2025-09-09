function [S_thresholded, Param_AT] = seal_automatic_thresholding(S_amp, G_adj, varargin)
% SEAL_AUTOMATIC_THRESHOLDING Implements the SISSY automatic thresholding.
%
%   Performs a three-stage process to delineate active source regions from a
%   source amplitude map, based on the method described in the SISSY paper
%   by Becker et al., 2017.
%
%   Usage:
%   [S_thresh, P] = seal_automatic_thresholding(S_amp, G_adj, 'final_threshold_ratio', 0.1);
%
%   Inputs:
%       S_amp (double vector): Source amplitudes (N_sources x 1). This is
%                              typically the L2-norm of the source time series
%                              at the time point of interest.
%       G_adj (sparse matrix): Adjacency matrix for the cortical mesh
%                              (N_sources x N_sources). G_adj(i, j) = 1 if
%                              sources i and j are neighbors.
%
%   Optional Name-Value Pair Inputs:
%       'final_threshold_ratio' (double): Relative threshold for the final
%                                         delineation step. Default: 0.1 (10%).
%       'max_merge_iter' (integer): A safeguard to prevent infinite loops
%                                   during region merging. Default: 5000.
%
%   Outputs:
%       S_thresholded (double vector): A vector of source amplitudes where
%                                      inactive sources are set to zero.
%       Param_AT (struct): A struct containing the intermediate results,
%                          including the initial and final region labels.
%
% Copyright 2025 SEAL (Source Electromagnetics & Analysis Lab)

    %% --- 1. Input Parsing ---
    p = inputParser;
    addRequired(p, 'S_amp', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'G_adj', @(x) issparse(x) && ismatrix(x));
    addParameter(p, 'final_threshold_ratio', 0.1, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    addParameter(p, 'max_merge_iter', 5000, @(x) isnumeric(x) && isscalar(x) && x>0);
    parse(p, S_amp, G_adj, varargin{:});

    final_threshold_ratio = p.Results.final_threshold_ratio;
    max_merge_iter = p.Results.max_merge_iter;
    
    fprintf('Starting Automatic Thresholding...\n');

    %% --- STEP 1: Watershed Segmentation ---
    fprintf('  Step 1: Performing graph-based watershed segmentation...\n');
    initial_labels = local_graph_watershed(S_amp, G_adj);
    Param_AT.InitialLabels = initial_labels;
    fprintf('    Found %d initial regions.\n', length(unique(initial_labels)));

    %% --- STEP 2: Iterative Region Merging ---
    fprintf('  Step 2: Performing iterative region merging...\n');
    final_labels = local_iterative_merging(initial_labels, S_amp, G_adj, max_merge_iter);
    Param_AT.FinalLabels = final_labels;
    fprintf('    Finished with %d final regions.\n', length(unique(final_labels)));
    
    %% --- STEP 3: Final Thresholding and Delineation ---
    fprintf('  Step 3: Delineating final active regions...\n');
    S_thresholded = local_final_thresholding(final_labels, S_amp, final_threshold_ratio);
    Param_AT.FinalThresholdRatio = final_threshold_ratio;
    
    fprintf('Automatic Thresholding complete.\n');
end

%--------------------------------------------------------------------------
% Local Helper Functions
%--------------------------------------------------------------------------

function labels = local_graph_watershed(amps, adj)
    % Implements a graph-based watershed algorithm.
    % Nodes are sorted by amplitude, and a Disjoint Set Union (DSU) data
    % structure is used to efficiently track the merging of catchment basins.
    
    num_nodes = length(amps);
    [sorted_amps, sorted_indices] = sort(amps, 'ascend');
    
    labels = zeros(num_nodes, 1);
    
    % DSU data structure: parent(i) points to the parent of element i.
    parent = 1:num_nodes;
    
    % Process nodes in increasing order of amplitude
    for i = 1:num_nodes
        node_idx = sorted_indices(i);
        
        % Find neighbors of the current node
        neighbors = find(adj(node_idx, :));
        
        % Check neighbors that have already been processed (and thus labeled)
        processed_neighbors = neighbors(labels(neighbors) > 0);
        
        if isempty(processed_neighbors)
            % This node is a new local minimum, start a new basin
            labels(node_idx) = node_idx;
        else
            % Find the unique labels of the neighbors
            neighbor_roots = unique(arrayfun(@(n) find_root(n, parent), processed_neighbors));
            
            if length(neighbor_roots) == 1
                % All processed neighbors belong to the same basin. Join it.
                root = neighbor_roots(1);
                labels(node_idx) = labels(root);
                parent(node_idx) = root;
            else
                % This node is a ridge between multiple basins.
                % For this simplified watershed, we don't create explicit
                % watershed lines; we just join the first one we find.
                % A more complex implementation could mark this as a border.
                root = neighbor_roots(1);
                labels(node_idx) = labels(root);
                parent(node_idx) = root;
                % And merge the other basins into this one
                for j = 2:length(neighbor_roots)
                    parent(neighbor_roots(j)) = root;
                end
            end
        end
    end
    
    % Finalize labels by ensuring all nodes point to their final root
    unique_roots = unique(arrayfun(@(n) find_root(n, parent), 1:num_nodes));
    root_map = containers.Map(unique_roots, 1:length(unique_roots));
    
    final_labels = zeros(num_nodes, 1);
    for i = 1:num_nodes
        final_labels(i) = root_map(find_root(i, parent));
    end
    labels = final_labels;

    % Helper for DSU
    function root = find_root(i, p)
        while i ~= p(i)
            p(i) = p(p(i)); % Path compression
            i = p(i);
        end
        root = i;
    end
end

function final_labels = local_iterative_merging(initial_labels, amps, adj, max_iter)
    % Iteratively merges regions based on the paper's criteria.
    current_labels = initial_labels;
    
    for iter = 1:max_iter
        % Get properties of all current regions
        region_props = local_get_region_properties(current_labels, amps, adj);
        
        % Get a sorted list of all potential merges
        merge_candidates = local_get_all_merge_candidates(region_props);
        
        if isempty(merge_candidates)
            fprintf('No more adjacent regions to merge. Stopping.\n');
            break;
        end
        
        merged_this_pass = false;
        for i = 1:length(merge_candidates)
            candidate = merge_candidates(i);
            
            % Check if this merge is allowed
            is_allowed = local_check_merge_validity(candidate.region1_idx, candidate.region2_idx, region_props, amps, adj);
            
            if is_allowed
                % Perform the merge
                label1 = str2double(region_props(candidate.region1_idx).label);
                label2 = str2double(region_props(candidate.region2_idx).label);
                current_labels(current_labels == label2) = label1;
                
                fprintf('Iter %d: Merged region %d into %d (diff=%.2e).\n', iter, label2, label1, candidate.diff);
                merged_this_pass = true;
                break; % Restart with the new set of labels
            end
        end
        
        if ~merged_this_pass
            fprintf('No valid merges found in this pass. Stopping.\n');
            break; % Finished if a full pass results in no valid merges
        end
    end
    
    if iter == max_iter
        warning('AT:MaxIterReached', 'Region merging reached max iterations.');
    end
    
    final_labels = current_labels;
end

function region_props = local_get_region_properties(labels, amps, adj)
    % Calculates properties for each unique region/label.
    unique_labels = unique(labels);
    num_regions = length(unique_labels);
    
    % Pre-allocate a struct array
    region_props(num_regions) = struct('label', [], 'member_indices', [], 'avg_amp', [], 'adj_regions', []);

    for i = 1:num_regions
        l = unique_labels(i);
        member_indices = find(labels == l);
        
        region_props(i).label = num2str(l);
        region_props(i).member_indices = member_indices;
        region_props(i).avg_amp = mean(amps(member_indices));
        
        % Find adjacent regions
        region_adj_mask = any(adj(member_indices, :), 1);
        adj_labels = unique(labels(region_adj_mask));
        region_props(i).adj_regions = setdiff(adj_labels, l);
    end
end

function candidates = local_get_all_merge_candidates(region_props)
    % Finds all pairs of adjacent regions and sorts them by amplitude difference.
    candidates = [];
    num_regions = length(region_props);
    
    for i = 1:num_regions
        for j = 1:length(region_props(i).adj_regions)
            adj_label = region_props(i).adj_regions(j);
            % Find the index in region_props for this adjacent label
            adj_idx = find(strcmp({region_props.label}, num2str(adj_label)));

            % Avoid duplicates (e.g., A-B and B-A)
            if i < adj_idx
                diff = abs(region_props(i).avg_amp - region_props(adj_idx).avg_amp);
                candidates = [candidates; struct('region1_idx', i, 'region2_idx', adj_idx, 'diff', diff)];
            end
        end
    end
    
    if ~isempty(candidates)
        [~, sort_idx] = sort([candidates.diff]);
        candidates = candidates(sort_idx);
    end
end

function is_allowed = local_check_merge_validity(r1_idx, r2_idx, region_props, amps, adj)
    % Checks if merging r1 and r2 would reduce the number of local extrema.
    
    % 1. Count extrema BEFORE the merge
    [extrema_count_before, ~] = local_count_extrema(region_props);
    
    % 2. Create hypothetical new labels after merge
    labels_hypothetical = zeros(length(amps), 1);
    num_regions = length(region_props);
    for i = 1:num_regions
        l_val = str2double(region_props(i).label);
        labels_hypothetical(region_props(i).member_indices) = l_val;
    end
    label1 = str2double(region_props(r1_idx).label);
    label2 = str2double(region_props(r2_idx).label);
    labels_hypothetical(labels_hypothetical == label2) = label1;
    
    % 3. Get properties and count extrema AFTER hypothetical merge
    props_hypothetical = local_get_region_properties(labels_hypothetical, amps, adj);
    [extrema_count_after, ~] = local_count_extrema(props_hypothetical);
    
    % 4. The merge is allowed ONLY if the number of extrema does not decrease
    is_allowed = (extrema_count_after >= extrema_count_before);
end

function [extrema_count, is_extrema] = local_count_extrema(region_props)
    % Counts the number of local maxima and minima among the regions.
    num_regions = length(region_props);
    is_max = false(num_regions, 1);
    is_min = false(num_regions, 1);
    
    for i = 1:num_regions
        current_amp = region_props(i).avg_amp;
        adj_region_labels = region_props(i).adj_regions;
        
        if isempty(adj_region_labels)
            continue; % Cannot be an extremum if it has no neighbors
        end
        
        neighbor_amps = [];
        for j = 1:length(adj_region_labels)
             adj_idx = find(strcmp({region_props.label}, num2str(adj_region_labels(j))));
             neighbor_amps = [neighbor_amps, region_props(adj_idx).avg_amp];
        end
        
        if all(current_amp > neighbor_amps)
            is_max(i) = true;
        end
        if all(current_amp < neighbor_amps)
            is_min(i) = true;
        end
    end
    
    is_extrema = is_max | is_min;
    extrema_count = sum(is_extrema);
end

function S_out = local_final_thresholding(labels, amps, ratio)
    % Applies the final relative amplitude threshold to the merged regions.
    region_props = local_get_region_properties(labels, amps, []); % No need for adjacency here
    
    if isempty(region_props)
        S_out = zeros(size(amps));
        return;
    end
    
    avg_amps = [region_props.avg_amp];
    max_avg_amp = max(avg_amps);
    
    if max_avg_amp == 0
        S_out = zeros(size(amps));
        return;
    end
    
    threshold_val = ratio * max_avg_amp;
    
    final_active_mask = false(size(amps));
    for i = 1:length(region_props)
        if region_props(i).avg_amp >= threshold_val
            final_active_mask(region_props(i).member_indices) = true;
        end
    end
    
    S_out = zeros(size(amps));
    S_out(final_active_mask) = amps(final_active_mask);
end
