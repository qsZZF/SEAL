function Results = seal_evaluate_spatial_metrics(estimated_source_activity, true_active_indices, active_seedVox, sourceSpaceInfo, varargin)
%SEAL_EVALUATE_SPATIAL_METRICS Computes a suite of spatial performance metrics for ESI.
%   Results = seal_evaluate_spatial_metrics(estimated_source_activity, ...
%                                   true_active_indices, sourceSpaceInfo, 'ParameterName', ParameterValue, ...)
%
%   This function calculates various spatial metrics including DLE, SD, ROC AUC,
%   Precision-Recall AUC, F1-score, Precision, Recall, and APrime.
%   Helper functions for individual metrics are embedded as local functions.
%
%   Inputs:
%       estimated_source_activity (vector): Nsources x 1 vector of estimated source
%                                   activity magnitudes or energy (can be signed, abs is taken).
%       true_active_indices (vector): Indices of the true active source
%                                     locations in sourceSpaceInfo.Vertices.
%       active_seedVox (vector): Vector of indices of seedvox to active
%                                   source locations in sourceSpaceInfo.Vertices.
%       sourceSpaceInfo (struct): Structure containing source space information.
%           .Vertices (Nsources x 3 matrix): Coordinates of all source locations.
%           .VertConn (Nsources x Nsources sparse matrix, optional): Vertex connectivity.
%                                   Required for 'geodesic' distance if GeodesicDistanceMatrix
%                                   is not provided, and for ROC close field.
%           .GeodesicDistanceMatrix (Nsources x Nsources matrix, optional):
%                                   Precomputed geodesic distances. Takes precedence for 'geodesic' DistanceType.
%
%   Optional Name-Value Pair Inputs:
%       'MetricsToCompute' (cell array of strings): Which metrics to compute.
%           Options: {'DLE', 'SD', 'ROC_AUC', 'PR_AUC', 'F1', 'Precision', 'Recall', 'APrime'}.
%           Default: {'DLE', 'SD', 'ROC_AUC', 'PR_AUC', 'F1'}.
%       'DistanceType' (char/string): For DLE and SD. 'euclidean' (default), 'geodesic'.
%                                     If 'geodesic', uses GeodesicDistanceMatrix if available,
%                                     otherwise calculates from VertConn.
%       'ThresholdBinary' (scalar or char/string): Threshold to binarize
%                                     estimated_source_activity for F1, Precision, Recall, APrime.
%                                     - If numeric < 1: relative to max activity (e.g., 0.1 for 10%).
%                                     - If numeric >= 1: invalid, use Otsu method.
%                                     - 'otsu': Uses Otsu's method.
%                                     Default: 0.1 (relative).
%       'NumStepsROC' (integer): Number of threshold steps for ROC/PR AUC. Default: 100.
%       'RepsROC' (integer): Number of repetitions for ROC AUC bias correction. Default: 50.
%       'OrderROC' (integer): Max neighbor order for ROC AUC close field. Default: 6.
%
%   Outputs:
%       Results (struct): A structure containing the calculated metrics.
%           (Fields depend on 'MetricsToCompute')
%           .ParamsUsed (struct): All options used for the evaluation.
%
%   Example:
%       Nsrc = 200;
%       sourceSpaceInfo.Vertices = rand(Nsrc, 3) * 100;
%       adj = spdiags(ones(Nsrc,2), [-1 1], Nsrc, Nsrc);
%       sourceSpaceInfo.VertConn = adj; % For geodesic and ROC
%       true_idx = {sort(randperm(Nsrc, 20)')};
%       s_est = randn(Nsrc,1); s_est(true_idx(1:10)) = s_est(true_idx(1:10)) + 2;
%
%       evalResults = seal_evaluate_spatial_metrics(s_est, true_idx, sourceSpaceInfo, ...
%                       'MetricsToCompute', {'DLE', 'SD', 'ROC_AUC', 'F1'}, ...
%                       'DistanceType', 'euclidean', ...
%                       'ThresholdBinary', 0.2);
%       disp(evalResults);


    %% Input Parsing
    p = inputParser;
    addRequired(p, 'estimated_source_activity', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'true_active_indices', @(x) isnumeric(x) && isvector(x) && all(x > 0));
    addRequired(p, 'sourceSpaceInfo', @isstruct);
    addRequired(p, 'active_seedVox', @(x) isnumeric(x) && isvector(x) && all(x > 0))

    defaultMetrics = {'dle', 'sd', 'roc_auc', 'pr_auc', 'f1'};
    validMetrics = {'dle', 'sd', 'roc_auc', 'pr_auc', 'f1', 'precision', 'recall', 'aprime'};

    addParameter(p, 'MetricsToCompute', defaultMetrics, @(x) (iscellstr(x) || ischar(x) || isstring(x)) && all(ismember(lower(cellstr(x)), validMetrics)) );
    addParameter(p, 'DistanceType', 'euclidean', @(x) ismember(lower(x), {'euclidean', 'geodesic'}));
    addParameter(p, 'ThresholdBinary', 0.1, @(x) (isnumeric(x) && isscalar(x)) || (ischar(x) || isstring(x)));
    addParameter(p, 'NumStepsROC', 100, @(x) isnumeric(x) && isscalar(x) && x > 1);
    addParameter(p, 'RepsROC', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'OrderROC', 6, @(x) isnumeric(x) && isscalar(x) && x >= 0);

    parse(p, estimated_source_activity, true_active_indices,  active_seedVox, sourceSpaceInfo, varargin{:});

    s_est_input = p.Results.estimated_source_activity(:);
    true_active_idx_input = unique(p.Results.true_active_indices(:));
    active_seedVox = unique(p.Results.active_seedVox(:));
    
    srcInfo_main = p.Results.sourceSpaceInfo; % Renamed to avoid conflict in local functions
    
    metricsToCompute = lower(cellstr(p.Results.MetricsToCompute)); % Ensure cell array of lowercase strings

    Results.ParamsUsed = p.Results;
    Results.MetricsComputed = metricsToCompute;

    if ~isfield(srcInfo_main, 'Vertices')
        error('sourceSpaceInfo must contain a .Vertices field.');
    end
    Nsources = size(srcInfo_main.Vertices, 1);
    if length(s_est_input) ~= Nsources
        error('Length of estimated_source_activity must match Nsources.');
    end
    if any(true_active_idx_input > Nsources)
        error('true_active_indices contains out-of-bounds indices.');
    end

    s_est_abs = abs(s_est_input);
%     s_est_abs_for_mag_metrics = s_est_abs_for_mag_metrics/max(s_est_abs_for_mag_metrics);

    % Resolve internal distance type for helper functions
    distanceType_internal = p.Results.DistanceType;
    if strcmpi(p.Results.DistanceType, 'geodesic')
        if isfield(srcInfo_main, 'GeodesicDistanceMatrix') && ~isempty(srcInfo_main.GeodesicDistanceMatrix)
            distanceType_internal = 'geodesic_matrix';
        elseif isfield(srcInfo_main, 'VertConn') && ~isempty(srcInfo_main.VertConn)
            distanceType_internal = 'geodesic_graph';
        else
            warning('seal_evaluate_spatial_metrics:GeodesicUnavailable', ...
                    '''geodesic'' DistanceType selected, but neither GeodesicDistanceMatrix nor VertConn provided. Defaulting to Euclidean for DLE/SD.');
            distanceType_internal = 'euclidean'; % Fallback
            Results.ParamsUsed.DistanceTypeResolved = 'euclidean (fallback)';
        end
         Results.ParamsUsed.DistanceTypeResolved = distanceType_internal;
    else
         Results.ParamsUsed.DistanceTypeResolved = 'euclidean';
    end


    %% DLE & SD Calculation
    if any(strcmpi('dle', metricsToCompute)) || any(strcmpi('sd', metricsToCompute))
        try
            [dles, mean_dle, sd_val, peaks_idx] = compute_DLE_SD(true_active_idx_input, s_est_abs, srcInfo_main, distanceType_internal);
            if any(strcmpi('dle', metricsToCompute))
                Results.DLE_Values = dles;            
                Results.MeanDLE = mean_dle;
                Results.EstimatedPeakIndicesForDLE = peaks_idx;
            end
            if any(strcmpi('sd', metricsToCompute))
                Results.SD = sd_val;
            end
        catch ME_dle
            warning('seal_evaluate_spatial_metrics:DLEFailed', 'DLE/SD calculation failed: %s', ME_dle.message);
            Results.DLE_Values = NaN; Results.MeanDLE = NaN; Results.EstimatedPeakIndicesForDLE = NaN;
            Results.SD = NaN;
        end
    end

    %% ROC AUC Calculation
    if any(strcmpi('roc_auc', metricsToCompute))
        if ~isfield(srcInfo_main, 'VertConn') && strcmpi(distanceType_internal, 'geodesic_graph')
            warning('seal_evaluate_spatial_metrics:ROCSkipVertConn', 'ROC_AUC with bias correction (and geodesic_graph for close field) requires .VertConn. Skipping ROC_AUC if VertConn is missing.');
            Results.ROC_AUC = NaN; Results.ROC_Curve = [NaN NaN; NaN NaN];
        else
            try
                [roc_auc_val, roc_curve_val, param_roc] = compute_roc_auc(s_est_abs, true_active_idx_input, active_seedVox, srcInfo_main, Nsources, ...
                    p.Results.NumStepsROC, p.Results.RepsROC, p.Results.OrderROC);
                Results.ROC_AUC = roc_auc_val;
                Results.ROC_Curve = roc_curve_val;
                Results.ParamROC_details = param_roc;
            catch ME_roc
                warning('seal_evaluate_spatial_metrics:ROCFailed', 'ROC AUC calculation failed: %s', ME_roc.message);
                Results.ROC_AUC = NaN; Results.ROC_Curve = [NaN NaN; NaN NaN];
            end
        end
    end

    %% PR AUC Calculation
    if any(strcmpi('pr_auc', metricsToCompute))
        try
            [pr_auc_val, pr_curve_val, param_pr] = compute_pr_auc(s_est_abs, true_active_idx_input, Nsources, p.Results.NumStepsROC);
            Results.PR_AUC = pr_auc_val;
            Results.PR_Curve = pr_curve_val;
            Results.ParamPR_details = param_pr;
        catch ME_pr
            warning('seal_evaluate_spatial_metrics:PRFailed', 'PR AUC calculation failed: %s', ME_pr.message);
            Results.PR_AUC = NaN; Results.PR_Curve = [NaN NaN; NaN NaN];
        end
    end

    %% Binary Metrics (F1, Precision, Recall, APrime)
    needsBinaryMetrics = any(ismember({'f1', 'precision', 'recall', 'aprime'}, metricsToCompute));
    if needsBinaryMetrics
        s_est_binary = []; threshold_val_used = NaN; max_s_est = max(s_est_abs);
        if max_s_est < eps
            s_est_binary = false(Nsources, 1);
            warning('seal_evaluate_spatial_metrics:ZeroActivityForBinary', 'Max estimated activity is zero for binary metrics.');
        else
            if isnumeric(p.Results.ThresholdBinary)
                thresh_param = p.Results.ThresholdBinary;
                if thresh_param < 1 && thresh_param >=0
                    threshold_val_used = thresh_param * max_s_est;
                else
                    threshold_val_used = thresh_param;
                end
                s_est_binary = (s_est_abs >= threshold_val_used);
            elseif strcmpi(p.Results.ThresholdBinary, 'otsu')
                % Otsu logic (simplified, assuming custom otsu or graythresh)
                try
                    s_norm_for_otsu = s_est_abs / max_s_est;
                    if exist('graythresh','file'), otsu_thresh_norm = graythresh(s_norm_for_otsu);
                    elseif exist('otsu','file'), otsu_thresh_norm = otsu(s_norm_for_otsu); % User's function
                    else, error('Otsu method not found.'); 
                    end
                    threshold_val_used = otsu_thresh_norm * max_s_est;
                    s_est_binary = (s_est_abs >= threshold_val_used);
                catch ME_otsu_bin
                    warning('seal_evaluate_spatial_metrics:OtsuBinaryFailed', 'Otsu thresholding failed: %s. Using 0.1 relative.', ME_otsu_bin.message);
                    threshold_val_used = 0.1 * max_s_est; s_est_binary = (s_est_abs >= threshold_val_used);
                end
            else, error('Invalid ThresholdBinary value.');
            end
        end
        Results.ParamsUsed.ActualThresholdForBinaryMetrics = threshold_val_used;

        true_pos_map = false(Nsources, 1); true_pos_map(true_active_idx_input) = true;
        TP = sum(s_est_binary & true_pos_map); FP = sum(s_est_binary & ~true_pos_map);
        FN = sum(~s_est_binary & true_pos_map);

        if any(strcmpi('precision', metricsToCompute)) || any(strcmpi('f1', metricsToCompute))
            if (TP + FP) == 0, Results.Precision = 1; else, Results.Precision = TP / (TP + FP); end
        end
        if any(strcmpi('recall', metricsToCompute)) || any(strcmpi('f1', metricsToCompute))
            if (TP + FN) == 0, Results.Recall = 1; else, Results.Recall = TP / (TP + FN); end
        end
        if any(strcmpi('f1', metricsToCompute))
            if (Results.Precision + Results.Recall) == 0, Results.F1_Score = 0;
            else, Results.F1_Score = 2 * (Results.Precision * Results.Recall) / (Results.Precision + Results.Recall); end
        end
        if any(strcmpi('aprime', metricsToCompute))
            Hr = Results.Recall;
            actual_negatives_idx = setdiff((1:Nsources)', true_active_idx_input);
            N_actual_negatives = length(actual_negatives_idx);
            if N_actual_negatives == 0, Fr = 0; else, Fr = FP / N_actual_negatives; end
            Results.APrime = (Hr - Fr)/2 + 0.5;
            if isnan(Results.APrime), Results.APrime = 0.5; end
        end
    end
end

%% --- Functions for Metrics ---

function [dle_values, mean_dle, sd_value, estimated_peak_indices_for_true] = compute_DLE_SD(true_loc_idx_unique, s_est_abs, srcInfo_local, distanceType_resolved)
    % DLE and SD
    est_source_coords = srcInfo_local.Vertices;
    true_source_coords = srcInfo_local.Vertices(true_loc_idx_unique, :);
    Nsources_local = size(est_source_coords, 1);

    switch distanceType_resolved
        case 'euclidean'
            if exist('pdist2', 'file')
                dist_matrix_all_to_true = pdist2(est_source_coords, true_source_coords, 'euclidean');
            else
                dist_matrix_all_to_true = zeros(Nsources_local, length(true_loc_idx_unique));
                for i_s = 1:Nsources_local
                        dist_matrix_all_to_true(i_s,:) = sqrt(sum((repmat(est_source_coords(i_s,:),...
                            length(true_loc_idx_unique),1) - true_source_coords).^2,2))';
                end
            end
        case 'geodesic_matrix'
            dist_matrix_all_to_true = srcInfo_local.GeodesicDistanceMatrix(:, true_loc_idx_unique);
        case 'geodesic_graph'
            adj = srcInfo_local.VertConn; [row, col] = find(adj);
            weights = sqrt(sum((srcInfo_local.Vertices(row,:) - srcInfo_local.Vertices(col,:)).^2, 2));
            weighted_adj = sparse(row, col, weights, Nsources_local, Nsources_local);
            dist_matrix_all_to_true = zeros(Nsources_local, length(true_loc_idx_unique));
            for k_true_ptr = 1:length(true_loc_idx_unique)
                dist_matrix_all_to_true(:, k_true_ptr) = graphshortestpath(weighted_adj, true_loc_idx_unique(k_true_ptr))';
            end
        otherwise, error('Internal error: Unresolved distance type for DLE.');
    end
    
    % Distance to the nearest true source and the corresponding index for each estimated source
    [nearest_dist, nearest_true_indices] = min(dist_matrix_all_to_true, [], 2);
    nearest_true_indices = true_loc_idx_unique(nearest_true_indices);

    dle_values = NaN(length(true_loc_idx_unique), 1);
    estimated_peak_indices_for_true = NaN(length(true_loc_idx_unique), 1);

    for i_true = 1:length(true_loc_idx_unique)
        current_true_source_actual_idx = true_loc_idx_unique(i_true);
        % Find activated sources in estimated map nearest to the i_true source in the simulated map 
        partition_indices = find(nearest_true_indices == current_true_source_actual_idx);
        if isempty(partition_indices), continue; end
        activity_in_partition = s_est_abs(partition_indices);
        if ~any(activity_in_partition), continue; end
        % Find the estimated source with the largest energy
        [~, max_idx_within_partition] = max(activity_in_partition);
        peak_est_idx = partition_indices(max_idx_within_partition);        
        dle_values(i_true) = dist_matrix_all_to_true(peak_est_idx, i_true);       
        estimated_peak_indices_for_true(i_true) = peak_est_idx;
    end
    mean_dle = mean(dle_values, 'omitnan');
    % compute SD
    activity_weights_sq = abs(s_est_signed).^2;
    numerator = sum( (nearest_dist.^2) .* activity_weights_sq );
    denominator = sum(activity_weights_sq); 
    if denominator < Nsources_local * eps, sd_value = Inf; else, sd_value = sqrt(numerator / denominator); end
end

function [auc_roc, tpr_fpr_curve, ParamROC_local] = compute_roc_auc(s_est_abs, true_active_idx, seed_vox, srcInfo_local, num_steps, num_reps, close_order)
    ParamROC_local = struct();    
    max_s_est = max(s_est_abs);
    Nsources_local = size(s_est_abs, 1);
    
    if max_s_est < eps, auc_roc=0.5; tpr_fpr_curve=[0 0;1 1]; ParamROC_local.Mean_TPR=zeros(num_steps,1); ParamROC_local.Mean_FPR=zeros(num_steps,1); return; end
    thresholds = linspace(0, max_s_est, num_steps);
    
    all_tpr = zeros(num_steps, num_reps); all_fpr = zeros(num_steps, num_reps);
    
    adj = srcInfo_local.VertConn; close_field_global_idx = [];
    if close_order > 0 && ~isempty(adj) && exist('graphshortestpath', 'file')
        for k_true = 1:length(seed_vox)
            try
                dists_from_true_k = graphshortestpath(adj, seed_vox(k_true));
                close_field_global_idx = union(close_field_global_idx, find(dists_from_true_k <= close_order & dists_from_true_k >=0));
            catch 
            
            end % Ignore errors in graph path for robustness here
        end
    end
    if isempty(close_field_global_idx), close_field_global_idx = true_active_idx; end
    close_field_global_idx = unique(close_field_global_idx(:));

    for r = 1:num_reps
        close_inactive_cand = setdiff(close_field_global_idx, true_active_idx);
        far_inactive_cand = setdiff((1:Nsources_local)', close_field_global_idx);
        num_true = length(true_active_idx);
        n_close = floor(num_true/2); n_far = num_true - n_close;
        
        samp_close_inact = []; if ~isempty(close_inactive_cand) && n_close>0, idx_c = randperm(length(close_inactive_cand)); samp_close_inact = close_inactive_cand(idx_c(1:min(n_close,end))); end
        samp_far_inact = []; if ~isempty(far_inactive_cand) && n_far>0, idx_f = randperm(length(far_inactive_cand)); samp_far_inact = far_inactive_cand(idx_f(1:min(n_far,end))); end
        
        % Fill up if one set was too small
        if length(samp_close_inact) < n_close && ~isempty(far_inactive_cand)
            needed = n_close - length(samp_close_inact);
            rem_far = setdiff(far_inactive_cand, samp_far_inact);
            if ~isempty(rem_far), idx_rf = randperm(length(rem_far)); samp_close_inact = [samp_close_inact; rem_far(idx_rf(1:min(needed,end)))]; end
        end
        if length(samp_far_inact) < n_far && ~isempty(close_inactive_cand)
             needed = n_far - length(samp_far_inact);
             rem_close = setdiff(close_inactive_cand, samp_close_inact);
             if ~isempty(rem_close), idx_rc = randperm(length(rem_close)); samp_far_inact = [samp_far_inact; rem_close(idx_rc(1:min(needed,end)))]; end
        end
        inactive_sample_idx = unique([samp_close_inact; samp_far_inact]);
        if isempty(inactive_sample_idx) && num_true < Nsources_local, inactive_sample_idx = setdiff((1:Nsources_local)', true_active_idx); if length(inactive_sample_idx)>num_true, inactive_sample_idx = inactive_sample_idx(randperm(length(inactive_sample_idx),num_true)); end; end


        for t_idx = 1:num_steps
            est_pos = find(s_est_abs >= thresholds(t_idx));
            TP = length(intersect(est_pos, true_active_idx)); FN = length(true_active_idx) - TP;
            FP = length(intersect(est_pos, inactive_sample_idx)); TN = length(inactive_sample_idx) - FP;
            if (TP+FN)==0, all_tpr(t_idx,r)=0; else, all_tpr(t_idx,r) = TP/(TP+FN); end
            if (FP+TN)==0, all_fpr(t_idx,r)=0; else, all_fpr(t_idx,r) = FP/(FP+TN); end
        end
    end
    mean_tpr = mean(all_tpr,2); mean_fpr = mean(all_fpr,2);
    [unique_fpr_tpr, ~] = unique([mean_fpr, mean_tpr], 'rows', 'stable'); % Keep stable for sort
    [sorted_fpr, sort_idx] = sort(unique_fpr_tpr(:,1));
    sorted_tpr = unique_fpr_tpr(sort_idx,2);
    final_fpr = [0; sorted_fpr; 1]; final_tpr = [0; sorted_tpr; 1];
    [unique_final_pts, ia] = unique([final_fpr, final_tpr], 'rows', 'first'); % Use 'first' to keep original ordering for trapz
    tpr_fpr_curve = unique_final_pts;
    auc_roc = trapz(tpr_fpr_curve(:,1), tpr_fpr_curve(:,2));
    ParamROC_local.Mean_TPR = mean_tpr; ParamROC_local.Mean_FPR = mean_fpr;
end

function [auc_pr, precision_recall_curve, ParamPR_local] = compute_pr_auc(s_est_abs, true_active_idx,  num_steps)
    % Simplified PR AUC - inputs assumed valid
    ParamPR_local = struct();
    max_s_est = max(s_est_abs);
    if max_s_est < eps, auc_pr=0; precision_recall_curve=[0 1;1 0]; ParamPR_local.PrecisionValues=zeros(num_steps,1); ParamPR_local.RecallValues=zeros(num_steps,1); return; end
    thresholds = linspace(0, max_s_est, num_steps);
    precision_vals = zeros(num_steps,1); recall_vals = zeros(num_steps,1);
    num_actual_pos = length(true_active_idx);
    if num_actual_pos == 0, auc_pr=NaN; precision_recall_curve=[0 1;1 0]; ParamPR_local.PrecisionValues=NaN; ParamPR_local.RecallValues=NaN; return; end

    for t_idx = 1:num_steps
        est_pos = find(s_est_abs >= thresholds(t_idx));
        TP = length(intersect(est_pos, true_active_idx));
        FP = length(setdiff(est_pos, true_active_idx));
        if (TP+FP)==0, precision_vals(t_idx)=1; else, precision_vals(t_idx)=TP/(TP+FP); end
        recall_vals(t_idx) = TP/num_actual_pos;
    end
    ParamPR_local.PrecisionValues = precision_vals; ParamPR_local.RecallValues = recall_vals;
    
    % Sort by recall, then by precision descending for ties
    [sorted_pr_pts, ~] = sortrows([recall_vals, precision_vals], [1, -2]);
    
    % Add standard anchor points (0,1) and (1,0)
    final_recall = [0; sorted_pr_pts(:,1); 1];
    final_precision = [1; sorted_pr_pts(:,2); 0]; % P=1 at R=0, P=0 at R=1 (common convention)
    
    [unique_pr_pts, ~] = unique([final_recall, final_precision], 'rows', 'stable');
    
    % Ensure precision is non-increasing for AUC calculation
    unique_precision_mono = unique_pr_pts(:,2);
    for i = length(unique_precision_mono)-1:-1:1
        unique_precision_mono(i) = max(unique_precision_mono(i), unique_precision_mono(i+1));
    end
    precision_recall_curve = [unique_pr_pts(:,1), unique_precision_mono];
    auc_pr = trapz(precision_recall_curve(:,1), precision_recall_curve(:,2));
end
