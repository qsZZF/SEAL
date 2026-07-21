function Results = seal_evaluate_spatial_metrics(estimated_source_activity, true_active_indices, active_seedVox, sourceSpaceInfo, varargin)
%SEAL_EVALUATE_SPATIAL_METRICS Computes a suite of spatial performance metrics for ESI.
%   Results = seal_evaluate_spatial_metrics(estimated_source_activity, ...
%                                   true_active_indices, sourceSpaceInfo, 'ParameterName', ParameterValue, ...)
%
%   This function calculates various spatial metrics including DLE, SD, EMD,
%   ROC AUC, Precision-Recall AUC, F1-score, Precision, Recall, and APrime.
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
%           Options: {'DLE', 'SD', 'EMD', 'ROC_AUC', 'PR_AUC', 'F1',
%                     'Precision', 'Recall', 'APrime', 'MSE'}.
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
%       'NumStepsROC' (integer): Number of threshold steps for the auxiliary
%                                     spatially balanced ROC. Default: 100.
%       'NumStepsPR' (integer): Retained for backward compatibility. PR AUC
%                                     uses every distinct source score.
%       'RepsROC' (integer): Number of repetitions for the auxiliary
%                                     spatially balanced ROC. Default: 50.
%       'OrderROC' (integer): Max neighbor order for its close field. Default: 6.
%       'EMDLambda' (scalar): Entropic regularization for Sinkhorn. Default: 1.
%       'EMDMaxIter' (integer): Max Sinkhorn iterations. Default: 500.
%       'EMDTolerance' (scalar): Sinkhorn convergence tolerance. Default: 1e-5.
%       'EMDMaxSources' (integer/Inf): Max strongest estimated sources used for
%                                     EMD approximation. True active sources are
%                                     always included. Default: 1000.
%       'EMDCostMatrix' (matrix): Optional full or selected-source cost matrix.
%       'EMDTrueDistribution' (vector): Optional Nsources x 1 true distribution.
%                                     Defaults to uniform mass over true_active_indices.
%       'EstimatedSourceTimeSeries' (matrix): Nsources x Ntime estimated source
%                                     time series. Required for 'mse'.
%       'TrueSourceTimeSeries' (matrix): Nsources x Ntime ground-truth source
%                                     time series. Required for 'mse'.
%       'TemporalRoiIndices' (cell/vector): ROI/source indices used for temporal
%                                     correlation summaries. Defaults to active_seedVox.
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
    addRequired(p, 'active_seedVox', @(x) isnumeric(x) && isvector(x) && all(x > 0))
    addRequired(p, 'sourceSpaceInfo', @(x) isstruct(x));
    
    defaultMetrics = {'dle', 'sd', 'roc_auc', 'pr_auc', 'f1'};
    validMetrics = {'dle', 'sd', 'emd', 'roc_auc', 'pr_auc', 'f1', 'precision', 'recall', 'aprime','mse'};

    addParameter(p, 'MetricsToCompute', defaultMetrics, @(x) (iscellstr(x) || ischar(x) || isstring(x)) && all(ismember(lower(cellstr(x)), validMetrics)) );
    addParameter(p, 'DistanceType', 'euclidean', @(x) ismember(lower(x), {'euclidean', 'geodesic'}));
    addParameter(p, 'ThresholdBinary', 0.1, @(x) (isnumeric(x) && isscalar(x)) || (ischar(x) || isstring(x)));
    addParameter(p, 'NumStepsROC', 100, @(x) isnumeric(x) && isscalar(x) && x > 1);
    addParameter(p, 'NumStepsPR', 100, @(x) isnumeric(x) && isscalar(x) && x > 1);
    addParameter(p, 'RepsROC', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'OrderROC', 6, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'EMDLambda', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'EMDMaxIter', 500, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'EMDTolerance', 1e-5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'EMDMaxSources', 1000, @(x) isnumeric(x) && isscalar(x) && (isinf(x) || x >= 1));
    addParameter(p, 'EMDCostMatrix', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'EMDTrueDistribution', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
    addParameter(p, 'EstimatedSourceTimeSeries', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'TrueSourceTimeSeries', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'TemporalRoiIndices', [], @(x) isempty(x) || isnumeric(x) || iscell(x));

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
    s_est_time = p.Results.EstimatedSourceTimeSeries;
    s_real_time = p.Results.TrueSourceTimeSeries;
    if ~isempty(s_est_time) && size(s_est_time, 1) ~= Nsources
        error('Rows of EstimatedSourceTimeSeries must match Nsources.');
    end
    if ~isempty(s_real_time) && size(s_real_time, 1) ~= Nsources
        error('Rows of TrueSourceTimeSeries must match Nsources.');
    end
    if ~isempty(s_est_time) && ~isempty(s_real_time) && size(s_est_time, 2) ~= size(s_real_time, 2)
        error('EstimatedSourceTimeSeries and TrueSourceTimeSeries must have the same number of time points.');
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

    %% EMD
    if any(strcmpi('emd', metricsToCompute))
        try
            [emd_val, emd_details] = compute_emd_metric( ...
                s_est_abs, true_active_idx_input, srcInfo_main, distanceType_internal, ...
                p.Results.EMDLambda, p.Results.EMDMaxIter, p.Results.EMDTolerance, ...
                p.Results.EMDMaxSources, p.Results.EMDCostMatrix, p.Results.EMDTrueDistribution);
            Results.EMD = emd_val;
            Results.EMD_Details = emd_details;
        catch ME_emd
            warning('seal_evaluate_spatial_metrics:EMDFailed', ...
                'EMD calculation failed: %s', ME_emd.message);
            Results.EMD = NaN;
            Results.EMD_Details = struct();
        end
    end

    %% MSE
    if any(strcmpi('mse', metricsToCompute))
        if isempty(s_est_time) || isempty(s_real_time)
            warning('seal_evaluate_spatial_metrics:MSEMissingTimeSeries', ...
                'MSE requires EstimatedSourceTimeSeries and TrueSourceTimeSeries. Skipping MSE.');
            Results.MSE_Values = NaN;
            Results.CORRall = NaN;
            Results.CORRroi = NaN;
            Results.CORRavg = NaN;
        else
            try
                temporal_roi_indices = p.Results.TemporalRoiIndices;
                if isempty(temporal_roi_indices)
                    temporal_roi_indices = active_seedVox;
                end
                [MSE,CORRall,CORRroi,CORRavg] = ESItemporalMetric(s_est_time,s_real_time,temporal_roi_indices);
                Results.MSE_Values = MSE;
                Results.CORRall = CORRall;
                Results.CORRroi = CORRroi;
                Results.CORRavg = CORRavg;
            catch MSE_dle
                warning('seal_evaluate_spatial_metrics:MSEFailed', 'MSE calculation failed: %s', MSE_dle.message);
                Results.MSE_Values = NaN;
                Results.CORRall = NaN;
                Results.CORRroi = NaN;
                Results.CORRavg = NaN;
            end
        end
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
        try
            [roc_auc_val, roc_curve_val, param_roc] = compute_standard_roc_auc( ...
                s_est_abs, true_active_idx_input);
            Results.ROC_AUC = roc_auc_val;
            Results.ROC_Curve = roc_curve_val;
            Results.ParamROC_details = param_roc;

            % Preserve the previous close/far balanced calculation as a
            % separately named diagnostic. It is not a standard ROC-AUC.
            try
                [balanced_auc, balanced_curve, balanced_param] = compute_spatially_balanced_roc_auc( ...
                    s_est_abs, true_active_idx_input, active_seedVox, srcInfo_main, ...
                    p.Results.NumStepsROC, p.Results.RepsROC, p.Results.OrderROC);
                Results.ROC_AUC_SpatialBalanced = balanced_auc;
                Results.ROC_Curve_SpatialBalanced = balanced_curve;
                Results.ParamROC_SpatialBalanced = balanced_param;
            catch ME_balanced_roc
                warning('seal_evaluate_spatial_metrics:SpatialBalancedROCFailed', ...
                    'Spatially balanced ROC calculation failed: %s', ME_balanced_roc.message);
                Results.ROC_AUC_SpatialBalanced = NaN;
                Results.ROC_Curve_SpatialBalanced = [NaN NaN; NaN NaN];
            end
        catch ME_roc
            warning('seal_evaluate_spatial_metrics:ROCFailed', 'ROC AUC calculation failed: %s', ME_roc.message);
            Results.ROC_AUC = NaN; Results.ROC_Curve = [NaN NaN; NaN NaN];
        end
    end

    %% PR AUC Calculation
    if any(strcmpi('pr_auc', metricsToCompute))
        try
            [pr_auc_val, pr_curve_val, param_pr] = compute_pr_auc( ...
                s_est_abs, true_active_idx_input, p.Results.NumStepsPR);
            Results.PR_AUC = pr_auc_val;
            Results.AveragePrecision = pr_auc_val;
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
            Results.F1 = Results.F1_Score;
        end
        if any(strcmpi('aprime', metricsToCompute))
            if (TP + FN) == 0, Hr = 1; else, Hr = TP / (TP + FN); end
            actual_negatives_idx = setdiff((1:Nsources)', true_active_idx_input);
            N_actual_negatives = length(actual_negatives_idx);
            if N_actual_negatives == 0, Fr = 0; else, Fr = FP / N_actual_negatives; end
            Results.APrime = compute_aprime(Hr, Fr);
            if isnan(Results.APrime), Results.APrime = 0.5; end
        end
    end
end

%% --- Functions for Metrics ---

function [emd_val, emd_details] = compute_emd_metric(s_est_abs, true_active_idx, srcInfo_local, distanceType_resolved, lambda, max_iter, tol, max_sources, cost_matrix_input, true_distribution_input)
    Nsources_local = size(srcInfo_local.Vertices, 1);

    if isempty(true_distribution_input)
        true_distribution = zeros(Nsources_local, 1);
        true_distribution(true_active_idx) = 1;
    else
        true_distribution = true_distribution_input(:);
        if length(true_distribution) ~= Nsources_local
            error('EMDTrueDistribution length must match Nsources.');
        end
        if any(true_distribution < 0)
            error('EMDTrueDistribution must be nonnegative.');
        end
    end

    [dist_est, dist_true, emd_source_idx, is_approx] = prepare_emd_distributions( ...
        s_est_abs, true_distribution, max_sources);
    cost_matrix = resolve_emd_cost_matrix( ...
        cost_matrix_input, srcInfo_local, distanceType_resolved, emd_source_idx, Nsources_local);

    emd_val = calculate_emd_sinkhorn(dist_est, dist_true, cost_matrix, lambda, max_iter, tol);
    emd_details = struct();
    emd_details.SourceIndices = emd_source_idx;
    emd_details.NumSourcesUsed = numel(emd_source_idx);
    emd_details.IsApproximate = is_approx;
    emd_details.RequestedMaxSources = max_sources;
    emd_details.DistanceTypeResolved = distanceType_resolved;
    emd_details.Lambda = lambda;
    emd_details.MaxIter = max_iter;
    emd_details.Tolerance = tol;
end

function [dist_est, dist_true, source_idx, is_approx] = prepare_emd_distributions(s_est_abs, true_distribution, max_sources)
    s_est_abs = s_est_abs(:);
    true_distribution = true_distribution(:);
    Nsources_local = length(s_est_abs);

    if length(true_distribution) ~= Nsources_local
        error('Estimated and true EMD distributions must have the same length.');
    end
    if any(s_est_abs < 0) || any(true_distribution < 0)
        error('EMD distributions must be nonnegative.');
    end
    if sum(s_est_abs) <= eps
        error('Estimated source activity has zero total mass.');
    end
    if sum(true_distribution) <= eps
        error('True EMD distribution has zero total mass.');
    end

    if isinf(max_sources) || Nsources_local <= max_sources
        source_idx = (1:Nsources_local)';
    else
        est_scores = s_est_abs;
        est_scores(~isfinite(est_scores)) = 0;
        [~, sorted_idx] = sort(est_scores, 'descend');
        n_est_keep = min(floor(max_sources), Nsources_local);
        est_support = sorted_idx(1:n_est_keep);
        true_support = find(true_distribution > 0);
        source_idx = unique([true_support(:); est_support(:)]);
    end

    dist_est = s_est_abs(source_idx);
    dist_true = true_distribution(source_idx);
    dist_est = dist_est / sum(dist_est);
    dist_true = dist_true / sum(dist_true);
    is_approx = numel(source_idx) < Nsources_local;
end

function cost_matrix = resolve_emd_cost_matrix(cost_matrix_input, srcInfo_local, distanceType_resolved, source_idx, Nsources_local)
    n_used = numel(source_idx);

    if ~isempty(cost_matrix_input)
        if isequal(size(cost_matrix_input), [Nsources_local, Nsources_local])
            cost_matrix = cost_matrix_input(source_idx, source_idx);
        elseif isequal(size(cost_matrix_input), [n_used, n_used])
            cost_matrix = cost_matrix_input;
        else
            error('EMDCostMatrix must be either Nsources x Nsources or NumSourcesUsed x NumSourcesUsed.');
        end
    else
        switch distanceType_resolved
            case 'geodesic_matrix'
                cost_matrix = srcInfo_local.GeodesicDistanceMatrix(source_idx, source_idx);
            case 'geodesic_graph'
                cost_matrix = compute_geodesic_emd_cost_matrix(srcInfo_local, source_idx);
            otherwise
                cost_matrix = compute_pairwise_euclidean(srcInfo_local.Vertices(source_idx, :));
        end
    end

    cost_matrix = double(cost_matrix);
    bad_cost = ~isfinite(cost_matrix);
    if any(bad_cost(:))
        finite_cost = cost_matrix(isfinite(cost_matrix));
        finite_cost = finite_cost(finite_cost > 0);
        if isempty(finite_cost)
            replacement_cost = 0;
        else
            replacement_cost = 2 * max(finite_cost);
        end
        cost_matrix(bad_cost) = replacement_cost;
    end
    cost_matrix(1:n_used+1:end) = 0;
end

function cost_matrix = compute_pairwise_euclidean(coords)
    if exist('pdist2', 'file') == 2
        cost_matrix = pdist2(coords, coords, 'euclidean');
    else
        sq_norm = sum(coords.^2, 2);
        sq_dist = bsxfun(@plus, sq_norm, sq_norm') - 2 * (coords * coords');
        cost_matrix = sqrt(max(sq_dist, 0));
    end
end

function cost_matrix = compute_geodesic_emd_cost_matrix(srcInfo_local, source_idx)
    if ~isfield(srcInfo_local, 'VertConn') || isempty(srcInfo_local.VertConn) || exist('graphshortestpath', 'file') ~= 2
        warning('seal_evaluate_spatial_metrics:EMDGeodesicFallback', ...
            'Geodesic EMD cost matrix is unavailable. Falling back to Euclidean distances.');
        cost_matrix = compute_pairwise_euclidean(srcInfo_local.Vertices(source_idx, :));
        return;
    end

    Nsources_local = size(srcInfo_local.Vertices, 1);
    [row, col] = find(srcInfo_local.VertConn);
    weights = sqrt(sum((srcInfo_local.Vertices(row,:) - srcInfo_local.Vertices(col,:)).^2, 2));
    weighted_adj = sparse(row, col, weights, Nsources_local, Nsources_local);

    n_used = numel(source_idx);
    cost_matrix = zeros(n_used, n_used);
    for k = 1:n_used
        dists = graphshortestpath(weighted_adj, source_idx(k));
        cost_matrix(k, :) = dists(source_idx);
    end
end

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
    activity_weights_sq = abs(s_est_abs);
    numerator = sum( (nearest_dist.^2) .* activity_weights_sq );
    denominator = sum(activity_weights_sq); 
    if denominator < Nsources_local * eps, sd_value = Inf; else, sd_value = sqrt(numerator / denominator); end
end

function [auc_roc, tpr_fpr_curve, ParamROC_local] = compute_standard_roc_auc(s_est_abs, true_active_idx)
    % Exact empirical ROC using every cortical vertex and all unique scores.
    num_sources = numel(s_est_abs);
    labels = false(num_sources, 1);
    labels(true_active_idx) = true;
    num_positive = sum(labels);
    num_negative = num_sources - num_positive;

    ParamROC_local = struct('Method', 'exact_all_vertices', ...
        'NumPositive', num_positive, 'NumNegative', num_negative);
    if num_positive == 0 || num_negative == 0
        auc_roc = NaN;
        tpr_fpr_curve = [NaN NaN; NaN NaN];
        return;
    end

    [scores_descending, order] = sort(s_est_abs, 'descend');
    labels_descending = labels(order);
    group_end = [find(diff(scores_descending) ~= 0); num_sources];
    group_start = [1; group_end(1:end-1) + 1];
    positives_per_group = zeros(numel(group_end), 1);
    negatives_per_group = zeros(numel(group_end), 1);
    for i_group = 1:numel(group_end)
        group_labels = labels_descending(group_start(i_group):group_end(i_group));
        positives_per_group(i_group) = sum(group_labels);
        negatives_per_group(i_group) = numel(group_labels) - positives_per_group(i_group);
    end

    true_positive_rate = [0; cumsum(positives_per_group) / num_positive];
    false_positive_rate = [0; cumsum(negatives_per_group) / num_negative];
    tpr_fpr_curve = [false_positive_rate, true_positive_rate];
    auc_roc = trapz(false_positive_rate, true_positive_rate);
    ParamROC_local.Thresholds = [Inf; scores_descending(group_start)];
end

function [auc_roc, tpr_fpr_curve, ParamROC_local] = compute_spatially_balanced_roc_auc(s_est_abs, true_active_idx, seed_vox, srcInfo_local, num_steps, num_reps, close_order)
    ParamROC_local = struct();    
    max_s_est = max(s_est_abs);
    Nsources_local = size(s_est_abs, 1);
    
    if max_s_est < eps, auc_roc=0.5; tpr_fpr_curve=[0 0;1 1]; ParamROC_local.Mean_TPR=zeros(num_steps,1); ParamROC_local.Mean_FPR=zeros(num_steps,1); return; end
    thresholds = linspace(0, max_s_est, num_steps);
    
    all_tpr = zeros(num_steps, num_reps); all_fpr = zeros(num_steps, num_reps);
    
    if isfield(srcInfo_local, 'VertConn') && ~isempty(srcInfo_local.VertConn)
        adj = srcInfo_local.VertConn;
    else
        adj = [];
    end
    close_field_global_idx = [];
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

    use_exact_all_negatives = (close_order <= 0) || isempty(adj);

    for r = 1:num_reps
        close_inactive_cand = setdiff(close_field_global_idx, true_active_idx);
        far_inactive_cand = setdiff((1:Nsources_local)', close_field_global_idx);
        num_true = length(true_active_idx);
        n_close = floor(num_true/2); n_far = num_true - n_close;

        if use_exact_all_negatives
            inactive_sample_idx = setdiff((1:Nsources_local)', true_active_idx);
        else
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
            if isempty(inactive_sample_idx) && num_true < Nsources_local
                inactive_sample_idx = setdiff((1:Nsources_local)', true_active_idx);
                if length(inactive_sample_idx)>num_true, inactive_sample_idx = inactive_sample_idx(randperm(length(inactive_sample_idx),num_true)); end
            end
        end


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
    [unique_final_pts, ~] = unique([final_fpr, final_tpr], 'rows', 'first'); % Use 'first' to keep original ordering for trapz
    tpr_fpr_curve = unique_final_pts;
    auc_roc = trapz(tpr_fpr_curve(:,1), tpr_fpr_curve(:,2));
    ParamROC_local.Mean_TPR = mean_tpr; ParamROC_local.Mean_FPR = mean_fpr;
end

function [auc_pr, precision_recall_curve, ParamPR_local] = compute_pr_auc(s_est_abs, true_active_idx, ~)
    % Average precision over every distinct source score.  A fixed grid of
    % relative thresholds fails for a focal map: its recall can be constant
    % over the entire grid, spuriously returning a PR AUC of zero.
    scores = s_est_abs(:);
    num_sources = numel(scores);
    labels = false(num_sources, 1);
    labels(true_active_idx) = true;
    num_positive = sum(labels);
    ParamPR_local = struct('Method', 'ranked_score_average_precision', ...
        'NumPositive', num_positive, ...
        'NumNegative', num_sources - num_positive, ...
        'NumDistinctScores', 0, ...
        'NormalizedByRecallRange', false);
    if num_positive == 0
        auc_pr = NaN;
        precision_recall_curve = [NaN NaN; NaN NaN];
        return;
    end

    [scores_sorted, sort_index] = sort(scores, 'descend');
    labels_sorted = labels(sort_index);
    group_end = [find(diff(scores_sorted) ~= 0); num_sources];
    true_positive = cumsum(labels_sorted);
    predicted_positive = (1:num_sources).';

    recall = true_positive(group_end) / num_positive;
    precision = true_positive(group_end) ./ predicted_positive(group_end);
    thresholds = scores_sorted(group_end);

    % The initial point is a conventional empty selection.  Average
    % precision is a step integral, not a trapezoidal interpolation.
    recall = [0; recall];
    precision = [1; precision];
    thresholds = [Inf; thresholds];
    auc_pr = sum(diff(recall) .* precision(2:end));

    precision_recall_curve = [recall, precision];
    ParamPR_local.NumDistinctScores = numel(group_end);
    ParamPR_local.Thresholds = thresholds;
    ParamPR_local.PrecisionValues = precision;
    ParamPR_local.RecallValues = recall;
end

function aprime_val = compute_aprime(hit_rate, false_alarm_rate)
    if hit_rate == false_alarm_rate
        aprime_val = 0.5;
    elseif hit_rate > false_alarm_rate
        aprime_val = 0.5 + ((hit_rate - false_alarm_rate) * (1 + hit_rate - false_alarm_rate)) / ...
            (4 * hit_rate * (1 - false_alarm_rate));
    else
        aprime_val = 0.5 - ((false_alarm_rate - hit_rate) * (1 + false_alarm_rate - hit_rate)) / ...
            (4 * false_alarm_rate * (1 - hit_rate));
    end
end

function [MSE,CORRall,CORRroi,CORRavg] = ESItemporalMetric(s_est,s_real,activevoxseed)

    real_norm = norm(s_real, 'fro');
    est_norm = norm(s_est, 'fro');
    s_real_norm = s_real/real_norm;
    s_est_norm = s_est/est_norm;
    if real_norm < eps || est_norm < eps
        MSE = NaN;
    else
        MSE = norm(s_real_norm - s_est_norm, 'fro');
    end


    activevox = [];
    if isempty(activevoxseed)
        CORRavg = NaN;
        CORRroi = NaN;
        activevoxseed = {};
    elseif isnumeric(activevoxseed)
        activevoxseed = num2cell(activevoxseed(:)');
    end
    numK = numel(activevoxseed);
    Sgt = zeros(numK,size(s_est_norm,2));
    Src = Sgt;
    for k = 1:numK
        tseed = activevoxseed{k};
        activevox = union(activevox,tseed);
        Sgt(k,:) = mean(s_real_norm(tseed,:));
        Src(k,:) = mean(s_est_norm(tseed,:));
    end
    Sgt = Sgt-mean(Sgt,2); Src = Src-mean(Src,2);
    CORRavg = mean(abs(sum(Sgt.*Src,2))./sqrt(sum(Sgt.^2,2).*sum(Src.^2,2)+1e-16));

    S_centered = s_real_norm - mean(s_real_norm,2);
    S_est_centerd = s_est_norm - mean(s_est_norm,2);
    CORRvox = abs(sum(S_centered.*S_est_centerd,2))./sqrt(sum(S_centered.^2,2).*sum(S_est_centerd.^2,2)+1e-16);
    CORRall = mean(CORRvox);
    if ~isempty(activevox)
        CORRroi = mean(CORRvox(activevox));
    end


end
