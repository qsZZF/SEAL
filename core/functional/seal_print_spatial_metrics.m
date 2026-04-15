function print_spatial_metrics(spatial_metric)
% PRINT_SPATIAL_METRICS 打印指定的空间定位评估指标
%   print_spatial_metrics(spatial_metric) 接收结构体，仅打印以下指标：
%   PR_AUC, SD, DLE, ROC_AUC, APrime, F1, Precision, Recall
%
%   输入:
%       spatial_metric - 结构体，包含所需字段

    if nargin < 1 || ~isstruct(spatial_metric)
        error('输入必须是一个结构体 (spatial_metric)。');
    end

    fprintf('\n========== 空间定位评估指标 ==========\n');

    % 定义需要打印的字段（顺序可调）
    fields_to_print = { ...
        'PR_AUC',    'PR AUC',     ''; ...
        'SD',        'SD',         'm'; ...
        'MeanDLE',       'DLE',        'm'; ...
        'ROC_AUC',   'ROC AUC',    ''; ...
        'APrime',    "A'",         ''; ...
        'F1',        'F1 Score',   ''; ...
        'Precision', 'Precision',  ''; ...
        'Recall',    'Recall',     ''; ...
        };

    for i = 1:size(fields_to_print, 1)
        field_name = fields_to_print{i, 1};
        display_name = fields_to_print{i, 2};
        unit = fields_to_print{i, 3};

        if isfield(spatial_metric, field_name)
            value = spatial_metric.(field_name);
            if isnumeric(value) && isscalar(value)
                if ~isempty(unit)
                    fprintf('%-12s : %10.3f %s\n', display_name, value, unit);
                else
                    fprintf('%-12s : %10.3f\n', display_name, value);
                end
            elseif isnumeric(value) && ~isscalar(value)
                % 如果是向量或矩阵，打印统计信息
                fprintf('%-12s : 向量长度 %d, 范围 [%.4f, %.4f]\n', ...
                    display_name, numel(value), min(value(:)), max(value(:)));
            else
                fprintf('%-12s : %s (非数值)\n', display_name, class(value));
            end

        end
    end

    fprintf('=======================================\n\n');
end