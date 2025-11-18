classdef Validators
    properties(Constant)
        % 验证函数句柄
        isText = @(x) (ischar(x) && isrow(x)) || (isstring(x) && isscalar(x));
        isNonEmptyText = @(x) SEALConfig.IsText(x) && ~isempty(strtrim(x));
        isFilePath = @(x) SEALConfig.IsText(x) && (isfolder(x) || ~isempty(fileparts(x)));
    end
end

