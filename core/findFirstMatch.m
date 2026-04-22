% 辅助函数：在候选字段中找第一个匹配的别名
function matched = findFirstMatch(aliases, availableFields)
    matched = '';
    for i = 1:length(aliases)
        if ismember(aliases{i}, availableFields)
            matched = aliases{i};
            return;
        end
    end
end