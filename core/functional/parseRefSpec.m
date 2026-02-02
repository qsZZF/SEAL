function refSpec = parseRefSpec(refText, nCh, chanLabels)
% refText: string（用户输入）
% nCh:     通道总数
% chanLabels: 1×nCh cellstr，允许为空（则不能按名字解析）

    refText = strtrim(string(refText));
    refLower = lower(refText);

    % 1) 全通道平均
    if refLower == "average"
        refSpec.type = 'average-all';
        return;
    end

    % 2) 只有一个名字或一个索引（单电极参考）
    %    支持 A1 / 'A1' / 5 / '5'
    oneTok = stripQuotes(refText);
    asNum  = str2double(oneTok);
    if ~isnan(asNum) && isscalar(asNum)
        % 单索引
        checkIndexRange(asNum, nCh);
        refSpec.type = 'channels';
        refSpec.idx  = asNum;
        return;
    end
    if ~isempty(chanLabels)
        % 单名字
        idx = name2index(oneTok, chanLabels);
        if ~isempty(idx)
            refSpec.type = 'channels';
            refSpec.idx  = idx;
            return;
        end
    end

    % 3) 形如 [1,3,5] 或 {'A1','A2'} （多通道平均参考）
    if startsWith(refText, "[") && endsWith(refText, "]")
        vec = str2num(char(refText)); %#ok<ST2NM>  % 只解析数字向量
        validateattributes(vec, {'numeric'}, {'vector','integer','positive'});
        vec = unique(vec);
        vec(vec<1 | vec>nCh) = [];
        assert(~isempty(vec), '索引越界或为空。');
        refSpec.type = 'channels';
        refSpec.idx  = vec;
        return;
    end

    if startsWith(refText, "{") && endsWith(refText, "}")
        % 尝试两种形式：
        %   a) {'A1','A2','T7'}  -> channels
        %   b) {[1,3,5],[2,4,6]} or {{'A1','A2'},{'T7','T8'}} -> groups
        % 先判定是否有多个子块（用 [] 或 {} 分组）
        % i) 提取 [..] 数字组
        numGroups = regexp(char(refText), '\[([0-9,\s]+)\]', 'match');
        if ~isempty(numGroups)
            groups = cell(1,numel(numGroups));
            for k=1:numel(numGroups)
                vec = str2num(numGroups{k}); %#ok<ST2NM>
                vec = unique(vec);
                vec(vec<1 | vec>nCh) = [];
                groups{k} = vec;
            end
            % 若仅一个组，视为 channels；若多个组，视为 groups
            if numel(groups)==1
                refSpec.type = 'channels';
                refSpec.idx  = groups{1};
            else
                refSpec.type  = 'groups';
                refSpec.groups = groups;
            end
            return;
        end

        % ii) 提取 {'A1','A2'} 或 {{'A1','A2'},{'T7','T8'}}
        if ~isempty(chanLabels)
            % 外层去掉大括号，按 '}' ',' '{' 粗分子块，再从每个子块里抓“引号中的名字”
            inner = char(refText(2:end-1));
            % 分割子块（逗号分隔的大块可能还是整个列表；下面按引号抓名字）
            names = regexp(inner, '''([^'']+)''', 'tokens');
            if ~isempty(names)
                nm = cellfun(@(c)c{1}, names, 'uni', 0);
                % 试图判断是否是分组形式：出现 `}{` 或 `}},{{` 视作多组
                isGrouped = contains(inner, '}{') || contains(inner, '},{');
                if ~isGrouped
                    % 视作 channels
                    idx = names2indices(nm, chanLabels);
                    assert(~isempty(idx), '未匹配到任何通道名。');
                    refSpec.type = 'channels';
                    refSpec.idx  = idx;
                else
                    % 分组：按 `}` 与 `{` 重新切块
                    rawGroups = regexp(inner, '\{([^}]+)\}', 'tokens');
                    groups = cell(1, numel(rawGroups));
                    for g = 1:numel(rawGroups)
                        nn = regexp(rawGroups{g}{1}, '''([^'']+)''', 'tokens');
                        nm_g = cellfun(@(c)c{1}, nn, 'uni', 0);
                        groups{g} = names2indices(nm_g, chanLabels);
                    end
                    refSpec.type   = 'groups';
                    refSpec.groups = groups;
                end
                return;
            end
        end
    end

    error('无法解析参考参数。请参照示例输入：''A1''、{''A1'',''A2''}、''Average''、[1,3,5]、{[1,3,5],[2,4]}');

end


function s = stripQuotes(s)
    s = string(s);
    s = strtrim(s);
    if (startsWith(s,"'") && endsWith(s,"'")) || (startsWith(s,"""") && endsWith(s,""""))
        s = extractBetween(s, 2, strlength(s)-1);
        s = string(s);
    end
end

function checkIndexRange(idx, nCh)
    assert(all(idx>=1 & idx<=nCh), '索引越界（1..%d）。', nCh);
end

function idx = name2index(name, chanLabels)
    idx = [];
    if isempty(chanLabels), return; end
    hit = find(strcmpi(strtrim(name), strtrim(chanLabels)));
    if ~isempty(hit), idx = hit(1); end
end

function idx = names2indices(namesCell, chanLabels)
    idx = [];
    for i=1:numel(namesCell)
        j = name2index(namesCell{i}, chanLabels);
        if ~isempty(j), idx(end+1)=j; end %#ok<AGROW>
    end
    idx = unique(idx);
end

