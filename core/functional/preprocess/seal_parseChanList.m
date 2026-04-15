function idx = seal_parseChanList(text, nCh, labels)
    % 确保输入为字符串类型
    if ~(isstring(text) || ischar(text)), text = string(text); end
    text = strtrim(text);
    if strlength(text) == 0, idx = []; return; end
    
    % 中文标点替换为英文标点
    text = replace(text, ["，","；","、"], ",");
    % 统一破折号/减号（注意：此处统一后，既可能表示区间，也可能是通道名的一部分）
    text = replace(text, ["—","–","－"], "-");
    
    tokens = strtrim(split(text, ','));
    idx = [];
    
    for i = 1:numel(tokens)
        tk = tokens(i);
        
        % 1. 优先按标签名称全字匹配（解决带连字符的通道名，如 'CP-1'）
        hit = find(strcmpi(tk, labels)); 
        if ~isempty(hit)
            % 如果找到了对应的通道名，直接取第一个匹配项的索引
            idx = [idx, hit(1)];
            continue; % 成功匹配名称后，直接跳过后续判断，处理下一个 token
        end
        
        % 2. 只有在名字没匹配上的情况下，才去解析数字逻辑
        if contains(tk, "-")
            % 尝试解析数字区间，例如 '1-5'
            ab = str2double(split(tk, "-"));
            if numel(ab) == 2 && all(isfinite(ab))
                % 保证正序生成索引（防止用户输入反向区间如 '5-1' 导致生成空数组）
                startIdx = min(ab(1), ab(2));
                endIdx = max(ab(1), ab(2));
                idx = [idx, startIdx:endIdx]; 
            end
        else
            % 3. 尝试解析纯单个数字，例如 '12'
            v = str2double(tk);
            if isfinite(v)
                idx = [idx, v];
            end
        end
        % ==============================================
    end
    
    % 去重并剔除超出合法范围的索引
    idx = unique(idx);
    idx = idx(idx >= 1 & idx <= nCh);
end