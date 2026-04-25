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

    unrecognized = {};
    for i = 1:numel(tokens)
        tk = tokens(i);
        matched = false;

        % 名称匹配
        hit = find(strcmpi(tk, labels));
        if ~isempty(hit)
            idx = [idx, hit(1)];
            matched = true;
        elseif contains(tk, "-")
            ab = str2double(split(tk, "-"));
            if numel(ab) == 2 && all(isfinite(ab))
                idx = [idx, min(ab):max(ab)];
                matched = true;
            end
        else
            v = str2double(tk);
            if isfinite(v)
                idx = [idx, v];
                matched = true;
            end
        end

        if ~matched
            unrecognized{end+1} = char(tk);
        end
    end

    if ~isempty(unrecognized)
        warning('SEAL:ParseChan:Unrecognized', ...
            '以下 token 无法识别(已忽略): %s', strjoin(unrecognized, ', '));
    end

    % 去重并剔除超出合法范围的索引
    idx = unique(idx);
    idx = idx(idx >= 1 & idx <= nCh);
end