function idx = parseChanList(text, nCh, labels)
    if ~(isstring(text)||ischar(text)), text = string(text); end
    text = strtrim(text);
    if strlength(text)==0, idx = []; return; end

    % 中文标点替换
    text = replace(text, ["，","；","、"], ",");
    text = replace(text, ["—","–","－"], "-");

    tokens = strtrim(split(text, ','));
    idx = [];
    for i=1:numel(tokens)
        tk = tokens(i);
        if contains(tk,"-")
            ab = str2double(split(tk,"-"));
            if numel(ab)==2 && all(isfinite(ab)), idx = [idx, ab(1):ab(2)]; end
        else
            v = str2double(tk);
            if isfinite(v)
                idx = [idx, v];
            else
                % 按标签匹配
                hit = find(strcmpi(strtrim(tk), labels));
                if ~isempty(hit), idx = [idx, hit(1)]; end
            end
        end
    end
    idx = unique(idx);
    idx = idx(idx>=1 & idx<=nCh);
end
