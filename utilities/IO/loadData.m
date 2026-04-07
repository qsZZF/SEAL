function data = loadData(filename)
    % 加载MAT文件中的所有变量到结构体
    loaded = load(filename);
    
    fields = fieldnames(loaded);
    numFields = length(fields);
    
    if numFields == 1
        data = loaded.(fields{1});
    else
        data = loaded;
    end
end