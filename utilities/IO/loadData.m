function data = loadData(filename)
    % 加载MAT文件中的所有变量到结构体
    loaded = load(filename);
    
    fields = fieldnames(loaded);
    numFields = length(fields);
    
    if numFields == 1
        data = loaded.(fields{1});
    else
        error('MAT文件包含 %d 个字段，只允许包含1个字段\n包含的字段: %s', ...
              numFields, strjoin(fields, ', '));
    end
end