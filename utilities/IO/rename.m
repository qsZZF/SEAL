function success = rename(src, tgt)
    % RENAME 重命名文件或文件夹
    %   如果源路径不存在，不进行任何操作
    %   不会修改扩展名

    success = false;
    
    src = strtrim(src);
    tgt = strtrim(tgt);
    
    if isempty(src) || isempty(tgt)
        error('Neither src nor tgt should be empty.');
    end
    
    % 检查源路径是否存在
    if ~exist(src, 'file')
        % 源不存在，静默返回
        return;
    end
    
    % 检查目标路径是否已存在
    if exist(tgt, 'file')
        error('Target path/file exists.');
    end
    
    % 执行重命名
    try
        [status, ~] = movefile(src, tgt);
        success = (status == 1);
    catch
        success = false;
    end
end