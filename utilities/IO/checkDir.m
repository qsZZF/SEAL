function checkDir(path)
    if path == ""
        error('SEAL:InvalidPath', ...
            'No valid path specified.');
    end
    
    path = strtrim(path);
    if ~isempty(path) && (path(end) == '/' || path(end) == '\')
        path = path(1:end-1);
    end
    
    if isempty(path)
        error('SEAL:InvalidPath', ...
            'No valid path specified.');
    end
    
    [folderPath, ~, ext] = fileparts(path);
    
    if ~isempty(ext)
        if ~isempty(folderPath) && ~isfolder(folderPath)
            mkdir(folderPath);
        end
    else
        if ~isfolder(path)
            mkdir(path);
        end
    end
end