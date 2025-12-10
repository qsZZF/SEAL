function filesrcPath = getmat(srcPath)
    if isfolder(srcPath)
        % 如果是文件夹，查找项目文件
        [~, projectName, ~] = fileparts( srcPath);
         Files = dir(fullfile( srcPath, strcat(projectName, '.mat')));
        if isempty( Files)
            error('SEAL:IO:FileNotFound', ...
                'No such file in path: %s',  srcPath);
        end
    else
         [ srcPath, projectName, ~] = fileparts( srcPath);
    end
    filesrcPath = fullfile( srcPath, strcat(projectName, ".mat"));
end