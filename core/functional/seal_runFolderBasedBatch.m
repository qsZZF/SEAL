function seal_runFolderBasedBatch(app,  parentFolder, taskRecipe)
% ===================================================
% 1. 扫描总文件夹下的所有子文件夹 (实验单元)
% ===================================================
allItems = dir(parentFolder);
% 过滤掉 '.' 和 '..' 以及非文件夹的项
projectNode = app.projectNode;
subDirs = allItems([allItems.isdir] & ~startsWith({allItems.name}, '.'));

if isempty(subDirs)
    uialert(app.UIFigure, '选中的文件夹内没有发现任何有效的子文件夹！', '批处理失败', 'Icon', 'error');
    return;
end

totalTasks = length(subDirs);
wb = waitbar(0, '正在初始化批处理引擎...');

for i = 1:totalTasks
    currentFolderName = subDirs(i).name;
    currentFolderPath = fullfile(subDirs(i).folder, currentFolderName);

    try
        waitbar(i/totalTasks, wb, sprintf('正在构建与计算 (%d/%d): %s', i, totalTasks, currentFolderName));

        dataFile = dir(fullfile(currentFolderPath, '*data*.mat'));
        cortexFile = dir(fullfile(currentFolderPath, '*cortex*.mat'));
        lfFile = dir(fullfile(currentFolderPath, '*leadfield*.mat'));
        noise_cov_file = dir(fullfile(currentFolderPath, '*noise_cov*.mat'));

        if isempty(dataFile) || isempty(cortexFile) || isempty(lfFile)
            warning('文件夹跳过 [%s]: 未凑齐 Data, Cortex 和 Leadfield 三个必需文件。', currentFolderName);
            continue;
        end


        cortexFullPath = fullfile(cortexFile(1).folder, cortexFile(1).name);
        lfFullPath = fullfile(lfFile(1).folder, lfFile(1).name);

        newProtocol = projectNode.createNewProtocol(currentFolderName);

        for j=1:length(dataFile)
            dataFullPath = fullfile(dataFile(j).folder, dataFile(j).name);
            dNode = newProtocol.openDataFromData(dataFullPath);
            cNode = newProtocol.openCortexFromData(cortexFullPath);
            lNode = newProtocol.openLeadfieldFromData(lfFullPath);

            app.flushFileTree();

            signalMatrix = dNode.data;
            L_matrix = lNode.data;
            CortexModel = cNode.data;
            currentSession = dNode.parent;

            if isempty(noise_cov_file)
                % 1. 如果没有找到文件，默认使用单位矩阵
                NoiseCovariance = eye(size(signalMatrix, 1));
            else

                noise_covFullPath = fullfile(noise_cov_file(1).folder, noise_cov_file(1).name);

                % 3. 安全加载：load 返回的是一个结构体 (struct)
                tempData = load(noise_covFullPath);

                % 4. 获取结构体内的所有变量名
                fns = fieldnames(tempData);

                if ~isempty(fns)
                    % 5. 动态提取文件里的第一个变量（获取真正的矩阵）
                    NoiseCovariance = tempData.(fns{1});

                    % 6. 现在对真正的矩阵进行维度校验
                    if size(NoiseCovariance, 1) ~= size(signalMatrix, 1) || size(NoiseCovariance, 2) ~= size(signalMatrix, 1)
                        warning('加载的协方差矩阵维度 (%d x %d) 与通道数 (%d) 不匹配！将回退使用单位矩阵。', ...
                            size(NoiseCovariance,1), size(NoiseCovariance,2), size(signalMatrix,1));
                        NoiseCovariance = eye(size(signalMatrix, 1));
                    end
                else
                    % 防御：如果用户存了一个空的 .mat 文件
                    warning('noise_cov 文件内没有任何变量！将回退使用单位矩阵。');
                    NoiseCovariance = eye(size(signalMatrix, 1));
                end
            end


            if ndims(signalMatrix)==3
                signalMatrix= concatenate_trials(signalMatrix);
            end


            for n=1:length(taskRecipe)
                taskRecipe(n).Parameters.NoiseCovariance = NoiseCovariance;

                paramsStruct = taskRecipe(n).Parameters;
                currentStruct = paramsStruct;
                if isempty(paramsStruct)
                    paramCell = {}; % 如果没有参数，就给个空元胞
                else
                    % 去除字段名和对应的值
                    fieldsToRemove = {'srate', 'Algorithm','name','ResultFolderPath'};
                    if any(isfield(paramsStruct, fieldsToRemove))
                        paramsStruct = rmfield(paramsStruct, intersect(fieldnames(paramsStruct), fieldsToRemove));
                    end

                    pNames = fieldnames(paramsStruct);
                    pVals = struct2cell(paramsStruct);

                    % 交替拼装成 {'Name1', Value1, 'Name2', Value2...} 的格式
                    paramCell = cell(1, 2 * length(pNames));
                    paramCell(1:2:end) = pNames;
                    paramCell(2:2:end) = pVals;
                end

                numInputs = abs(nargin(taskRecipe(n).AlgorithmName));

                try
                    if numInputs == 4
                        [sourceEstimate, ~] = feval(taskRecipe(n).AlgorithmName, ...
                            signalMatrix, L_matrix, CortexModel, paramCell{:});

                    elseif numInputs == 3
                        [sourceEstimate, ~] = feval(taskRecipe(n).AlgorithmName, ...
                            signalMatrix, L_matrix, paramCell{:});

                    else
                        error('算法 [%s] 的输入参数数量 (%d) 不符合平台标准规范。', ...
                            taskRecipe(n).AlgorithmName, numInputs);
                    end
                catch ME
                    error('算法 [%s] 执行失败: %s', taskRecipe(n).AlgorithmName, ME.message);
                end

                if paramsStruct.NumOrientations == 3
                    temp = reshape(sourceEstimate, 3, [], size(sourceEstimate, 2));
                    sourceEstimate = squeeze(sum(temp.^2, 1));    % 沿第一维（3行）计算平方和
                end

                savepath = fullfile(currentSession.path, "Results");
                resultFileName = sprintf('%s_Result.mat',  taskRecipe(n).AlgorithmName);
                resultFilePath = fullfile(savepath, resultFileName);
                currentStruct.name = resultFileName;
                currentSession.openResultFromData(resultFilePath,currentStruct);
                app.flushFileTree();
                resultDir = fileparts(resultFilePath); % 瞬间剥离出纯文件夹路径
                if ~isfolder(resultDir)
                    mkdir(resultDir); % 如果不存在，MATLAB 会连环创建沿途所有缺失的文件夹！
                end
                % 结果存盘
                save(resultFilePath, 'sourceEstimate');
            end
            %%

            dNode.unload();
            cNode.unload();
            lNode.unload();
            clear signalMatrix L_matrix sourceEstimate extraInfo;
        end
    catch ME
        warning('文件夹 [%s] 处理失败跳过: %s', currentFolderName, ME.message);
        continue; % 报错不中断整个循环，继续处理下一个文件夹
    end
end

close(wb);
uialert(app.UIFigure, '所有的批处理任务已执行完毕！', '批处理完成', 'Icon', 'success');
end


function concatenated_data = concatenate_trials(data_3d)

% 获取数组尺寸
[n_channels, n_points_per_trial, n_trials] = size(data_3d);

% 预分配拼接后的二维矩阵
total_points = n_points_per_trial * n_trials;
concatenated_data = zeros(n_channels, total_points);

% 循环处理每个试次
for trial = 1:n_trials
    % 计算当前试次的起始和结束索引
    start_idx = (trial - 1) * n_points_per_trial + 1;
    end_idx = trial * n_points_per_trial;

    % 提取当前试次数据
    current_trial_data = data_3d(:, :, trial);

    % 将数据放入拼接矩阵
    concatenated_data(:, start_idx:end_idx) = current_trial_data;

end
end