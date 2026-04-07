classdef DataNode < BaseNode
    %DATANODE 节点
    % 管理协议的信息，包含名称、位置等元数据
    
    properties (Constant)
        type = "DataNode"  % 节点类型标识符
    end
    
    properties
        % 特定属性
        dataInfo DataInfo
        data_cache
        res_cache
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        data
        result
        infoFile strings
        size
        srate
    end
    
    methods
        function obj = DataNode()
            %DATANODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, path)
            %OPEN 打开数据
            
            file = path;
            [path, name, ~] = fileparts(path);
            
            % 验证协议文件存在
            if ~isfile(file)
                error('SEAL:NoFile', ...
                    'No such file in path: %s', path);
            end

            % 设置路径
            obj.path = path;
            
            try
                obj.dataInfo = DataInfo.openExisting(file);
            catch ME
                error('SEAL:DataNode:LoadFailed', ...
                    'Failed to open data data: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存数据
            obj.createDirectoryStructure();
            obj.dataInfo.save(obj.infoFile);
        end

        function load(obj)
            %LOAD 加载数据
            if obj.isLoaded
                return;
            end

            % 加载元数据
            if ~isempty(obj.dataInfo)
                if isfile(obj.dataInfo.dataPath)
                    % 先把文件读进临时变量
                    rawData = load(obj.dataInfo.dataPath, '-mat'); % 或者用 load()

                    % --- 自动拆开“俄罗斯套娃” ---
                    outerFields = fieldnames(rawData);
                    if length(outerFields) == 1 && isstruct(rawData.(outerFields{1}))
                        rawData = rawData.(outerFields{1});
                    end

                    % --- 智能字段匹配 ---
                    availableFields = fieldnames(rawData);
                    dataAliases = {'data', 'Data', 'DataFile', 'matrix', 'signal', 'EEG'};

                    matchedField = '';
                    for i = 1:length(dataAliases)
                        if ismember(dataAliases{i}, availableFields)
                            matchedField = dataAliases{i};
                            break;
                        end
                    end

                    % --- 最终赋值给缓存 ---
                    if ~isempty(matchedField)
                        tempData = rawData.(matchedField);
                    else
                        % 如果实在找不到，弹窗让用户选（调用我们上次说的防御机制）
                        [indx, tf] = listdlg('ListString', availableFields, ...
                            'SelectionMode', 'single', ...
                            'PromptString', {'数据加载失败！未识别到标准矩阵字段。', '请手动指定代表脑电信号的变量:'}, ...
                            'Name', '手动匹配字段');
                        if tf
                            obj.data_cache = rawData.(availableFields{indx});
                        else
                            error('SEAL:DataNode:LoadFailed', '用户取消了字段匹配，数据加载中止。');
                        end
                    end
                    [~, ~, fileExt] = fileparts(obj.dataInfo.dataPath);

                    if (ischar(tempData) || isstring(tempData)) && strcmpi(fileExt, '.set')
                        % 鉴定结果：它只是个字符串！说明真实数据在外部的 .fdt 文件里
                        fdtFilename = char(tempData);
                        [setDir, ~, ~] = fileparts(obj.dataInfo.dataPath);
                        fdtPath = fullfile(setDir, fdtFilename);

                        if isfile(fdtPath)
                            % 打开并读取 IEEE Little-Endian 的 32位浮点数二进制文件
                            fileID = fopen(fdtPath, 'r', 'ieee-le');
                            if fileID == -1
                                error('SEAL:DataNode:FDTError', '无法打开 .fdt 文件读取权限: %s', fdtPath);
                            end

                            try
                                rawBin = fread(fileID, Inf, 'float32');
                                fclose(fileID);

                                % 尝试利用 rawData (即 EEG 结构体) 里的参数重构 3D 矩阵
                                if isfield(rawData, 'nbchan') && isfield(rawData, 'pnts') && isfield(rawData, 'trials')
                                    obj.data_cache = reshape(rawBin, [rawData.nbchan, rawData.pnts, rawData.trials]);
                                elseif isprop(obj.dataInfo, 'size') && ~isempty(obj.dataInfo.size)
                                    % 备用方案：使用 DataInfo 里提前存好的 size
                                    obj.data_cache = reshape(rawBin, obj.dataInfo.size);
                                else
                                    warning('SEAL:DataNode:MissingDimensions', '找不到重构矩阵的维度参数，数据可能变形。');
                                    obj.data_cache = rawBin;
                                end

                            catch ME
                                if exist('fileID', 'var') && fileID ~= -1, fclose(fileID); end
                                error('读取 .fdt 二进制文件失败: %s', ME.message);
                            end
                        else
                            error('SEAL:DataNode:MissingFDT', '找不到依赖的外部数据文件: %s', fdtPath);
                        end
                    else
                        % 鉴定结果：是真正的数值矩阵 (单文件 .set 或 普通 .mat)
                        % 直接安全赋值
                        obj.data_cache = tempData;
                    end
                end
                obj.isLoaded = true;
            end
        end

        function unload(obj)
            %UNLOAD 卸载数据，释放内存
            obj.data_cache = [];
            obj.isLoaded = false;
        end

        function deleteFromDisk(obj)
            delete(obj.infoFile);
        end
        
        %% 依赖属性get方法
        function dataName = get.name(obj)
            dataName = obj.dataInfo.name;
        end

        function rsize = get.size(obj)
            rsize = obj.dataInfo.size;
        end

        function srate = get.srate(obj)
            srate = obj.dataInfo.srate;
        end
        
        function data = get.data(obj)
            if ~obj.isLoaded
                obj.load();
            end
            data = obj.data_cache;
        end

        function result = get.result(obj)
            if ~obj.isLoaded
                obj.load();
            end
            result = obj.res_cache;
        end

        function setResultPath(obj, path)
            obj.dataInfo.resultPath = path;
        end

        function path = get.infoFile(obj)
            path = fullfile(obj.path, strcat(obj.name, ".mat"));
        end

        function metadata = getMetadata(obj)
            metadata =  obj.getMetadata@BaseNode();
            metadata = [metadata;
                "Size", mat2str(obj.size);
                "Sampling Rate", string(obj.srate)];
        end
    end

    methods (Static)
        function data = fromInfo(dataInfo)
            data = DataNode();
            data.dataInfo = dataInfo;
        end

        function data = fromData(path)
            data = DataNode();
            data.dataInfo = DataInfo.fromData(path);
        end

        function data = openExisting(srcPath)
            %OPENDATA 打开现有数据
            % 输入:
            %   srcPath - 文件路径或文件夹路径
            %   parentProtocol - 父协议节点
            
            % 创建节点
            data = DataNode();

            % 加载数据
            data.open(getmat(srcPath));
        end
    end

    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建目录结构
            
            if obj.path == ""
                return;
            end
            
            % 创建主目录
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
        end
    end
end