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

        infoFile string

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
                    obj.data_cache = loadData(obj.dataInfo.dataPath);
                end
                if isfile(obj.dataInfo.resultPath)
                    obj.res_cache = loadData(obj.dataInfo.resultPath);
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
            data.dataInfo = dataInfo();
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