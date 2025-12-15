classdef CortexNode < BaseNode
    %CORTEXNODE 皮层节点
    % 管理协议的皮层信息，包含皮层名称、位置等元数据
    
    properties (Constant)
        type = "CortexNode"  % 节点类型标识符
    end
    
    properties
        % 皮层特定属性
        cortexInfo CortexInfo
        cache
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        data
        infoFile string
    end
    
    methods
        function obj = CortexNode()
            %CORTEXNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, path)
            %OPEN 打开数据
            
            % 检查路径是文件夹还是文件
            if isfile(path)
                [path, name, ~] = fileparts(path);
            end

            % 设置路径
            obj.path = path;
            
            % 验证协议文件存在
            if ~isfile(obj.infoFile)
                error('SEAL:NoFile', ...
                    'No such file in path: %s', obj.infoFile);
            end
            
            try
                obj.cortexInfo = CortexInfo.openExisting(obj.infoFile);
            catch ME
                error('SEAL:ChannelNode:LoadFailed', ...
                    'Failed to open channel data: %s', ME.message);
            end
        end


        function save(obj)
            %SAVE 保存皮层数据
            obj.createDirectoryStructure();
            obj.cortexInfo.save(obj.infoFile);
        end

        function load(obj)
            %LOAD 加载皮层数据
            if obj.isLoaded
                return;
            end
            
            % 加载皮层元数据
            if ~isempty(obj.cortexInfo)
                [~, name, ~] = fileparts(obj.cortexInfo.dataPath);
                if isfile(obj.cortexInfo.dataPath)
                    obj.cache = loadData(obj.cortexInfo.dataPath);
                    obj.isLoaded = true;
                end
            end
        end

        function unload(obj)
            %UNLOAD 卸载皮层数据，释放内存
            obj.cache = [];
            obj.isLoaded = false;
        end
        
        %% 依赖属性get方法
        function cortexName = get.name(obj)
            cortexName = obj.cortexInfo.name;
        end
        
        function data = get.data(obj)
            if ~obj.isLoaded
                obj.load();
            end
            data = obj.cache;
        end

        function path = get.infoFile(obj)
            path = fullfile(obj.path, "cortex_info.mat");
        end
    
    end

    methods (Static)

        function cortex = fromInfo(cortexInfo)
            cortex = CortexNode();
            cortex.cortexInfo = cortexInfo();
        end

        function cortex = fromData(path)
            cortex = CortexNode();
            cortex.cortexInfo = CortexInfo.fromData(path);
        end

        function cortex = openExisting(srcPath)
            %OPENCORTEX 打开现有皮层数据
            % 输入:
            %   srcPath - 皮层文件路径或皮层文件夹路径
            %   parentProtocol - 父协议节点
            
            % 创建皮层节点
            cortex = CortexNode();

            % 加载皮层数据
            cortex.open(getmat(srcPath));
        end
    end

    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建皮层目录结构
            
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