classdef LeadfieldNode < BaseNode
    %LEADFIELDNODE 导联场节点
    % 管理协议的导联场信息，包含导联场名称、位置等元数据
    
    properties (Constant)
        type = "LeadfieldNode"  % 节点类型标识符
    end
    
    properties
        % 导联场特定属性
        leadfieldInfo LeadfieldInfo
        cache
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        data
        infoFile string
    end
    
    methods
        function obj = LeadfieldNode()
            %LEADFIELDNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, path)
            %OPEN 打开数据
            
            % 检查路径是文件夹还是文件
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
                obj.leadfieldInfo = LeadfieldInfo.openExisting(obj.infoFile);
            catch ME
                error('SEAL:LeadfieldNode:LoadFailed', ...
                    'Failed to open leadfield data: %s', ME.message);
            end
        end


        function save(obj)
            %SAVE 保存导联场数据
            obj.createDirectoryStructure();
            obj.leadfieldInfo.save(obj.infoFile);
        end

        function load(obj)
            %LOAD 加载导联场数据
            if obj.isLoaded
                return;
            end
            
            % 加载导联场元数据
            if ~isempty(obj.leadfieldInfo)
                [~, name, ~] = fileparts(obj.leadfieldInfo.dataPath);
                if isfile(obj.leadfieldInfo.dataPath)
                    obj.cache = loadData(obj.leadfieldInfo.dataPath);
                    obj.isLoaded = true;
                end
            end
        end

        function unload(obj)
            %UNLOAD 卸载导联场数据，释放内存
            obj.cache = [];
            obj.isLoaded = false;
        end
        
        %% 依赖属性get方法
        function leadfieldName = get.name(obj)
            leadfieldName = obj.leadfieldInfo.name;
        end
        
        function data = get.data(obj)
            if ~obj.isLoaded
                obj.load();
            end
            data = obj.cache;
        end

        function path = get.infoFile(obj)
            path = fullfile(obj.path, "leadfield_info.mat");
        end
    
    end

    methods (Static)
        
        function leadfield = fromInfo(leadfieldInfo)
            leadfield = LeadfieldNode();
            leadfield.leadfieldInfo = leadfieldInfo();
        end

        function leadfield = openExisting(srcPath)
            %OPENLEADFIELD 打开现有导联场数据
            % 输入:
            %   srcPath - 导联场文件路径或导联场文件夹路径
            %   parentProtocol - 父协议节点
            
            % 创建导联场节点
            leadfield = LeadfieldNode();

            % 加载导联场数据
            leadfield.open(getmat(srcPath));
        end
    end

    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建导联场目录结构
            
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