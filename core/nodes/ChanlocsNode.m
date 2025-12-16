classdef ChanlocsNode < BaseNode
    %CHANLOCSNODE 通道节点
    % 管理协议的通道信息，包含通道名称、位置等元数据
    
    properties (Constant)
        type = "ChanlocsNode"  % 节点类型标识符
    end
    
    properties
        % 通道特定属性
        chanlocsInfo ChanlocsInfo
        cache
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        data
        infoFile string
    end
    
    methods
        function obj = ChanlocsNode()
            %CHANLOCSNODE 构造函数
            
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
                obj.chanlocsInfo = ChanlocsInfo.openExisting(obj.infoFile);
            catch ME
                error('SEAL:ChanlocsNode:LoadFailed', ...
                    'Failed to open chanlocs data: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存通道数据
            obj.createDirectoryStructure();
            obj.chanlocsInfo.save(obj.infoFile);
        end

        function load(obj)
            %LOAD 加载通道数据
            if obj.isLoaded
                return;
            end
            
            % 加载通道元数据
            if ~isempty(obj.chanlocsInfo)
                [~, name, ~] = fileparts(obj.chanlocsInfo.dataPath);
                if isfile(obj.chanlocsInfo.dataPath)
                    obj.cache = loadData(obj.chanlocsInfo.dataPath);
                    obj.isLoaded = true;
                end
            end
        end

        function unload(obj)
            %UNLOAD 卸载通道数据，释放内存
            obj.cache = [];
            obj.isLoaded = false;
        end

        function deleteFromDisk(obj)
            delete(obj.infoFile);
        end
        
        %% 依赖属性get方法
        function chanlocsName = get.name(obj)
            chanlocsName = obj.chanlocsInfo.name;
        end
        
        function data = get.data(obj)
            if ~obj.isLoaded
                obj.load();
            end
            data = obj.cache;
        end

        function path = get.infoFile(obj)
            path = fullfile(obj.path, "chanlocs_info.mat");
        end
    
    end

    methods (Static)
        function chanlocs = fromInfo(chanlocsInfo)
            chanlocs = ChanlocsNode();
            chanlocs.chanlocsInfo = chanlocsInfo;
        end

        function chanlocs = fromData(path)
            chanlocs = ChanlocsNode();
            chanlocs.chanlocsInfo = ChanlocsInfo.fromData(path);
        end

        function chanlocs = openExisting(srcPath)
            %OPENCHANLOCS 打开现有通道数据
            % 输入:
            %   srcPath - 通道文件路径或通道文件夹路径
            %   parentProtocol - 父协议节点
            
            % 创建通道节点
            chanlocs = ChanlocsNode();

            % 加载通道数据
            chanlocs.open(getmat(srcPath));
        end
    end

    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建通道目录结构
            
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