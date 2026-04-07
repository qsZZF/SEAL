classdef ResultNode < BaseNode
    %RESULTNODE 结果节点
    % 管理计算结果在节点树中的状态，包含算法信息、参数和结果矩阵
    
    properties (Constant)
        type = "ResultNode"  % 节点类型标识符
    end
    
    properties
        % 特定属性
        resultInfo ResultInfo  % 关联的核心结果数据对象
        res_cache              % 用于缓存庞大的计算结果矩阵（如源激活矩阵）
        sessionpath string
    end
    
    properties (Dependent)
        % 依赖属性 (对外暴露的安全接口)
        name string
        algorithm string
        sourceActivity
        params struct
        infoFile string
       
    end
    
    methods
        function obj = ResultNode()
            %RESULTNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, path)
            %OPEN 打开结果数据
            
            file = path;
            [path, ~, ~] = fileparts(path);
            
            % 验证文件存在
            if ~isfile(file)
                error('SEAL:NoFile', ...
                    'No such file in path: %s', path);
            end
            % 设置路径
            obj.path = path;
            %obj.sessionpath = obj.parent.path;
            try
                
                obj.resultInfo = ResultInfo.openExisting(file);
            catch ME
                error('SEAL:ResultNode:LoadFailed', ...
                    'Failed to open result data: %s', ME.message);
            end
        end
        
        function save(obj)
            %SAVE 保存结果数据
            %obj.createDirectoryStructure();
            obj.resultInfo.save(obj.infoFile);
        end
        
        function load(obj)
            %LOAD 加载数据到内存
            if obj.isLoaded
                return;
            end
            
            % 加载缓存矩阵
            if ~isempty(obj.resultInfo)
                % 将 ResultData 中的源激活矩阵加载到缓存中
                % 这样可以支持在 unload 时清空大矩阵，而保留基础参数
                obj.res_cache = loadData(obj.resultInfo.datapath); 
                obj.isLoaded = true;
            end
        end
        
        function unload(obj)
            %UNLOAD 卸载数据，释放内存
            % 在做源成像时，源空间矩阵动辄几百MB，不看的时候一定要卸载
            obj.res_cache = [];
            obj.isLoaded = false;
        end
        
        function deleteFromDisk(obj)
            % 从硬盘物理删除
            if isfile(obj.infoFile)
                delete(obj.infoFile);
            end
            if isfile(obj.path)
                delete(obj.path);
            end
        end
        
        %% 依赖属性 get 方法 (实现属性与封装对象的联动)
        function resName = get.name(obj)
            resName = obj.resultInfo.name;
        end
        
        
        function algo = get.algorithm(obj)
            algo = obj.resultInfo.algorithm;
        end
        
        function p = get.params(obj)
            p = obj.resultInfo.params;
        end

        function act = get.sourceActivity(obj)
            if ~obj.isLoaded
                obj.load();
            end
            act = obj.res_cache;
        end
        
        function path = get.infoFile(obj)
            % 自动拼凑该节点对应的 .mat 文件路径
            path = fullfile(obj.sessionpath, strcat(obj.name, ".mat"));
        end
        
        function metadata = getMetadata(obj)
            %GETMETADATA 获取属性面板要展示的信息
            metadata =  obj.getMetadata@BaseNode();
            
            % 追加 Result 特有的元数据
            metadata = [metadata;
                "Algorithm", obj.algorithm];
               % "Matrix Size", mat2str(size(obj.sourceActivity))];
        end

    end
    
    %% 静态方法 (工厂模式创建节点)
    methods (Static)
        function node = fromResultData(resultInfo)
            % 从内存中直接生成一个新节点 (常用于刚计算完时)
            node = ResultNode();
            node.resultInfo = resultInfo;
        end
        
        function node = openExisting(srcPath)
            %OPENDATA 打开现有的结果文件并生成节点
            node = ResultNode();
            % 兼容你现有的 getmat 辅助函数机制
            node.open(srcPath); 
        end
    end
    
    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建目录结构
            if obj.path == ""
                return;
            end
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
        end
    end
end