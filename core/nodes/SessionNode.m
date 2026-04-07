classdef SessionNode < BaseNode
    %SESSIONNODE 会话节点
    % 管理单个会话，包含原始数据和预处理数据
    
    properties (Constant)
        type = "SessionNode"  % 节点类型标识符
    end
    
    properties
        % 会话特定属性
        sessionInfo SessionInfo
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        sessionType string
        desc string
        dataCount int32             % 总数据节点数量
        
        infoFile string
    end
    
    methods
        function obj = SessionNode()
            %SESSIONNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, sessionPath)
            %OPEN 打开会话，加载轻量化的Info数据
            % 输入:
            %   sessionPath - 会话路径
            
            if obj.isLoaded
                return;
            end
            
            try
                % 检查路径是文件夹还是文件
                if isfile(sessionPath)
                    [sessionPath, name, ~] = fileparts(sessionPath);
                end

                % 设置路径
                obj.path = sessionPath;
                
                % 验证会话文件存在
                if ~isfile(obj.infoFile)
                    error('SEAL:SessionNode:NoSessionFile', ...
                        'No such session in path: %s', obj.infoFile);
                end
                
                % 加载会话元数据
                obj.sessionInfo = SessionInfo.openExisting(obj.infoFile);
                obj.isLoaded = false;
                
                % 递归打开子节点
                obj.openChildNodes();
                
            catch ME
                error('SEAL:SessionNode:OpenFailed', ...
                    'Failed to open session: %s', ME.message);
            end
        end

        function load(obj)
            %LOAD 加载会话完整数据，用于后续操作
            
            if obj.isLoaded
                return;
            end
            
            try
                % 加载所有子节点的数据
                for i = 1:obj.childCount
                    if ~obj.children(i).isLoaded
                        obj.children(i).load();
                    end
                end
                
                obj.isLoaded = true;
                
            catch ME
                error('SEAL:SessionNode:LoadFailed', ...
                    'Failed to load session data: %s', ME.message);
            end
        end

        function unload(obj)
            %UNLOAD 卸载会话数据，释放内存
            
            try
                % 递归卸载所有子节点的数据
                for i = 1:obj.childCount
                    if obj.children(i).isLoaded
                        obj.children(i).unload();
                    end
                end
                
                % 标记为未加载
                obj.isLoaded = false;
                
            catch ME
                warning('SEAL:SessionNode:UnloadFailed', ...
                    'Failed to unload session data: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存会话数据
            
            try
                obj.createDirectoryStructure();
                obj.sessionInfo.save(obj.infoFile);
                
                % 递归保存所有子节点
                for i = 1:obj.childCount
                    obj.children(i).save();
                end
                
            catch ME
                error('SEAL:SessionNode:SaveFailed', ...
                    'Failed to save session: %s', ME.message);
            end
        end

        function deleteFromDisk(obj)
            rmdir(obj.path, 's');
        end
        
        %% 数据节点管理
        function addDataNode(obj, dataNode)
            %ADDDATANODE 添加数据节点
            % 输入:
            %   dataNode - DataNode对象
            
            % 验证节点类型
            if ~isa(dataNode, 'DataNode')
                error('SEAL:SessionNode:InvalidDataNode', ...
                    '只能添加DataNode类型的子节点');
            end
            
            % 添加到子节点列表
            obj.addChild(dataNode);
            
%             % 如果是原始数据，更新会话信息
%             if dataNode.isOriginal
%                 obj.sessionInfo.metadata.hasOriginalData = true;
%             end
        end
        
        function removeDataNode(obj, dataNode)
            %REMOVEDATANODE 移除数据节点
            % 输入:
            %   dataNode - 要移除的DataNode对象或名称
            
            if ischar(dataNode) || isstring(dataNode)
                % 根据名称查找节点
                dataNode = obj.getDataNodeByName(dataNode);
            end
            
            if isempty(dataNode)
                warning('SEAL:SessionNode:DataNodeNotFound', ...
                    '未找到指定的数据节点');
                return;
            end
            
            % 如果是原始数据，更新会话信息
            if dataNode.isOriginal
                obj.sessionInfo.metadata.hasOriginalData = false;
            end
            
            % 从子节点列表中移除
            obj.removeChild(dataNode);
        end
        
        function dataNode = getDataNodeByName(obj, name)
            %GETDATANODEBYNAME 根据名称获取数据节点
            % 输入:
            %   name - 数据节点名称
            % 输出:
            %   dataNode - DataNode对象（如果存在）
            
            dataNode = [];
            for i = 1:obj.childCount
                if strcmp(obj.children(i).name, name) && isa(obj.children(i), 'DataNode')
                    dataNode = obj.children(i);
                    break;
                end
            end
        end
        
        %% 依赖属性get方法
        function name = get.name(obj)
            if ~isempty(obj.sessionInfo)
                name = obj.sessionInfo.name;
            else
                name = "Unnamed Session";
            end
        end

        function path = get.infoFile(obj)
            [~, name, ~] = fileparts(obj.path);
            path = fullfile(obj.path, "session_info.mat");
        end
        
        function sessionType = get.sessionType(obj)
            if ~isempty(obj.sessionInfo)
                sessionType = obj.sessionInfo.type;
            else
                sessionType = "";
            end
        end
        
        function desc = get.desc(obj)
            if ~isempty(obj.sessionInfo)
                desc = obj.sessionInfo.desc;
            else
                desc = "";
            end
        end
        
        function count = get.dataCount(obj)
            % 统计数据节点数量
            count = 0;
            for i = 1:obj.childCount
                if isa(obj.children(i), 'DataNode')
                    count = count + 1;
                end
            end
        end

        function node = openDataFromData(obj, dataPath)
            dataInfo = DataInfo.fromData(dataPath);
            dataNode = DataNode.fromInfo(dataInfo);
            dataNode.path = obj.path;
            obj.addChild(dataNode)
            dataNode.save();
            node = dataNode;
        end

        function node = openResultFromData(obj, datapath,paramStruct)
            resultInfo = ResultInfo.fromData(datapath,paramStruct);
            resultNode = ResultNode.fromResultData(resultInfo);
            resultNode.path = datapath;
            resultNode.sessionpath = obj.path;
            obj.addChild(resultNode);
            resultNode.save();
            node = resultNode;
        end

        function node = openDataFromEEGLAB(obj, dataPath)
            dataInfo = DataInfo.fromEEGLAB(dataPath);
            dataNode = DataNode.fromInfo(dataInfo);
            dataNode.path = obj.path;
            obj.addChild(dataNode)
            dataNode.save();
            node = dataNode;
        end
    end

    methods (Static)
        function session = openExisting(srcPath)
            %OPENEXISTING 打开现有会话
            % 输入:
            %   srcPath - 会话文件路径或会话文件夹路径
            % 输出:
            %   session - SessionNode对象
            
            % 创建会话节点
            session = SessionNode();
            
            % 打开会话
            session.open(srcPath);
        end
        
        function session = createNew(name, path, varargin)
            %CREATENEW 创建新会话
            % 输入:
            %   name - 会话名称
            %   path - 存储路径（通常是Protocol的Sessions文件夹）
            %   varargin - 可选参数对
            
            p = inputParser;
            addRequired(p, 'name', Validators.isText);
            addRequired(p, 'path', Validators.isText);
            addParameter(p, 'Type', '', Validators.isText);
            addParameter(p, 'Description', '', Validators.isText);
            
            parse(p, name, path, varargin{:});
            
            % 构建会话路径
            sessionPath = fullfile(path, name);
            
            % 创建会话节点
            session = SessionNode();
            session.path = sessionPath;
            
            % 创建会话信息
            session.sessionInfo = SessionInfo.createNew(...
                p.Results.name, ...
                p.Results.Type, ...
                p.Results.Description);
            
            % 保存会话基本信息
            session.save();
            
            % 打开会话
            session.open(sessionPath);
        end
    end
    
    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建会话目录结构
            
            if obj.path == ""
                return;
            end

            % 创建主目录
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
        end

        function openChildNodes(obj)
            %OPENCHILDNODES 智能打开子节点
            % 这个方法会在打开会话时自动扫描，并区分数据节点和结果节点

            files = dir(fullfile(obj.path, '*.mat'));
            % 过滤掉 session_info 等配置文件
            dataFiles = files(~cellfun(@isempty, regexp({files.name}, '^((?!_info).)*\.mat$')));

            for i = 1:length(dataFiles)
                filePath = fullfile(obj.path, dataFiles(i).name);

                try
                    % 1. 【雷达探测】：不加载数据，先偷看文件里面装了什么类
                    fileVars = whos('-file', filePath);
                    varClasses = {fileVars.class};

                    % 2. 【智能分拣】：根据类名调用不同的打开函数
                    if ismember('ResultInfo', varClasses) || ismember('ResultNode', varClasses)
                        % 如果文件里包含 ResultInfo，说明它是结果节点
                        childNode = ResultNode.openExisting(filePath);
                    else
                        % 否则，默认作为原始脑电数据节点打开
                        % (假设 DataInfo 封装在里面，或者只是个普通的包含 data 的结构体)
                        childNode = DataNode.openExisting(filePath);
                    end

                    % 3. 挂载到子节点并建立父子关系
                    childNode.parent = obj;
                    obj.addChild(childNode);

                catch ME
                    % 降级为 warning，保证某一个文件损坏不会导致整个 App 崩溃启动不了
                    warning('SEAL:SessionNode:OpenChildFailed', ...
                        'Failed to open node file: %s.\nError: %s', ...
                        filePath, ME.message);
                end
            end
        end
    end
end