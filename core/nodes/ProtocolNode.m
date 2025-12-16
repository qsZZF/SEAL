classdef ProtocolNode < BaseNode
    %PROTOCOLNODE 协议节点
    % 管理研究协议，包含会话、皮层模型、导联场等
    
    properties (Constant)
        type = "ProtocolNode"  % 节点类型标识符
    end
    
    properties
        % 协议特定属性
        protocolInfo ProtocolInfo
        chanlocsNode ChanlocsNode           % 通道节点
        cortexNode CortexNode             % 皮层节点
        leadfieldNode LeadfieldNode       % 导联场节点
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        protocolType string
        desc string
        sessionCount int32                 % 会话数量

        chanlocsInfoFile string
        cortexInfoFile string
        leadfieldInfoFile string

        infoFile string
        sessionFolder string
    end
    
    methods
        function obj = ProtocolNode()
            %PROTOCOLNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
            
            % 初始化属性
            obj.chanlocsNode = ChanlocsNode.empty;
            obj.cortexNode = CortexNode.empty;
            obj.leadfieldNode = LeadfieldNode.empty;
        end
        
        function open(obj, protocolPath)
            %OPEN 打开协议，加载轻量化的Info数据
            % 输入:
            %   protocolPath - 协议路径
            
            if obj.isLoaded
                return;
            end
            
            try
                % 检查路径是文件夹还是文件
                if isfile(protocolPath)
                    % 如果是文件路径，提取文件夹路径
                    [protocolPath, ~, ~] = fileparts(protocolPath);
                end

                % 设置路径
                obj.path = protocolPath;
                
                % 验证协议文件存在
                if ~isfile(obj.infoFile)
                    error('SEAL:ProtocolNode:NoProtocolFile', ...
                        'No such protocol in path: %s', protocolPath);
                end
                
                % 加载协议元数据
                obj.protocolInfo = ProtocolInfo.openExisting(obj.infoFile);
                obj.isLoaded = false;
                
                % 递归打开子节点
                obj.openChildNodes();
                
            catch ME
                error('SEAL:ProtocolNode:OpenFailed', ...
                    'Failed to open protocol: %s', ME.message);
            end
        end

        function load(obj)
            %LOAD 加载协议完整数据，用于后续操作
            
            if obj.isLoaded
                return;
            end
            
            try
                if ~isempty(obj.chanlocsNode) && ~obj.chanlocsNode.isLoaded
                    obj.chanlocsNode.load();
                end
                
                if ~isempty(obj.cortexNode) && ~obj.cortexNode.isLoaded
                    obj.cortexNode.load();
                end
                
                if ~isempty(obj.leadfieldNode) && ~obj.leadfieldNode.isLoaded
                    obj.leadfieldNode.load();
                end

                for i = 1:obj.childCount
                    if ~obj.children(i).isLoaded
                        obj.children(i).load();
                    end
                end
                
                obj.isLoaded = true;
                
            catch ME
                error('SEAL:ProtocolNode:LoadFailed', ...
                    'Failed to load protocol data: %s', ME.message);
            end
        end

        function unload(obj)
            %UNLOAD 卸载协议数据，释放内存
            
            try
                
                for i = 1:obj.childCount
                    if obj.children(i).isLoaded
                        obj.children(i).unload();
                    end
                end
                
                if ~isempty(obj.chanlocsNode) && obj.chanlocsNode.isLoaded
                    obj.chanlocsNode.unload();
                end
                
                if ~isempty(obj.cortexNode) && obj.cortexNode.isLoaded
                    obj.cortexNode.unload();
                end
                
                if ~isempty(obj.leadfieldNode) && obj.leadfieldNode.isLoaded
                    obj.leadfieldNode.unload();
                end
                
                obj.isLoaded = false;
                
            catch ME
                warning('SEAL:ProtocolNode:UnloadFailed', ...
                    'Failed to unload protocol data: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存协议数据
            try
                
                obj.createDirectoryStructure();
                % 递归保存所有子节点
                for i = 1:obj.childCount
                    obj.children(i).save();
                end
                
                % 保存单例节点
                if ~isempty(obj.chanlocsNode)
                    obj.chanlocsNode.save();
                end
                
                if ~isempty(obj.cortexNode)
                    obj.cortexNode.save();
                end
                
                if ~isempty(obj.leadfieldNode)
                    obj.leadfieldNode.save();
                end
                
                % 保存更新后的protocolInfo
                obj.protocolInfo.save(obj.infoFile);
                
            catch ME
                error('SEAL:ProtocolNode:SaveFailed', ...
                    'Failed to save protocol: %s', ME.message);
            end
        end

        function deleteFromDisk(obj)
            rmdir(obj.path, 's');
        end

        
        %% 子节点管理方法
        function addSession(obj, sessionNode)
            %ADDSESSION 添加会话节点
            obj.addChild(sessionNode);
        end
        
        function removeSession(obj, sessionNode)
            %REMOVESESSION 移除会话节点
            obj.removeChild(sessionNode);
        end
        
        function setChanlocsNode(obj, chanlocsNode)
            %SETCHANLOCSNODE 设置通道节点
            % 移除现有的通道节点
            obj.chanlocsNode = chanlocsNode;
            chanlocsNode.parent = obj;
        end
        
        function setCortexNode(obj, cortexNode)
            %SETCORTEXNODE 设置皮层节点
            % 移除现有的皮层节点
            obj.cortexNode = cortexNode;
            cortexNode.parent = obj;
        end
        
        function setLeadfieldNode(obj, leadfieldNode)
            %SETLEADFIELDNODE 设置导联场节点
            % 移除现有的导联场节点
            obj.leadfieldNode = leadfieldNode;
            leadfieldNode.parent = obj;
        end

        function node = openChanlocsFromData(obj, dataPath)
            chanlocsInfo = ChanlocsInfo.fromData(dataPath);
            obj.chanlocsNode = ChanlocsNode.fromInfo(chanlocsInfo);
            obj.chanlocsNode.path = obj.path;
            obj.chanlocsNode.parent = obj;
            obj.chanlocsNode.save();
            node = obj.chanlocsNode;
        end

        function node = openCortexFromData(obj, dataPath)
            cortexInfo = CortexInfo.fromData(dataPath);
            obj.cortexNode = CortexNode.fromInfo(cortexInfo);
            obj.cortexNode.path = obj.path;
            obj.cortexNode.parent = obj;
            obj.cortexNode.save();
            node = obj.cortexNode;
        end

       function node = openLeadfieldFromData(obj, dataPath)
            leadfieldInfo = LeadfieldInfo.fromData(dataPath);
            obj.leadfieldNode = LeadfieldNode.fromInfo(leadfieldInfo);
            obj.leadfieldNode.path = obj.path;
            obj.leadfieldNode.parent = obj;
            obj.leadfieldNode.save();
            node = obj.leadfieldNode;
       end

       function node = openDataFromData(obj, dataPath)
            sessionNode = SessionNode.createNew(findAvailableName(obj.sessionFolder, "Session #"), obj.sessionFolder);
            node = sessionNode.openDataFromData(dataPath);
            obj.addChild(sessionNode);
            sessionNode.save();
       end
        
        %% 依赖属性get方法
        function count = get.sessionCount(obj)
            % 统计会话节点数量
            count = 0;
            for i = 1:obj.childCount
                if isa(obj.children(i), 'SessionNode')
                    count = count + 1;
                end
            end
        end
        
        function name = get.name(obj)
            if ~isempty(obj.protocolInfo)
                name = obj.protocolInfo.name;
            else
                name = "Unnamed Protocol";
            end
        end

        function path = get.infoFile(obj)
            [~, name, ~] = fileparts(obj.path);
            path = fullfile(obj.path, strcat(name, ".mat"));
        end

        function path = get.sessionFolder(obj)
            path = fullfile(obj.path, "Sessions");
        end

        function path = get.chanlocsInfoFile(obj)
            path = fullfile(obj.path, "chanlocs_info.mat");
        end

        function path = get.cortexInfoFile(obj)
            path = fullfile(obj.path, "cortex_info.mat");
        end

        function path = get.leadfieldInfoFile(obj)
            path = fullfile(obj.path, "leadfield_info.mat");
        end
        
        function protocolType = get.protocolType(obj)
            if ~isempty(obj.protocolInfo)
                protocolType = obj.protocolInfo.type;
            else
                protocolType = "";
            end
        end
        
        function desc = get.desc(obj)
            if ~isempty(obj.protocolInfo)
                desc = obj.protocolInfo.desc;
            else
                desc = "";
            end
        end
        
        %% 重写父类方法
        function addChild(obj, childNode)
            %ADDCHILD 重写添加子节点方法，确保只添加允许的类型
            
            allowedTypes = {'SessionNode'};  % ProtocolNode只能添加SessionNode子节点
            if ~any(cellfun(@(type) isa(childNode, type), allowedTypes))
                error('SEAL:ProtocolNode:InvalidChildType', ...
                    'ProtocolNode只能添加SessionNode子节点');
            end
            
            % 调用父类方法
            addChild@BaseNode(obj, childNode);
        end
        
        function sessions = getAllSessions(obj)
            %GETALLSESSIONS 获取所有会话节点
            sessions = [];
            for i = 1:obj.childCount
                if isa(obj.children(i), 'SessionNode')
                    if isempty(sessions)
                        sessions = obj.children(i);
                    else
                        sessions(end+1) = obj.children(i);
                    end
                end
            end
        end
    end

    methods (Static)
        function protocol = openExisting(srcPath)
            %OPENEXISTING 打开现有协议
            % 输入:
            %   srcPath - 协议文件路径或协议文件夹路径
            % 输出:
            %   protocol - ProtocolNode对象
            
            % 创建协议节点
            protocol = ProtocolNode();
            
            % 打开协议
            protocol.open(srcPath);
        end
        
        function protocol = createNew(protocolName, path, varargin)
            %CREATENEW 创建新协议
            % 输入:
            %   protocolName - 协议名称
            %   path - 存储路径
            %   varargin - 可选参数对
            
            p = inputParser;
            addRequired(p, 'protocolName', Validators.isText);
            addRequired(p, 'path', Validators.isText);
            addParameter(p, 'Type', '', Validators.isText);
            addParameter(p, 'Description', '', Validators.isText);
            
            parse(p, protocolName, path, varargin{:});
            
            % 构建协议路径
            protocolPath = fullfile(path, protocolName);
            
            % 创建协议节点
            protocol = ProtocolNode();
            protocol.path = protocolPath;
            protocol.protocolInfo = ProtocolInfo.createNew(...
                p.Results.protocolName, p.Results.Type, p.Results.Description);
            
            % 保存协议基本信息
            protocol.save();
        end
    end

    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建协议目录结构
            
            if obj.path == ""
                return;
            end
            
            % 创建主目录
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
            
            % 创建会话目录
            if ~isfolder(obj.sessionFolder)
                mkdir(obj.sessionFolder);
            end
        end
        
        function openChildNodes(obj)
            %OPENCHILDNODES 打开子节点
            % 这个方法会在打开协议时自动扫描并打开现有的子节点
            
            % 打开单例节点
            obj.openSingletonNodes();
            
            % 打开会话节点
            obj.openSessionNodes();
        end
        
        function openSingletonNodes(obj)
            %OPENSINGLETONNODES 打开单例节点
            
            % 打开通道节点
            if ~isempty(obj.chanlocsInfoFile) && isfile(obj.chanlocsInfoFile)
                try
                    obj.chanlocsNode = ChanlocsNode.openExisting(obj.chanlocsInfoFile);
                    obj.chanlocsNode.parent = obj;
                catch ME
                    warning('SEAL:ProtocolNode:OpenChanlocsFailed', ...
                        'Failed to open chanlocs node: %s. Error: %s', ...
                        obj.chanlocsInfoFile, ME.message);
                end
            end
            
            % 打开皮层节点
            if ~isempty(obj.cortexInfoFile) && isfile(obj.cortexInfoFile)
                try
                    obj.cortexNode = CortexNode.openExisting(obj.cortexInfoFile);
                    obj.cortexNode.parent = obj;
                catch ME
                    warning('SEAL:ProtocolNode:OpenCortexFailed', ...
                        'Failed to open cortex node: %s. Error: %s', ...
                        obj.cortexInfoFile, ME.message);
                end
            end
            
            % 打开导联场节点
            if ~isempty(obj.leadfieldInfoFile) && isfile(obj.leadfieldInfoFile)
                try
                    fprintf('  Opening leadfield node from: %s\n', obj.leadfieldInfoFile);
                    obj.leadfieldNode = LeadfieldNode.openExisting(obj.leadfieldInfoFile);
                    obj.leadfieldNode.parent = obj;
                catch ME
                    warning('SEAL:ProtocolNode:OpenLeadfieldFailed', ...
                        'Failed to open leadfield node: %s. Error: %s', ...
                        obj.leadfieldInfoFile, ME.message);
                end
            end
        end
        
        function openSessionNodes(obj)
            %OPENSESSIONNODES 打开会话节点
            
            % 检查会话目录是否存在
            if ~isfolder(obj.sessionFolder)
                return;
            end
            
            % 扫描会话目录
            sessionDirs = dir(obj.sessionFolder);
            sessionCount = 0;
            
            for i = 1:length(sessionDirs)
                sessionDir = sessionDirs(i);
                
                % 跳过 . 和 .. 目录以及非目录项
                if ~sessionDir.isdir || startsWith(sessionDir.name, '.')
                    continue;
                end
                
                sessionFile = fullfile(obj.sessionFolder, sessionDir.name, "session_info.mat");
                
                % 检查是否存在会话文件
                if ~isfile(sessionFile)
                    continue;
                end
                
                try
                    if exist('SessionNode', 'class')
                        sessionNode = SessionNode.openExisting(sessionFile);
                        sessionNode.parent = obj;
                        obj.addChild(sessionNode);
                        sessionCount = sessionCount + 1;
                    end
                    
                catch ME
                    warning('SEAL:ProtocolNode:OpenSessionFailed', ...
                        'Failed to open session: %s. Error: %s', ...
                        sessionPath, ME.message);
                end
            end
        end
        
        function node = getSingletonNode(obj, nodeType)
            %GETSINGLETONNODE 获取单例类型节点
            node = [];
            switch nodeType
                case 'ChanlocsNode'
                    node = obj.chanlocsNode;
                case 'CortexNode'
                    node = obj.cortexNode;
                case 'LeadfieldNode'
                    node = obj.leadfieldNode;
            end
        end
    end
end