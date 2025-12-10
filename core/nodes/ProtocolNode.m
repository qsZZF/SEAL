classdef ProtocolNode < BaseNode
    %PROTOCOLNODE 协议节点
    % 管理研究协议，包含会话、皮层模型、导联场等
    
    properties (Constant)
        type = "ProtocolNode"  % 节点类型标识符
    end
    
    properties
        % 协议特定属性
        protocolInfo ProtocolInfo
        channelNode ChannelNode           % 通道节点
%         cortexNode CortexNode             % 皮层节点
%         leadfieldNode LeadfieldNode       % 导联场节点
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        protocolType string
        desc string
        sessionCount int32                 % 会话数量
    end
    
    methods
        function obj = ProtocolNode()
            %PROTOCOLNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, protocolPath)
            %OPEN 打开协议
            if obj.isLoaded
                return;
            end
            
            try
                % 检查协议文件是否存在
                if isfolder(protocolPath)
                    % 如果是文件夹，查找项目文件
                    [~, projectName, ~] = fileparts(protocolPath);
                    protocolFiles = dir(fullfile(protocolPath, strcat(projectName, '.mat')));
                    if isempty(protocolFiles)
                        error('SEAL:ProtocolNode:NoProtocolFile', ...
                            'No such protocol in path: %s', protocolPath);
                    end
                else
                     [protocolPath, projectName, ~] = fileparts(protocolPath);
                end
                
                protocolFile = fullfile(protocolPath, strcat(projectName, ".mat"));

                % 加载协议元数据
                obj.protocolInfo = ProtocolInfo.openExisting(protocolFile);
                
                % 标记为已加载
                obj.isLoaded = true;
                
                % 加载子节点
                obj.loadChildNodes();
                
            catch ME
                error('SEAL:ProtocolNode:LoadFailed', ...
                    'Failed to open protocol: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存协议数据
            obj.createDirectoryStructure();
            obj.protocolInfo.save(obj.path);
            
            % 递归保存所有子节点
            for i = 1:obj.childCount
                obj.children(i).save();
            end
        end

        function load(obj)
            %LOAD 加载协议数据
            if obj.isLoaded
                return;
            end
            
            % 加载协议元数据
            if isempty(obj.protocolInfo)
                [~, name, ~] = fileparts(obj.path);
                protocolFile = fullfile(obj.path, strcat(name, '.mat'));
                if isfile(protocolFile)
                    obj.protocolInfo = ProtocolInfo.openExisting(protocolFile);
                end
            end
            
            % 递归加载所有子节点
            for i = 1:obj.childCount
                obj.children(i).load();
            end
            
            obj.isLoaded = true;
        end

        function unload(obj)
            %UNLOAD 卸载协议数据，释放内存
            % 递归卸载所有子节点
            for i = 1:obj.childCount
                obj.children(i).unload();
            end
            
            obj.isLoaded = false;
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
        
        function setChannelNode(obj, channelNode)
            %SETCHANNELNODE 设置通道节点
            % 移除现有的通道节点
            obj.channelNode = channelNode;
            obj.protocolInfo.channelPath = obj.channelNode.path;
            channelNode.parent = obj;
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
            
            allowedprotocolTypes = {'SessionNode', 'ChannelNode', 'CortexNode', 'LeadfieldNode'};
            if ~any(cellfun(@(type) isa(childNode, type), allowedTypes))
                error('SEAL:ProtocolNode:InvalidChildType', ...
                    'ProtocolNode只能添加SessionNode、ChannelNode、CortexNode或LeadfieldNode子节点');
            end
            
            % 对于单例节点类型，确保只有一个实例
            singletonTypes = {'ChannelNode', 'CortexNode', 'LeadfieldNode'};
            childType = class(childNode);
            if any(strcmp(childType, singletonTypes))
                existingNode = obj.getSingletonNode(childType);
                if ~isempty(existingNode)
                    error('SEAL:ProtocolNode:DuplicateSingleton', ...
                        'ProtocolNode只能有一个%s实例', childType);
                end
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
            %   parentProject - 父项目节点
            
            if isfolder(srcPath)
                % 如果是文件夹，查找协议文件
                protocolFiles = dir(fullfile(srcPath, '*.mat'));
                if isempty(protocolFiles)
                    error('SEAL:ProtocolNode:NoProtocolFile', ...
                        'No such protocol in path: %s', srcPath);
                end
            end
            
            % 创建协议节点
            protocol = ProtocolNode();

            % 加载协议数据
            protocol.open(srcPath);
        end
        
        function protocol = createNew(protocolName, path, varargin)
            %CREATENEW 创建新协议
            % 输入:
            %   protocolName - 协议名称
            %   parentProject - 父项目节点
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
            
            % 保存协议
            protocol.save();
            protocol.load();
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
            sessionsPath = fullfile(obj.path, 'Sessions');
            if ~isfolder(sessionsPath)
                mkdir(sessionsPath);
            end
        end
        
        function loadChildNodes(obj)
            %LOADCHILDNODES 加载子节点
            % 这个方法会在打开协议时自动扫描并加载现有的子节点
            
            % 加载会话节点
            sessionsPath = fullfile(obj.path, 'Sessions');
            if isfolder(sessionsPath)
                sessionDirs = dir(sessionsPath);
                for i = 1:length(sessionDirs)
                    if sessionDirs(i).isdir && ~startsWith(sessionDirs(i).name, '.')
                        sessionPath = fullfile(sessionsPath, sessionDirs(i).name);
                        try
                            sessionNode = SessionNode.openSession(sessionPath, obj);
                            obj.addChild(sessionNode);
                        catch ME
                            warning('SEAL:ProtocolNode:LoadSessionFailed', ...
                                'Failed to load session: %s. Error: %s', ...
                                sessionPath, ME.message);
                        end
                    end
                end
            end
            
            % 加载其他单例节点（通道、皮层、导联场）
            % 这里需要根据实际文件结构来实现
            % 暂时留空，待相关节点类实现后补充
        end
        
        function node = getSingletonNode(obj, nodeType)
            %GETSINGLETONNODE 获取单例类型节点
            node = [];
            for i = 1:obj.childCount
                if isa(obj.children(i), nodeType)
                    node = obj.children(i);
                    return;
                end
            end
        end
    end
end
