classdef ChannelNode < BaseNode
    %CHANNELNODE 通道节点
    % 管理协议的通道信息，包含通道名称、位置等元数据
    
    properties (Constant)
        type = "ChannelNode"  % 节点类型标识符
    end
    
    properties
        % 通道特定属性
        channelInfo ChannelInfo
        cache
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        data
    end
    
    methods
        function obj = ChannelNode()
            %CHANNELNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, channelPath)
            %OPEN 打开通道数据
            obj.path = channelPath;
            if obj.isLoaded
                return;
            end
            
            try
                obj.channelInfo = ChannelInfo.openExisting(channelPath);
                
            catch ME
                error('SEAL:ChannelNode:LoadFailed', ...
                    'Failed to open channel data: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存通道数据
            obj.createDirectoryStructure();
            obj.channelInfo.save(obj.path);
        end

        function load(obj)
            %LOAD 加载通道数据
            if obj.isLoaded
                return;
            end
            
            % 加载通道元数据
            if ~isempty(obj.channelInfo)
                [~, name, ~] = fileparts(obj.channelInfo.dataPath);
                if isfile(obj.channelInfo.dataPath)
                    obj.cache = loadData(obj.channelInfo.dataPath);
                    obj.isLoaded = true;
                end
            end
        end

        function unload(obj)
            %UNLOAD 卸载通道数据，释放内存
            obj.cache = [];
            obj.isLoaded = false;
        end
        
        %% 依赖属性get方法
        function channelName = get.name(obj)
            channelName = obj.channelInfo.name;
        end
        
        function data = get.data(obj)
            if ~obj.isLoaded
                obj.load();
            end
            data = obj.cache;
        end
    
    end

    methods (Static)
        function channel = openExisting(srcPath)
            %OPENCHANNEL 打开现有通道数据
            % 输入:
            %   srcPath - 通道文件路径或通道文件夹路径
            %   parentProtocol - 父协议节点
            
            % 创建通道节点
            channel = ChannelNode();

            % 加载通道数据
            channel.open(getmat(srcPath));
        end

        function channel = fromData(dataPath, varargin)
            %FROMDATA 从通道数据文件直接创建ChannelNode
            % 输入:
            %   dataPath - 通道数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 通道名称（可选，默认为文件名）
            %     'Description' - 通道描述（可选）
            %     'Metadata' - 额外元数据（可选）
            %     'NodePath' - 节点保存路径（可选）
            % 输出:
            %   channel - ChannelNode对象
            
            % 参数解析
            p = inputParser;
            addRequired(p, 'dataPath', @(x) isfile(x) && endsWith(x, '.mat'));
            addParameter(p, 'Name', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Description', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Metadata', struct(), @isstruct);
            addParameter(p, 'NodePath', '', @(x) ischar(x) || isstring(x));
            
            parse(p, dataPath, varargin{:});
            
            % 如果未提供名称，使用文件名
            if isempty(p.Results.Name)
                [~, name, ~] = fileparts(dataPath);
                channelName = name;
            else
                channelName = p.Results.Name;
            end
            
            % 创建ChannelInfo
            channelInfo = ChannelInfo.fromData(...
                dataPath, ...
                'Name', channelName, ...
                'Description', p.Results.Description, ...
                'Metadata', p.Results.Metadata);
            
            % 创建ChannelNode
            channel = ChannelNode();
            channel.channelInfo = channelInfo;
            
            % 设置节点路径
            if ~isempty(p.Results.NodePath)
                channel.path = p.Results.NodePath;
            else
                % 如果没有提供节点路径，使用数据文件所在目录
                [dataDir, ~, ~] = fileparts(dataPath);
                channel.path = fullfile(dataDir, channelName);
            end
            
            % 尝试加载数据
            try
                channel.load();
            catch ME
                warning('SEAL:ChannelNode:AutoLoadFailed', ...
                    '自动加载数据失败: %s', ME.message);
            end
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