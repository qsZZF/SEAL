classdef ChannelNode < BaseNode
    %CHANNELNODE 通道节点
    % 管理协议的通道信息，包含通道名称、位置等元数据
    
    properties (Constant)
        type = "ChannelNode"  % 节点类型标识符
    end
    
    properties
        % 通道特定属性
        channel
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        metadata struct            % 元数据
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
                obj.load(obj.path);
                
            catch ME
                error('SEAL:ChannelNode:LoadFailed', ...
                    'Failed to open channel data: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存通道数据
            obj.createDirectoryStructure();
            save(obj.path, obj.channel, '-mat');
        end

        function load(obj)
            %LOAD 加载通道数据
            if obj.isLoaded
                return;
            end
            
            % 加载通道元数据
            if isempty(obj.channel)
                [~, name, ~] = fileparts(obj.path);
                channelFile = fullfile(obj.path, strcat(name, '.mat'));
                if isfile(channelFile)
                    obj.channel = open(channelFile);
                end
            end
            
            obj.isLoaded = true;
        end

        function unload(obj)
            %UNLOAD 卸载通道数据，释放内存
            clear(obj.channel);
            obj.isLoaded = false;
        end
        
        %% 依赖属性get方法
        function channelName = get.name(obj)
            channelName = "Channel";
        end
        
        function meta = get.metadata(obj)
            if ~isempty(obj.channel)
                meta = obj.channel.metadata;
            else
                meta = struct();
            end
        end
    end

    methods (Static)
        function channel = openExisting(srcPath)
            %OPENCHANNEL 打开现有通道数据
            % 输入:
            %   srcPath - 通道文件路径或通道文件夹路径
            %   parentProtocol - 父协议节点
            
            if isfolder(srcPath)
                % 如果是文件夹，查找通道文件
                channelFiles = dir(fullfile(srcPath, '*.mat'));
                if isempty(channelFiles)
                    error('SEAL:ChannelNode:NoChannelFile', ...
                        'No such channel data in path: %s', srcPath);
                end
            end
            
            % 创建通道节点
            channel = ChannelNode();

            % 加载通道数据
            channel.open(srcPath);
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
        
        function updateModifiedDate(obj)
            %UPDATEMODIFIEDDATE 更新修改时间
            obj.modifiedDate = datetime('now');
        end
    end
end