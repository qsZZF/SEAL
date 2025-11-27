classdef ChannelInfo < handle
    %CHANNELINFO 通道信息管理类
    % 负责通道数据的元数据管理和存储

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        
        % 通道数据
        dataPath string            % 通道数据存储
        data struct
        
        % 元数据
        metadata struct        % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime
    end
    
    properties (Dependent)
        channelCount int32     % 通道数量
    end
    
    methods
        function obj = ChannelInfo(name, dataPath, desc, metadata)
            %CHANNELINFO 构造函数
            % 输入:
            %   name - 通道集合名称
            %   dataPath - 原始数据路径
            %   desc - 描述
            %   metadata - 元数据

            if nargin >= 1
                obj.name = name;
            end
            if nargin >= 2
                obj.names = channelNames;
            else
                obj.names = {};
            end
            if nargin >= 3
                obj.desc = desc;
            end
            if nargin >= 4
                obj.metadata = metadata;
            else
                obj.metadata = struct();
            end
            
            % 初始化时间戳
            obj.createdDate = datetime('now');
            obj.modifiedDate = obj.createdDate;
            obj.data = struct();
            obj.dataPath = dataPath;
        end
        
        %% 文件操作方法
        function save(obj, path)
            %SAVE 保存通道信息
            if path == ""
                error('SEAL:ChannelInfo:InvalidPath', ...
                    'No valid channel path specified.');
            end
            
            if ~isfolder(path)
                mkdir(path);
            end
            
            % 保存通道元数据
            obj.updateModifiedDate();
            save(fullfile(path, strcat(obj.name, '.mat')), 'obj');
        end
    end
    
    methods (Access = private)
        function updateModifiedDate(obj)
            %UPDATEMODIFIEDDATE 更新修改时间
            obj.modifiedDate = datetime('now');
        end
    end
    
    %% 静态方法
    methods (Static)
        function channel = createNew(name, dataPath, desc, metadata)
            %CREATENEW 创建新通道信息
            % 输入:
            %   name - 通道集合名称
            %   dataPath - 通道名称列表
            %   desc - 描述
            %   metadata - 元数据
            % 输出:
            %   channel - 新创建的ChannelInfo对象
            channel = ChannelInfo(name, dataPath, desc, metadata);
        end
        
        function channel = openExisting(channelPath)
            %OPENEXISTING 打开现有通道信息
            % 输入:
            %   channelPath - 通道文件路径
            % 输出:
            %   channel - 加载的ChannelInfo对象

            channelFile = channelPath;
            
            if ~isfile(channelFile)
                error('SEAL:ChannelInfo:FileNotFound', ...
                    'No such file exists: %s', channelFile);
            end
            
            % 加载通道数据
            channel = load(channelFile, 'obj').obj;
            
            % 验证加载的对象类型
            if ~isa(channel, 'ChannelInfo')
                error('SEAL:ChannelInfo:InvalidFileFormat', ...
                    'File does not contain a valid ChannelInfo object: %s', channelFile);
            end
        end
    end
end