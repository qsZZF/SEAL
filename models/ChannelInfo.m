classdef ChannelInfo < handle
    %CHANNELINFO 通道信息管理类
    % 负责通道数据的元数据管理和存储

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        
        % 通道数据地址
        dataPath string            % 通道数据存储
        
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
                obj.dataPath = dataPath;
            else
                obj.dataPath = {};
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
            obj.dataPath = dataPath;
        end
        
        %% 文件操作方法
        function save(obj, path)
            %SAVE 保存通道信息

            checkDir(path);
            obj.updateModifiedDate();
            saveData(path, obj);
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
            channel = loadData(channelFile);
            
            % 验证加载的对象类型
            if ~isa(channel, 'ChannelInfo')
                error('SEAL:ChannelInfo:InvalidFileFormat', ...
                    'File does not contain a valid ChannelInfo object: %s', channelFile);
            end
        end

        function channel = fromData(dataPath, varargin)
            %FROMDATA 从通道数据文件直接创建ChannelInfo
            % 输入:
            %   dataPath - 通道数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 通道名称（可选，默认为文件名）
            %     'Description' - 通道描述
            %     'Metadata' - 额外元数据
            % 输出:
            %   channel - ChannelInfo对象
            
            % 参数解析
            p = inputParser;
            addRequired(p, 'dataPath', @(x) isfile(x) && endsWith(x, '.mat'));
            addParameter(p, 'Name', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Description', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Metadata', struct(), @isstruct);
            addParameter(p, 'AutoExtract', true, @islogical);
            
            parse(p, dataPath, varargin{:});
            
            % 如果未提供名称，使用文件名
            if isempty(p.Results.Name)
                [~, name, ~] = fileparts(dataPath);
                channelName = name;
            else
                channelName = p.Results.Name;
            end
            
            % 从数据文件中提取元数据
            metadata = p.Results.Metadata;
            
            % 创建ChannelInfo对象
            channel = ChannelInfo(...
                channelName, ...
                dataPath, ...
                p.Results.Description, ...
                metadata);
        end

        function meta = extractMetadata(dataPath)
        end
    end
end