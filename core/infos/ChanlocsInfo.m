classdef ChanlocsInfo < handle
    %CHANLOCSINFO 通道信息管理类
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
        chanlocsCount int32     % 通道数量
    end
    
    methods
        function obj = ChanlocsInfo(name, dataPath, desc, metadata)
            %CHANLOCSINFO 构造函数
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
        function chanlocs = openExisting(chanlocsPath)
            %OPENEXISTING 打开现有通道信息
            % 输入:
            %   chanlocsPath - 通道文件路径
            % 输出:
            %   chanlocs - 加载的ChanlocsInfo对象

            chanlocsFile = chanlocsPath;
            
            if ~isfile(chanlocsFile)
                error('SEAL:ChanlocsInfo:FileNotFound', ...
                    'No such file exists: %s', chanlocsFile);
            end
            
            % 加载通道数据
            chanlocs = loadData(chanlocsFile);
            
            % 验证加载的对象类型
            if ~isa(chanlocs, 'ChanlocsInfo')
                error('SEAL:ChanlocsInfo:InvalidFileFormat', ...
                    'File does not contain a valid ChanlocsInfo object: %s', chanlocsFile);
            end
        end

        function chanlocs = fromData(dataPath, varargin)
            %FROMDATA 从通道数据文件直接创建ChanlocsInfo
            % 输入:
            %   dataPath - 通道数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 通道名称（可选，默认为文件名）
            %     'Description' - 通道描述
            %     'Metadata' - 额外元数据
            % 输出:
            %   chanlocs - ChanlocsInfo对象
            
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
                chanlocsName = name;
            else
                chanlocsName = p.Results.Name;
            end
            
            % 从数据文件中提取元数据
            metadata = p.Results.Metadata;
            
            % 创建ChanlocsInfo对象
            chanlocs = ChanlocsInfo(...
                chanlocsName, ...
                dataPath, ...
                p.Results.Description, ...
                metadata);
        end

        function meta = extractMetadata(dataPath)
        end
    end
end