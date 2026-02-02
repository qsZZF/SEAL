classdef DataInfo < handle
    %DATAINFO 信息管理类
    % 负责数据的元数据管理和存储

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        
        % 数据地址
        dataPath string
        resultPath string
        
        % 元数据
        metadata struct        % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime
    end
    
    properties (Dependent)
        dataCount int32     % 数量
    end
    
    methods
        function obj = DataInfo(name, dataPath, resultPath, desc, metadata)
            %DATAINFO 构造函数
            % 输入:
            %   name - 集合名称
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
                obj.resultPath = resultPath;
            else
                obj.resultPath = {};
            end
            if nargin >= 4
                obj.desc = desc;
            end
            if nargin >= 5
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
            %SAVE 保存信息

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
        function data = openExisting(dataPath)
            %OPENEXISTING 打开现有信息
            % 输入:
            %   dataPath - 文件路径
            % 输出:
            %   data - 加载的DataInfo对象

            dataFile = dataPath;
            
            if ~isfile(dataFile)
                error('SEAL:DataInfo:FileNotFound', ...
                    'No such file exists: %s', dataFile);
            end
            
            % 加载数据
            data = loadData(dataFile);
            
            % 验证加载的对象类型
            if ~isa(data, 'DataInfo')
                error('SEAL:DataInfo:InvalidFileFormat', ...
                    'File does not contain a valid DataInfo object: %s', dataFile);
            end
        end

        function data = fromData(dataPath, varargin)
            %FROMDATA 从数据文件直接创建DataInfo
            % 输入:
            %   dataPath - 数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 名称（可选，默认为文件名）
            %     'Description' - 描述
            %     'Metadata' - 额外元数据
            % 输出:
            %   data - DataInfo对象
            
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
                dataName = name;
            else
                dataName = p.Results.Name;
            end
            
            % 从数据文件中提取元数据
            metadata = p.Results.Metadata;
            
            % 创建DataInfo对象
            data = DataInfo(...
                dataName, ...
                dataPath, ...
                "", ...
                p.Results.Description, ...
                metadata);
        end

        function meta = extractMetadata(dataPath)
        end
    end
end