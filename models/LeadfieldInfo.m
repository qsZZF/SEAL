classdef LeadfieldInfo < handle
    %LEADFIELDINFO 导联场信息管理类
    % 负责导联场数据的元数据管理和存储

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        headModelType string
        orientation = []
        
        % 导联场数据地址
        dataPath string            % 导联场数据存储
        
        % 元数据
        metadata struct        % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime
    end
    
    properties (Dependent)
        leadfieldCount int32     % 导联场数量
    end
    
    methods
        function obj = LeadfieldInfo(name, dataPath, desc, metadata)
            %LEADFIELDINFO 构造函数
            % 输入:
            %   name - 导联场集合名称
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
            %SAVE 保存导联场信息

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
        function leadfield = openExisting(leadfieldPath)
            %OPENEXISTING 打开现有导联场信息
            % 输入:
            %   leadfieldPath - 导联场文件路径
            % 输出:
            %   leadfield - 加载的LeadfieldInfo对象

            leadfieldFile = leadfieldPath;
            
            if ~isfile(leadfieldFile)
                error('SEAL:LeadfieldInfo:FileNotFound', ...
                    'No such file exists: %s', leadfieldFile);
            end
            
            % 加载导联场数据
            leadfield = loadData(leadfieldFile);
            
            % 验证加载的对象类型
            if ~isa(leadfield, 'LeadfieldInfo')
                error('SEAL:LeadfieldInfo:InvalidFileFormat', ...
                    'File does not contain a valid LeadfieldInfo object: %s', leadfieldFile);
            end
        end

        function leadfield = fromData(dataPath, varargin)
            %FROMDATA 从导联场数据文件直接创建LeadfieldInfo
            % 输入:
            %   dataPath - 导联场数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 导联场名称（可选，默认为文件名）
            %     'Description' - 导联场描述
            %     'Metadata' - 额外元数据
            % 输出:
            %   leadfield - LeadfieldInfo对象
            
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
                leadfieldName = name;
            else
                leadfieldName = p.Results.Name;
            end
            
            % 从数据文件中提取元数据
            metadata = p.Results.Metadata;
%             if p.Results.AutoExtract
%                 metadata = LeadfieldInfo.extractMetadata(dataPath, metadata);
%             end
            
            % 创建LeadfieldInfo对象
            try
                leadfield = LeadfieldInfo(...
                    leadfieldName, ...
                    dataPath, ...
                    p.Results.Description, ...
                    metadata);

                loadedData = loadData(dataPath);

                leadfield.headModelType = string(loadedData.HeadModelType);
                if isfield(loadedData, "Orientation")
                    leadfield.orientation = loadedData.Orientation;
                else
                    leadfield.orientation = [];
                end

            catch ME
               error('SEAL:LeadfieldInfo:FileCorrupted', ...
                    'File corrupted: %s', dataPath);
            end

        end

        function meta = extractMetadata(dataPath)
        end
    end
end