classdef CortexInfo < handle
    %CORTEXINFO 皮层信息管理类
    % 负责皮层数据的元数据管理和存储

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        
        % 皮层数据地址
        dataPath string            % 皮层数据存储
        
        % 元数据
        metadata struct        % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime
    end
    
    properties (Dependent)
        cortexCount int32     % 皮层数量
    end
    
    methods
        function obj = CortexInfo(name, dataPath, desc, metadata)
            %CORTEXINFO 构造函数
            % 输入:
            %   name - 皮层集合名称
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
            %SAVE 保存皮层信息
            if path == ""
                error('SEAL:CortexInfo:InvalidPath', ...
                    'No valid cortex path specified.');
            end
            
            % 保存皮层元数据
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
        function cortex = openExisting(cortexPath)
            %OPENEXISTING 打开现有皮层信息
            % 输入:
            %   cortexPath - 皮层文件路径
            % 输出:
            %   cortex - 加载的CortexInfo对象

            cortexFile = cortexPath;
            
            if ~isfile(cortexFile)
                error('SEAL:CortexInfo:FileNotFound', ...
                    'No such file exists: %s', cortexFile);
            end
            
            % 加载皮层数据
            cortex = loadData(cortexFile);
            
            % 验证加载的对象类型
            if ~isa(cortex, 'CortexInfo')
                error('SEAL:CortexInfo:InvalidFileFormat', ...
                    'File does not contain a valid CortexInfo object: %s', cortexFile);
            end
        end

        function cortex = fromData(dataPath, varargin)
            %FROMDATA 从皮层数据文件直接创建CortexInfo
            % 输入:
            %   dataPath - 皮层数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 皮层名称（可选，默认为文件名）
            %     'Description' - 皮层描述
            %     'Metadata' - 额外元数据
            % 输出:
            %   cortex - CortexInfo对象
            
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
                cortexName = name;
            else
                cortexName = p.Results.Name;
            end
            
            % 从数据文件中提取元数据
            metadata = p.Results.Metadata;
%             if p.Results.AutoExtract
%                 metadata = CortexInfo.extractMetadata(dataPath, metadata);
%             end
            
            % 创建CortexInfo对象
            cortex = CortexInfo(...
                cortexName, ...
                dataPath, ...
                p.Results.Description, ...
                metadata);
        end

        function meta = extractMetadata(dataPath)
        end
    end
end