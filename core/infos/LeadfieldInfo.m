classdef LeadfieldInfo < handle
    %LEADFIELDINFO 导联场信息管理类
    % 负责导联场数据的元数据管理和存储

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        headModelType string
        orientation = []
        normals = [] 
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
                
                % 1. 读取文件
                % 1. 读取文件
                loadedData = loadData(dataPath);
                
                matchedGainField = ''; % 初始化为空
                
                % ==========================================
                % 2. 类型防御 与 极简模式/拆盒机制
                % ==========================================
                if isnumeric(loadedData)
                    % 情况 A：直接读出来就是一个纯数字矩阵 (double)
                    % 强行给它穿件马甲，伪装成标准结构体，名字就叫 'Gain_Matrix'
                    loadedData = struct('Gain_Matrix', loadedData);
                    matchedGainField = 'Gain_Matrix';
                    availableFields = fieldnames(loadedData);
                    
                elseif isstruct(loadedData)
                    % 情况 B：读出来是正常的结构体 (mat文件的默认行为)
                    outerFields = fieldnames(loadedData);
                    
                    if length(outerFields) == 1
                        loneVar = loadedData.(outerFields{1});
                        
                        if isstruct(loneVar)
                            % 发现俄罗斯套娃！剥开它
                            loadedData = loneVar;
                            availableFields = fieldnames(loadedData);
                        elseif isnumeric(loneVar)
                            % 结构体里只有唯一一个变量，且是数字矩阵！直接认领
                            matchedGainField = outerFields{1};
                            availableFields = outerFields;
                        else
                            availableFields = outerFields;
                        end
                    else
                        % 结构体里有多个变量，老老实实进入后面的字典匹配
                        availableFields = outerFields;
                    end
                else
                    error('不支持的数据类型: 读入的数据既不是结构体也不是数值矩阵。');
                end
                
                % ==========================================
                % 3. 刚需防线：寻找导联场矩阵 (Gain)
                % ==========================================
                if isempty(matchedGainField)
                    gainAliases = {'Gain', 'L', 'Leadfield', 'LF', 'Forward', 'G', 'matrix'};
                    for i = 1:length(gainAliases)
                        if ismember(gainAliases{i}, availableFields)
                            matchedGainField = gainAliases{i};
                            break;
                        end
                    end
                end
                

                % leadfield.size = size(loadedData.(matchedGainField));
                
                % ==========================================
                % 3. 选填项宽容读取 (非刚需，绝不报错)
                % ==========================================
                if isfield(loadedData, 'HeadModelType')
                    leadfield.headModelType = string(loadedData.HeadModelType);
                else
                    leadfield.headModelType = "Unknown"; % 给个默认值
                end
                
                % 注意你之前的代码用的是双引号字符串 "Orientation"，这里统一改为单引号查字段名更规范
                if isfield(loadedData, 'Orientation')
                    leadfield.orientation = loadedData.Orientation;
                else
                    leadfield.orientation = []; % 没有就不强求
                end

                normalAliases = {'GridOrient', 'Normals', 'Orientations', 'CorticalNormals', 'Norm'};
                matchedNormalField = '';
                for i = 1:length(normalAliases)
                    if ismember(normalAliases{i}, availableFields)
                        matchedNormalField = normalAliases{i};
                        break;
                    end
                end

                if ~isempty(matchedNormalField)
                    leadfield.normals = loadedData.(matchedNormalField);
                else
                    leadfield.normals = [];
                end

            catch ME
                % 精准抛出底层错误
                error('SEAL:LeadfieldInfo:ImportFailed', ...
                    '导联场数据提取失败: %s\n详细原因: %s', dataPath, ME.message);
            end
        end

        function meta = extractMetadata(dataPath)
        end
    end
end