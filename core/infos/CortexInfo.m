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

            % 从数据文件中提取元数据前缀
            metadata = p.Results.Metadata;
            if ~isstruct(metadata)
                metadata = struct();
            end

            % ==========================================
            % 1. 海关安检：读取并拆盒
            % ==========================================
            try
                loadedData = loadData(dataPath);
                if isstruct(loadedData)
                    outerFields = fieldnames(loadedData);
                    if length(outerFields) == 1 && isstruct(loadedData.(outerFields{1}))
                        loadedData = loadedData.(outerFields{1});
                    end
                    availableFields = fieldnames(loadedData);
                else
                    error('皮层模型必须是结构体或 .mat 文件！');
                end

                % ==========================================
                % 2. 字段智能匹配 (查字典)
                % ==========================================
                vertAliases = {'Vertices', 'pos', 'Position', 'vertices', 'nodes'};
                facesAliases = {'Faces', 'tri', 'Triangles', 'faces', 'polygons'};
                atlasAliases = {'Atlas', 'Scout', 'atlas', 'regions', 'parcellation'};
                connAliases = {'VertConn'};
                matchedVert = ''; matchedFaces = ''; matchedAtlas = '';matchedConn='';

                for i = 1:length(vertAliases), if ismember(vertAliases{i}, availableFields), matchedVert = vertAliases{i}; break; end, end
                for i = 1:length(facesAliases), if ismember(facesAliases{i}, availableFields), matchedFaces = facesAliases{i}; break; end, end
                for i = 1:length(atlasAliases), if ismember(atlasAliases{i}, availableFields), matchedAtlas = atlasAliases{i}; break; end, end
                for i = 1:length(connAliases), if ismember(connAliases{i}, availableFields), matchedConn = connAliases{i}; break; end, end
                % ==========================================
                % 3. 刚需防线：人工兜底 (只在导入这一步问一次！)
                % ==========================================
                if isempty(matchedVert)
                    [indx, tf] = listdlg('ListString', availableFields, 'SelectionMode', 'single', ...
                        'PromptString', {'未识别到【顶点坐标 (Vertices)】', '请指定:'}, 'Name', '皮层匹配 (1/4)');
                    if tf, matchedVert = availableFields{indx}; else, error('用户取消导入。'); end
                end

                if isempty(matchedFaces)
                    [indx, tf] = listdlg('ListString', availableFields, 'SelectionMode', 'single', ...
                        'PromptString', {'未识别到【拓扑连接/面 (Faces)】', '请指定:'}, 'Name', '皮层匹配 (2/4)');
                    if tf, matchedFaces = availableFields{indx}; else, error('用户取消导入。'); end
                end

                if isempty(matchedAtlas)
                    [indx, tf] = listdlg('ListString', availableFields, 'SelectionMode', 'single', ...
                        'PromptString', {'未识别到【脑区图谱 (Atlas)】', '请指定 (若无图谱可直接取消):'}, 'Name', '皮层匹配 (3/4)');
                    if tf, matchedAtlas = availableFields{indx}; end
                end

                if isempty(matchedConn)
                    [indx, tf] = listdlg('ListString', availableFields, 'SelectionMode', 'single', ...
                        'PromptString', {'未识别到【链接（VertConn）】', '请指定 :'}, 'Name', '皮层匹配 (4/4)');
                    if tf, matchedConn = availableFields{indx}; end
                end
                % ==========================================
                % 4. 提取轻量级信息 & 登记暗号
                % ==========================================
                % 数一数点和面的数量，供 Metadata UI 面板装X用
                metadata.VertexCount = size(loadedData.(matchedVert), 1);
                metadata.FaceCount = size(loadedData.(matchedFaces), 1);
                metadata.HasAtlas = ~isempty(matchedAtlas);
                
                % 【精髓】：把千辛万苦匹配到的真名，作为“暗号”写进签证里！
                metadata.Key_Vertices = matchedVert;
                metadata.Key_VertConn = matchedFaces;
                metadata.Key_Atlas = matchedAtlas;
                metadata.Key_VertConn = matchedConn;
                % ==========================================
                % 5. 准许入境：创建对象
                % ==========================================
                cortex = CortexInfo(...
                    cortexName, ...
                    dataPath, ...
                    p.Results.Description, ...
                    metadata);

            catch ME
                error('SEAL:CortexInfo:ImportFailed', ...
                    '皮层数据导入失败: %s\n详细原因: %s', dataPath, ME.message);
            end
        end

        function meta = extractMetadata(dataPath)
        end
    end
end