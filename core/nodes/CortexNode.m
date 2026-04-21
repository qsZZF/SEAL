classdef CortexNode < BaseNode
    %CORTEXNODE 皮层节点
    % 管理协议的皮层信息，包含皮层名称、位置等元数据
    
    properties (Constant)
        type = "CortexNode"  % 节点类型标识符
    end
    
    properties
        % 皮层特定属性
        cortexInfo CortexInfo
        cache
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        data
        infoFile string
    end
    
    methods
        function obj = CortexNode()
            %CORTEXNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, path)
            %OPEN 打开数据
            
            % 检查路径是文件夹还是文件
            if isfile(path)
                [path, name, ~] = fileparts(path);
            end

            % 设置路径
            obj.path = path;
            
            % 验证协议文件存在
            if ~isfile(obj.infoFile)
                error('SEAL:NoFile', ...
                    'No such file in path: %s', obj.infoFile);
            end
            
            try
                obj.cortexInfo = CortexInfo.openExisting(obj.infoFile);
            catch ME
                error('SEAL:ChannelNode:LoadFailed', ...
                    'Failed to open channel data: %s', ME.message);
            end
        end


        function save(obj)
            %SAVE 保存皮层数据
            obj.createDirectoryStructure();
            obj.cortexInfo.save(obj.infoFile);
        end

        function load(obj)
            %LOAD 加载皮层模型，按需提取并组装标准结构体
            if obj.isLoaded
                return;
            end
            
            if ~isempty(obj.cortexInfo) && isfile(obj.cortexInfo.dataPath)
                
                % 1. 临时读取原始数据
                rawData = loadData(obj.cortexInfo.dataPath); 
                
                % ==========================================
                % 2. 自动拆开“俄罗斯套娃”
                % ==========================================
                if isstruct(rawData)
                    outerFields = fieldnames(rawData);
                    if length(outerFields) == 1 && isstruct(rawData.(outerFields{1}))
                        rawData = rawData.(outerFields{1}); % 剥去外壳
                    end
                    availableFields = fieldnames(rawData);
                else
                    error('皮层模型必须是结构体或 .mat 文件！');
                end
                
                % ==========================================
                % 3. 定义三大刚需的别名查找字典
                % ==========================================
                % 顶点坐标 (N x 3)
                vertAliases = {'Vertices', 'pos', 'Position', 'vertices', 'nodes'};
                % 顶点连接/面 (M x 3)
                FacesAliases = { 'Faces', 'tri', 'Triangles', 'faces', 'polygons'};
                % 脑区图谱 (结构体数组或索引)
                atlasAliases = {'Atlas', 'Scout', 'atlas', 'regions', 'parcellation'};
                connAliases = {'VertConn'};
                % 初始化找到的真实字段名
                matchedVert = '';
                matchedFaces = '';
                matchedAtlas = '';
                matchedConn ='';
                % 自动查字典匹配
                for i = 1:length(vertAliases), if ismember(vertAliases{i}, availableFields), matchedVert = vertAliases{i}; break; end, end
                for i = 1:length(FacesAliases), if ismember(FacesAliases{i}, availableFields), matchedFaces = FacesAliases{i}; break; end, end
                for i = 1:length(atlasAliases), if ismember(atlasAliases{i}, availableFields), matchedAtlas = atlasAliases{i}; break; end, end
                for i = 1:length(connAliases), if ismember(connAliases{i}, availableFields), matchedConn = connAliases{i}; break; end, end
                % ==========================================
                % 4. 触发防御机制：缺哪个，就弹窗问哪个
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
                % 5. 终极净化与组装：装入纯净的 data_cache
                % ==========================================
                % 不管原始字段叫什么鬼名字，在这里统一“换发身份证”！
                obj.cache = struct();
                obj.cache.Vertices = rawData.(matchedVert);
                obj.cache.Faces = rawData.(matchedFaces);
                obj.cache.VertConn = rawData.(matchedConn);

                %Atlas校验
                if ~isempty(matchedAtlas)
                    A = rawData.(matchedAtlas);
                    ok = isstruct(A) && all(isfield(A, {'Name','Scouts'}));
                    if ok
                        s1 = A(1).Scouts;
                        ok = isstruct(s1) && all(isfield(s1, {'Vertices','Region','Color'}));
                    end
                    if ~ok
                        warning('SEAL:CortexNode:AtlasFormat', ...
                            'Atlas 缺少必需字段 (Name / Scouts / Scouts.Vertices / Region / Color)，已置空。');
                        obj.cache.Atlas = [];
                    else
                        obj.cache.Atlas = A;
                    end
                else
                    obj.cache.Atlas = [];
                end

                % (函数结束，庞杂的 rawData 被内存自动回收，只留下纯净的 obj.data_cache)
                obj.isLoaded = true;
            end
        end

        function unload(obj)
            %UNLOAD 卸载皮层数据，释放内存
            obj.cache = [];
            obj.isLoaded = false;
        end

        function deleteFromDisk(obj)
            delete(obj.infoFile);
        end
        
        %% 依赖属性get方法
        function cortexName = get.name(obj)
            cortexName = obj.cortexInfo.name;
        end
        
        function data = get.data(obj)
            if ~obj.isLoaded
                obj.load();
            end
            data = obj.cache;
        end

        function path = get.infoFile(obj)
            path = fullfile(obj.path, "cortex_info.mat");
        end
    
    end

    methods (Static)

        function cortex = fromInfo(cortexInfo)
            cortex = CortexNode();
            cortex.cortexInfo = cortexInfo();
        end

        function cortex = fromData(path)
            cortex = CortexNode();
            cortex.cortexInfo = CortexInfo.fromData(path);
        end

        function cortex = openExisting(srcPath)
            %OPENCORTEX 打开现有皮层数据
            % 输入:
            %   srcPath - 皮层文件路径或皮层文件夹路径
            %   parentProtocol - 父协议节点
            
            % 创建皮层节点
            cortex = CortexNode();

            % 加载皮层数据
            cortex.open(getmat(srcPath));
        end
    end

    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建皮层目录结构
            
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