classdef LeadfieldNode < BaseNode
    %LEADFIELDNODE 导联场节点
    % 管理协议的导联场信息，包含导联场名称、位置等元数据
    
    properties (Constant)
        type = "LeadfieldNode"  % 节点类型标识符
    end
    
    properties
        % 导联场特定属性
        leadfieldInfo LeadfieldInfo
        cache
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        data
        infoFile string
        sizee
        headModelType string
        orientation
        normals
    end
    
    methods
        function obj = LeadfieldNode()
            %LEADFIELDNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode();
        end
        
        function open(obj, path)
            %OPEN 打开数据
            
            % 检查路径是文件夹还是文件
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
                obj.leadfieldInfo = LeadfieldInfo.openExisting(obj.infoFile);
            catch ME
                error('SEAL:LeadfieldNode:LoadFailed', ...
                    'Failed to open leadfield data: %s', ME.message);
            end
        end


        function save(obj)
            %SAVE 保存导联场数据
            obj.createDirectoryStructure();
            obj.leadfieldInfo.save(obj.infoFile);
        end

        function load(obj)
            %LOAD 加载导联场数据
            if obj.isLoaded
                return;
            end
            if ~isempty(obj.leadfieldInfo) && isfile(obj.leadfieldInfo.dataPath)

                % 1. 临时读取原始数据 (在这个函数结束时，多余的内存会被自动回收)
                rawData = loadData(obj.leadfieldInfo.dataPath);

                matchedGainField = '';

                % ==========================================
                % 2. 类型防御 与 拆盒机制
                % ==========================================
                if isnumeric(rawData)
                    % 极其痛快：如果直接就是矩阵，直接塞进缓存
                    obj.cache = rawData;
                    matchedGainField = 'DirectNumeric';
                    availableFields = '';
                    
                elseif isstruct(rawData)
                    outerFields = fieldnames(rawData);

                    if length(outerFields) == 1
                        loneVar = rawData.(outerFields{1});
                        if isstruct(loneVar)
                            % 剥开俄罗斯套娃
                            rawData = loneVar;
                            availableFields = fieldnames(rawData);
                        elseif isnumeric(loneVar)
                            % 结构体里只有一个纯矩阵，直接提取！
                            obj.cache = loneVar;
                            matchedGainField = outerFields{1};
                        else
                            availableFields = outerFields;
                        end
                    else
                        % 多个变量，准备查字典
                        availableFields = outerFields;
                    end
                else
                    error('不支持的导联场数据类型。');
                end

                % ==========================================
                % 3. 字典匹配与精确提取 (按需取用)
                % ==========================================
                if isempty(matchedGainField)
                    gainAliases = {'Gain', 'L', 'Leadfield', 'LF', 'Forward', 'G', 'matrix'};
                    for i = 1:length(gainAliases)
                        if ismember(gainAliases{i}, availableFields)
                            matchedGainField = gainAliases{i};
                            % 【核心】：只把矩阵掏出来放进 cache，丢弃整个结构体！
                            obj.cache = rawData.(matchedGainField);
                            break;
                        end
                    end
                end


                normalAliases = {'GridOrient', 'Normals', 'Orientations', 'CorticalNormals', 'Norm'};
                for i = 1:length(normalAliases)
                    if ismember(normalAliases{i}, availableFields)
                        % 提取法向量并存入属性
                        obj.leadfieldInfo.normals = rawData.(normalAliases{i});
                        break; % 找到了就立刻停止搜寻
                    end
                end

                % ==========================================
                % 4. 人工兜底防御
                % ==========================================
                if isempty(matchedGainField)
                    [indx, tf] = listdlg('ListString', availableFields, ...
                        'SelectionMode', 'single', ...
                        'PromptString', {'导联场加载: 未识别到标准矩阵字段。', '请手动指定:'}, ...
                        'Name', '导联场匹配向导', ...
                        'ListSize', [250, 150]);
                    if tf
                        % 用户选啥，我们就只提取啥
                        obj.cache = rawData.(availableFields{indx});
                    else
                        error('SEAL:LeadfieldNode:LoadFailed', '用户取消了导联场匹配。');
                    end
                end
                % 加载导联场元数据

                obj.isLoaded = true;
            end
        end

        function unload(obj)
            %UNLOAD 卸载导联场数据，释放内存
            obj.cache = [];
            obj.isLoaded = false;
        end

        function deleteFromDisk(obj)
            delete(obj.infoFile);
        end


        function str = computeOrientationString(obj)
            %COMPUTEORIENTATIONSTRING 基于 orientation 数值生成可读描述
            ratio = obj.orientation;   % 只读一次，不再触发计算链

            if isempty(ratio) || ~isnumeric(ratio) || ~isscalar(ratio)
                str = mat2str(ratio);
                return;
            end

            switch ratio
                case 1
                    str = "Constrained (1 - Normal to cortex)";
                case 3
                    str = "Unconstrained (3 - Free orientation)";
                case 2
                    str = "Loose (2)";
                otherwise
                    str = sprintf("Unknown (ratio = %.3f)", ratio);
            end
        end

        function cortexNode = findSiblingCortex(obj)
            %FINDSIBLINGCORTEX 在父节点下查找同级的 CortexNode
            cortexNode = [];

            if isempty(obj.parent) || ~isprop(obj.parent, 'children')
                return;
            end

            siblings = obj.parent.children;
            for i = 1:numel(siblings)
                if isa(siblings(i), 'CortexNode')
                    cortexNode = siblings(i);
                    return;
                end
            end
        end

        
        %% 依赖属性get方法
        function leadfieldName = get.name(obj)
            leadfieldName = obj.leadfieldInfo.name;
        end

        function headModelType = get.headModelType(obj)
            headModelType = obj.leadfieldInfo.headModelType;
        end

        function ori = get.orientation(obj)
            % 默认值：无法判断时回退到存储值
            ori = obj.leadfieldInfo.orientation;

            try
                cortexNode = obj.findSiblingCortex();
                if isempty(cortexNode)
                    return;
                end

                lfSize = obj.sizee;         % [nChannels, nCols]
                cxSize = cortexNode.sizee;   % [nVertices, 3]

                if numel(lfSize) < 2 || numel(cxSize) < 2
                    return;
                end

                nCols = lfSize(2);
                nVertices = cxSize(1);

                if nVertices == 0
                    return;
                end

                ori = nCols / nVertices;
            catch
                % 出错就保持默认值
            end
        end

        function data = get.data(obj)
            if ~obj.isLoaded
                obj.load();
            end
            data = obj.cache;
        end

        function path = get.infoFile(obj)
            path = fullfile(obj.path, "leadfield_info.mat");
        end

        function rsize = get.sizee(obj)
            rsize = size(obj.data);
        end

        function metadata = getMetadata(obj)
             %orientationStr = obj.computeOrientationString();
            metadata =  obj.getMetadata@BaseNode();
            metadata = [metadata;
                "Head Model Type", string(obj.headModelType);
                "Orientation", obj.orientation;
                 "Size", mat2str(obj.sizee)];
        end
    
        function n = get.normals(obj)
            n = obj.leadfieldInfo.normals;
        end
    end

    methods (Static)
        
        function leadfield = fromInfo(leadfieldInfo)
            leadfield = LeadfieldNode();
            leadfield.leadfieldInfo = leadfieldInfo;
        end

        function leadfield = fromData(path)
            leadfield = LeadfieldNode();
            leadfield.leadfieldInfo = LeadfieldInfo.fromData(path);
        end

        function leadfield = openExisting(srcPath)
            %OPENLEADFIELD 打开现有导联场数据
            % 输入:
            %   srcPath - 导联场文件路径或导联场文件夹路径
            %   parentProtocol - 父协议节点
            
            % 创建导联场节点
            leadfield = LeadfieldNode();

            % 加载导联场数据
            leadfield.open(getmat(srcPath));
        end
    end

    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建导联场目录结构
            
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