classdef SEAL_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        FileMenu             matlab.ui.container.Menu
        NewProjectMenu       matlab.ui.container.Menu
        LoadProjectMenu      matlab.ui.container.Menu
        NewprotocolMenu      matlab.ui.container.Menu
        ImportprotocolMenu   matlab.ui.container.Menu
        HelpMenu             matlab.ui.container.Menu
        OperationsPanel      matlab.ui.container.Panel
        RunBatchButton       matlab.ui.control.Button
        EvaluationButton     matlab.ui.control.Button
        SimulationButton     matlab.ui.control.Button
        ESIAnalysisButton    matlab.ui.control.Button
        VisualizationButton  matlab.ui.control.Button
        PreprocessButton     matlab.ui.control.Button
        StatusPanel          matlab.ui.container.Panel
        TextArea             matlab.ui.control.TextArea
        ProjectPanel         matlab.ui.container.Panel
        FileTree             matlab.ui.container.Tree
        MetadataPanel        matlab.ui.container.Panel
        GridLayout           matlab.ui.container.GridLayout
        ContextMenu          matlab.ui.container.ContextMenu
        CopyMenu             matlab.ui.container.Menu
        ContextMenu2         matlab.ui.container.ContextMenu
        ClearMenu            matlab.ui.container.Menu
    end


    properties (Access = private)
        DialogApp_showEEG
        DialogApp_showcortex
        Dialogapp_selectAlgorithm
        DialogApp_newprotocol
        DialogApp_datacheck
        DialogApp_preprocess
        DialogApp_simulation
        datanode={}
        Cortex={}               %导入后的Cortex数据
        L={}                    %导入后的lead-field Matrix数据
        Chanlocs={}             %导入后的channel locations数据
        OutputHistory           string  % 保存所有输出内容
        ExpandedNodeNames = string.empty;       %保存节点展开状态
        SeenNodeNames = string.empty;
        CurrentTaskRecipe
        importchanlocsmenu
        importcortexmenu
        importleadfieldmenu
        visTopoMenu
    end

    properties (Access = public)
        projectNode ProjectNode % root of project tree
        ProjectPath             % path to save project folders
        DataMapping             % For load .mat data
        wasCancelled
        newprotocol             string
        protocolnode={}         %树节点数组
        inode = 0               %树节点索引(自动累计)
        protocolmenu            %树节点上下文菜单
        cortexMenu
        eegDataFileMenu
        leadfieldMenu
        chanlocsMenu
        sessionMenu
        resultFileMenu
        PreprocessedDataMenu
        SimuParamMenu
        firsttime=1
        DialogApp_showResult = {}
    end

    methods (Access = private)
        % 核心函数：向文本区域添加内容

        function ImportDataSelected(app, dataType)
            % Get the currently selected protocol node
            selectedNode = app.FileTree.SelectedNodes;

            if numel(selectedNode) > 1
                uialert(app.UIFigure, ['Please select only one protocol to...' ...
                    ' import data into.'], 'Multiple Protocols Selected');
                return;
            end
            if  ~isa(selectedNode.Parent,'matlab.ui.container.Tree')
                uialert(app.UIFigure, 'Click the protocol first', 'Protocol not Selected');
                return
            end

            protocolNode = selectedNode.NodeData;

            % ---Overwrite check for Cortex/Leadfield, one for each protocol---
            if strcmp(dataType, 'Cortex') || strcmp(dataType, 'Leadfield')|| strcmp(dataType, 'Chanlocs')
                if isvalid(protocolNode.(strcat(lower(dataType), "Node")))
                    question = ['This protocol already has a ' strrep(dataType, '_', ' ') ' file. Overwrite it?'];
                    selection = uiconfirm(app.UIFigure, question, 'Confirm Overwrite', 'Icon', 'warning');
                    if ~strcmp(selection, 'OK')
                        app.appendOutput('Import cancelled by user.');
                        return;
                    end

                end
            end

          % --- Build filter based on dataType ---
switch dataType
    case 'Data'
        % Data 支持 .mat 和 .set
        filterPattern = '*.mat;*.set';
    case {'Cortex', 'Leadfield', 'Chanlocs'}
        % 其他类型只支持 .mat
        filterPattern = '*.mat';
    otherwise
        filterPattern = '*.*';
end

[fileName, filePath, ~] = uigetfile(filterPattern, ...
    ['Select ', strrep(dataType, '_', ' '), ' File']);

if isequal(fileName, 0)
    app.appendOutput('Import data cancelled');
    return;
end

originalPath = fullfile(filePath, fileName);
[~, ~, ext] = fileparts(fileName);

% --- Route by file extension ---
try
    switch lower(ext)
        case '.mat'
            switch dataType
                case 'Chanlocs'
                    protocolNode.openChanlocsFromData(originalPath);
                case 'Cortex'
                    protocolNode.openCortexFromData(originalPath);
                case 'Leadfield'
                    protocolNode.openLeadfieldFromData(originalPath);
                case 'Data'
                    protocolNode.openDataFromData(originalPath);
                otherwise
                    errordlg(['Unsupported dataType: ', char(dataType)], 'Import Error');
                    return;
            end
        case '.set'
            if strcmp(dataType, 'Data')
                protocolNode.openDataFromEEGLAB(originalPath);
            else
                errordlg(['.set files can only be imported as Data, not ', char(dataType)], ...
                         'Import Error');
                return;
            end
        otherwise
            errordlg(['Unsupported file extension: ', ext], 'Import Error');
            return;
    end
catch ME
    errordlg(['Failed to import file. Error: ', ME.message], 'Import Error');
    return;
end
       
            app.flushFileTree();
            text = strcat(dataType, ' imported successfully for protocol ''', protocolNode.name, '''.');
            app.appendOutput(text);
            app.bringToFront();
        end

        function DeleteNodeSelected(app, ~)
            selectedNodes = app.FileTree.SelectedNodes;
            if isempty(selectedNodes)
                uialert(app.UIFigure, 'Please select one or more protocols to delete.', 'No Selection');
                return;
            end
            numNodes = numel(selectedNodes);
            if numNodes == 1
                question = ['Are you sure you want to permanently delete the "', selectedNodes(1).Text, '"and all its subnodes? This operation will not affect the original file.'];
            else
                question = ['Are you sure you want to permanently delete ', num2str(numNodes), ' selected protocols and all their subnodes? This operation will not affect the original file.'];
            end

            selection = uiconfirm(app.UIFigure, question, 'Confirm Deletion', 'Icon', 'warning');
            if strcmp(selection, 'OK')
                for k = 1:numNodes
                    node = selectedNodes(k);
                    nodeData = node.NodeData;
                    parentNode = nodeData.parent;
                    if isvalid(parentNode)
                        parentNode.removeChild(nodeData);
                        % 单例节点删除后须清空 ProtocolNode 的具名引用，
                        % 否则 isempty 对已销毁 handle 仍返回 false
                        if isa(parentNode, 'ProtocolNode')
                            if isa(nodeData, 'CortexNode')
                                parentNode.cortexNode = CortexNode.empty;
                            elseif isa(nodeData, 'ChanlocsNode')
                                parentNode.chanlocsNode = ChanlocsNode.empty;
                            elseif isa(nodeData, 'LeadfieldNode')
                                parentNode.leadfieldNode = LeadfieldNode.empty;
                            end
                        end
                    end
                    if isvalid(node)
                        nodeData.deleteFromDisk();
                        nodeData.delete;
                    end
                end
                app.flushFileTree();
            end
        end

        function VisualizeDataSelected(app, vizType)
            selectedNode = app.FileTree.SelectedNodes;
            if isempty(selectedNode), return; end

            switch vizType
                case 'EEG_TimeSeries'
                    dataNode=selectedNode.NodeData;
                    app.DialogApp_showEEG=SEAL_showEEG(dataNode,app.UIFigure.Position);

                case 'EEG_Topo'
                    dataNode = selectedNode.NodeData;
                    topoplotX([],dataNode.dataInfo.chanlocs);
                case 'Cortex'
                    app.DialogApp_showcortex=SEAL_showCortex(app);
                case 'Results'
                    app.DialogApp_showResult =SEAL_showCortex(app);
                otherwise
                    uialert(app.UIFigure, ['Unknown visualization type: ' vizType], 'Error');
            end
        end

        function RunSourceAnalysisSelected(app, ~)
            selectedNode = app.FileTree.SelectedNodes;
            if isempty(selectedNode), return; end
            app.Dialogapp_selectAlgorithm=SEAL_selectAlgorithm(app,selectedNode);


        end

        function [Cortex, metadata] = import_from_bst_surface(~, filePath)
            % Converts a Brainstorm surface .mat file to the SEAL Cortex structure
            s = load(filePath);
            s = struct2cell(s); s = s{1};
            Cortex = struct();
            Cortex.Vertices = s.Vertices;
            Cortex.Faces = s.Faces;
            tIdx = strcmp('Structures',string(struct2table(s.Atlas).Name));
            Cortex.Structure = s.Atlas(tIdx);
            als_idx = setdiff(1:numel(s.Atlas), tIdx);
            if ~isempty(als_idx)
                Cortex.Atlas = s.Atlas(als_idx);
            end

            if isfield(s, 'VertConn'), Cortex.VertConn = s.VertConn; end
            if isfield(s, 'Comment'), Cortex.Comment = s.Comment; end
            if isfield(s, 'VertNormals'), Cortex.VertNormals = s.VertNormals; end
            if isfield(s, 'Curvature'), Cortex.Curvature = s.Curvature; end
            if isfield(s, 'SulciMap'), Cortex.SulciMap = s.SulciMap; end
            if isfield(s, 'tess2mri_interp'), Cortex.tess2mri_interp = s.tess2mri_interp; end
            if isfield(s, 'Reg'), Cortex.Reg = s.Reg; end
            metadata.Cortex = 'Yes';
        end

        function ExtractRecipeMenuSelected(app, ~)
            % 从 Result 的 .mat 文件中提取批处理配方 (ResultNode 右键菜单回调)
            
            % 1. 获取当前在树状图中右键选中的节点
            selectedUITreeNode = app.FileTree.SelectedNodes;
            if isempty(selectedUITreeNode)
                return;
            end
            
            targetNode = selectedUITreeNode.NodeData;
            
            % 2. 严格防御：确保点中的确实是 ResultNode
            if ~isa(targetNode, 'ResultNode')
                uialert(app.UIFigure, '只能从【结果节点 (ResultNode)】中提取计算配方！', ...
                    '节点类型错误', 'Icon', 'warning');
                return;
            end
            
            % ===================================================
            % 3. 核心逻辑：定位并加载物理 .mat 文件
            % ===================================================
            % 获取该结果文件在硬盘上的绝对路径
            matFilePath = targetNode.infoFile; 
            
            if ~isfile(matFilePath)
                uialert(app.UIFigure, '找不到该结果对应的物理 .mat 文件！可能已被移动或删除。', ...
                    '文件丢失', 'Icon', 'error');
                return;
            end
            
            % 默默加载这个文件中的所有变量
            try
                loadedData = load(matFilePath);
            catch ME
                uialert(app.UIFigure, sprintf('读取结果文件失败:\n%s', ME.message), ...
                    '读取错误', 'Icon', 'error');
                return;
            end
            resultObj = loadedData.resultObj;
            % ===================================================
            % 4. 提取参数与算法名 (需要根据你的 save 习惯微调)
            % ===================================================
            % 假设你在算完保存时，是直接把变量存进去的：
            % save(path, 'sourceEstimate', 'AlgorithmName', 'Parameters')
            % 那么它们现在就是 loadedData.AlgorithmName 和 loadedData.Parameters
            
            % 我们加一个智能判断，看你到底是直接存的变量，还是包在了 extraInfo 结构体里
            if isprop(resultObj, 'algorithm') && isprop(resultObj, 'params')
                algoName = resultObj.algorithm; 
                algoParams = resultObj.params;
                            
            else
                % 兜底防御：文件里没有存这些信息
                uialert(app.UIFigure, '在该 .mat 文件中未找到算法名称(AlgorithmName)或参数(Parameters)，无法提取配方！', ...
                    '缺少记录', 'Icon', 'error');
                return;
            end
            
            % ===================================================
            % 5. 组装 taskRecipe 并暂存在 App 中
            % ===================================================
            taskRecipe = struct();
            taskRecipe.AlgorithmName = algoName;
            taskRecipe.Parameters = algoParams;
            
            app.CurrentTaskRecipe = taskRecipe; % 保存到主 App 属性供后续批处理使用
            app.RunBatchButton.Enable = "on";
            % 6. 友好的用户反馈
            successMsg = sprintf('批处理配方已成功从文件中提取！\n\n目标算法: %s\n\n您现在可以点击 [Run Batch] 按钮，将此配方应用到其他被试文件夹了。', taskRecipe.AlgorithmName);
            uialert(app.UIFigure, successMsg, '配方提取成功', 'Icon', 'success');
            
        end
    
        function ExtractSessionRecipeMenuSelected(app, ~)
            % 1. 获取当前选取的 SessionNode 对象
            % (假设你通过 FileTree 的 SelectedNodes 获取自定义对象)
            sessionNode = app.FileTree.SelectedNodes.NodeData;

            % 2. 获取该 Session 下的所有子节点，并严格筛选出 ResultNode
            allChildren = sessionNode.children;
            resultNodes = []; % 初始化一个空数组用于存放真正的 ResultNode

            % 遍历所有的子节点，使用 isa 进行类型判断
            for k = 1:length(allChildren)
                % 注意：这里的 'ResultNode' 必须是你定义的类名（严格区分大小写）
                % 如果 allChildren 是元胞数组 (cell)，请改为 isa(allChildren{k}, 'ResultNode')
                if isa(allChildren(k), 'ResultNode')
                    resultNodes = [resultNodes, allChildren(k)];
                end
            end

            if isempty(resultNodes)
                uialert(app.UIFigure, '该 Session 下没有任何处理结果，无法提取配方。', '提示');
                return;
            end

            % 3. 初始化管线结构体数组 (Pipeline)
            % 字段设计：算法类型，具体参数
            pipelineRecipe = struct('AlgorithmName', {}, 'Parameters', {});

            for i = 1:length(resultNodes)
                currentNode = resultNodes(i);
                infoPath = currentNode.infoFile; % 使用你提到的 Infofile 路径

                if exist(infoPath, 'file') == 2
                    % 加载 ResultNode 的 .mat 文件
                    tempData = load(infoPath);
                    tempData = tempData.resultObj;
                    stepInfo.AlgorithmName = tempData.algorithm;
                    stepInfo.Parameters = tempData.params;

                    % 将当前步骤追加到管线中
                    pipelineRecipe(end+1) = stepInfo;
                else
                    warning('跳过节点 %s：未找到对应的 Infofile (%s)', currentNode.name, infoPath);
                end
            end

            % 5. 校验并保存
            if isempty(pipelineRecipe)
                uialert(app.UIFigure, '未能提取到任何有效的步骤参数！', '提取失败');
                return;
            end

             app.CurrentTaskRecipe = pipelineRecipe; % 保存到主 App 属性供后续批处理使用
            app.RunBatchButton.Enable = "on";
            % 6. 友好的用户反馈
            successMsg = sprintf('批处理配方已成功从文件中提取！\n\n您现在可以点击 [Run Batch] 按钮，将此配方应用到其他被试文件夹了。');
            uialert(app.UIFigure, successMsg, '配方提取成功', 'Icon', 'success');
        end

        
        function ImportChanlocs(app, ~)
            %ADDCHANLOCS 为当前选中的 DataNode 添加通道位置

            % 1. 获取当前选中的节点
            selectedNode = app.FileTree.SelectedNodes;
            if isempty(selectedNode)
                uialert(app.UIFigure, '请先选择一个数据节点！', '未选中节点');
                return;
            end

            % 从树节点拿到对应的 DataNode（取决于你的存储方式）
            dataNode = selectedNode.NodeData;
            if ~isa(dataNode, 'DataNode')
                uialert(app.UIFigure, '请选择一个数据节点！', '节点类型错误');
                return;
            end

            % 2. 选择文件
            [file, path] = uigetfile({
                '*.mat',  'MAT 文件 (*.mat)'; ...
                '*.locs', 'EEGLAB locs (*.locs)'; ...
                '*.ced',  'EEGLAB ced (*.ced)'; ...
                '*.xyz',  'XYZ 坐标 (*.xyz)'; ...
                '*.sfp',  'SFP 坐标 (*.sfp)'; ...
                '*.*',    '所有文件 (*.*)'}, ...
                '选择通道位置文件');
            if isequal(file, 0)
                return;
            end

            filePath = fullfile(path, file);
            [~, ~, ext] = fileparts(filePath);

            % 3. 解析文件
            try
                switch lower(ext)
                    case '.mat'
                        chanlocs = app.parseChanlocsMat(filePath);
                    case {'.locs', '.ced', '.xyz', '.sfp'}
                        chanlocs = app.parseChanlocsText(filePath, ext);
                    otherwise
                        uialert(app.UIFigure, ...
                            sprintf('不支持的文件格式: %s', ext), '格式错误');
                        return;
                end

                % 4. 通道数校验
                if ~isempty(dataNode.size) && dataNode.size(1) > 0
                    nChan = dataNode.size(1);
                    nLocs = length(chanlocs);
                    if nChan ~= nLocs
                        answer = uiconfirm(app.UIFigure, ...
                            sprintf('通道数不匹配！\n数据: %d 通道\n位置文件: %d 通道\n是否仍然导入？', nChan, nLocs), ...
                            '通道数警告', ...
                            'Options', {'导入', '取消'}, ...
                            'DefaultOption', '取消');
                        if strcmp(answer, '取消')
                            return;
                        end
                    end
                end

                % 5. 赋值并保存
                dataNode.dataInfo.chanlocs = chanlocs;
                dataNode.save();
                app.bringToFront();
                uialert(app.UIFigure, ...
                    sprintf('成功导入 %d 个通道位置。', length(chanlocs)), ...
                    '导入成功', 'Icon', 'success');

            catch ME
                uialert(app.UIFigure, ...
                    sprintf('通道位置导入失败:\n%s', ME.message), ...
                    '导入错误', 'Icon', 'error');
            end
        end

        function chanlocs = parseChanlocsMat(app, filePath)
            rawData = load(filePath);

            % 拆盒
            fields = fieldnames(rawData);
            if length(fields) == 1 && isstruct(rawData.(fields{1}))
                rawData = rawData.(fields{1});
                fields = fieldnames(rawData);
            end

            % 字典匹配
            aliases = {'chanlocs', 'Channel', 'channels', 'chaninfo', ...
                'electrodes', 'locs', 'labels', 'elec'};
            matchedField = '';
            for i = 1:length(aliases)
                if ismember(aliases{i}, fields)
                    matchedField = aliases{i};
                    break;
                end
            end

            if isempty(matchedField)
                [indx, tf] = listdlg('ListString', fields, ...
                    'SelectionMode', 'single', ...
                    'PromptString', '未识别到通道位置字段，请手动选择:', ...
                    'Name', '字段匹配', ...
                    'ListSize', [250, 150]);
                if tf
                    matchedField = fields{indx};
                else
                    error('用户取消了字段选择。');
                end
            end

            raw = rawData.(matchedField);

            % 已经是结构体 → 补全缺失字段
            if isstruct(raw)
                chanlocs = app.normalizeChanlocsStruct(raw);

                % 纯坐标矩阵 [N×3] → 从笛卡尔坐标构建
            elseif isnumeric(raw) && size(raw, 2) >= 3
                chanlocs = app.buildChanlocsFromXYZ(raw, {});
            else
                error('无法识别的通道位置数据格式。');
            end
        end

        function chanlocs = parseChanlocsText(app, filePath, ext)
            fid = fopen(filePath, 'r');
            if fid == -1
                error('无法打开文件: %s', filePath);
            end
            cleanup = onCleanup(@() fclose(fid));

            switch lower(ext)
                case '.xyz'
                    % 格式: index X Y Z label
                    data = textscan(fid, '%d %f %f %f %s');
                    coords = [data{2}, data{3}, data{4}];
                    labels = data{5};
                    chanlocs = app.buildChanlocsFromXYZ(coords, labels);

                case '.sfp'
                    % 格式: label X Y Z
                    data = textscan(fid, '%s %f %f %f');
                    coords = [data{2}, data{3}, data{4}];
                    labels = data{1};
                    chanlocs = app.buildChanlocsFromXYZ(coords, labels);

                case {'.locs', '.ced'}
                    if exist('readlocs', 'file')
                        raw = readlocs(filePath);
                        chanlocs = app.normalizeChanlocsStruct(raw);
                    else
                        % .locs 格式: index theta radius label
                        data = textscan(fid, '%d %f %f %s');
                        chanlocs = app.buildChanlocsFromPolar(data{2}, data{3}, data{4});
                    end
            end
        end

        function chanlocs = buildChanlocsFromXYZ(~, coords, labels)
            %BUILDCHANLOCS 从笛卡尔坐标 [N×3] 构建完整 chanlocs

            n = size(coords, 1);
            chanlocs = struct();

            for i = 1:n
                x = coords(i, 1);
                y = coords(i, 2);
                z = coords(i, 3);

                % 标签
                if ~isempty(labels) && i <= length(labels)
                    chanlocs(i).labels = char(labels{i});
                else
                    chanlocs(i).labels = sprintf('Ch%d', i);
                end

                chanlocs(i).type = '';

                % 笛卡尔 → 球面坐标
                [sph_theta, sph_phi, sph_radius] = cart2sph(x, y, z);
                sph_theta_deg = rad2deg(sph_theta);
                sph_phi_deg   = rad2deg(sph_phi);

                % 球面 → EEGLAB 极坐标
                % EEGLAB 约定: theta = -sph_theta, radius 投影到 2D
                theta  = -sph_theta_deg;
                radius = 0.5 - sph_phi_deg / 180;

                chanlocs(i).theta      = theta;
                chanlocs(i).radius     = radius;
                chanlocs(i).X          = x;
                chanlocs(i).Y          = y;
                chanlocs(i).Z          = z;
                chanlocs(i).sph_theta  = sph_theta_deg;
                chanlocs(i).sph_phi    = sph_phi_deg;
                chanlocs(i).sph_radius = sph_radius;
                chanlocs(i).urchan     = i;
                chanlocs(i).ref        = '';
            end
        end

        function chanlocs = buildChanlocsFromPolar(~, theta, radius, labels)
            %BUILDCHANLOCSFROMPOLAR 从 EEGLAB 极坐标构建完整 chanlocs

            n = length(theta);
            chanlocs = struct();

            for i = 1:n
                th  = theta(i);
                rad = radius(i);

                % 标签
                if ~isempty(labels) && i <= length(labels)
                    chanlocs(i).labels = char(labels{i});
                else
                    chanlocs(i).labels = sprintf('Ch%d', i);
                end

                chanlocs(i).type   = '';
                chanlocs(i).theta  = th;
                chanlocs(i).radius = rad;

                % EEGLAB 极坐标 → 球面坐标
                sph_theta_deg = -th;
                sph_phi_deg   = (0.5 - rad) * 180;
                sph_radius    = 85;  % EEGLAB 默认头半径

                % 球面 → 笛卡尔
                [x, y, z] = sph2cart(deg2rad(sph_theta_deg), ...
                    deg2rad(sph_phi_deg), ...
                    sph_radius);

                chanlocs(i).X          = x;
                chanlocs(i).Y          = y;
                chanlocs(i).Z          = z;
                chanlocs(i).sph_theta  = sph_theta_deg;
                chanlocs(i).sph_phi    = sph_phi_deg;
                chanlocs(i).sph_radius = sph_radius;
                chanlocs(i).urchan     = i;
                chanlocs(i).ref        = '';
            end
        end

        function chanlocs = normalizeChanlocsStruct(app, raw)
            %NORMALIZECHANLOCS 补全已有结构体中缺失的字段

            n = length(raw);

            % 检测有哪些坐标可用
            hasXYZ   = isfield(raw, 'X') && isfield(raw, 'Y') && isfield(raw, 'Z');
            hasSph   = isfield(raw, 'sph_theta') && isfield(raw, 'sph_phi');
            hasPolar = isfield(raw, 'theta') && isfield(raw, 'radius');

            if hasXYZ
                % 从 XYZ 反算其他坐标
                coords = [[raw.X]', [raw.Y]', [raw.Z]'];
                labels = {};
                for i = 1:n
                    if isfield(raw, 'labels')
                        labels{i} = raw(i).labels;
                    else
                        labels{i} = sprintf('Ch%d', i);
                    end
                end
                chanlocs = app.buildChanlocsFromXYZ(coords, labels);

                % 保留原始结构体中有但标准字段没覆盖到的额外字段
                extraFields = setdiff(fieldnames(raw), fieldnames(chanlocs));
                for i = 1:n
                    for j = 1:length(extraFields)
                        chanlocs(i).(extraFields{j}) = raw(i).(extraFields{j});
                    end
                end

            elseif hasPolar
                theta  = [raw.theta];
                radius = [raw.radius];
                labels = {};
                for i = 1:n
                    if isfield(raw, 'labels')
                        labels{i} = raw(i).labels;
                    else
                        labels{i} = sprintf('Ch%d', i);
                    end
                end
                chanlocs = app.buildChanlocsFromPolar(theta, radius, labels);

            elseif hasSph
                % 球面 → 笛卡尔 → 走 XYZ 流程
                sph_r = 85 * ones(n, 1);
                if isfield(raw, 'sph_radius')
                    sph_r = [raw.sph_radius]';
                end
                [x, y, z] = sph2cart(deg2rad([raw.sph_theta]'), ...
                    deg2rad([raw.sph_phi]'), ...
                    sph_r);
                coords = [x, y, z];
                labels = {};
                for i = 1:n
                    if isfield(raw, 'labels')
                        labels{i} = raw(i).labels;
                    else
                        labels{i} = sprintf('Ch%d', i);
                    end
                end
                chanlocs = app.buildChanlocsFromXYZ(coords, labels);

            else
                % 只有标签没有坐标，保留原样，补上空字段
                chanlocs = raw;
                emptyFields = {'type','theta','radius','X','Y','Z', ...
                    'sph_theta','sph_phi','sph_radius','urchan','ref'};
                for i = 1:n
                    for j = 1:length(emptyFields)
                        if ~isfield(chanlocs, emptyFields{j})
                            chanlocs(i).(emptyFields{j}) = [];
                        end
                    end
                    if ~isfield(raw, 'urchan') || isempty(chanlocs(i).urchan)
                        chanlocs(i).urchan = i;
                    end
                end
            end
        end

        function EEGMenuOpening(app, ~)
            % 获取当前选中节点
            selectedNode = app.FileTree.SelectedNodes;
            if isempty(selectedNode)
                return;
            end

            dataNode = selectedNode.NodeData;
            if ~isa(dataNode, 'DataNode')
                return;
            end

            % 检测 chanlocs
            hasChanlocs = ~isempty(dataNode.dataInfo) ...
                && isprop(dataNode.dataInfo, 'chanlocs') ...
                && ~isempty(dataNode.dataInfo.chanlocs);

            if hasChanlocs

                % 已有 chanlocs，改菜单文字提示用户可以重新导入
                app.importchanlocsmenu.Enable = 'off';
                 app.visTopoMenu.Enable = 'on';
            else
                app.importchanlocsmenu.Enable = 'on';
                 app.visTopoMenu.Enable = 'off';
            end
        end


        function ProtocolMenuOpening(app, ~)
            selectedNode = app.FileTree.SelectedNodes;
            if isempty(selectedNode)
                return;
            end
            dataNode = selectedNode.NodeData;
            if ~isa(dataNode, 'ProtocolNode')
                return;
            end

            % 检测 chanlocs
            hasCortex = ~isempty(dataNode.cortexNode) && isvalid(dataNode.cortexNode);
            hasLeadfield = ~isempty(dataNode.leadfieldNode) && isvalid(dataNode.leadfieldNode);

            if hasCortex

                % 已有 chanlocs，改菜单文字提示用户可以重新导入
                app.importcortexmenu.Enable = 'off';
            else
                app.importcortexmenu.Enable = 'on';
            end


            if hasLeadfield

                % 已有 chanlocs，改菜单文字提示用户可以重新导入
                app.importleadfieldmenu.Enable = 'off';
            else
                app.importleadfieldmenu.Enable = 'on';
            end
        end


        function showSplash(app, imagePath, duration)
            % 显示开屏图片
            % imagePath: 图片路径（如 'docImg/splash.png'）
            % duration:  停留秒数（如 2）

            img = imread(imagePath);
            [h, w, ~] = size(img);

            % 获取屏幕尺寸
            screenSize = get(0, 'ScreenSize');
            screenW = screenSize(3);
            screenH = screenSize(4);

            % 目标尺寸：占屏幕的 40%（按需调整）
            targetRatio = 0.4;
            maxW = screenW * targetRatio;
            maxH = screenH * targetRatio;

            % 按比例缩放，保持长宽比
            scale = min(maxW / w, maxH / h);
            if scale < 1
                w = round(w * scale);
                h = round(h * scale);
            end

            % 屏幕居中
            posX = round((screenW - w) / 2);
            posY = round((screenH - h) / 2);

            % 创建无边框 figure
            splashFig = figure( ...
                'Position',    [posX, posY, w, h], ...
                'MenuBar',     'none', ...
                'ToolBar',     'none', ...
                'NumberTitle', 'off', ...
                'Name',        'SEAL', ...
                'Resize',      'off', ...
                'WindowStyle', 'modal', ...
                'Color',       'w');

            % 铺满整个 figure 显示图片
            ax = axes('Parent', splashFig, ...
                'Units', 'normalized', ...
                'Position', [0 0 1 1]);
            imshow(img, 'Parent', ax);
            axis(ax, 'off');
            figure(splashFig);
            drawnow;
            try
                warnState = warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
                warning('off', 'MATLAB:ui:javaframe:PropertyToBeRemoved');

                jFrame = get(splashFig, 'JavaFrame');
                jWindow = jFrame.fHG2Client.getWindow();
                jWindow.setAlwaysOnTop(true);

                % 恢复 warning 状态
    warning(warnState);
            catch
                % 新版本 MATLAB 没有 JavaFrame，忽略即可
            end
            pause(duration);

            if isvalid(splashFig)
                close(splashFig);
            end
        end

    end

    methods (Access = public)
        function appendOutput(app, text)
            % 拼接新内容（添加换行符模拟命令行）
            app.OutputHistory = app.OutputHistory + text + newline;

            % 更新文本区域显示
            app.TextArea.Value = app.OutputHistory;

            scroll(app.TextArea, 'bottom');
        end

        function flushFileTree(app)
            app.FileTree.Children.delete;

            for i=1:length(app.projectNode.children)
                protocolNode = app.projectNode.children(i);
                if ~isvalid(protocolNode)
                    continue;
                end
                protocolView = uitreenode(app.FileTree);
                protocolView.Text = string(protocolNode.name);
                protocolView.ContextMenu = app.protocolmenu;
                protocolView.NodeData = protocolNode;
                protocolView.Icon = 'icon_protocol.svg';
                if isvalid(protocolNode.chanlocsNode)
                    chanlocsView = uitreenode(protocolView);
                    chanlocsView.Text = string(protocolNode.chanlocsNode.name);
                    chanlocsView.ContextMenu = app.chanlocsMenu;
                    chanlocsView.NodeData = protocolNode.chanlocsNode;
                end
                if isvalid(protocolNode.cortexNode)
                    cortexView = uitreenode(protocolView);
                    cortexView.Text = string(protocolNode.cortexNode.name);
                    cortexView.ContextMenu = app.cortexMenu;
                    cortexView.NodeData = protocolNode.cortexNode;
                    cortexView.Icon = 'icon_cortex.svg';
                end
                if isvalid(protocolNode.leadfieldNode)
                    leadfieldView = uitreenode(protocolView);
                    leadfieldView.Text = string(protocolNode.leadfieldNode.name);
                    leadfieldView.ContextMenu = app.leadfieldMenu;
                    leadfieldView.NodeData = protocolNode.leadfieldNode;
                    leadfieldView.Icon = 'icon_lead_field.svg';
                end
                for j=1:length(protocolNode.children)
                    sessionNode = protocolNode.children(j);

                    if ~isvalid(sessionNode) ||~isa(sessionNode,'SessionNode')
                        continue;
                    end
                    sessionView = uitreenode(protocolView);
                    sessionView.Text = string(sessionNode.name);
                    sessionView.ContextMenu = app.sessionMenu;
                    sessionView.NodeData = sessionNode;
                    sessionView.Icon = 'icon_session.svg';
                    for k = 1:length(sessionNode.children)
                        childNode = sessionNode.children(k);

                        if ~isvalid(childNode)
                            continue;
                        end

                        childView = uitreenode(sessionView);
                        childView.Text = string(childNode.name);
                        childView.NodeData = childNode;

                        if isa(childNode, 'DataNode')
                            % 如果是原始/预处理的脑电数据
                            childView.ContextMenu = app.eegDataFileMenu;
                            childView.Icon = 'icon_data.svg';
                        elseif isa(childNode, 'ResultNode')
                            % 如果是源成像/计算结果
                            childView.ContextMenu = app.resultFileMenu;
                            childView.Icon = 'icon_result.svg';
                        else
                            
                        end
                    end
                    sessionNameStr = string(sessionNode.name);
                    % 如果是第一次见到这个 Session（刚启动加载，或刚新建的）
                    if ~ismember(sessionNameStr, app.SeenNodeNames)
                        app.SeenNodeNames(end+1) = sessionNameStr;      
                        app.ExpandedNodeNames(end+1) = sessionNameStr;   % 默认给它展开特权
                    end
                    if ismember(string(sessionNode.name), app.ExpandedNodeNames)
                        expand(sessionView);
                    end
                end
                protocolNameStr = string(protocolNode.name);
                % 如果是第一次见到这个 Protocol
                if ~ismember(protocolNameStr, app.SeenNodeNames)
                    app.SeenNodeNames(end+1) = protocolNameStr;      % 登记到花名册
                    app.ExpandedNodeNames(end+1) = protocolNameStr;  % 默认给它展开特权
                end
                if ismember(string(protocolNode.name), app.ExpandedNodeNames)
                    expand(protocolView);
                end
            end
           
            app.SeenNodeNames = unique(app.SeenNodeNames);
            app.ExpandedNodeNames = unique(app.ExpandedNodeNames);
        end

        function clearMetadataGrid(app)
            app.GridLayout.Children.delete;
        end

        function flushMetadataGrid(app, node)
            app.clearMetadataGrid();
            % 【核心防御】清理之前残留的右键菜单，防止内存泄漏！
            oldMenus = findall(app.UIFigure, 'Type', 'uicontextmenu', 'Tag', 'PathCopyMenu');
            delete(oldMenus);

            metadata = node.NodeData.getMetadata();
            
            for i = 1:size(metadata, 1)
                keyText = string(metadata(i, 1));
                valText = string(metadata(i, 2));
                
                % 1. 创建 Key 标签
                key = uilabel(app.GridLayout, ...
                    "Text", keyText, ...
                    "FontWeight", "bold", ...
                    "Tooltip", keyText);
                
                % 2. 智能判断：如果是路径，挂载“右键菜单”！
                if contains(lower(keyText), "path")
                    % 视觉上绝对纯净，和普通文本一模一样
                    value = uilabel(app.GridLayout, ...
                        "Text", valText, ...
                        "Tooltip", valText + newline + "右键点击可复制路径");
                    
                    % 动态创建专属右键菜单，打上 Tag 方便下次清理
                    cm = uicontextmenu(app.UIFigure, 'Tag', 'PathCopyMenu');
                    uimenu(cm, 'Text', 'Copy', ...
                           'MenuSelectedFcn', @(src, event) copyPathToClipboard(app, valText));
                    
                    % 将右键菜单挂载到这个纯净的 Label 上
                    value.ContextMenu = cm;
                else
                    % 不是路径，正常显示标签
                    value = uilabel(app.GridLayout, ...
                        "Text", valText, ...
                        "Tooltip", valText);
                end
            end

        end
        
        function copyPathToClipboard(app, textToCopy)
            % 塞入系统剪贴板 (注意要把 string 转成 char，兼容性更好)
            clipboard('copy', char(textToCopy));
            app.appendOutput('路径已成功复制到剪贴板！');
            % 弹出一个不会中断程序的提示框，1.5秒后用户点不点它都不影响操作
            % uialert(app.UIFigure, '路径已成功复制到剪贴板！', '复制成功', 'Icon', 'success');
        end


        function bringToFront(app)
            % 把主窗口重新拉到最前面，防止被系统对话框踢到后台
            if isvalid(app.UIFigure)
                figure(app.UIFigure);
            end
        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            screenSize = get(0, 'ScreenSize');
            screenWidth = screenSize(3);
            screenHeight = screenSize(4);
            
            app.UIFigure.Position(1) = screenWidth*0.4;
            app.UIFigure.Position(2) = screenHeight*0.3;
            app.UIFigure.Resize = 'off';   
            

           

            addpath(genpath('algorithms'));
            addpath(genpath('core'));
            addpath(genpath('external'));
            addpath(genpath('docImg'));
            addpath(genpath('docs'));
            addpath(genpath('utilities'));
            addpath(genpath('app'));
            app.UIFigure.Name='SEAL';
           

            app.protocolmenu = uicontextmenu(app.UIFigure);
            uimenu(app.protocolmenu, 'Text', 'Import EEG Data...', 'MenuSelectedFcn', @(s,e) ImportDataSelected(app, 'Data'));
            app.importcortexmenu = uimenu(app.protocolmenu, 'Text', 'Import Cortex...', 'MenuSelectedFcn', @(s,e) ImportDataSelected(app, 'Cortex'));
            app.importleadfieldmenu = uimenu(app.protocolmenu, 'Text', 'Import Lead Field...', 'MenuSelectedFcn', @(s,e) ImportDataSelected(app, 'Leadfield'));
            uimenu(app.protocolmenu, 'Text', 'Delete Protocol(s)', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));
            app.protocolmenu.ContextMenuOpeningFcn = @(s,e) ProtocolMenuOpening(app, e);

            app.cortexMenu = uicontextmenu(app.UIFigure);
            uimenu(app.cortexMenu, 'Text', 'Visualize...', 'MenuSelectedFcn', @(s,e) VisualizeDataSelected(app, 'Cortex'));
            uimenu(app.cortexMenu, 'Text', 'Delete', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));

            app.leadfieldMenu = uicontextmenu(app.UIFigure);
            uimenu(app.leadfieldMenu, 'Text', 'Delete', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));

            app.chanlocsMenu = uicontextmenu(app.UIFigure);
            uimenu(app.chanlocsMenu, 'Text', 'Delete', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));

            app.sessionMenu = uicontextmenu(app.UIFigure);
            uimenu(app.sessionMenu, 'Text', 'Extract Session Recipe...', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) ExtractSessionRecipeMenuSelected(app, e));
            uimenu(app.sessionMenu, 'Text', 'Delete', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));

            app.eegDataFileMenu = uicontextmenu(app.UIFigure);
            uimenu(app.eegDataFileMenu, 'Text', 'Visualize Time Series...', 'MenuSelectedFcn', @(s,e) VisualizeDataSelected(app, 'EEG_TimeSeries'));
            app.importchanlocsmenu = uimenu(app.eegDataFileMenu, 'Text', 'Import Channel Location...', 'MenuSelectedFcn', @(s,e) ImportChanlocs(app, e));
            app.visTopoMenu = uimenu(app.eegDataFileMenu, 'Text', 'Visualize Topography...', 'MenuSelectedFcn', @(s,e) VisualizeDataSelected(app, 'EEG_Topo'));
            uimenu(app.eegDataFileMenu, 'Text', 'Run Source Analysis...', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) RunSourceAnalysisSelected(app, e));
            uimenu(app.eegDataFileMenu, 'Text', 'Delete', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));
            app.eegDataFileMenu.ContextMenuOpeningFcn = @(s,e) EEGMenuOpening(app, e);

            app.PreprocessedDataMenu = uicontextmenu(app.UIFigure);
            uimenu(app.PreprocessedDataMenu, 'Text', 'Visualize Time Series...', 'MenuSelectedFcn', @(s,e) VisualizeDataSelected(app, 'EEG_TimeSeries'));
            uimenu(app.PreprocessedDataMenu, 'Text', 'Run Source Analysis...', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) RunSourceAnalysisSelected(app, e));
            uimenu(app.PreprocessedDataMenu, 'Text', 'Delete', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));

            app.SimuParamMenu = uicontextmenu(app.UIFigure);
            uimenu(app.SimuParamMenu, 'Text', 'Delete', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e))

            app.resultFileMenu = uicontextmenu(app.UIFigure);
            uimenu(app.resultFileMenu, 'Text', 'Visualize Results...', 'Separator', 'on','MenuSelectedFcn', @(s,e) VisualizeDataSelected(app, 'Results'));
            uimenu(app.resultFileMenu, 'Text', 'Extract Recipe...', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) ExtractRecipeMenuSelected(app, e));
            uimenu(app.resultFileMenu, 'Text', 'Delete', 'Separator', 'on', 'MenuSelectedFcn', @(s,e) DeleteNodeSelected(app, e));



            app.OutputHistory = "";
            if ispref('SEAL_Toolbox', 'ProjectPath')
                app.ProjectPath = getpref('SEAL_Toolbox', 'ProjectPath');
                try
                    app.projectNode = ProjectNode.openExisting(app.ProjectPath);
                    app.ProjectPanel.Title = strcat('Project:',app.projectNode.name);
                    app.appendOutput('Click the protocolnode to update protocol information!!');
                    app.flushFileTree();
                catch ME
                    app.appendOutput(strcat("Preference """, app.ProjectPath, " "" seems to be illegal."));
                    rmpref('SEAL_Toolbox', 'ProjectPath');
                    app.ProjectPath = "";
                    app.projectNode = ProjectNode.empty();
                    return;
                end

            else
                app.ProjectPath = "";
            end
            app.PreprocessButton.Enable = "off";
            app.VisualizationButton.Enable = "off";
            app.ESIAnalysisButton.Enable = "off";
            app.SimulationButton.Enable = "off";
            app.RunBatchButton.Enable = "off";
            app.EvaluationButton.Enable = "off";

            splashPath = fullfile('docImg', 'background.jpg');  % 你的图片路径
            if exist(splashPath, 'file')
                app.showSplash(splashPath, 3);  % 停留 2 秒
            end
        end

        % Menu selected function: NewprotocolMenu
        function NewprotocolMenuSelected(app, event)
            app.DialogApp_newprotocol=SEAL_Newprotocol(app);
        end

        % Selection changed function: FileTree
        function FileTreeSelectionChanged(app, event)
            %选择的节点变化，信息需要更新
            %需要限制只能选择Protocol节点
            oldNodes = event.PreviousSelectedNodes;
            newNodes = event.SelectedNodes;

            % 取消之前选中的所有节点
            if ~isempty(oldNodes)
                for i = 1:length(oldNodes)
                    node = oldNodes(i);
                    nodeData = node.NodeData;
                    if ~isempty(nodeData) && ismethod(nodeData, 'deselect')
                        try
                            nodeData.deselect();
                            nodeData.unload();
                        catch ME
                            warning('Failed to deselect node: %s', ME.message);
                        end
                    end
                end
            end

            % 选中新节点
            if ~isempty(newNodes)
                for i = 1:length(newNodes)
                    node = newNodes(i);
                    nodeData = node.NodeData;
                    if ~isempty(nodeData) && ismethod(nodeData, 'select')
                        try
                            if isa(nodeData,'DataNode')
                                app.PreprocessButton.Enable = "on";
                                app.VisualizationButton.Enable = "on";
                                app.ESIAnalysisButton.Enable = "on";
                                app.SimulationButton.Enable = "off";
                                
                                app.EvaluationButton.Enable = "off";
                            elseif isa(nodeData,'CortexNode')
                                app.PreprocessButton.Enable = "off";
                                app.VisualizationButton.Enable = "on";
                                app.ESIAnalysisButton.Enable = "off";
                                app.SimulationButton.Enable = "off";
                                
                                app.EvaluationButton.Enable = "off";
                            elseif isa(nodeData,'LeadfieldNode')
                                app.PreprocessButton.Enable = "off";
                                app.VisualizationButton.Enable = "off";
                                app.ESIAnalysisButton.Enable = "off";
                                app.SimulationButton.Enable = "off";
                                
                                app.EvaluationButton.Enable = "off";
                            elseif isa(nodeData,'SessionNode')
                                app.PreprocessButton.Enable = "off";
                                app.VisualizationButton.Enable = "off";
                                app.ESIAnalysisButton.Enable = "off";
                                app.SimulationButton.Enable = "off";
                                
                                app.EvaluationButton.Enable = "off";
                            elseif isa(nodeData,'ProtocolNode')
                                app.PreprocessButton.Enable = "off";
                                app.VisualizationButton.Enable = "off";
                                app.ESIAnalysisButton.Enable = "off";
                                app.SimulationButton.Enable = "on";
                                
                                app.EvaluationButton.Enable = "off";
                            elseif isa(nodeData,'ResultNode')
                                if nodeData.resultInfo.isSimu == 1
                                    app.EvaluationButton.Enable = "on";
                                end
                                app.PreprocessButton.Enable = "off";
                                app.VisualizationButton.Enable = "on";
                                app.ESIAnalysisButton.Enable = "off";
                                app.SimulationButton.Enable = "off";
                               
                            end

                            nodeData.select();
                        catch ME
                            warning('Failed to select node: %s', ME.message);
                        end
                    end
                end
            end

            if numel(newNodes) == 1
                app.flushMetadataGrid(newNodes);
            else
                app.clearMetadataGrid();
            end
        end

        % Menu selected function: NewProjectMenu
        function NewProjectMenuSelected(app, event)
            if isempty(app.ProjectPath)
                defaultPath = pwd; % If no project path exists, use current path as default.
            else
                [defaultPath, ~, ~] = fileparts(app.ProjectPath);

                % If on root path, set it as default.
                if isempty(defaultPath)
                    defaultPath = app.ProjectPath; % Otherwise, set location of last project as default.
                end
            end
            selectedPath = uigetdir(defaultPath, 'Select a location to create the new Project');

            if ~isequal(selectedPath, 0) && ~isempty(strtrim(selectedPath))

                defaultProjectName = findAvailableName(selectedPath);
                projectName = inputdlg('Enter a name for the new project:', 'New Project', [1 50], {defaultProjectName});

                if isempty(projectName)
                    app.appendOutput('Create Project cancelled.');
                    return;
                end

                projectName = strtrim(projectName{1});

                if isempty(projectName)
                    errordlg('Project name cannot be empty.', 'Invalid Project Name');
                    return;
                end

                newProjectPath = fullfile(selectedPath, projectName);

                if exist(newProjectPath, 'dir')
                    errordlg('A project with this name already exists in that location.', 'Project Exists');
                    return;
                end

                app.projectNode = ProjectNode.createNew(projectName, selectedPath);

                app.ProjectPanel.Title = ['Project: ', projectName];


                app.ProjectPath = newProjectPath;
                setpref('SEAL_Toolbox', 'ProjectPath', newProjectPath);
                app.bringToFront();
                app.appendOutput(['Create directory for Project: ', newProjectPath]);
            else
                app.appendOutput('Create Project cancelled.');
            end
        end

        % Menu selected function: ImportprotocolMenu
        function ImportprotocolMenuSelected(app, event)

            % Let user select a protocol folder from anywhere
            sourceProtocolPath = uigetdir(app.ProjectPath, 'Select Protocol Folder to Import');
            if sourceProtocolPath == 0, return; end

            % Validate it's a real protocol
            if ~exist(fullfile(sourceProtocolPath, 'protocol_info.mat'), 'file')
                uialert(app.UIFigure, 'The selected folder is not a valid SEAL protocol.', 'Invalid Folder');
                return;
            end

            % Check for naming conflicts in the current project
            [~, protocolName] = fileparts(sourceProtocolPath);
            targetProtocolPath = fullfile(app.ProjectPath, protocolName);

            protocol = app.projectNode.openProtocol(sourceProtocolPath);

            app.flushFileTree();
            app.appendOutput(['Successfully import protocol: ', protocolName]);
            app.bringToFront();
        end

        % Menu selected function: LoadProjectMenu
        function LoadProjectMenuSelected(app, event)
            % Let user select an existing project folder
            selectedPath = uigetdir(app.ProjectPath, 'Select Project Folder');
            if selectedPath == 0
                app.appendOutput(['Load Project cancelled']);
                return;
            end

            app.projectNode = ProjectNode.openExisting(selectedPath);
            %
            %             [~, projectName] = fileparts(selectedPath);
            app.ProjectPanel.Title = strcat("Project: ", app.projectNode.name);
            setpref('SEAL_Toolbox', 'ProjectPath', selectedPath);
            app.ProjectPath = app.projectNode.path;
            %app = scanAndLoadProtocols(app);
            %             app.scanAndLoadProtocols();
            app.appendOutput(strcat('Project loaded: ', selectedPath));
            app.bringToFront();
        end

        % Button pushed function: PreprocessButton
        function PreprocessButtonPushed(app, event)
            app.DialogApp_preprocess=SEAL_Preprocess(app);
        end

        % Menu selected function: CopyMenu
        function CopyMenuSelected(app, event)
            textToCopy = app.ProtocolPathValue.Text;

            % 复制到剪贴板
            clipboard('copy', textToCopy);

            % 提示用户
            app.appendOutput(['Copied text：' textToCopy]);
        end

        % Node expanded function: FileTree
        function FileTreeNodeExpanded(app, event)
            nodeName = string(event.Node.Text);
            app.ExpandedNodeNames(end+1) = nodeName;
            app.ExpandedNodeNames = unique(app.ExpandedNodeNames); % 去重，防止重复记录
        end

        % Node collapsed function: FileTree
        function FileTreeNodeCollapsed(app, event)
            nodeName = string(event.Node.Text);
            app.ExpandedNodeNames(app.ExpandedNodeNames == nodeName) = [];
        end

        % Button pushed function: ESIAnalysisButton
        function ESIAnalysisButtonPushed(app, event)
            selectedNode = app.FileTree.SelectedNodes;
            if isempty(selectedNode)
                uialert(app.UIFigure, '运行失败！请先在左侧树状图中选中一个【脑电数据节点】。', ...
                    '未选择数据', 'Icon', 'warning');
                return; % 直接终止后续代码，子窗口不会弹出
            end
            
            % 检查二：用户选错节点类型了吗？(使用 isa 函数检查是否为 DataNode)
            if ~isa(selectedNode.NodeData, 'DataNode')
                % 如果选中的是 ProtocolNode 或其他，给出友好的中文提示
                uialert(app.UIFigure, sprintf('运行源分析需要脑电数据！\n您当前选中的是 [%s] 类型的节点。\n\n请点选带有波形图标的【DataNode】后再重试。', class(selectedNode)), ...
                    '节点类型错误', 'Icon', 'error');
                return; % 拦截！
            end
            app.RunSourceAnalysisSelected();
        end

        % Button pushed function: RunBatchButton
        function RunBatchButtonPushed(app, event)
            if isempty(app.CurrentTaskRecipe)
                uialert(app.UIFigure, 'Select the result node and extract recipe.！', ...
                    '缺少批处理模板', 'Icon', 'warning');
                return
            end

            parentFolder = uigetdir('', '请选择包含所有被试(Session)的批处理总文件夹');

            if isequal(parentFolder, 0)
                return;
            end

            numSteps = length(app.CurrentTaskRecipe);

            if numSteps > 1 
                % =====================================================
                % 情况 A：多步管线配方 (循环调用)
                % =====================================================
                disp(['[批处理准备] 检测到：多步管线配方，共需执行 ', num2str(numSteps), ' 步算法。']);

                % 禁用按钮，防止用户在批处理期间乱点
                app.RunBatchButton.Enable = 'off';

                try
                    disp(['[批处理准备] 共检测到 ', num2str(length(app.CurrentTaskRecipe)), ' 步算法。']);

                    % 把整个配方数组直接传给底层函数！
                    seal_runFolderBasedBatch(app, parentFolder, app.CurrentTaskRecipe);

                    uialert(app.UIFigure, '所有被试已全部批处理完成！', '批处理成功', 'Icon', 'success');
                catch ME
                    errormsg = sprintf('批处理出错:\n%s', ME.message);
                    uialert(app.UIFigure, errormsg, '批处理中断', 'Icon', 'error');
                end

                % 恢复按钮状态
                app.RunBatchButton.Enable = 'on';
                app.bringToFront();
            else
                % =====================================================
                % 情况 B：单步算法配方 (直接调用)
                % =====================================================
                disp('[批处理准备] 检测到：单步算法配方。');
                seal_runFolderBasedBatch(app, parentFolder, app.CurrentTaskRecipe);
                uialert(app.UIFigure, '批处理完成！', '成功', 'Icon', 'success');
            end

        end

        % Button pushed function: SimulationButton
        function SimulationButtonPushed(app, event)
            selectedNode = app.FileTree.SelectedNodes;
            selectedNode = selectedNode.NodeData;
            if length(selectedNode)~=1
                app.appendOutput('Please choose a ProtocolNode');
                return
            else
                if ~isa(selectedNode,'ProtocolNode')
                    app.appendOutput('Please choose a ProtocolNode');
                    return
                else
                    childrenNodes = selectedNode.children;
               
                    % 2. 查找目标节点的位置索引
                    idxCortex = find(arrayfun(@(x) isa(x, 'CortexNode'), childrenNodes), 1);
                    idxLeadfield = find(arrayfun(@(x) isa(x, 'LeadfieldNode'), childrenNodes), 1);

                    % 3. 判断是否找到，并提取节点句柄
                    if isempty(idxCortex)
                        app.appendOutput('Cannot find a CortexNode');
                        return
                    end

                    if isempty(idxLeadfield)
                        app.appendOutput('Cannot find a LeadfieldNode');
                        return
                    end

                    app.DialogApp_simulation = SEAL_Simulation(app, selectedNode);
                end
            end

        end

        % Menu selected function: ClearMenu
        function ClearMenuSelected(app, event)
            app.TextArea.Value = '';
        end

        % Button pushed function: VisualizationButton
        function VisualizationButtonPushed(app, event)
            selectedNode = app.FileTree.SelectedNodes;
            selectedNode = selectedNode.NodeData;
            if length(selectedNode)~=1
                app.appendOutput('Please choose ONE Node');
                return
            else
                if isa(selectedNode,'ProtocolNode') || isa(selectedNode,'LeadfieldNode') || isa(selectedNode,'SessionNode')
                    app.appendOutput('Please choose a ResultNode\DataNode\CortexNode');
                    return
                elseif isa(selectedNode,'ResultNode')
                    vizType ='Results';
                elseif isa(selectedNode,'DataNode')
                    vizType ='EEG_TimeSeries';
                elseif isa(selectedNode,'CortexNode')
                    vizType='Cortex';
                end
            end
            app.VisualizeDataSelected(vizType);
        end

        % Button pushed function: EvaluationButton
        function EvaluationButtonPushed(app, event)
            node = app.FileTree.SelectedNodes;
            ResultNode = node.NodeData;
            protocolNode = ResultNode.parent.parent;
            cortexNode = protocolNode.cortexNode;
            data = load(ResultNode.resultInfo.SimuSourcePath);
            SimuInfo = data.data;

            s = ResultNode.sourceActivity;
            active_indices= SimuInfo.activeIndices;
            seedvox = SimuInfo.SourceData.SeedVertices;
            s_energy = sqrt(sum(s.^2, 2));   % 用 RMS 而非能量
            true_idx = unique(cell2mat(active_indices(:)));  % 显式去重并保证列向量
            disttype = 'euclidean'; % 'geodesic'
            threshold = 0.1;
            metrics_to_compute = {'DLE','SD','EMD','ROC_AUC','PR_AUC','F1','Precision','Recall','APrime'};
            metric_args = {'EMDMaxSources', 1000};
            has_true_source = isfield(SimuInfo, 's_real') && ~isempty(SimuInfo.s_real);

            if has_true_source
                metric_args = [metric_args, {'EMDTrueDistribution', sqrt(sum(SimuInfo.s_real.^2, 2))}];
            end

            if has_true_source && ~isvector(s)
                metrics_to_compute{end+1} = 'MSE';
                metric_args = [metric_args, {'EstimatedSourceTimeSeries', s, ...
                    'TrueSourceTimeSeries', SimuInfo.s_real, ...
                    'TemporalRoiIndices', active_indices}];
            else
                warning('SEAL:Evaluation:MissingTrueSourceTimeSeries', ...
                    '未找到 s_real 或估计结果不是时序矩阵，跳过 temporal MSE。');
            end

            spatial_metric = seal_evaluate_spatial_metrics( ...
                s_energy, true_idx, seedvox, cortexNode.data , ...
                'DistanceType', disttype, ...
                'ThresholdBinary', threshold, ...   % 建议 0.1~0.5 之间或 'otsu'
                'MetricsToCompute', metrics_to_compute, ...
                metric_args{:});
            seal_print_spatial_metrics(spatial_metric);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 719 558];
            app.UIFigure.Name = 'MATLAB App';

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create NewProjectMenu
            app.NewProjectMenu = uimenu(app.FileMenu);
            app.NewProjectMenu.MenuSelectedFcn = createCallbackFcn(app, @NewProjectMenuSelected, true);
            app.NewProjectMenu.Text = 'New Project';

            % Create LoadProjectMenu
            app.LoadProjectMenu = uimenu(app.FileMenu);
            app.LoadProjectMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadProjectMenuSelected, true);
            app.LoadProjectMenu.Text = 'Load Project';

            % Create NewprotocolMenu
            app.NewprotocolMenu = uimenu(app.FileMenu);
            app.NewprotocolMenu.MenuSelectedFcn = createCallbackFcn(app, @NewprotocolMenuSelected, true);
            app.NewprotocolMenu.Text = 'New protocol';

            % Create ImportprotocolMenu
            app.ImportprotocolMenu = uimenu(app.FileMenu);
            app.ImportprotocolMenu.MenuSelectedFcn = createCallbackFcn(app, @ImportprotocolMenuSelected, true);
            app.ImportprotocolMenu.Text = 'Import protocol';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.UIFigure);
            app.HelpMenu.Text = 'Help';

            % Create MetadataPanel
            app.MetadataPanel = uipanel(app.UIFigure);
            app.MetadataPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.MetadataPanel.Title = 'Metadata';
            app.MetadataPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.MetadataPanel.FontWeight = 'bold';
            app.MetadataPanel.Position = [285 291 435 268];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.MetadataPanel);
            app.GridLayout.ColumnWidth = {'1x', '3x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create ProjectPanel
            app.ProjectPanel = uipanel(app.UIFigure);
            app.ProjectPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.ProjectPanel.Title = 'Project';
            app.ProjectPanel.FontWeight = 'bold';
            app.ProjectPanel.Position = [1 191 285 368];

            % Create FileTree
            app.FileTree = uitree(app.ProjectPanel);
            app.FileTree.Multiselect = 'on';
            app.FileTree.SelectionChangedFcn = createCallbackFcn(app, @FileTreeSelectionChanged, true);
            app.FileTree.NodeExpandedFcn = createCallbackFcn(app, @FileTreeNodeExpanded, true);
            app.FileTree.NodeCollapsedFcn = createCallbackFcn(app, @FileTreeNodeCollapsed, true);
            app.FileTree.Position = [0 0 284 346];

            % Create StatusPanel
            app.StatusPanel = uipanel(app.UIFigure);
            app.StatusPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.StatusPanel.Title = 'Status';
            app.StatusPanel.FontWeight = 'bold';
            app.StatusPanel.Position = [1 1 719 191];

            % Create TextArea
            app.TextArea = uitextarea(app.StatusPanel);
            app.TextArea.Position = [1 -3 719 174];

            % Create OperationsPanel
            app.OperationsPanel = uipanel(app.UIFigure);
            app.OperationsPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.OperationsPanel.Title = 'Operations';
            app.OperationsPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.OperationsPanel.FontWeight = 'bold';
            app.OperationsPanel.Position = [285 191 435 101];

            % Create PreprocessButton
            app.PreprocessButton = uibutton(app.OperationsPanel, 'push');
            app.PreprocessButton.ButtonPushedFcn = createCallbackFcn(app, @PreprocessButtonPushed, true);
            app.PreprocessButton.Position = [37 47 110 24];
            app.PreprocessButton.Text = 'Preprocess';

            % Create VisualizationButton
            app.VisualizationButton = uibutton(app.OperationsPanel, 'push');
            app.VisualizationButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizationButtonPushed, true);
            app.VisualizationButton.Position = [169 47 115 24];
            app.VisualizationButton.Text = 'Visualization';

            % Create ESIAnalysisButton
            app.ESIAnalysisButton = uibutton(app.OperationsPanel, 'push');
            app.ESIAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @ESIAnalysisButtonPushed, true);
            app.ESIAnalysisButton.FontWeight = 'bold';
            app.ESIAnalysisButton.Position = [305 47 109 24];
            app.ESIAnalysisButton.Text = 'ESI Analysis';

            % Create SimulationButton
            app.SimulationButton = uibutton(app.OperationsPanel, 'push');
            app.SimulationButton.ButtonPushedFcn = createCallbackFcn(app, @SimulationButtonPushed, true);
            app.SimulationButton.Position = [38 14 110 24];
            app.SimulationButton.Text = 'Simulation';

            % Create EvaluationButton
            app.EvaluationButton = uibutton(app.OperationsPanel, 'push');
            app.EvaluationButton.ButtonPushedFcn = createCallbackFcn(app, @EvaluationButtonPushed, true);
            app.EvaluationButton.Position = [170 14 114 24];
            app.EvaluationButton.Text = 'Evaluation';

            % Create RunBatchButton
            app.RunBatchButton = uibutton(app.OperationsPanel, 'push');
            app.RunBatchButton.ButtonPushedFcn = createCallbackFcn(app, @RunBatchButtonPushed, true);
            app.RunBatchButton.Position = [306 14 108 24];
            app.RunBatchButton.Text = 'Run Batch ';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create CopyMenu
            app.CopyMenu = uimenu(app.ContextMenu);
            app.CopyMenu.MenuSelectedFcn = createCallbackFcn(app, @CopyMenuSelected, true);
            app.CopyMenu.Text = 'Copy';

            % Create ContextMenu2
            app.ContextMenu2 = uicontextmenu(app.UIFigure);

            % Create ClearMenu
            app.ClearMenu = uimenu(app.ContextMenu2);
            app.ClearMenu.MenuSelectedFcn = createCallbackFcn(app, @ClearMenuSelected, true);
            app.ClearMenu.Text = 'Clear';
            
            % Assign app.ContextMenu2
            app.TextArea.ContextMenu = app.ContextMenu2;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SEAL_GUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
