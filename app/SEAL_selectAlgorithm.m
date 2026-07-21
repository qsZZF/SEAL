classdef SEAL_selectAlgorithm < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        LoadPriorButton                matlab.ui.control.Button
        SelectPriorButton              matlab.ui.control.Button
        GridLayout                     matlab.ui.container.GridLayout
        GeneralSettingsPanel           matlab.ui.container.Panel
        NoiseWhiteningDropDown         matlab.ui.control.DropDown
        NoiseWhiteningDropDown_2Label  matlab.ui.control.Label
        RunSelectedAlgorithmsButton    matlab.ui.control.Button
        DepthWeightingEditField        matlab.ui.control.NumericEditField
        DepthWeightingEditFieldLabel   matlab.ui.control.Label
        AlgorithmsPanel                matlab.ui.container.Panel
        AlgorithmsTree                 matlab.ui.container.CheckBoxTree
        MinimumL2_NormNode             matlab.ui.container.TreeNode
        MNENode                        matlab.ui.container.TreeNode
        LOGURANode                     matlab.ui.container.TreeNode
        dSPMNode                       matlab.ui.container.TreeNode
        sLORETANode                    matlab.ui.container.TreeNode
        eLORETANode                    matlab.ui.container.TreeNode
        LORETANode                     matlab.ui.container.TreeNode
        LAURANode                      matlab.ui.container.TreeNode
        BeamformerNode                 matlab.ui.container.TreeNode
        LCMVNode                       matlab.ui.container.TreeNode
        MCMVNode                       matlab.ui.container.TreeNode
        SparsemodelNode                matlab.ui.container.TreeNode
        Lq_ADMMNode                    matlab.ui.container.TreeNode
        ChampagneNode                  matlab.ui.container.TreeNode
        BlkChamNode                    matlab.ui.container.TreeNode
        TSChamNode                     matlab.ui.container.TreeNode
        SmoothChamNode                 matlab.ui.container.TreeNode
        IRESNode                       matlab.ui.container.TreeNode
        FOCUSSNode                     matlab.ui.container.TreeNode
        SpatioTemporalNode             matlab.ui.container.TreeNode
        MxNENode                       matlab.ui.container.TreeNode
        irMxNENode                     matlab.ui.container.TreeNode
        SISSY_L21Node                  matlab.ui.container.TreeNode
        VSSI_LpNode                    matlab.ui.container.TreeNode
        TFMxNENode                     matlab.ui.container.TreeNode
        DebiasedGroupLassoNode         matlab.ui.container.TreeNode
        dSPNNode                       matlab.ui.container.TreeNode
        FastIRESNode                   matlab.ui.container.TreeNode
        STOUTNode                      matlab.ui.container.TreeNode
        dMAPNode                       matlab.ui.container.TreeNode
        BESTIESNode                    matlab.ui.container.TreeNode
        SISTBFNode                     matlab.ui.container.TreeNode
        uSTARNode                      matlab.ui.container.TreeNode
        STARTSNode                     matlab.ui.container.TreeNode
        sSTARTSNode                    matlab.ui.container.TreeNode
        VSSI_GGDNode                   matlab.ui.container.TreeNode
        ContextMenu                    matlab.ui.container.ContextMenu
        ImportDataMenu                 matlab.ui.container.Menu
    end


    properties (Access = private)
        protocolinfo = {}
        Data
        Cortex = {}
        LeadField = {}
        LeadField_fixed = {}
        DialogApp_datacheck
        AdditionalDataNode = {}
        AdditionalDataNodeNum = 0
        AdditionalData={}
        verticalSpacing = 9
        horizontalSpacing = 9
        allPanels = {}
        originalheight
        originalwidth
        TrialIndices
        MinimumL2_NormUIFigure
        BeamformerUIFigure
        SpatioUIFigure
        NumOrientations=3
        DepthWeightExp=0.5
        Callingapp
        CurrentID
        SaveIndex
        protocolNode ProtocolNode
        sessionNode SessionNode
        dataNode DataNode
        resultdata=struct();
        counters=zeros(50,1);
        trials = 1
        SavePath
        srate
        PanelCache            % containers.Map: algName -> panel handle
        AlgParams             % containers.Map: algName -> struct of params
        pendingAlgNames = {}
        NoiseStrategy = 'brainstorm'
    end

    properties (Access = public)
        MNESettingsPanel %自定义控件
        MNEsourceCovarianceEditField
        MNEsourceCovarianceLabel
        MNERegularizationParameterEditField
        MNERegularizationParameterLabel
        RegularizationParameter
        sourceCovariance
        dSPMSettingsPanel
        sLORETASettingsPanel
        eLORETASettingsPanel
        LORETASettingsPanel
        LAURASettingsPanel
        STARTSSettingsPanel
        STARTSNumTBFsEditField
        STARTSNumTBFsLabel
        NumTBFs = 4
        STARTSNeiborOrderEditField
        STARTSNeiborOrderLabel
        NeighborOrder = 2
        STARTLearnNoiseDropDown
        STARTSLearnNoiseLabel
        LearnNoise = false
        STARTNoiseObservationDropDown
        STARTSNoiseObservationLabel
        NoiseObservationMatrix
        SparseModeUIFigure
        blockChampagneSettingsPanel
        blockSNeiborOrderEditField
        blockSNeiborOrderLabel
        blockLearnNoiseDropDown
        blockSLearnNoiseLabel
        blockNoiseObservationDropDown
        blockNoiseObservationLabel
        blockSAtlasConstraintsEditField
        blockSAtlasConstraintsLabel
        uSTARSettingsPanel
        uSTARNumTBFsEditField
        uSTARNumTBFsLabel
        uSTARMicrostateMapsLabel
        uSTARMicrostateMapDropDown
        MicrostateMaps = []
        uSTARMicrostateOptionsLabel
        uSTARRunAlgorithmDropDown
        uSTARRunAlgorithmLabel
        uSTARRunAlgorithm= 'STAR'
        uSTARPriorTypeDropDown
        uSTARPriorTypeLabel
        importMicrostateOptionsButton
        importMicrostateMapsButton
        firsttime = 1

    end


    methods (Access = private)
      function applyLayout(app, panel, allvisiblePanels)
    visiblePanels = allvisiblePanels;   % 直接用，不需要再 for 循环复制
    if isempty(visiblePanels), return; end
    
    panelWidth = visiblePanels{1}.Position(3);
    panelHeights = cellfun(@(p) p.Position(4), visiblePanels);
    
    % --- 关键改动 1：figure 的可用高度 = max(原始高度, 最高 panel 高度 + 上下留白) ---
    % 这样单个超高 panel 也能塞进一列而不会溢出
    minFigureHeight = app.originalheight;
    tallestPanel    = max(panelHeights);
    availableHeight = max(minFigureHeight, tallestPanel + 2*app.verticalSpacing);
    
    % 按高度降序排序（大的在左）
    [sortedHeights, sortIdx] = sort(panelHeights, 'descend');
    visiblePanels = visiblePanels(sortIdx);
    
    % 分列
    columnHeights = zeros(1, length(visiblePanels));
    columns = {};
    currentColumn = 1;
    columns{currentColumn} = [];
    
    for i = 1:length(visiblePanels)
        % --- 关键改动 2：用动态 availableHeight 而不是 originalheight ---
        if columnHeights(currentColumn) + sortedHeights(i) <= availableHeight - 2*app.verticalSpacing
            columns{currentColumn}(end+1) = i;
            columnHeights(currentColumn) = columnHeights(currentColumn) + sortedHeights(i) + app.verticalSpacing;
        else
            currentColumn = currentColumn + 1;
            columns{currentColumn} = i;
            columnHeights(currentColumn) = sortedHeights(i) + app.verticalSpacing;
        end
    end
    
    columnHeights = columnHeights + app.verticalSpacing;
    
    % 剔除空列
    for i = length(columnHeights):-1:1
        if columnHeights(i) == app.verticalSpacing
            columnHeights(i) = [];
        end
    end
    
    % 按列高度降序（高列在左）
    [~, colOrder] = sort(columnHeights, 'descend');
    columns = columns(colOrder);
    
    newWidth = (panelWidth + app.horizontalSpacing) * length(columns) + app.horizontalSpacing;
    
    % --- 关键改动 3：同时设置 figure 的宽度和高度 ---
    % 保持 figure 左下角位置不变，只改尺寸
    oldPos = panel.Position;
    panel.Position = [oldPos(1) oldPos(2) newWidth availableHeight];
    
    % 布局 panel
    currentX = app.horizontalSpacing;
    for colIdx = 1:length(columns)
        colPanels = columns{colIdx};
        % --- 关键改动 4：从 figure 顶部往下排，基于 availableHeight ---
        currentY = availableHeight - app.verticalSpacing;
        
        for panelIdx = colPanels
            p = visiblePanels{panelIdx};
            pHeight = p.Position(4);
            p.Position = [currentX currentY-pHeight panelWidth pHeight];
            currentY = currentY - pHeight - app.verticalSpacing;
        end
        
        currentX = currentX + panelWidth + app.horizontalSpacing;
    end
    
    panel.Visible = 'on';
end

        function MNEsourceCovarianceValueChanged(app, ~)
            app.sourceCovariance = app.MNEsourceCovarianceEditField.value;
        end

        function MNERegularizationParameterEditFieldValueChanged(app, ~)
            app.RegularizationParameter = app.MNERegularizationParameterEditField.value;
        end

        function closeMinL2FigureCallback(app, src, ~)
            delete(src);
            % 4. 重置属性
            app.MinimumL2_NormUIFigure = [];
        end

        function NodeDatas = getOnlyLeafNodes(app)
            % 获取所有被勾选的节点
            checkedNodes = app.AlgorithmsTree.CheckedNodes;

            % 如果没有任何节点被勾选，直接返回空字符串数组
            if isempty(checkedNodes)
                NodeDatas = strings(0);
                return;
            end

            % 1. 核心判断逻辑：找出所有“叶子节点”
            % isempty(n.Children) 如果为 true，说明它没有下级，是最底层的节点
            isLeaf = arrayfun(@(n) isempty(n.Children), checkedNodes);

            % 2. 过滤掉所有根节点和中间节点，只保留叶子节点
            leafNodes = checkedNodes(isLeaf);

            % 3. 提取保留下来的 NodeData
            NodeDatas = strings(0); % 初始化空字符串数组
            count = 1;

            for i = 1:length(leafNodes)
                if ~isempty(leafNodes(i).NodeData)
                    % 强制转换为 string 存入数组，防止 char 数组拼接报错
                    NodeDatas(count) = string(leafNodes(i).NodeData);
                    count = count + 1;
                end
            end
        end
        function results = deal_results(app ,Source)%  3V*T 数据降维 -> V*T
            temp = reshape(Source, 3, [], size(Source, 2));
            results = squeeze(sqrt(sum(temp.^2, 1)));    % 沿第一维（3行）计算平方和
        end

        function STARTSNumTBFsEditFieldValueChanged(app,~)
            app.NumTBFs=app.STARTSNumTBFsEditField.Value;
        end

        function STARTSNeiborOrderEditFieldValueChanged(app,~)
            app.NeighborOrder=app.STARTSNeiborOrderEditField.Value;
        end

        function STARTLearnNoiseDropDownValueChanged(app,~)
            app.LearnNoise=app.STARTLearnNoiseDropDown.Value;
        end

        function getIndex(app,filelist)
            for i=1:length(filelist)
                if contains(filelist(i).name, 'MNE')
                    app.counters(1)=app.counters(1)+1;
                elseif contains(filelist(i).name, 'STARTS')
                    app.counters(7)=app.counters(7)+1;
                elseif contains(filelist(i).name, 'BlockChampagne')
                    app.counters(8)=app.counters(8)+1;
                elseif contains(filelist(i).name, 'uSTAR')
                    app.counters(9)=app.counters(9)+1;
                end
                %为每一个都更新参数
            end
        end

        function uSTARRunAlgorithmDropDownValueChanged(app,~)
            app.uSTARRunAlgorithm=app.uSTARRunAlgorithmDropDown.Value;
        end

        function importMicrostateMapsButtonPushed(app,~)
            filterSpec = { ...
                '*.mat', 'MATLAB (.mat)'; ...
                '*.set', 'EEGLAB (.set)'; ...
                '*.*',   'All Files (*.*)'};

            [fileName, filePath,~] = uigetfile(filterSpec, ['Select ', strrep('Data', '_', ' '), ' File']);

            originalPath = fullfile(filePath, fileName);
            app.MicrostateMaps=load(originalPath);
            app.MicrostateMaps=struct2cell(app.MicrostateMaps);
            app.MicrostateMaps=app.MicrostateMaps{1};
            app.uSTARMicrostateMapDropDown.Items={fileName};
            app.uSTARMicrostateMapDropDown.Value=fileName;
        end
    
        function concatenated_data = concatenate_trials(app, data_3d)

            % 获取数组尺寸
            [n_channels, n_points_per_trial, n_trials] = size(data_3d);

            % 预分配拼接后的二维矩阵
            total_points = n_points_per_trial * n_trials;
            concatenated_data = zeros(n_channels, total_points);

            % 循环处理每个试次
            for trial = 1:n_trials
                % 计算当前试次的起始和结束索引
                start_idx = (trial - 1) * n_points_per_trial + 1;
                end_idx = trial * n_points_per_trial;

                % 提取当前试次数据
                current_trial_data = data_3d(:, :, trial);

                % 将数据放入拼接矩阵
                concatenated_data(:, start_idx:end_idx) = current_trial_data;


            end

        end


        function closeAllWindows(app)
            if isvalid(app.MinimumL2_NormUIFigure), delete(app.MinimumL2_NormUIFigure); end
            if isvalid(app.SpatioUIFigure), delete(app.SpatioUIFigure); end
            if isvalid(app.SparseModeUIFigure), delete(app.SparseModeUIFigure); end
            if isvalid(app.BeamformerUIFigure), delete(app.BeamformerUIFigure); end
            delete(app.UIFigure);
        end


        function p = buildPanelFromConfig(app, algName, cfg)
            switch cfg.parent
                case 'minL2',  parentFig = app.MinimumL2_NormUIFigure;
                case 'spatio', parentFig = app.SpatioUIFigure;
                case 'sparse', parentFig = app.SparseModeUIFigure;
                case 'beamformer', parentFig = app.BeamformerUIFigure;
            end

            nParams = numel(cfg.params);

            % --- 布局常量 ---
            rowH      = 50;   % 每行高度（label + ctrl + 间距）
            topPad    = 15;   % 距 panel 顶部（标题下方）的额外留白
            bottomPad = 15;   % 距 panel 底部的留白
            titleH    = 20;   % panel 标题占用高度

            % panel 总高度 = 标题 + 顶部留白 + 所有行 + 底部留白
            panelH = max(80, titleH + topPad + nParams*rowH + bottomPad);

            p = uipanel(parentFig, ...
                'Title', cfg.title, ...
                'FontWeight','bold', ...
                'Position',[10 10 174 panelH], ...
                'Visible','on', ...
                'Tag', algName);

            handles = struct();

            for i = 1:nParams
                spec = cfg.params{i};
                fieldName = spec{1};
                labelText = spec{2};
                ctrlType  = spec{3};
                default   = spec{4};

                currentVal = app.getAlgParam(algName, fieldName, default);

                % 关键：从顶部往下排，留出 topPad
                % panel 内部坐标系原点在左下角，所以 y 要从 (panelH - titleH - topPad) 开始减
                rowTop = panelH - titleH - topPad - (i-1)*rowH;
                yLabel = rowTop - 22;          % label 高 22
                yCtrl  = yLabel - 25;          % ctrl 在 label 下方

                uilabel(p,'Position',[5 yLabel 150 22],'Text',labelText);

                switch ctrlType
                    case 'numeric'
                        h = uieditfield(p,'numeric','Position',[10 yCtrl 150 22], ...
                            'Value',currentVal);
                    case 'edit'
                        h = uieditfield(p,'Position',[10 yCtrl 150 22], ...
                            'Value',currentVal);
                    case 'dropdown'
                        opts = spec{5};
                        h = uidropdown(p,'Position',[10 yCtrl 150 22], ...
                            'Items',opts,'Value',currentVal);
                end
                app.setAlgParam(algName, fieldName, currentVal);
                h.ValueChangedFcn = @(src,~) app.setAlgParam(algName, fieldName, src.Value);
                handles.(fieldName) = h;
            end

            if isfield(cfg,'customBuilder') && ~isempty(cfg.customBuilder)
                feval(cfg.customBuilder, app, p, handles);
            end

            p.UserData = handles;
        end

        % 通用 getter / setter
        function v = getAlgParam(app, algName, fieldName, defaultVal)
            if isKey(app.AlgParams, algName)
                s = app.AlgParams(algName);
                if isfield(s, fieldName), v = s.(fieldName); return; end
            end
            v = defaultVal;
        end

        function setAlgParam(app, algName, fieldName, value)
            if isKey(app.AlgParams, algName), s = app.AlgParams(algName); else, s = struct(); end
            s.(fieldName) = value;
            app.AlgParams(algName) = s;
        end

        function prm = fillAlgDefaults(app, algName, prm)
            if nargin < 3 || isempty(prm) || ~isstruct(prm)
                prm = struct();
            end
            reg = seal_getAlgorithmRegistry();
            registryName = algName;
            if strcmpi(registryName, 'Smooth-Cham')
                registryName = 'SmoothChampagne';
            end
            if ~isfield(reg, registryName) || ~isfield(reg.(registryName), 'params')
                return;
            end

            specs = reg.(registryName).params;
            for ii = 1:numel(specs)
                fieldName = specs{ii}{1};
                defaultValue = specs{ii}{4};
                if ~isfield(prm, fieldName) || isempty(prm.(fieldName))
                    prm.(fieldName) = defaultValue;
                end
            end
        end

        function relayoutAll(app)
            figs = {app.MinimumL2_NormUIFigure, app.SpatioUIFigure, app.SparseModeUIFigure, app.BeamformerUIFigure};
            for i = 1:numel(figs)
                f = figs{i};
                if ~isvalid(f), continue; end
                kids = f.Children;
                panels = num2cell(kids(arrayfun(@(c) isa(c,'matlab.ui.container.Panel'), kids)));
                if isempty(panels)
                    f.Visible = 'off';        % 没有 panel 就藏起来
                else
                    app.applyLayout(f, panels);  % 内部会设回 'on'
                end
            end
        end


        function ensureSettingFiguresExist(app)
            % 三个 settings figure 只在第一次（或被关闭后）才创建

            % 计算停靠位置：主窗口右侧
            anchorX = app.UIFigure.Position(1) + app.UIFigure.Position(3);
            anchorY = app.UIFigure.Position(2);
            figW = 563;
            figH = 264;

            % --- Minimum L2_Norm Settings ---
            if isempty(app.MinimumL2_NormUIFigure) || ~isvalid(app.MinimumL2_NormUIFigure)
                app.MinimumL2_NormUIFigure = uifigure( ...
                    'Visible','off', ...
                    'AutoResizeChildren','off', ...
                    'Position',[anchorX anchorY figW figH], ...
                    'Name','MinimumL2_Norm Settings', ...
                    'CloseRequestFcn', @(s,e) app.closeSettingFigureCallback(s,'minL2'));
            end

            % --- Spatio Temporal Settings ---
            if isempty(app.SpatioUIFigure) || ~isvalid(app.SpatioUIFigure)
                app.SpatioUIFigure = uifigure( ...
                    'Visible','off', ...
                    'AutoResizeChildren','off', ...
                    'Position',[anchorX anchorY figW figH], ...
                    'Name','Spatio Temporal Settings', ...
                    'CloseRequestFcn', @(s,e) app.closeSettingFigureCallback(s,'spatio'));
            end

            % --- SparseMode Settings ---
            if isempty(app.SparseModeUIFigure) || ~isvalid(app.SparseModeUIFigure)
                app.SparseModeUIFigure = uifigure( ...
                    'Visible','off', ...
                    'AutoResizeChildren','off', ...
                    'Position',[anchorX anchorY figW figH], ...
                    'Name','SparseMode Settings', ...
                    'CloseRequestFcn', @(s,e) app.closeSettingFigureCallback(s,'sparse'));
            end

            if isempty(app.BeamformerUIFigure) || ~isvalid(app.BeamformerUIFigure)
                app.BeamformerUIFigure = uifigure( ...
                    'Visible','off', ...
                    'AutoResizeChildren','off', ...
                    'Position',[anchorX anchorY figW figH], ...
                    'Name','Beamformer Settings', ...
                    'CloseRequestFcn', @(s,e) app.closeSettingFigureCallback(s,'beamformer'));
            end
        end


        function closeSettingFigureCallback(app, src, which)
            % 用户手动关闭 settings figure 时触发
            % 1) 销毁 figure；2) 清空对应属性；3) 清理 PanelCache 里指向已销毁 panel 的 key

            delete(src);
            switch which
                case 'minL2',  app.MinimumL2_NormUIFigure = [];
                case 'spatio', app.SpatioUIFigure = [];
                case 'sparse', app.SparseModeUIFigure = [];
                case 'beamformer', app.BeamformerUIFigure = [];
            end

            % 清理失效缓存
            if ~isempty(app.PanelCache) && app.PanelCache.Count > 0
                ks = keys(app.PanelCache);
                for i = 1:numel(ks)
                    if ~isvalid(app.PanelCache(ks{i}))
                        remove(app.PanelCache, ks{i});
                    end
                end
            end
        end


        function out = postProcess(app, S, nO)
            % 3 方向分量降维到 1
            if nO == 1
                out = S;
            else
                out = app.deal_results(S);
            end
        end

        function timeVector = getCurrentDataTimeVector(app, nTime)
            timeVector = [];
            if isempty(app.dataNode) || isempty(app.dataNode.dataInfo) || ~isprop(app.dataNode.dataInfo, 'dataPath')
                return;
            end

            dataPath = string(app.dataNode.dataInfo.dataPath);
            if strlength(dataPath) == 0 || ~isfile(dataPath)
                return;
            end

            try
                vars = whos('-file', char(dataPath));
                names = {vars.name};
                if any(strcmp(names, 'Time'))
                    tmp = load(char(dataPath), 'Time', '-mat');
                    candidate = double(tmp.Time(:))';
                    if numel(candidate) == nTime
                        timeVector = candidate;
                        return;
                    end
                end
            catch ME
                warning('SEAL:Result:TimeVectorLoadFailed', ...
                    '读取数据时间向量失败: %s', ME.message);
            end
        end

        function saveOneResult(app, algName, result, algoInfo, savePath, trialIdx, trialCount, batchID)
            % 保存单次（单算法 × 单 trial）结果
            % trialIdx   : 当前 trial 序号（1-based）；为空表示二维数据，走原行为
            % trialCount : 该算法总 trial 数（用于 showeeg 识别一组）
            % batchID    : 本次 run 生成的唯一批次 ID（同算法同批次所有 trial 共享）

            if nargin < 6 || isempty(trialIdx)
                trialTag = '';
                nameTag  = algName;
                isTrialSplit = false;
            else
                trialTag = sprintf('_trial%03d', trialIdx);
                nameTag  = sprintf('%s_trial%03d', algName, trialIdx);
                isTrialSplit = true;
            end

            timestamp = datestr(now,'yyyymmdd_HHMMSS_FFF');  % 加毫秒避免同秒冲突
            fileName  = sprintf('%s_%s%s_%s.mat', app.dataNode.name, algName, trialTag, timestamp);
            fullPath  = fullfile(savePath, fileName);

            resultInfo = struct();
            resultInfo.Algorithm        = ['seal_' algName];
            resultInfo.NumOrientations  = app.NumOrientations;
            resultInfo.srate            = app.srate;
            resultInfo.ResultFolderPath = savePath;
            resultInfo.name             = sprintf('%s_%s', app.dataNode.name, nameTag);
            resultInfo.parentDataName   = app.dataNode.name;

            timeVector = app.getCurrentDataTimeVector(size(result, 2));
            if ~isempty(timeVector)
                resultInfo.TimeVector = timeVector;
                resultInfo.TimeStart = timeVector(1);
                resultInfo.TimeEnd = timeVector(end);
            end

            % —— 关键：trial 联动所需的元数据 ——
            if isTrialSplit
                nTperTrial = size(result, 2);  % 单 trial 的时间点数
                tStart = (trialIdx - 1) * nTperTrial + 1;  % 在拼接长时序里的起始列（1-based）
                tEnd   = trialIdx * nTperTrial;

                resultInfo.IsTrialSplit = true;
                resultInfo.TrialIndex   = trialIdx;
                resultInfo.TrialCount   = trialCount;
                resultInfo.TrialRange   = [tStart, tEnd];  % 全局时间索引范围，闭区间
                resultInfo.BatchID      = batchID;         % 同批次共享
                resultInfo.AlgorithmKey = algName;         % 方便按算法分组
            else
                resultInfo.IsTrialSplit = false;
            end

            % 合并算法参数
            if isstruct(algoInfo)
                f = fieldnames(algoInfo);
                for k = 1:numel(f)
                    resultInfo.(f{k}) = algoInfo.(f{k});
                end
            end

            saveData(fullPath, result);
            app.dataNode.setResultPath(fullPath);
            node = app.sessionNode.openResultFromData(fullPath, resultInfo);

            if app.dataNode.dataInfo.isSimu == 1
                node.resultInfo.SimuSourcePath = app.dataNode.dataInfo.dataPath;
                node.resultInfo.isSimu = 1;
            end
            app.dataNode.save();
            node.save();
        end
        function est = estimateOutputSize(app, algNames, scaling)
            % 估算每个选中算法的输出大小和峰值内存占用
            % 返回结构体数组，每个元素对应一个算法的估计值

            nSrc  = size(app.protocolNode.cortexNode.data.Vertices, 1);
            nT    = size(app.Data, 2);
            nCh   = size(app.Data, 1);
            nO    = app.NumOrientations;

            nAlgs = numel(algNames);
            est(nAlgs) = struct( ...
                'algName','', ...
                'outputBytes',0, ...
                'peakMemBytes',0, ...
                'diskBytes',0, ...
                'warnings',{{}});

            for i = 1:nAlgs
                algName = char(algNames(i));

                % 输出 S 的大小（后处理前是 3*nSrc x nT，后处理降维后是 nSrc x nT）
                outputElems_raw  = nO * nSrc * nT;      % 算法原始输出
                outputElems_post = nSrc * nT;           % 后处理降维后
                outputBytes = outputElems_post * 8;     % double 8 字节

                % 算法特定的峰值内存估计
                peakMem = estimatePeakMemory(algName, nCh, nSrc, nT, nO);

                % 磁盘：.mat 文件大小大约等于 outputBytes（未压缩），压缩后 0.3~0.7 倍
                diskBytes = outputBytes * 0.7;

                % 收集警告
                warns = {};
                if outputBytes > 2e9
                    warns{end+1} = sprintf('输出 %.2f GB，超过 2GB', outputBytes/1e9);
                end
                if peakMem > 8e9
                    warns{end+1} = sprintf('峰值内存约 %.2f GB，超过 8GB', peakMem/1e9);
                end

                est(i).algName      = algName;
                est(i).outputBytes  = outputBytes;
                est(i).peakMemBytes = peakMem;
                est(i).diskBytes    = diskBytes;
                est(i).warnings     = warns;
            end
        end


       function sys = queryAvailableResources(app, savePath)
    % 返回结构体:
    %   availMem      - 总可用内存 (bytes), 用于"所有算法叠加"判断
    %   largestBlock  - 单次最大可分配连续块 (bytes), 用于"单个大矩阵"判断
    %                   仅 Windows 能精确获取,Mac/Linux 用 availMem 代替
    %   totalMem      - 物理总内存 (bytes)
    %   availDisk     - 目标路径所在卷的剩余空间 (bytes)
    %   platform      - 'windows' | 'mac' | 'linux'
    %   memSource     - 内存数据来源 (用于诊断)
    %   diskSource    - 磁盘数据来源

    sys = struct( ...
        'availMem',     [], ...
        'largestBlock', [], ...
        'totalMem',     [], ...
        'availDisk',    [], ...
        'platform',     '', ...
        'memSource',    'fallback', ...
        'diskSource',   'fallback');

    % ============================================================
    %  内存查询
    % ============================================================
    if ispc
        sys.platform = 'windows';
        try
            m = memory;

            % 总可用 (这才是真正的"可用",不是"最大单块")
            sys.availMem = m.MemAvailableAllArrays;

            % 最大连续块 (受碎片影响,可能远小于 availMem)
            largest = m.MaxPossibleArrayBytes;

            % R2022a+ 在某些情况下返回接近 intmax('int64') (9.2e18)
            % 这种"理论上限"对实际分配毫无意义,要降级成 availMem
            if isempty(largest) || largest <= 0 || largest > 1e15
                largest = sys.availMem;
            end
            sys.largestBlock = largest;

            % 物理总内存:用 wmic 查
            try
                [s, out] = system('wmic ComputerSystem get TotalPhysicalMemory /value');
                if s == 0
                    tok = regexp(out, 'TotalPhysicalMemory=(\d+)', 'tokens', 'once');
                    if ~isempty(tok)
                        sys.totalMem = str2double(tok{1});
                    end
                end
            catch
            end
            if isempty(sys.totalMem) || isnan(sys.totalMem)
                sys.totalMem = sys.availMem;
            end

            sys.memSource = 'memory_function';
        catch ME
            warning('queryAvailableResources:WindowsMemoryFailed', ...
                'memory() 调用失败: %s', ME.message);
        end

    elseif ismac
        sys.platform = 'mac';
        try
            % 物理总内存
            [s, out] = system('sysctl -n hw.memsize');
            if s == 0
                sys.totalMem = str2double(strtrim(out));
            end

            % page size (Intel Mac = 4096, Apple Silicon = 16384)
            [s, out] = system('sysctl -n hw.pagesize');
            if s == 0
                pageSize = str2double(strtrim(out));
            else
                pageSize = 4096;   % 兜底
            end

            % vm_stat: 用 free + inactive + speculative 估算"实际可用"
            [s, vm] = system('vm_stat');
            if s == 0
                getPages = @(name) parsePages(vm, name);
                freeP     = getPages('Pages free');
                inactiveP = getPages('Pages inactive');
                specP     = getPages('Pages speculative');

                pages = 0;
                if ~isnan(freeP),     pages = pages + freeP;     end
                if ~isnan(inactiveP), pages = pages + inactiveP; end
                if ~isnan(specP),     pages = pages + specP;     end

                if pages > 0
                    sys.availMem  = pages * pageSize;
                    sys.memSource = 'vm_stat';
                end
            end

            % vm_stat 解析失败时,退而求其次:总内存的 60%
            if isempty(sys.availMem) && ~isempty(sys.totalMem)
                sys.availMem  = sys.totalMem * 0.6;
                sys.memSource = 'estimate_60pct';
            end

            sys.largestBlock = sys.availMem;   % macOS 无法精确知道连续块上限
        catch ME
            warning('queryAvailableResources:MacMemoryFailed', ...
                'macOS 内存查询失败: %s', ME.message);
        end

    else
        sys.platform = 'linux';
        try
            % 优先读 /proc/meminfo (最可靠)
            [s, out] = system('cat /proc/meminfo');
            if s == 0
                tokTotal = regexp(out, 'MemTotal:\s+(\d+)\s*kB',     'tokens', 'once');
                tokAvail = regexp(out, 'MemAvailable:\s+(\d+)\s*kB', 'tokens', 'once');
                if ~isempty(tokTotal)
                    sys.totalMem = str2double(tokTotal{1}) * 1024;
                end
                if ~isempty(tokAvail)
                    sys.availMem  = str2double(tokAvail{1}) * 1024;
                    sys.memSource = 'meminfo';
                end
            end

            % fallback: free -b + awk 直接取列 (避免数 token 越界)
            if isempty(sys.availMem)
                [s, out] = system('free -b | awk ''/^Mem:/ {print $7}''');
                if s == 0
                    val = str2double(strtrim(out));
                    if ~isnan(val) && val > 0
                        sys.availMem  = val;
                        sys.memSource = 'free';
                    end
                end
            end
            if isempty(sys.totalMem)
                [s, out] = system('free -b | awk ''/^Mem:/ {print $2}''');
                if s == 0
                    val = str2double(strtrim(out));
                    if ~isnan(val) && val > 0
                        sys.totalMem = val;
                    end
                end
            end

            sys.largestBlock = sys.availMem;
        catch ME
            warning('queryAvailableResources:LinuxMemoryFailed', ...
                'Linux 内存查询失败: %s', ME.message);
        end
    end

    % ============================================================
    %  内存兜底 (上面任何分支没拿到值都会落到这里)
    % ============================================================
    if isempty(sys.availMem) || isnan(sys.availMem) || sys.availMem <= 0
        sys.availMem  = 8e9;
        sys.memSource = 'fallback_8GB';
    end
    if isempty(sys.totalMem) || isnan(sys.totalMem) || sys.totalMem <= 0
        sys.totalMem = max(sys.availMem, 8e9);
    end
    if isempty(sys.largestBlock) || isnan(sys.largestBlock) || sys.largestBlock <= 0
        sys.largestBlock = sys.availMem;
    end

    % ============================================================
    %  磁盘查询
    % ============================================================
    queryPath = savePath;
    if ~isfolder(queryPath)
        % 路径还没创建时,查父目录
        queryPath = fileparts(queryPath);
        if isempty(queryPath) || ~isfolder(queryPath)
            queryPath = pwd;
        end
    end

    try
        if ispc
            drive = queryPath(1:2);   % "C:"

            % 优先用 wmic (输出稳定,不受系统语言影响)
            [s, out] = system(['wmic logicaldisk where "DeviceID=''' drive '''" get FreeSpace /value']);
            if s == 0
                tok = regexp(out, 'FreeSpace=(\d+)', 'tokens', 'once');
                if ~isempty(tok)
                    sys.availDisk  = str2double(tok{1});
                    sys.diskSource = 'wmic';
                end
            end

            % fallback: dir /-C
            if isempty(sys.availDisk)
                [s, out] = system(['dir "' queryPath '" /-C']);
                if s == 0
                    tok = regexp(out, '(\d+)\s+bytes free', 'tokens');
                    if ~isempty(tok)
                        sys.availDisk  = str2double(tok{end}{1});
                        sys.diskSource = 'dir';
                    end
                end
            end
        else
            % macOS + Linux 通用: df -k (POSIX 标准,返回 1KB 块数)
            [s, out] = system(['df -k "' queryPath '"']);
            if s == 0
                lines = strsplit(out, newline);
                lines = lines(~cellfun(@isempty, strtrim(lines)));
                if numel(lines) >= 2
                    cols = strsplit(strtrim(lines{end}));
                    if numel(cols) >= 4
                        val = str2double(cols{4});
                        if ~isnan(val) && val > 0
                            sys.availDisk  = val * 1024;
                            sys.diskSource = 'df';
                        end
                    end
                end
            end
        end
    catch ME
        warning('queryAvailableResources:DiskFailed', ...
            '磁盘查询失败: %s', ME.message);
    end

    if isempty(sys.availDisk) || isnan(sys.availDisk) || sys.availDisk <= 0
        sys.availDisk  = 50e9;
        sys.diskSource = 'fallback_50GB';
    end

    % ============================================================
    %  内嵌辅助函数 (R2016b+ 支持函数内的局部函数)
    % ============================================================
    function val = parsePages(text, fieldName)
        pat = [fieldName, ':\s+(\d+)'];
        tok = regexp(text, pat, 'tokens', 'once');
        if isempty(tok)
            val = NaN;
        else
            val = str2double(tok{1});
        end
    end
       end

        function decision = preCheckAndConfirm(app, algNames, scaling, savePath)
            % 返回决策：'proceed' | 'abort' | 'proceed_with_subset'
            % 如果用户选择部分跑，modified algNames 通过 app.pendingAlgNames 传出

            decision = 'proceed';

            est = app.estimateOutputSize(algNames, scaling);
            sys = app.queryAvailableResources(savePath);

            totalOutput = sum([est.outputBytes]);
            totalDisk   = sum([est.diskBytes]);
            maxPeakMem  = max([est.peakMemBytes]);

            % --- 硬限制（超过就强制拦截）---
            HARD_TOTAL_LIMIT   = sys.availMem * 0.85;     % 总占用不超过可用 85%
            HARD_SINGLE_LIMIT  = sys.largestBlock * 0.95; % 单矩阵不超过最大块 95%
            HARD_DISK_LIMIT    = sys.availDisk * 0.95;

            hardBlockReasons = {};

            if maxPeakMem > HARD_SINGLE_LIMIT
                hardBlockReasons{end+1} = sprintf( ...
                    '单算法峰值 %.2f GB 超过最大可分配连续块 %.2f GB(可能因内存碎片)', ...
                    maxPeakMem/1e9, sys.largestBlock/1e9);
            end

            if maxPeakMem > HARD_TOTAL_LIMIT
                hardBlockReasons{end+1} = sprintf( ...
                    '单算法峰值 %.2f GB 超过总可用内存 %.2f GB', ...
                    maxPeakMem/1e9, sys.availMem/1e9);
            end

            if totalDisk > HARD_DISK_LIMIT
                hardBlockReasons{end+1} = sprintf( ...
                    '总输出 %.2f GB 超过可用磁盘 %.2f GB', ...
                    totalDisk/1e9, sys.availDisk/1e9);
            end

            if ~isempty(hardBlockReasons)
                msg = sprintf(['无法运行：\n\n%s\n\n' ...
                    '建议：\n' ...
                    '  • 减少勾选的算法数量\n' ...
                    '  • 减少数据时间长度（裁剪 epoch）\n' ...
                    '  • 使用较稀疏的 cortex\n' ...
                    '  • 释放磁盘空间或选择其他保存位置'], ...
                    strjoin(hardBlockReasons, sprintf('\n')));
                uialert(app.UIFigure, msg, '资源不足', 'Icon','error');
                decision = 'abort';
                return;
            end

            % --- 软警告（提示但允许继续）---
            softWarnings = {};
            if totalOutput > 5e9
                softWarnings{end+1} = sprintf('预计总输出: %.2f GB', totalOutput/1e9);
            end
            if maxPeakMem > sys.availMem * 0.7
                softWarnings{end+1} = sprintf( ...
                    '峰值内存接近上限: %.2f GB / %.2f GB 可用', ...
                    maxPeakMem/1e9, sys.availMem/1e9);
            end

            perAlgWarnings = {};
            for i = 1:numel(est)
                if ~isempty(est(i).warnings)
                    perAlgWarnings{end+1} = sprintf('  • %s: %s', ...
                        est(i).algName, strjoin(est(i).warnings, '; '));
                end
            end

            % --- 构建详情表 ---
            detailLines = {'算法预估详情:'};
            for i = 1:numel(est)
                detailLines{end+1} = sprintf( ...
                    '  %-20s  输出:%7.1f MB  峰值内存:%7.1f MB', ...
                    est(i).algName, ...
                    est(i).outputBytes/1e6, ...
                    est(i).peakMemBytes/1e6);
            end

            % 如果有警告，弹确认框
            if ~isempty(softWarnings) || ~isempty(perAlgWarnings)
                msg = sprintf('资源评估:\n\n%s\n\n%s\n\n%s\n\n是否继续？', ...
                    strjoin(softWarnings, sprintf('\n')), ...
                    strjoin(perAlgWarnings, sprintf('\n')), ...
                    strjoin(detailLines, sprintf('\n')));

                sel = uiconfirm(app.UIFigure, msg, '预检查结果', ...
                    'Options',{'继续运行','只运行安全算法','取消'}, ...
                    'DefaultOption',1, ...
                    'CancelOption',3, ...
                    'Icon','warning');

                switch sel
                    case '继续运行'
                        decision = 'proceed';
                    case '只运行安全算法'
                        % 过滤掉超标的算法
                        safeIdx = arrayfun(@(e) isempty(e.warnings), est);
                        app.pendingAlgNames = algNames(safeIdx);
                        if isempty(app.pendingAlgNames)
                            uialert(app.UIFigure,'没有安全算法可运行。','提示');
                            decision = 'abort';
                        else
                            decision = 'proceed_with_subset';
                        end
                    case '取消'
                        decision = 'abort';
                end
            else
                % 全部算法都在安全范围内，可直接继续（也可选择显示简要确认）
                totalMB = totalOutput/1e6;
                if totalMB > 200   % 超过 200 MB 显示确认
                    msg = sprintf(['资源评估通过:\n\n' ...
                        '  算法数量: %d\n' ...
                        '  预计总输出: %.1f MB\n' ...
                        '  最大峰值内存: %.1f MB\n\n' ...
                        '%s\n\n是否继续？'], ...
                        numel(algNames), totalMB, maxPeakMem/1e6, ...
                        strjoin(detailLines, sprintf('\n')));
                    sel = uiconfirm(app.UIFigure, msg, '预检查', ...
                        'Options',{'继续','取消'}, 'DefaultOption',1);
                    if strcmp(sel,'取消')
                        decision = 'abort';
                    end
                end
            end
        end

        function Cov_n = resolveNoiseCovariance(app, nCh)
            Cov_n = [];

            if ~isempty(app.dataNode) && ~isempty(app.dataNode.dataInfo) && ...
                    isprop(app.dataNode.dataInfo, 'NoiseCovMat') && ~isempty(app.dataNode.dataInfo.NoiseCovMat)
                Cov_n = app.selectNoiseCovChannels(app.dataNode.dataInfo.NoiseCovMat, nCh);
                if ~isempty(Cov_n)
                    return;
                end
            end

            dataPath = "";
            if ~isempty(app.dataNode) && ~isempty(app.dataNode.dataInfo) && isprop(app.dataNode.dataInfo, 'dataPath')
                dataPath = string(app.dataNode.dataInfo.dataPath);
            end
            if strlength(dataPath) == 0 || ~isfile(dataPath)
                return;
            end

            dataDir = fileparts(char(dataPath));
            candidates = {
                fullfile(dataDir, 'NoiseCovariance.mat')
                fullfile(dataDir, 'noiseCovariance.mat')
                fullfile(dataDir, 'noisecov_full.mat')
                };

            for ii = 1:numel(candidates)
                if ~isfile(candidates{ii})
                    continue;
                end
                try
                    s = load(candidates{ii}, '-mat');
                    Cov_full = app.extractNoiseCovariance(s);
                    Cov_n = app.selectNoiseCovChannels(Cov_full, nCh);
                    if ~isempty(Cov_n)
                        fprintf('[Whitening] 使用噪声协方差文件: %s\n', candidates{ii});
                        return;
                    end
                catch ME
                    warning('SEAL:Whitening:NoiseCovLoadFailed', ...
                        '读取噪声协方差文件失败 %s: %s', candidates{ii}, ME.message);
                end
            end
        end

        function Cov = extractNoiseCovariance(app, s)
            Cov = [];
            if isfield(s, 'NoiseCov')
                Cov = s.NoiseCov;
                return;
            end
            if isfield(s, 'noiseCovariance')
                obj = s.noiseCovariance;
                if isstruct(obj) && isfield(obj, 'NoiseCov')
                    Cov = obj.NoiseCov;
                    return;
                elseif isobject(obj) && isprop(obj, 'NoiseCov')
                    Cov = obj.NoiseCov;
                    return;
                end
            end

            fields = fieldnames(s);
            for ii = 1:numel(fields)
                obj = s.(fields{ii});
                if isnumeric(obj) && ismatrix(obj) && size(obj,1) == size(obj,2)
                    Cov = obj;
                    return;
                elseif isstruct(obj) && isfield(obj, 'NoiseCov')
                    Cov = obj.NoiseCov;
                    return;
                elseif isobject(obj) && isprop(obj, 'NoiseCov')
                    Cov = obj.NoiseCov;
                    return;
                end
            end
        end

        function Cov = selectNoiseCovChannels(app, Cov_full, nCh)
            Cov = [];
            if isempty(Cov_full) || ~isnumeric(Cov_full) || ~ismatrix(Cov_full)
                return;
            end
            Cov_full = double(Cov_full);
            if isequal(size(Cov_full), [nCh, nCh])
                Cov = Cov_full;
                return;
            end
            if size(Cov_full,1) ~= size(Cov_full,2) || size(Cov_full,1) < nCh
                return;
            end

            idx = [];
            dataPath = "";
            if ~isempty(app.dataNode) && ~isempty(app.dataNode.dataInfo) && isprop(app.dataNode.dataInfo, 'dataPath')
                dataPath = string(app.dataNode.dataInfo.dataPath);
            end
            if strlength(dataPath) > 0 && isfile(dataPath)
                try
                    vars = whos('-file', char(dataPath));
                    names = {vars.name};
                    if any(strcmp(names, 'iEEG'))
                        tmp = load(char(dataPath), 'iEEG', '-mat');
                        idx = double(tmp.iEEG(:));
                    elseif any(strcmp(names, 'GoodChannel'))
                        tmp = load(char(dataPath), 'GoodChannel', '-mat');
                        idx = double(tmp.GoodChannel(:));
                    end
                catch
                    idx = [];
                end
            end

            if numel(idx) == nCh
                if min(idx) == 0
                    idx = idx + 1;
                end
                if all(idx >= 1) && all(idx <= size(Cov_full,1))
                    Cov = Cov_full(idx, idx);
                end
            end
        end

        function [W, Cov_reg, rankNoise, condNumber] = brainstormRegWhitener(app, Cov_n, noiseReg)
            Cov_n = (double(Cov_n) + double(Cov_n)') / 2;
            [U_C, S_C] = svd(Cov_n, 'econ');
            eigVals = diag(S_C);
            if isempty(eigVals) || eigVals(1) <= 0
                error('SEAL:Whitening:InvalidNoiseCov', '噪声协方差矩阵没有正特征值。');
            end

            noiseStd = sqrt(eigVals);
            tol = length(noiseStd) * eps(single(noiseStd(1)));
            rankNoise = sum(noiseStd > tol);
            U_keep = U_C(:, 1:rankNoise);
            noiseStd = noiseStd(1:rankNoise);

            ridgeFactor = mean(eigVals) * noiseReg;
            Cov_reg = U_keep * diag(noiseStd.^2) * U_keep' + ridgeFactor * eye(size(Cov_n,1));
            W = U_keep * diag(1 ./ sqrt(noiseStd.^2 + ridgeFactor)) * U_keep';

            regEig = eig((Cov_reg + Cov_reg') / 2);
            posEig = regEig(regEig > 0);
            condNumber = max(posEig) / min(posEig);
        end

        function depthLeadfield = getMneDepthLeadfield(app, nCh, nSources)
            depthLeadfield = [];
            if isempty(app.protocolNode) || isempty(app.protocolNode.leadfieldNode) || ...
                    isempty(app.protocolNode.leadfieldNode.leadfieldInfo)
                return;
            end

            lfPath = string(app.protocolNode.leadfieldNode.leadfieldInfo.dataPath);
            if strlength(lfPath) == 0 || ~isfile(lfPath)
                return;
            end

            try
                vars = whos('-file', char(lfPath));
                names = {vars.name};
                candidates = {'L_free_avgref', 'Gain_free_avgref', 'Leadfield_free_avgref', ...
                    'L_avgref', 'Gain_avgref', 'Leadfield_avgref', 'L_free', 'Gain_raw'};
                for ii = 1:numel(candidates)
                    if ~any(strcmp(names, candidates{ii}))
                        continue;
                    end
                    tmp = load(char(lfPath), candidates{ii}, '-mat');
                    value = double(tmp.(candidates{ii}));
                    if isequal(size(value), [nCh, 3 * nSources])
                        depthLeadfield = value;
                        fprintf('[MNE] Brainstorm depth leadfield: %s\n', candidates{ii});
                        return;
                    end
                end
            catch ME
                warning('SEAL:MNE:DepthLeadfieldLoadFailed', ...
                    '读取 Brainstorm 深度权重导联场失败: %s', ME.message);
            end
        end

        function [W, R_white, info] = computeWhitener(app, X, strategy)
            % COMPUTEWHITENER  根据 NoiseStrategy 计算白化矩阵 W
            %   X        : 传感器数据 nCh x nT
            %   strategy : 'brainstorm' | 'identity' | 'baseline' | 'data_cov' | 'pass_through'
            %
            % 返回:
            %   W        : nCh x nCh 白化矩阵。满足 W * Cov_n * W' ≈ I
            %   R_white  : 白化后算法应传入的 NoiseCovariance (通常是 I)
            %   info     : 结构体,记录实际使用的策略、baseline 长度、条件数等

            nCh = size(X, 1);
            info = struct('strategy', strategy, 'baselineSamples', 0, ...
                'condNumber', NaN, 'rankDeficient', false);

            switch lower(strategy)

                case {'brainstorm', 'brainstorm_reg'}
                    Cov_n = app.resolveNoiseCovariance(nCh);
                    if isempty(Cov_n)
                        baselineLen = min(90, size(X, 2));
                        warning('SEAL:Whitening:NoBrainstormCov', ...
                            '未找到 Brainstorm 噪声协方差,已回退到 baseline 前 %d 个采样点', baselineLen);
                        if baselineLen < 2
                            warning('SEAL:Whitening:BaselineTooShort', ...
                                'Baseline length (%d) too short; falling back to identity.', baselineLen);
                            W       = eye(nCh);
                            R_white = eye(nCh);
                            info.strategy = 'brainstorm_fallback_identity';
                            info.source = 'fallback_identity';
                            return;
                        end
                        if baselineLen < nCh
                            warning('SEAL:Whitening:RankDeficient', ...
                                'Baseline length (%d) < channels (%d); covariance will be regularized.', ...
                                baselineLen, nCh);
                            info.rankDeficient = true;
                        end
                        Bpre  = X(:, 1:baselineLen);
                        Bpre  = Bpre - mean(Bpre, 2);
                        Cov_n = (Bpre * Bpre') / max(baselineLen - 1, 1);
                        info.strategy = 'brainstorm_fallback_baseline';
                        info.source = 'fallback_baseline_first90';
                        info.baselineSamples = baselineLen;
                    else
                        [W, ~, rankNoise, condNumber] = app.brainstormRegWhitener(Cov_n, 0.1);
                        R_white = eye(nCh);
                        info.source = 'brainstorm_reg_noise_cov';
                        info.noiseReg = 0.1;
                        info.rankDeficient = rankNoise < nCh;
                        info.rankNoise = rankNoise;
                        info.condNumber = condNumber;
                        return;
                    end

                case 'identity'
                    % 不做白化,等价于假设噪声已经是单位阵
                    W       = eye(nCh);
                    R_white = eye(nCh);
                    return;

                case 'pass_through'
                    % 不动数据,让算法内部自己处理噪声(传空矩阵)
                    W       = eye(nCh);
                    R_white = [];
                    return;

                case 'precomputed'
                    
                        Cov_n = app.resolveNoiseCovariance(nCh);
                    
                    if isempty(Cov_n)
                        warning('SEAL:Whitening:NoPrecomputedCov', ...
                            '未找到预计算的 NoiseCovMat/NoiseCovariance.mat,已回退到 identity');
                        W       = eye(nCh);
                        R_white = eye(nCh);
                        info.source = 'fallback_identity';
                        return;
                    end

                    if ~isequal(size(Cov_n), [nCh, nCh])
                        warning('SEAL:Whitening:SizeMismatch', ...
                            '预计算 NoiseCovMat 尺寸 (%dx%d) 与当前数据通道数 (%d) 不匹配,已回退到 identity', ...
                            size(Cov_n,1), size(Cov_n,2), nCh);
                        W       = eye(nCh);
                        R_white = eye(nCh);
                        info.source = 'fallback_identity';
                        return;
                    end

                    fprintf('[Whitening] 使用预计算的 NoiseCovMat\n');
                    info.source = 'precomputed';

                case 'baseline'
                    % 用刺激前 baseline 估计噪声协方差
                    % 默认用前 90 个采样点;若不足则用全部
                    baselineLen = min(90, size(X, 2));
                    if baselineLen < nCh
                        warning('SEAL:Whitening:RankDeficient', ...
                            'Baseline length (%d) < channels (%d); covariance will be regularized.', ...
                            baselineLen, nCh);
                        info.rankDeficient = true;
                    end
                    Bpre  = X(:, 1:baselineLen);
                    Bpre  = Bpre - mean(Bpre, 2);
                    Cov_n = (Bpre * Bpre') / max(baselineLen - 1, 1);
                    info.baselineSamples = baselineLen;

                case 'data_cov'
                    % 用整段数据估协方差(只在无 baseline 时用,会把信号当噪声)
                    Xc    = X - mean(X, 2);
                    Cov_n = (Xc * Xc') / max(size(X,2) - 1, 1);

                otherwise
                    warning('SEAL:Whitening:UnknownStrategy', ...
                        'Unknown NoiseStrategy "%s", falling back to identity.', strategy);
                    W       = eye(nCh);
                    R_white = eye(nCh);
                    return;
            end

            % --- 对称化 + SVD + 特征值下限截断 ---
            Cov_n = (Cov_n + Cov_n') / 2;
            [U_C, S_C] = svd(Cov_n);
            s = diag(S_C);

            if isempty(s) || ~isfinite(s(1)) || s(1) <= 0
                warning('SEAL:Whitening:InvalidCovariance', ...
                    'Noise covariance has no positive finite eigenvalue; falling back to identity.');
                W       = eye(nCh);
                R_white = eye(nCh);
                info.source = 'fallback_identity_invalid_cov';
                info.condNumber = NaN;
                return;
            end

            sMax       = s(1);
            rankTol    = max(size(Cov_n)) * eps(sMax);
            rankNoise  = sum(s > rankTol);
            lowerBound = sMax * 1e-6;
            s(s < lowerBound) = lowerBound;

            W        = diag(1 ./ sqrt(s)) * U_C';
            R_white  = eye(nCh);
            info.rankNoise = rankNoise;
            info.rankDeficient = info.rankDeficient || rankNoise < nCh;
            info.condNumber = sMax / s(end);
        end

function S = runOneAlgorithm(app, algName, X, L, C, R, Eye_N, nO, prm, toLog, cancelCheck)
    % 单次调用某个算法，返回原始 S（3*nSrc x nT 或 nSrc x nT）
    % 后处理（降维）和保存由调用方统一处理。
    if nargin < 11 || isempty(cancelCheck)
        cancelCheck = @() false;
    end
    S = [];
    switch algName

        % ---------- Minimum L2 Norm ----------
        case 'MNE'
            if ~isfield(prm, 'RegularizationParameter') || isempty(prm.RegularizationParameter)
                prm.RegularizationParameter = 1/3;
            end
            if ~isfield(prm, 'sourceCovariance') || isempty(prm.sourceCovariance)
                prm.sourceCovariance = 'brainstorm_depth';
            end
            depthLeadfield = [];
            useBrainstormDepth = (ischar(prm.sourceCovariance) || isstring(prm.sourceCovariance)) && ...
                strcmpi(char(string(prm.sourceCovariance)), 'brainstorm_depth');
            if useBrainstormDepth
                nSources = size(L, 2) / nO;
                depthLeadfield = app.getMneDepthLeadfield(size(L, 1), nSources);
            end
            [S,~] = seal_MNE(X, L, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'SourceCovariance',        prm.sourceCovariance, ...
                'DepthWeightExp',          app.DepthWeightExp, ...
                'DepthWeightLimit',        10, ...
                'DepthLeadfield',          depthLeadfield, ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'dSPM'
            [S,~] = seal_dSPM(X, L, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'sLORETA'
            [S,~] = seal_sLORETA(X, L, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'eLORETA'
            [S,~] = seal_eLORETA(X, L, ...
                'SourceCovariance', prm.sourceCovariance, ...
                'MaxIterations',    prm.MaxIterations, ...
                'Tolerance',        prm.Tolerance, ...
                'NoiseCovariance',  R, ...
                'NumOrientations',  nO);

        case 'LORETA'
            [S,~] = seal_LORETA(X, L, C, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'DepthWeightingExponent',  prm.DepthWeightingExponent, ...
                'LaplacianRegFactor',      prm.LaplacianRegFactor, ...
                'SourceCovariance',        prm.sourceCovariance, ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'LAURA'
            [S,~] = seal_LAURA(X, L, C, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'DepthWeightingExponent',  prm.DepthWeightingExponent, ...
                'Exponent_e_LAURA',        prm.Exponent_e_LAURA, ...
                'MaxNeighborsN_LAURA',     prm.MaxNeighborsN_LAURA, ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'LOGURA'
            [S,~] = seal_LOGURA(X, L, C, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'DepthWeightingExponent',  prm.DepthWeightingExponent, ...
                'Exponent_e_LOGURA',       prm.Exponent_e_LOGURA, ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        % ---------- Sparse Model ----------
        case 'Champagne'
            [S,~] = seal_Champagne(X, L, ...
                'PruneThreshold', prm.PruneThreshold, ...
                'UpdateNoise',    toLog(prm.UpdateNoise), ...
                'MaxIterations',  prm.MaxIterations, ...
                'Tolerance',      prm.Tolerance, ...
                'NoiseCovariance',R, ...
                'NumOrientations',nO);

        case 'BlockChampagne'
            [S,~] = seal_BlockChampagne(X, L, C, ...
                'NeighborOrder',          prm.NeighborOrder, ...
                'LearnNoise',             toLog(prm.LearnNoise), ...
                'NoiseObservationMatrix', Eye_N, ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case 'TS_Champagne'
            [S,~] = seal_TS_Champagne(X, L, C.VertConn,...
                'NumOrientations', nO);

        case {'SmoothChampagne', 'Smooth-Cham'}
            [S,~] = seal_SmoothChampagne(X, L, C, ...
                'KernelWidth',            prm.KernelWidth, ...
                'KernelPower',            prm.KernelPower, ...
                'KernelThreshold',        prm.KernelThreshold, ...
                'TileSize',               prm.TileSize, ...
                'MaxIterations',          prm.MaxIterations, ...
                'Tolerance',              prm.Tolerance, ...
                'UpdateNoiseCovariance',  toLog(prm.UpdateNoiseCovariance), ...
                'Verbose',                toLog(prm.Verbose), ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case 'IRES'
            [S,~] = seal_IRES(X, L, C, ...
                'Alpha',                  prm.Alpha, ...
                'MaxIterations',          prm.MaxIterations, ...
                'MaxADMMIterations',      prm.MaxADMMIterations, ...
                'Tolerance',              prm.Tolerance, ...
                'ADMMTolerance',          prm.ADMMTolerance, ...
                'ProgressInterval',       prm.ProgressInterval, ...
                'Rho',                    prm.Rho, ...
                'WeightEpsilon',          prm.WeightEpsilon, ...
                'Verbose',                toLog(prm.Verbose), ...
                'CancelCheck',            cancelCheck, ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case 'FOCUSS'
            [S,~] = seal_FOCUSS(X, L, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'RegularizationMode',      prm.RegularizationMode, ...
                'MaxIterations',           prm.MaxIterations, ...
                'Tolerance',               prm.Tolerance, ...
                'PruneThreshold',          prm.PruneThreshold, ...
                'DepthWeighting',          prm.DepthWeighting, ...
                'WeightMode',              prm.WeightMode, ...
                'Verbose',                 toLog(prm.Verbose), ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'Lq_ADMM'
            [S,~] = seal_Lq_ADMM(X, L, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'DataFidelity',            prm.DataFidelity, ...
                'ADMM_Rho',                prm.ADMM_Rho, ...
                'q_norm',                  prm.q_norm, ...
                'MaxIterations',           prm.MaxIterations, ...
                'Tolerance',               prm.Tolerance, ...
                'NumOrientations',         nO);

            % ---------- Spatio-Temporal ----------
        case 'STARTS'
            [S,~] = seal_STARTS(X, L, C, ...
                'NumTBFs',                prm.NumTBFs, ...
                'MaxIterations',          prm.MaxIterations, ...
                'Tolerance',              prm.Tolerance, ...
                'UpdateNoiseCovariance',  prm.UpdateNoiseCovariance, ...
                'NeighborOrder',          prm.NeighborOrder, ...
                'BlockCorrelation',       prm.BlockCorrelation, ...
                'PruneGammas',            prm.PruneGammas, ...
                'PruneAlphas',            prm.PruneAlphas, ...
                'CwkRegularization',      prm.CwkRegularization, ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case 'sSTARTS'
            [S,~] = seal_sSTARTS(X, L, C, ...
                'NumTBFs',                prm.NumTBFs, ...
                'LaplacianTau',           prm.LaplacianTau, ...
                'MaxIterations',          prm.MaxIterations, ...
                'Tolerance',              prm.Tolerance, ...
                'UpdateNoiseCovariance',  prm.UpdateNoiseCovariance, ...
                'PruneGammas',            prm.PruneGammas, ...
                'PruneAlphas',            prm.PruneAlphas, ...
                'CwkRegularization',      prm.CwkRegularization, ...
                'Verbose',                toLog(prm.Verbose), ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case 'uSTAR'
            [S,~] = seal_uSTAR(X, L, C, ...
                'NumTBFs',        prm.NumTBFs, ...
                'MaxIterations',  prm.MaxIterations, ...
                'Tolerance',      prm.Tolerance, ...
                'RunAlgorithm',   prm.uSTARRunAlgorithm, ...
                'PriorType',      prm.PriorType, ...
                'CwkRegularization', prm.CwkRegularization, ...
                'MicrostateMaps', app.MicrostateMaps, ...
                'NumOrientations',nO);

        case 'BESTIES'
            [S,~] = seal_BESTIES(X, L, C, ...
                'NumTBFs',                prm.NumTBFs, ...
                'TemporalBasisMode',      prm.TemporalBasisMode, ...
                'InitialTau',             prm.InitialTau, ...
                'LearnTau',               toLog(prm.LearnTau), ...
                'ComputeFreeEnergy',      toLog(prm.ComputeFreeEnergy), ...
                'NoiseObservationMatrix', Eye_N, ...
                'MaxIterations',          prm.MaxIterations, ...
                'Tolerance',              prm.Tolerance, ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case {'STBFSI', 'SI-STBF'}
            [S,~] = seal_STBFSI(X, L, C, ...
                'NumTBFs',                prm.NumTBFs, ...
                'MaxSpatialBases',        prm.MaxSpatialBases, ...
                'NeighborOrder',          prm.NeighborOrder, ...
                'DiffusionSigma',         prm.DiffusionSigma, ...
                'DiffusionOrder',         prm.DiffusionOrder, ...
                'SeedSelection',          prm.SeedSelection, ...
                'MaxIterations',          prm.MaxIterations, ...
                'Tolerance',              prm.Tolerance, ...
                'PruneAlphas',            prm.PruneAlphas, ...
                'PruneGammas',            prm.PruneGammas, ...
                'PruneSpatialBases',      toLog(prm.PruneSpatialBases), ...
                'CovarianceMode',         prm.CovarianceMode, ...
                'Verbose',                toLog(prm.Verbose), ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case 'MxNE'
            [S,~] = seal_MxNE(X, L, ...
                'Alpha',           prm.Alpha, ...
                'WeightsMin',      prm.WeightsMin, ...
                'SolverMaxIter',   prm.SolverMaxIter, ...
                'SolverTolerance', prm.SolverTolerance, ...
                'DebiasMaxIter',   prm.DebiasMaxIter, ...
                'DebiasTolerance', prm.DebiasTolerance, ...
                'Debias',          toLog(prm.Debias), ...
                'PCAWhitening',    toLog(prm.PCAWhitening), ...
                'NoiseCovariance', R, ...
                'NumOrientations', nO);

        case 'irMxNE'
            [S,~] = seal_irMxNE(X, L, ...
                'Alpha',           prm.Alpha, ...
                'NumIRMXNEIter',   prm.NumIRMXNEIter, ...
                'WeightsMin',      prm.WeightsMin, ...
                'SolverMaxIter',   prm.SolverMaxIter, ...
                'SolverTolerance', prm.SolverTolerance, ...
                'IRMXNETolerance', prm.IRMXNETolerance, ...
                'DebiasMaxIter',   prm.DebiasMaxIter, ...
                'DebiasTolerance', prm.DebiasTolerance, ...
                'Debias',          toLog(prm.Debias), ...
                'PCAWhitening',    toLog(prm.PCAWhitening), ...
                'NoiseCovariance', R, ...
                'NumOrientations', nO);

        case 'FASTIRES'
            [S,~] = seal_FASTIRES(X, L, C, ...
                'Alpha',                  prm.Alpha, ...
                'NumTBFs',                prm.NumTBFs, ...
                'Epsilon',                prm.Epsilon, ...
                'MaxIterations',          prm.MaxIterations, ...
                'MaxInnerIterationsX',    prm.MaxInnerIterationsX, ...
                'MaxInnerIterationsY',    prm.MaxInnerIterationsY, ...
                'MaxInnerIterationsProj', prm.MaxInnerIterationsProj, ...
                'Tolerance',              prm.Tolerance, ...
                'NoiseCovariance',        R, ...
                'NumOrientations',        nO);

        case 'VSSI_Lp'
            [S,~] = seal_VSSI_Lp(X, L, C, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'P',                       prm.P, ...
                'Rho',                     prm.Rho, ...
                'MaxIterations',           prm.MaxIterations, ...
                'Tolerance',               prm.Tolerance, ...
                'ProgressInterval',        prm.ProgressInterval, ...
                'AdaptiveRho',             toLog(prm.AdaptiveRho), ...
                'NormalizeProblem',        toLog(prm.NormalizeProblem), ...
                'Verbose',                 toLog(prm.Verbose), ...
                'CancelCheck',             cancelCheck, ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'VSSI_GGD'
            [S,~] = seal_VSSI_GGD(X, L, C, ...
                'SparseWeight',       prm.SparseWeight, ...
                'P',                  prm.P, ...
                'MaxOuterIterations', prm.MaxOuterIterations, ...
                'MaxADMMIterations',  prm.MaxADMMIterations, ...
                'MaxIRLSIterations',  prm.MaxIRLSIterations, ...
                'Tolerance',          prm.Tolerance, ...
                'ADMMTolerance',      prm.ADMMTolerance, ...
                'ProgressInterval',   prm.ProgressInterval, ...
                'Rho',                prm.Rho, ...
                'HutchinsonProbes',   prm.HutchinsonProbes, ...
                'AdaptiveRho',        toLog(prm.AdaptiveRho), ...
                'NormalizeProblem',   toLog(prm.NormalizeProblem), ...
                'Verbose',            toLog(prm.Verbose), ...
                'CancelCheck',        cancelCheck, ...
                'NoiseCovariance',    R, ...
                'NumOrientations',    nO);

        case 'STOUT'
            [S,~] = seal_STOUT(X, L, ...
                'ReferenceSpectrum',         prm.ReferenceSpectrum, ...
                'SamplingRate',              app.srate, ...
                'ConstraintMode',            prm.ConstraintMode, ...
                'HybridGamma',               prm.HybridGamma, ...
                'RegularizationParameter',   prm.RegularizationParameter, ...
                'NoiseCovariance',           R, ...
                'NumOrientations',           nO, ...
                'Verbose',                   toLog(prm.Verbose));

        case 'dMAP'
            G_adj = C.VertConn;
            [S,~] = seal_dMAP(X, L, G_adj, ...
                'lambda_F',        prm.lambda_F, ...
                'alpha_F',         prm.alpha_F, ...
                'max_iter',        prm.max_iter, ...
                'tol',             prm.tol, ...
                'NoiseCovariance', R, ...
                'NumOrientations', nO);

        case 'dSPN'
            initialSourceCov = prm.InitialSourceCovScale * ones(size(L, 2), 1);
            [S,~] = seal_dSPN(X, L, C.VertConn, ...
                'NoiseCovariance',  R, ...
                'InitialSourceCov', initialSourceCov, ...
                'NumOrientations',  nO, ...
                'max_iter',         prm.max_iter, ...
                'tol',              prm.tol, ...
                'rho',              prm.rho, ...
                'q_min',            prm.q_min, ...
                'sigma_min',        prm.sigma_min, ...
                'riccati_iter',     prm.riccati_iter, ...
                'riccati_tol',      prm.riccati_tol, ...
                'InnovationJitter', prm.InnovationJitter, ...
                'smooth_pass',      toLog(prm.smooth_pass));

        case 'debiased_GroupLasso'
            [S,~] = seal_debiased_GroupLasso(X, L, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'DynamicLambda',           toLog(prm.DynamicLambda), ...
                'JointEstNum',             prm.JointEstNum, ...
                'Debias',                  toLog(prm.Debias), ...
                'MaxIterations',           prm.MaxIterations, ...
                'Tolerance',               prm.Tolerance, ...
                'ToleranceNorm',           prm.ToleranceNorm, ...
                'JointTolerance',          prm.JointTolerance, ...
                'VaryingRho',              toLog(prm.VaryingRho), ...
                'ClearNotSelect',          toLog(prm.ClearNotSelect), ...
                'PrewhitenNoise',          toLog(prm.PrewhitenNoise), ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'SISSY_L21'
            [S,~] = seal_SISSY_L21(X, L, C, ...
                'lambda',          prm.lambda, ...
                'DynamicLambda',   toLog(prm.DynamicLambda), ...
                'JointEstNum',     prm.JointEstNum, ...
                'alpha',           prm.alpha, ...
                'rho',             prm.rho, ...
                'max_iter',        prm.max_iter, ...
                'tol',             prm.tol, ...
                'NoiseCovariance', R, ...
                'NumOrientations', nO);

        case 'TFMxNE'
            [S,~] = seal_TFMxNE(X, L, ...
                'SpatialReg',      prm.SpatialReg, ...
                'TemporalReg',     prm.TemporalReg, ...
                'TimeStep',        prm.TimeStep, ...
                'WindowSize',      prm.WindowSize, ...
                'MaxIterations',   prm.MaxIterations, ...
                'Tolerance',       prm.Tolerance, ...
                'DepthBiasExponent', prm.DepthBiasExponent, ...
                'DepthBiasLimit',    prm.DepthBiasLimit, ...
                'LooseOrientationWeight', prm.LooseOrientationWeight, ...
                'Debias',          toLog(prm.Debias), ...
                'NoiseCovariance', R, ...
                'NumOrientations', nO);

        % ---------- Beamformer ----------
        case 'LCMV'
            [S,~] = seal_LCMV(X, L, ...
                'Type',                    prm.Type, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'ScalarBeamformer',        toLog(prm.ScalarBeamformer), ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        case 'MCMV'
            [S,~] = seal_MCMV(X, L, ...
                'RegularizationParameter', prm.RegularizationParameter, ...
                'MaxIterations',           prm.MaxIterations, ...
                'AutoScaleNoise',          toLog(prm.AutoScaleNoise), ...
                'Verbose',                 double(toLog(prm.Verbose)), ...
                'NoiseCovariance',         R, ...
                'NumOrientations',         nO);

        otherwise
            warning('算法 "%s" 尚未实现，已跳过。', algName);
    end
end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp, selectednode)
          app.UIFigure.Name='Select Algorithms';
            app.Callingapp=mainapp;
%%调整app窗口位置
            if ~isempty(app.Callingapp) && isvalid(app.Callingapp)
                % 1. 获取主窗口的位置信息 [Left, Bottom, Width, Height]
                mainPos = app.Callingapp.UIFigure.Position;
                
                % 2. 获取当前子窗口的自身大小
                childWidth = app.UIFigure.Position(3);
                childHeight = app.UIFigure.Position(4);
                
                % 3. 核心计算：让子窗口的中心对准主窗口的中心
                % 新的左边距 = 主窗口的左边距 + (主窗口宽度 - 子窗口宽度) / 2
                newLeft = mainPos(1) + (mainPos(3) - childWidth) / 2;
                
                % 新的下边距 = 主窗口的下边距 + (主窗口高度 - 子窗口高度) / 2
                newBottom = mainPos(2) + (mainPos(4) - childHeight) / 2;
                
                % 4. 强制刷新子窗口的位置
                app.UIFigure.Position = [newLeft, newBottom, childWidth, childHeight];
                
            else
                % 兜底防御：如果直接单独运行了这个子 App（没有主窗口），那就让它在屏幕中央显示
                movegui(app.UIFigure, 'center');
            end
%%
            app.dataNode = selectednode.NodeData;
            app.srate = app.dataNode.srate;
            tempdata = app.dataNode.data; 

            if ndims(tempdata) == 3
                app.Data   = tempdata;            % 保留三维 nCh x nT x nTrials
                app.trials = size(tempdata, 3);
            else
                app.Data   = tempdata;            % 二维 nCh x nT
                app.trials = 1;
            end
            app.sessionNode = selectednode.NodeData.parent;
            app.protocolNode=selectednode.NodeData.parent.parent;
            %需要转存才能访问
            app.originalheight=app.UIFigure.Position(4);
            app.originalwidth=app.UIFigure.Position(3);
            app.MinimumL2_NormUIFigure = [];
            app.SpatioUIFigure = [];
            app.PanelCache = containers.Map('KeyType','char','ValueType','any');
            app.AlgParams  = containers.Map('KeyType','char','ValueType','any');
        end

        % Menu selected function: ImportDataMenu
        function ImportDataMenuSelected(app, event)
            filterSpec = { ...
                '*.mat', 'MATLAB (.mat)'; ...
                '*.set', 'EEGLAB (.set)'; ...
                '*.*',   'All Files (*.*)'};

            [fileName, filePath,~] = uigetfile(filterSpec, ['Select ', strrep('Data', '_', ' '), ' File']);

            originalPath = fullfile(filePath, fileName);

            app.AdditionalDataNodeNum = app.AdditionalDataNodeNum + 1;
            app.AdditionalData{app.AdditionalDataNodeNum}.filename=fileName;
            app.AdditionalData{app.AdditionalDataNodeNum}.filepath=filePath;

            app.AdditionalDataNode{app.AdditionalDataNodeNum}= uitreenode(app.DataTree, ...
                'Text', ['Data: ', fileName], 'NodeData',originalPath);

        end

        % Callback function: AlgorithmsTree
        function AlgorithmsTreeCheckedNodesChanged(app, event)
            app.ensureSettingFiguresExist();   
            reg = seal_getAlgorithmRegistry();
            currentSet = cellstr(app.getOnlyLeafNodes());

            % 1) 删除：缓存里有但本次未勾选的
            cached = keys(app.PanelCache);
            for i = 1:numel(cached)
                if ~ismember(cached{i}, currentSet)
                    p = app.PanelCache(cached{i});
                    if isvalid(p), delete(p); end
                    remove(app.PanelCache, cached{i});
                end
            end

            % 2) 新建：本次勾选但缓存里没有的（且在 registry 里有定义）
            for i = 1:numel(currentSet)
                name = currentSet{i};
                if isKey(app.PanelCache, name), continue; end
                if ~isfield(reg, name), continue; end       % 没配置的算法跳过
                app.PanelCache(name) = app.buildPanelFromConfig(name, reg.(name));
            end

            % 3) 重新布局每个 figure 上的 panel
            app.relayoutAll();

       
        end

        % Button pushed function: RunSelectedAlgorithmsButton
        function RunSelectedAlgorithmsButtonPushed(app, event)
            % ========== 1. 准备工作 ==========
            app.SavePath = fullfile(app.sessionNode.path, "Results");
            if ~isfolder(app.SavePath), mkdir(app.SavePath); end

            % 判断数据维度
            isEpoched = (ndims(app.Data) == 3);
            if isEpoched
                [nCh, nT, nTrials] = size(app.Data);
            else
                [nCh, nT] = size(app.Data);
                nTrials = 1;
            end

            % 能量自动缩放（基于全体数据计算）
            if isEpoched
                X_all = reshape(app.Data, nCh, nT * nTrials);
            else
                X_all = app.Data;
            end
            avg_power = (norm(X_all,'fro')^2) / numel(X_all);
            fprintf('单点平均能量: %g\n', avg_power);
            fprintf('[DataScale] 自动缩放已禁用，保持原始数据单位。\n');

            % ========== 2. 白化与常用引用 ==========
   
            L_raw   = app.protocolNode.leadfieldNode.data;

            % --- 应用白化策略 ---
            [W, R, whInfo] = app.computeWhitener(X_all, app.NoiseStrategy);

            L = W * L_raw;
        
            fprintf('[Whitening] strategy=%s  baselineSamples=%d  cond(Cov_n)=%.3g\n', ...
                whInfo.strategy, whInfo.baselineSamples, whInfo.condNumber);

  
            % 维度校验
            nColL   = size(app.protocolNode.leadfieldNode.data, 2);
            nPhysSrc= size(app.protocolNode.cortexNode.data.Vertices, 1);
            nd = nColL / nPhysSrc;
            if nd ~= 1 && nd ~= 3
                error('SEAL:Algorithm:DimensionMismatch','导联场维度异常！');
            end
            app.NumOrientations = nd;

            % ========== 2. 常用引用（减少后面代码长度）==========
           
            C       = app.protocolNode.cortexNode.data;  
            Eye_N   =  eye(nCh);      % NoiseObservationMatrix
            nO      = app.NumOrientations;
            toLog = @(s) (islogical(s) && isscalar(s) && s) || ...
                (isnumeric(s) && isscalar(s) && s ~= 0) || ...
                ((ischar(s) || isstring(s)) && strcmpi(char(string(s)), 'true'));

            NodeDatas = app.getOnlyLeafNodes();
            if isempty(NodeDatas)
                uialert(app.UIFigure,'请先勾选至少一个算法。','未选择算法');
                return;
            end

            savePath = uigetdir(app.SavePath,'选择结果保存文件夹');
            if isequal(savePath,0), return; end

            app.pendingAlgNames = NodeDatas;
            decision = app.preCheckAndConfirm(NodeDatas, [], savePath);

            switch decision
                case 'abort'
                    return;
                case 'proceed_with_subset'
                    NodeDatas = app.pendingAlgNames;
                case 'proceed'
                    % 继续
            end
            % ========== 3. 进度条 ==========
            totalAlgs = numel(NodeDatas);
            totalSteps = totalAlgs * nTrials;
            hWait = uiprogressdlg(app.UIFigure, ...
                'Title','正在运行源分析...','Message','准备中...','Cancelable','on');

            succeeded = 0;
            failed = {};
            cancelled = false;

            % ========== 4. 主循环 ==========
            for i = 1:totalAlgs
                if ~isvalid(hWait) || hWait.CancelRequested || cancelled, break; end

                algName = char(NodeDatas(i));

                if app.AlgParams.isKey(algName)
                    prm = app.AlgParams(algName);
                else
                    prm = struct();
                end
                prm = app.fillAlgDefaults(algName, prm);

                batchID = sprintf('%s_%s_%s_%s', ...
                    app.dataNode.name, algName, ...
                    datestr(now,'yyyymmdd_HHMMSS'), ...
                    char(java.util.UUID.randomUUID.toString));

                for t = 1:nTrials
                    if ~isvalid(hWait) || hWait.CancelRequested || cancelled, break; end

                    stepIdx = (i-1)*nTrials + t;
                    hWait.Value = (stepIdx-1)/totalSteps;
                    if isEpoched
                        hWait.Message = sprintf('%s  trial %d/%d  (算法 %d/%d)', ...
                            algName, t, nTrials, i, totalAlgs);
                    else
                        hWait.Message = sprintf('正在运行 %s (%d/%d)...', algName, i, totalAlgs);
                    end

                    % 取当前 trial 的数据并白化
                    if isEpoched
                        X_t = W * app.Data(:, :, t);
                    else
                        X_t = W * app.Data;
                    end

                    % 跑算法
                    try
                        S = app.runOneAlgorithm(algName, X_t, L, C, R, Eye_N, nO, prm, toLog, ...
                            @() ~isvalid(hWait) || hWait.CancelRequested);
                    catch ME
                        if strcmp(ME.identifier, 'SEAL:Cancelled')
                            cancelled = true;
                            hWait.Message = '正在取消当前求解...';
                            break;
                        end
                        warning('算法 %s (trial %d) 运行失败: %s', algName, t, ME.message);
                        failed{end+1} = sprintf('%s trial%d: %s', algName, t, ME.message); %#ok<AGROW>
                        continue;
                    end

                    if isempty(S)
                        failed{end+1} = sprintf('%s trial%d: 未实现或返回空', algName, t); %#ok<AGROW>
                        continue;
                    end

                    % 后处理 + 保存（每个 trial 一个文件、一个 ResultNode）
                    result = app.postProcess(S, nO);   % nSrc x nT
                    
                    try
                        prmToSave = prm;
                        prmToSave.NoiseWhitening = whInfo;
                        if isEpoched
                            app.saveOneResult(algName, result, prmToSave, savePath, t, nTrials, batchID);
                        else
                            app.saveOneResult(algName, result, prmToSave, savePath, [], [], '');
                        end
                        succeeded = succeeded + 1;
                    catch ME
                        warning('保存失败 %s trial %d: %s', algName, t, ME.message);
                        failed{end+1} = sprintf('%s trial%d (保存): %s', algName, t, ME.message); %#ok<AGROW>
                    end
                end
                if cancelled
                    break;
                end
            end

            % ========== 5. 收尾 ==========
            cancelled = cancelled || (~isvalid(hWait) || hWait.CancelRequested);
            if isvalid(hWait)
                close(hWait);
            end
            app.Callingapp.bringToFront();
            app.Callingapp.flushFileTree();
            nFailed = numel(failed);
            totalSteps_actual = totalAlgs * nTrials;   % 上面已经算过 totalSteps，复用即可

            if cancelled
                iconType  = 'info';
                titleText = '已取消';
                headLine  = '源分析已按请求取消。';
            elseif nFailed == 0
                iconType  = 'success';
                titleText = '运行完成';
                headLine  = '源分析全部完成！';
            elseif succeeded == 0
                iconType  = 'error';
                titleText = '运行失败';
                headLine  = '源分析全部失败！';
            else
                iconType  = 'warning';
                titleText = '部分失败';
                headLine  = '源分析部分完成（有失败项，请检查）';
            end

            msg = sprintf('%s\n成功: %d / %d\n结果目录:\n%s', ...
                headLine, succeeded, totalSteps_actual, savePath);
            if ~isempty(failed)
                msg = sprintf('%s\n\n失败算法:\n  • %s', msg, strjoin(failed, sprintf('\n  • ')));
            end
            uialert(app.UIFigure, msg, titleText, ...
                'Icon', iconType, 'CloseFcn', @(~,~) closeAllWindows(app));

        end

        % Value changed function: DepthWeightingEditField
        function DepthWeightingEditFieldValueChanged(app, event)
            value = app.DepthWeightingEditField.Value;
            app.DepthWeightExp=value;
        end

        % Value changed function: NoiseWhiteningDropDown
        function NoiseWhiteningDropDownValueChanged(app, event)
            app.NoiseStrategy = app.NoiseWhiteningDropDown.Value;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 398 300];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.RowHeight = {'1x'};

            % Create AlgorithmsPanel
            app.AlgorithmsPanel = uipanel(app.GridLayout);
            app.AlgorithmsPanel.AutoResizeChildren = 'off';
            app.AlgorithmsPanel.Title = 'Algorithms';
            app.AlgorithmsPanel.Layout.Row = 1;
            app.AlgorithmsPanel.Layout.Column = 1;
            app.AlgorithmsPanel.FontWeight = 'bold';

            % Create AlgorithmsTree
            app.AlgorithmsTree = uitree(app.AlgorithmsPanel, 'checkbox');
            app.AlgorithmsTree.Position = [6 1 173 253];

            % Create MinimumL2_NormNode
            app.MinimumL2_NormNode = uitreenode(app.AlgorithmsTree);
            app.MinimumL2_NormNode.Text = 'Minimum L2_Norm';

            % Create MNENode
            app.MNENode = uitreenode(app.MinimumL2_NormNode);
            app.MNENode.NodeData = 'MNE';
            app.MNENode.Text = 'MNE';

            % Create LOGURANode
            app.LOGURANode = uitreenode(app.MinimumL2_NormNode);
            app.LOGURANode.NodeData = 'LOGURA';
            app.LOGURANode.Text = 'LOGURA';

            % Create dSPMNode
            app.dSPMNode = uitreenode(app.MinimumL2_NormNode);
            app.dSPMNode.NodeData = 'dSPM';
            app.dSPMNode.Text = 'dSPM';

            % Create sLORETANode
            app.sLORETANode = uitreenode(app.MinimumL2_NormNode);
            app.sLORETANode.NodeData = 'sLORETA';
            app.sLORETANode.Text = 'sLORETA';

            % Create eLORETANode
            app.eLORETANode = uitreenode(app.MinimumL2_NormNode);
            app.eLORETANode.NodeData = 'eLORETA';
            app.eLORETANode.Text = 'eLORETA';

            % Create LORETANode
            app.LORETANode = uitreenode(app.MinimumL2_NormNode);
            app.LORETANode.NodeData = 'LORETA';
            app.LORETANode.Text = 'LORETA';

            % Create LAURANode
            app.LAURANode = uitreenode(app.MinimumL2_NormNode);
            app.LAURANode.NodeData = 'LAURA';
            app.LAURANode.Text = 'LAURA';

            % Create BeamformerNode
            app.BeamformerNode = uitreenode(app.AlgorithmsTree);
            app.BeamformerNode.Text = 'Beamformer';

            % Create LCMVNode
            app.LCMVNode = uitreenode(app.BeamformerNode);
            app.LCMVNode.NodeData = 'LCMV';
            app.LCMVNode.Text = 'LCMV';

            % Create MCMVNode
            app.MCMVNode = uitreenode(app.BeamformerNode);
            app.MCMVNode.NodeData = 'MCMV';
            app.MCMVNode.Text = 'MCMV';

            % Create SparsemodelNode
            app.SparsemodelNode = uitreenode(app.AlgorithmsTree);
            app.SparsemodelNode.Text = 'Sparse model';

            % Create Lq_ADMMNode
            app.Lq_ADMMNode = uitreenode(app.SparsemodelNode);
            app.Lq_ADMMNode.NodeData = 'Lq_ADMM';
            app.Lq_ADMMNode.Text = 'Lq_ADMM';

            % Create ChampagneNode
            app.ChampagneNode = uitreenode(app.SparsemodelNode);
            app.ChampagneNode.NodeData = 'Champagne';
            app.ChampagneNode.Text = 'Champagne';

            % Create BlkChamNode
            app.BlkChamNode = uitreenode(app.SparsemodelNode);
            app.BlkChamNode.NodeData = 'BlockChampagne';
            app.BlkChamNode.Text = 'Blk-Cham';

            % Create TSChamNode
            app.TSChamNode = uitreenode(app.SparsemodelNode);
            app.TSChamNode.NodeData = 'TS_Champagne';
            app.TSChamNode.Text = 'TS-Cham';

            % Create SmoothChamNode
            app.SmoothChamNode = uitreenode(app.SparsemodelNode);
            app.SmoothChamNode.NodeData = 'SmoothChampagne';
            app.SmoothChamNode.Text = 'Smooth-Cham';

            % Create IRESNode
            app.IRESNode = uitreenode(app.SparsemodelNode);
            app.IRESNode.NodeData = 'IRES';
            app.IRESNode.Text = 'IRES';

            % Create FOCUSSNode
            app.FOCUSSNode = uitreenode(app.SparsemodelNode);
            app.FOCUSSNode.NodeData = 'FOCUSS';
            app.FOCUSSNode.Text = 'FOCUSS';

            % Create SpatioTemporalNode
            app.SpatioTemporalNode = uitreenode(app.AlgorithmsTree);
            app.SpatioTemporalNode.Text = 'SpatioTemporal';

            % Create MxNENode
            app.MxNENode = uitreenode(app.SpatioTemporalNode);
            app.MxNENode.NodeData = 'MxNE';
            app.MxNENode.Text = 'MxNE';

            % Create irMxNENode
            app.irMxNENode = uitreenode(app.SpatioTemporalNode);
            app.irMxNENode.NodeData = 'irMxNE';
            app.irMxNENode.Text = 'irMxNE';

            % Create SISSY_L21Node
            app.SISSY_L21Node = uitreenode(app.SpatioTemporalNode);
            app.SISSY_L21Node.NodeData = 'SISSY_L21';
            app.SISSY_L21Node.Text = 'SISSY_L21';

            % Create VSSI_LpNode
            app.VSSI_LpNode = uitreenode(app.SpatioTemporalNode);
            app.VSSI_LpNode.NodeData = 'VSSI_Lp';
            app.VSSI_LpNode.Text = 'VSSI_Lp';

            % Create TFMxNENode
            app.TFMxNENode = uitreenode(app.SpatioTemporalNode);
            app.TFMxNENode.NodeData = 'TFMxNE';
            app.TFMxNENode.Text = 'TFMxNE';

            % Create DebiasedGroupLassoNode
            app.DebiasedGroupLassoNode = uitreenode(app.SpatioTemporalNode);
            app.DebiasedGroupLassoNode.NodeData = 'debiased_GroupLasso';
            app.DebiasedGroupLassoNode.Text = 'DeESI';

            % Create dSPNNode
            app.dSPNNode = uitreenode(app.SpatioTemporalNode);
            app.dSPNNode.NodeData = 'dSPN';
            app.dSPNNode.Text = 'dSPN';

            % Create FastIRESNode
            app.FastIRESNode = uitreenode(app.SpatioTemporalNode);
            app.FastIRESNode.NodeData = 'FASTIRES';
            app.FastIRESNode.Text = 'FASTIRES';

            % Create STOUTNode
            app.STOUTNode = uitreenode(app.SpatioTemporalNode);
            app.STOUTNode.NodeData = 'STOUT';
            app.STOUTNode.Text = 'STOUT';

            % Create dMAPNode
            app.dMAPNode = uitreenode(app.SpatioTemporalNode);
            app.dMAPNode.NodeData = 'dMAP';
            app.dMAPNode.Text = 'dMAP';

            % Create BESTIESNode
            app.BESTIESNode = uitreenode(app.SpatioTemporalNode);
            app.BESTIESNode.NodeData = 'BESTIES';
            app.BESTIESNode.Text = 'BESTIES';

            % Create SISTBFNode
            app.SISTBFNode = uitreenode(app.SpatioTemporalNode);
            app.SISTBFNode.NodeData = 'STBFSI';
            app.SISTBFNode.Text = 'SI-STBF';

            % Create uSTARNode
            app.uSTARNode = uitreenode(app.SpatioTemporalNode);
            app.uSTARNode.NodeData = 'uSTAR';
            app.uSTARNode.Text = 'uSTAR';

            % Create STARTSNode
            app.STARTSNode = uitreenode(app.SpatioTemporalNode);
            app.STARTSNode.NodeData = 'STARTS';
            app.STARTSNode.Text = 'STARTS';

            % Create sSTARTSNode
            app.sSTARTSNode = uitreenode(app.SpatioTemporalNode);
            app.sSTARTSNode.NodeData = 'sSTARTS';
            app.sSTARTSNode.Text = 'sSTARTS';

            % Create VSSI_GGDNode
            app.VSSI_GGDNode = uitreenode(app.SpatioTemporalNode);
            app.VSSI_GGDNode.NodeData = 'VSSI_GGD';
            app.VSSI_GGDNode.Text = 'VSSI_GGD';

            % Assign Checked Nodes
            app.AlgorithmsTree.CheckedNodesChangedFcn = createCallbackFcn(app, @AlgorithmsTreeCheckedNodesChanged, true);

            % Create GeneralSettingsPanel
            app.GeneralSettingsPanel = uipanel(app.GridLayout);
            app.GeneralSettingsPanel.AutoResizeChildren = 'off';
            app.GeneralSettingsPanel.Title = 'General Settings';
            app.GeneralSettingsPanel.Layout.Row = 1;
            app.GeneralSettingsPanel.Layout.Column = 2;
            app.GeneralSettingsPanel.FontWeight = 'bold';

            % Create DepthWeightingEditFieldLabel
            app.DepthWeightingEditFieldLabel = uilabel(app.GeneralSettingsPanel);
            app.DepthWeightingEditFieldLabel.FontSize = 10;
            app.DepthWeightingEditFieldLabel.Position = [10 225 94 22];
            app.DepthWeightingEditFieldLabel.Text = 'Depth Weighting';

            % Create DepthWeightingEditField
            app.DepthWeightingEditField = uieditfield(app.GeneralSettingsPanel, 'numeric');
            app.DepthWeightingEditField.ValueChangedFcn = createCallbackFcn(app, @DepthWeightingEditFieldValueChanged, true);
            app.DepthWeightingEditField.Position = [110 225 57 22];
            app.DepthWeightingEditField.Value = 0.5;

            % Create RunSelectedAlgorithmsButton
            app.RunSelectedAlgorithmsButton = uibutton(app.GeneralSettingsPanel, 'push');
            app.RunSelectedAlgorithmsButton.ButtonPushedFcn = createCallbackFcn(app, @RunSelectedAlgorithmsButtonPushed, true);
            app.RunSelectedAlgorithmsButton.Position = [9 111 158 29];
            app.RunSelectedAlgorithmsButton.Text = 'Run Selected Algorithms';

            % Create NoiseWhiteningDropDown_2Label
            app.NoiseWhiteningDropDown_2Label = uilabel(app.GeneralSettingsPanel);
            app.NoiseWhiteningDropDown_2Label.FontSize = 10;
            app.NoiseWhiteningDropDown_2Label.Position = [9 196 93 22];
            app.NoiseWhiteningDropDown_2Label.Text = 'Noise Whitening';

            % Create NoiseWhiteningDropDown
            app.NoiseWhiteningDropDown = uidropdown(app.GeneralSettingsPanel);
            app.NoiseWhiteningDropDown.Items = {'brainstorm', 'identity', 'Precomputed', 'data_cov', 'baseline', 'pass_through'};
            app.NoiseWhiteningDropDown.ValueChangedFcn = createCallbackFcn(app, @NoiseWhiteningDropDownValueChanged, true);
            app.NoiseWhiteningDropDown.FontSize = 10;
            app.NoiseWhiteningDropDown.Position = [107 196 60 22];
            app.NoiseWhiteningDropDown.Value = 'brainstorm';

            % Create SelectPriorButton
            app.SelectPriorButton = uibutton(app.UIFigure, 'push');
            app.SelectPriorButton.Position = [460 150 78 22];
            app.SelectPriorButton.Text = 'Select Prior';

            % Create LoadPriorButton
            app.LoadPriorButton = uibutton(app.UIFigure, 'push');
            app.LoadPriorButton.Position = [462 109 76 22];
            app.LoadPriorButton.Text = 'Load Prior';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create ImportDataMenu
            app.ImportDataMenu = uimenu(app.ContextMenu);
            app.ImportDataMenu.MenuSelectedFcn = createCallbackFcn(app, @ImportDataMenuSelected, true);
            app.ImportDataMenu.Text = 'Import Data';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SEAL_selectAlgorithm(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

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
