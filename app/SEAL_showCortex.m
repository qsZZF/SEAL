classdef SEAL_showCortex < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        TimePointEditField   matlab.ui.control.NumericEditField
        EditFieldLabel       matlab.ui.control.Label
        RegionTree           matlab.ui.container.CheckBoxTree
        RegionsNode          matlab.ui.container.TreeNode
        Label_2              matlab.ui.control.Label
        Label                matlab.ui.control.Label
        ResectXYZPanel       matlab.ui.container.Panel
        SliderZ              matlab.ui.control.Slider
        SliderY              matlab.ui.control.Slider
        SliderX              matlab.ui.control.Slider
        ResetButton          matlab.ui.control.Button
        StructureButton      matlab.ui.control.Button
        RightButton          matlab.ui.control.Button
        LeftButton           matlab.ui.control.Button
        RegionsPanel         matlab.ui.container.Panel
        ButtonGroup_3        matlab.ui.container.ButtonGroup
        jump2signalButton    matlab.ui.control.Button
        SelectpointButton    matlab.ui.control.StateButton
        ButtonGroup_2        matlab.ui.container.ButtonGroup
        HideFrameButton      matlab.ui.control.ToggleButton
        ShowFrameButten      matlab.ui.control.ToggleButton
        ButtonGroup          matlab.ui.container.ButtonGroup
        AtlasDrop            matlab.ui.control.DropDown
        NoColorButton        matlab.ui.control.ToggleButton
        HalfColorButton      matlab.ui.control.ToggleButton
        FullColorButten      matlab.ui.control.ToggleButton
        SurfaceoprionsPanel  matlab.ui.container.Panel
        EdgeButton           matlab.ui.control.StateButton
        ColorButton          matlab.ui.control.Button
        TranspSlider         matlab.ui.control.Slider
        TranspLabel          matlab.ui.control.Label
        ThresholdSlider      matlab.ui.control.Slider
        ThresholdLabel       matlab.ui.control.Label
        ThresholdEdit        matlab.ui.control.NumericEditField
        SmoothSlider         matlab.ui.control.Slider
        SmoothLabel          matlab.ui.control.Label
        Transpedit           matlab.ui.control.NumericEditField
        smoothedit           matlab.ui.control.NumericEditField
    end


    properties (Access = private)
        CallingApp
        cortexData
        smoothvalue=30
        smoothedVC
        transp = 1          %透明度
        handles1            %弹窗图像句柄
        isClosing = false;  % 标记是否正在关闭中
        hemisphere=''
        plothemIdx
        originalVertexColors
        originalVertices
        originalAxisLimits
        OffsetVertices ={}   %每个region的vertices
        sepVertices          %存储strcture后的vertices
        sepregionVertices={}  %存储strcture后的每个region的vertices
        alreadystructure = 0 % 0 未structure,1 已structure
        selectedVertices = []
        selectedVerticesNum = 0
        highlightedVertices=struct('vertexIndex', 0, 'originalColors', []);
        regioncolorset
        AtlasIdx = 1
        alphasignal=0.4;
        edgealphasignal=1;
        resultData
        PoriodData
        srate
        TimeVector = []
        SourceDisplayMode = 'signed'
        SourceDataThreshold = []
        SourceColorLimit = []
        SourceColorMap = []
        SourceDisplayScale = 1
        SourceDisplayUnits = ''
        CurrentSourceData = []
        CurrentSourceRGB = []
        CurrentUseGlobalLimit = true

        figOffset       % 图像窗口相对 App 的偏移 [dx, dy]
        lastFigPos      % 上次图像窗口位置
        lastAppPos      % 上次 App 位置
        posTimer        % 位置监听定时器
        isSyncing = false  % 防止递归

        IsTrialSplit   logical = false     % 是否分 trial 结果
        TrialRange     double  = []        % 采样点闭区间 [s1 s2]（全局时间轴上）
        TrialIndex     double  = []        % 本窗口所属 trial 编号
        TrialCount     double  = []        % 总 trial 数
        BatchID        char    = ''        % 同批次标识（备用）
        BaseTitle      char    = ''        % 初始标题，复原用
        OutOfTrialMsg  char    = ''        % 当前外显的"out of trial"文字，避免重复刷新

    end


    methods (Access = private)

        function timeVector = resolveResultTimeVector(app, node, params, nTime)
            timeVector = [];

            if ~isempty(params) && isstruct(params)
                if isfield(params, 'TimeVector') && numel(params.TimeVector) == nTime
                    timeVector = double(params.TimeVector(:))';
                    return;
                elseif isfield(params, 'Time') && numel(params.Time) == nTime
                    timeVector = double(params.Time(:))';
                    return;
                elseif isfield(params, 'TimeStart') && isfield(params, 'srate')
                    timeVector = double(params.TimeStart) + (0:nTime-1) ./ double(params.srate);
                    return;
                end
            end

            try
                session = node.parent;
                if ~isempty(session) && isprop(session, 'children')
                    for ii = 1:session.childCount
                        child = session.children(ii);
                        if ~isa(child, 'DataNode') || isempty(child.dataInfo) || ~isprop(child.dataInfo, 'dataPath')
                            continue;
                        end
                        dataPath = string(child.dataInfo.dataPath);
                        if strlength(dataPath) == 0 || ~isfile(dataPath)
                            continue;
                        end
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
                    end
                end
            catch
                timeVector = [];
            end

            timeVector = (0:nTime-1) ./ app.srate;
        end

        function colIdx = timeToColumn(app, selectTime)
            if ~isempty(app.TimeVector) && numel(app.TimeVector) == size(app.resultData, 2)
                [~, colIdx] = min(abs(app.TimeVector - selectTime));
            else
                colIdx = round(selectTime * app.srate) + 1;
            end
            colIdx = max(1, min(colIdx, size(app.resultData, 2)));
        end

        function configureSourceDisplay(app, resultNode)
            if nargin < 2
                resultNode = [];
            end
            app.SourceDisplayMode = 'brainstorm';
            app.SourceDataThreshold = 0.5;
            app.SourceDisplayScale = 1e12;
            app.SourceDisplayUnits = 'pA.m';
            maxVal = app.brainstormSourceDisplayMax(resultNode);
            app.SourceColorLimit = [0 maxVal * app.SourceDisplayScale];
            app.SourceColorMap = seal_brainstorm_colormap(256);
        end

        function maxVal = brainstormSourceDisplayMax(app, resultNode)
            maxVal = [];

            dataMatrix = app.findMatchingDataMatrix(resultNode);
            if ~isempty(dataMatrix) && size(dataMatrix, 2) == size(app.resultData, 2)
                try
                    gfp = sum(double(dataMatrix).^2, 1);
                    [~, iMax] = max(gfp);
                    maxVal = max(abs(app.resultData(:, iMax)));
                catch
                    maxVal = [];
                end
            end

            if isempty(maxVal) || ~isfinite(maxVal) || maxVal <= 0
                finiteValues = abs(app.resultData(isfinite(app.resultData)));
                if isempty(finiteValues)
                    maxVal = eps;
                else
                    maxVal = max(finiteValues(:));
                    if maxVal <= 0
                        maxVal = eps;
                    end
                end
            end
        end

        function dataMatrix = findMatchingDataMatrix(app, resultNode)
            dataMatrix = [];
            if isempty(resultNode) || ~isprop(resultNode, 'parent') || isempty(resultNode.parent)
                return;
            end

            params = [];
            try
                params = resultNode.params;
            catch
                params = [];
            end

            parentName = "";
            if isstruct(params) && isfield(params, 'parentDataName')
                parentName = string(params.parentDataName);
            end

            session = resultNode.parent;
            try
                for ii = 1:session.childCount
                    child = session.children(ii);
                    if ~isa(child, 'DataNode')
                        continue;
                    end
                    if strlength(parentName) > 0 && isprop(child, 'name') && string(child.name) ~= parentName
                        continue;
                    end

                    candidate = double(child.data);
                    if ndims(candidate) == 3
                        [nChan, nTime, nTrial] = size(candidate);
                        candidate = reshape(candidate, nChan, nTime * nTrial);
                    end
                    if ismatrix(candidate) && size(candidate, 2) == size(app.resultData, 2)
                        dataMatrix = candidate;
                        return;
                    end
                end
            catch
                dataMatrix = [];
            end
        end

        function [sourceData, colorMap, colorLimit, threshold] = sourceDisplayArgs(app, sourceData, useGlobalLimit)
            if nargin < 3
                useGlobalLimit = true;
            end

            colorMap = [];
            colorLimit = [];
            threshold = [];

            if strcmpi(app.SourceDisplayMode, 'brainstorm')
                sourceData = abs(sourceData) * app.SourceDisplayScale;
                colorMap = app.SourceColorMap;
                if isempty(colorMap)
                    colorMap = seal_brainstorm_colormap(256);
                end
                threshold = app.SourceDataThreshold;
                if isempty(threshold)
                    threshold = 0.5;
                end

                if useGlobalLimit && ~isempty(app.SourceColorLimit)
                    colorLimit = app.SourceColorLimit;
                else
                    finiteValues = sourceData(isfinite(sourceData));
                    if isempty(finiteValues)
                        maxVal = eps;
                    else
                        maxVal = max(abs(finiteValues(:)));
                        if maxVal <= 0
                            maxVal = eps;
                        end
                    end
                    colorLimit = [0 maxVal];
                end
            end
        end

        function updateSourceColorbarLabel(app)
            try
                if isfield(app.handles1, 'hColorbar') && isvalid(app.handles1.hColorbar) && ...
                        ~isempty(app.SourceDisplayUnits)
                    hText = xlabel(app.handles1.hColorbar, app.SourceDisplayUnits);
                    set(hText, 'Color', [1 1 1]);
                end
            catch
            end
        end

        function redrawSourceData(app, sourceData, rgb, useGlobalLimit)
            if nargin < 3
                rgb = [];
            end
            if nargin < 4
                useGlobalLimit = true;
            end

            app.CurrentSourceData = sourceData;
            app.CurrentSourceRGB = rgb;
            app.CurrentUseGlobalLimit = useGlobalLimit;

            [sourceData, colorMap, colorLimit, threshold] = app.sourceDisplayArgs(sourceData, useGlobalLimit);
            plotArgs = {};
            if ~isempty(colorLimit)
                plotArgs = [plotArgs, {'clim', colorLimit}];
            end
            if ~isempty(threshold)
                plotArgs = [plotArgs, {'threshold', threshold}];
            end

            [Valpha, cdata] = seal_bulidsource(sourceData, ...
                app.cortexData.Vertices, rgb, colorMap, app.handles1.axes, plotArgs{:});
            set(app.handles1.hp, 'FaceAlpha','flat', 'AlphaDataMapping','none', ...
                'FaceVertexAlphaData', Valpha, ...
                'FaceVertexCData',     cdata, ...
                'facecolor', 'interp');
            app.updateSourceColorbarLabel();
            app.handles1.h.UserData = 0;
            set(app.handles1.h, 'CloseRequestFcn', @(src, evt) seal_showCortexClose(src, evt));
        end

        function redrawSourceAtColumn(app, colIdx, rgb)
            if nargin < 3
                rgb = [];
            end
            app.redrawSourceData(app.resultData(:, colIdx), rgb);
        end

        function value = clampSourceThreshold(app, value)
            if isempty(value) || ~isnumeric(value) || ~all(isfinite(value(:)))
                value = app.SourceDataThreshold;
            end
            if isempty(value) || ~isnumeric(value) || ~all(isfinite(value(:)))
                value = 0.5;
            end
            value = max(0, min(1, double(value(1))));
        end

        function syncSourceThresholdControls(app)
            value = app.clampSourceThreshold(app.SourceDataThreshold);
            if ~isempty(app.ThresholdSlider) && isvalid(app.ThresholdSlider)
                app.ThresholdSlider.Value = value;
            end
            if ~isempty(app.ThresholdEdit) && isvalid(app.ThresholdEdit)
                app.ThresholdEdit.Value = value;
            end
        end

        function setSourceThresholdControlsEnabled(app, isEnabled)
            if isEnabled
                enableState = 'on';
            else
                enableState = 'off';
            end
            if ~isempty(app.ThresholdSlider) && isvalid(app.ThresholdSlider)
                app.ThresholdSlider.Enable = enableState;
            end
            if ~isempty(app.ThresholdEdit) && isvalid(app.ThresholdEdit)
                app.ThresholdEdit.Enable = enableState;
            end
        end

        function refreshCurrentSourceDisplay(app)
            if isempty(app.CurrentSourceData) || isempty(app.handles1) || ...
                    ~isfield(app.handles1, 'hp') || ~isvalid(app.handles1.hp)
                return;
            end
            app.redrawSourceData(app.CurrentSourceData, ...
                app.CurrentSourceRGB, app.CurrentUseGlobalLimit);
        end

        function setSourceThreshold(app, value, sourceControl)
            if nargin < 3
                sourceControl = '';
            end
            app.SourceDataThreshold = app.clampSourceThreshold(value);
            if ~strcmp(sourceControl, 'slider') && ...
                    ~isempty(app.ThresholdSlider) && isvalid(app.ThresholdSlider)
                app.ThresholdSlider.Value = app.SourceDataThreshold;
            end
            if ~strcmp(sourceControl, 'edit') && ...
                    ~isempty(app.ThresholdEdit) && isvalid(app.ThresholdEdit)
                app.ThresholdEdit.Value = app.SourceDataThreshold;
            end
            app.refreshCurrentSourceDisplay();
        end

        function stopPositionTimer(app)
            try
                if ~isempty(app.posTimer) && isvalid(app.posTimer)
                    stop(app.posTimer);
                    delete(app.posTimer);
                end
            catch
            end
            app.posTimer = [];
        end

        function deleteSourceFigure(app)
            try
                if isstruct(app.handles1) && isfield(app.handles1, 'h') && ...
                        ~isempty(app.handles1.h) && isvalid(app.handles1.h)
                    set(app.handles1.h, 'CloseRequestFcn', '');
                    delete(app.handles1.h);
                end
            catch
            end
        end

        function deleteControlFigure(app)
            try
                if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                    set(app.UIFigure, 'CloseRequestFcn', '');
                    delete(app.UIFigure);
                end
            catch
            end
        end

        function RightorLeft(app , h )
            if isfield(app.handles1, 'Multimarker') && ~isempty(app.handles1.Multimarker)
                for i = 1:length(app.handles1.Multimarker)
                    if isvalid(app.handles1.Multimarker(i))
                        delete(app.handles1.Multimarker(i));
                    end
                end
                app.handles1.Multimarker = [];
            end
            % 清除旧标记
            if isfield(app.handles1, 'singlemarker')
                delete(app.handles1.singlemarker);
            end
            vc=app.cortexData.Vertices;
            app.originalVertices = vc;
            alphaValue = app.transp;
            tIdx = strcmp('Structures',string(struct2table(app.cortexData.Atlas).Name));
            switch h
                case 'reset'
                    Valpha = ones(size(vc,1),1)*alphaValue;
                    app.plothemIdx = 1:size(vc,1);
                    if isfield(app.handles1, 'regionPatches') && ~isempty(app.handles1.regionPatches)
                        for i = 1:length(app.handles1.regionPatches)
                            if isvalid(app.handles1.regionPatches(i))
                                delete(app.handles1.regionPatches(i));
                            end
                        end
                        app.handles1.regionPatches = [];
                    end
                    app.RegionTree.CheckedNodes=[];
                    set(app.handles1.hp, 'FaceVertexAlphaData', Valpha, 'Vertices', app.smoothedVC, ...
                        'FaceColor',[0.9,0.7,0.7]);
                    app.alreadystructure = 0 ;
                case'structrue'
                    % ---- 1) 取左右半球索引（已存在的代码） ----
                    leftIdx  = app.cortexData.Atlas(tIdx).Scouts(1).Vertices;
                    rightIdx = app.cortexData.Atlas(tIdx).Scouts(2).Vertices;

                    % 用“未分离”的顶点做判定
                    V0 = app.smoothedVC;              % 或 app.originalVertices (两者一致即可)
                    Lc = mean(V0(leftIdx,  :), 1);    % 左半球质心 [x y z]
                    Rc = mean(V0(rightIdx, :), 1);    % 右半球质心 [x y z]
                    d  = Lc - Rc;                     % 质心差向量

                    % ---- 2) 自动选择“左右轴”：哪个坐标分量差异最大 ----
                    [~, axisIdx] = max(abs(d));       % axisIdx = 1/2/3 → X/Y/Z

                    % ---- 3) 分离距离与方向（让左右沿本来就分开的方向继续拉开）----
                    % 以选中轴上的范围做基准（0.5 可按需调整）
                    sepRange = max(V0(:, axisIdx))-min(V0(:, axisIdx));
                    separationDistance = max(sepRange, eps) * 0.5;

                    % 方向：让“质心更大的一侧”朝正向移动
                    sgn = sign(d(axisIdx));           % d>0 说明左在正向一侧
                    if sgn==0, sgn = 1; end

                    leftMove  = [0 0 0];  rightMove = [0 0 0];
                    leftMove(axisIdx)  = +sgn * separationDistance/2;
                    rightMove(axisIdx) = -sgn * separationDistance/2;

                    % ---- 4) 应用到主表面 ----
                    separatedVertices = V0;  % 从未分离版本出发，避免累计偏移
                    separatedVertices(leftIdx,  :) = separatedVertices(leftIdx,  :) + leftMove;
                    separatedVertices(rightIdx, :) = separatedVertices(rightIdx, :) + rightMove;

                    Valpha = ones(size(V0,1),1) * alphaValue;
                    set(app.handles1.hp, 'Vertices', separatedVertices, 'FaceVertexAlphaData', Valpha);
                    app.sepVertices = separatedVertices;

                   % ---- 5) 若有已勾选的 Region，同步更新（恢复完整双半球 Faces）----
                   if ~isempty(app.RegionTree.CheckedNodes)
                       checkednodes = app.popregionnode(app.RegionTree.CheckedNodes);

                       Fall = app.cortexData.Faces;
                       for i = 1:numel(checkednodes)
                           nodeIndex  = checkednodes(i).NodeData;
                           scoutVerts = app.cortexData.Atlas(app.AtlasIdx).Scouts(nodeIndex).Vertices;

                           % 恢复该 scout 的完整 Faces（双半球）
                           facesMask   = any(ismember(Fall, scoutVerts), 2);
                           regionFaces = Fall(facesMask, :);

                           % 从"未分离 + 带法向偏移"的全脑顶点出发，再对左右半球做 structure 位移
                           if isempty(app.OffsetVertices) || numel(app.OffsetVertices) < nodeIndex ...
                                   || isempty(app.OffsetVertices{nodeIndex})
                               % 保险：如果没缓存，直接用 smoothedVC
                               RV = app.smoothedVC;
                           else
                               RV = app.OffsetVertices{nodeIndex};
                           end
                           RV(leftIdx,  :) = RV(leftIdx,  :) + leftMove;
                           RV(rightIdx, :) = RV(rightIdx, :) + rightMove;

                           h = app.handles1.regionPatches(nodeIndex);
                           if isvalid(h)
                               set(h, 'Vertices', RV, 'Faces', regionFaces);   % ← 关键：同时恢复 Faces
                           end
                           app.sepregionVertices{nodeIndex} = RV;
                       end
                   end

                    % ---- 6) 计算新坐标轴范围（带 10% 边距）----
                    % 计算新顶点范围
                    minVals = min(separatedVertices);
                    maxVals = max(separatedVertices);

                    % 计算扩展范围（增加10%的边界）
                    rangeX = maxVals(1) - minVals(1);
                    rangeY = maxVals(2) - minVals(2);
                    rangeZ = maxVals(3) - minVals(3);

                    paddingX = rangeX * 0.1;
                    paddingY = rangeY * 0.1;
                    paddingZ = rangeZ * 0.1;

                    newLimits = [minVals(1)-paddingX, maxVals(1)+paddingX, ...
                        minVals(2)-paddingY, maxVals(2)+paddingY, ...
                        minVals(3)-paddingZ, maxVals(3)+paddingZ];

                    % 设置新坐标轴范围
                    axis(app.handles1.axes, newLimits);
                case 'left'

                    app.plothemIdx = app.cortexData.Atlas(tIdx).Scouts(1).Vertices;

                    Fall  = app.cortexData.Faces;
                    keepF = all(ismember(Fall, app.plothemIdx), 2);
                    Fhem  = Fall(keepF, :);

                    if app.alreadystructure == 1
                        V = app.sepVertices;
                    else
                        V = app.smoothedVC;
                    end

                    % 只改几何，不碰 alpha / FaceColor / AlphaDataMapping
                    set(app.handles1.hp, 'Vertices', V, 'Faces', Fhem);

                    % region 同步也用 Faces 过滤
                    if ~isempty(app.RegionTree.CheckedNodes)
                        checkednodes = app.popregionnode(app.RegionTree.CheckedNodes);
                        for i = 1:length(checkednodes)
                            nodeIndex   = checkednodes(i).NodeData;
                            scoutVerts  = app.cortexData.Atlas(app.AtlasIdx).Scouts(nodeIndex).Vertices;

                            facesMask   = any(ismember(Fall, scoutVerts), 2);
                            regionFaces = Fall(facesMask, :);
                            regionFaces = regionFaces(all(ismember(regionFaces, app.plothemIdx), 2), :);

                            if app.alreadystructure == 1
                                RV = app.sepregionVertices{nodeIndex};
                            else
                                RV = app.OffsetVertices{nodeIndex};
                            end
                            set(app.handles1.regionPatches(nodeIndex), ...
                                'Vertices', RV, 'Faces', regionFaces);
                        end
                    end


                case 'right'
                    app.plothemIdx = app.cortexData.Atlas(tIdx).Scouts(2).Vertices;

                    Fall  = app.cortexData.Faces;
                    keepF = all(ismember(Fall, app.plothemIdx), 2);
                    Fhem  = Fall(keepF, :);

                    if app.alreadystructure == 1
                        V = app.sepVertices;
                    else
                        V = app.smoothedVC;
                    end

                    % 只改几何，不碰 alpha / FaceColor / AlphaDataMapping
                    set(app.handles1.hp, 'Vertices', V, 'Faces', Fhem);

                    % region 同步也用 Faces 过滤
                    if ~isempty(app.RegionTree.CheckedNodes)
                        checkednodes = app.popregionnode(app.RegionTree.CheckedNodes);
                        for i = 1:length(checkednodes)
                            nodeIndex   = checkednodes(i).NodeData;
                            scoutVerts  = app.cortexData.Atlas(app.AtlasIdx).Scouts(nodeIndex).Vertices;

                            facesMask   = any(ismember(Fall, scoutVerts), 2);
                            regionFaces = Fall(facesMask, :);
                            regionFaces = regionFaces(all(ismember(regionFaces, app.plothemIdx), 2), :);

                            if app.alreadystructure == 1
                                RV = app.sepregionVertices{nodeIndex};
                            else
                                RV = app.OffsetVertices{nodeIndex};
                            end
                            set(app.handles1.regionPatches(nodeIndex), ...
                                'Vertices', RV, 'Faces', regionFaces);
                        end
                    end
            end
        end

        function vertexClicked(app, src, event)
            % 获取点击位置
            clickedPoint = event.IntersectionPoint;

            %  获取所有顶点坐标
            vertices = get(src, 'Vertices');

            %  找到最近的顶点
            [vertexIndex, ~] = app.findClosestVertex(vertices, clickedPoint);

            isCtrlPressed = ismember('control', get(app.handles1.h, 'CurrentModifier'));
            isAltPressed = ismember('alt', get(app.handles1.h, 'CurrentModifier'));

            if isAltPressed
                % 如果当前这个点在已选集合中，就移除
                if ismember(vertexIndex, app.selectedVertices)
                    app.selectedVertices = setdiff(app.selectedVertices, vertexIndex);
                end

                if isempty(app.selectedVertices)
                    % 没有任何点被选中了，清空所有 marker
                    if isfield(app.handles1, 'Multimarker') && ~isempty(app.handles1.Multimarker)
                        for i = 1:length(app.handles1.Multimarker)
                            if isvalid(app.handles1.Multimarker(i))
                                delete(app.handles1.Multimarker(i));
                            end
                        end
                        app.handles1.Multimarker = [];
                    end

                    if isfield(app.handles1, 'singlemarker')
                        delete(app.handles1.singlemarker);
                        app.handles1.singlemarker = [];
                    end
                else
                    % 还有其他点，重新高亮剩余的多选点
                    app.selectedVerticesNum = 0;
                    app.highlightMultipleVertices(app.selectedVertices);
                end
                return;  % Alt 模式处理完就返回，不走后面的 Ctrl/单选逻辑
            end

            if isCtrlPressed
                % Ctrl+点击：多选模式

                % 检查是否已经选中了该顶点
                if ~ismember(vertexIndex, app.selectedVertices)
                    app.selectedVertices = union(app.selectedVertices, vertexIndex);
                end


                % 高亮所有选中的顶点
                app.selectedVerticesNum = 0;
                app.highlightMultipleVertices(app.selectedVertices);
            else
                % 清除现有的多选标记
                if isfield(app.handles1, 'Multimarker') && ~isempty(app.handles1.Multimarker)
                    for i = 1:length(app.handles1.Multimarker)
                        if isvalid(app.handles1.Multimarker(i))
                            delete(app.handles1.Multimarker(i));
                        end
                    end
                    app.handles1.Multimarker = [];
                end
                % 清除旧标记
                if isfield(app.handles1, 'singlemarker')
                    delete(app.handles1.singlemarker);
                end

                app.selectedVertices=vertexIndex;
                app.addVisualMarker(vertexIndex);
                disp('vertice selected');
            end
            %  更新信息显示
            % updateSelectionInfo(app, vertexIndex, clickedPoint);
        end

        function [vertexIndex, minDist] = findClosestVertex(app, vertices, point)
            % 计算所有顶点到点击点的距离
            dists = sum((vertices - point).^2, 2);

            % 找到最近的顶点
            [minDist, vertexIndex] = min(dists);
        end

        function highlightMultipleVertices(app,vertexIndex)
            % 获取顶点位置
            vertices = get(app.handles1.hp, 'Vertices');
            % 清除现有的多选标记
            if isfield(app.handles1, 'Multimarker') && ~isempty(app.handles1.Multimarker)
                for i = 1:length(app.handles1.Multimarker)
                    if isvalid(app.handles1.Multimarker(i))
                        delete(app.handles1.Multimarker(i));
                    end
                end
                app.handles1.Multimarker = [];
            end

            if isfield(app.handles1, 'singlemarker') && ~isempty(app.handles1.singlemarker)
                if isvalid(app.handles1.singlemarker)
                    delete(app.handles1.singlemarker);
                end
                app.handles1.singlemarker = [];
            end

            % 创建球体标记
            [x, y, z] = sphere(20);
            % 4. 动态计算标记大小（基于大脑尺寸）
            brainSize = max(max(vertices) - min(vertices));
            r = brainSize * 0.01; % 大脑尺寸的2%

            hold(app.handles1.axes, 'on');
            for i=1:length(vertexIndex)
                vertexPos = vertices(vertexIndex(i), :);
                app.selectedVerticesNum = app.selectedVerticesNum + 1 ;
                marker(app.selectedVerticesNum) = surf(app.handles1.axes, ...
                    x*r + vertexPos(1), ...
                    y*r + vertexPos(2), ...
                    z*r + vertexPos(3), ...
                    'FaceColor', 'r', ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.9, ...
                    'Tag', 'SelectionMarker', ...
                    'HitTest', 'off', ...         
                    'PickableParts', 'none');       % 对 uiaxes 生效，彻底透明
            end
            hold(app.handles1.axes, 'off');
            app.handles1.Multimarker = marker;
        end

        function addVisualMarker(app, vertexIndex)
            % 获取顶点位置
            vertices = get(app.handles1.hp, 'Vertices');
            vertexPos = vertices(vertexIndex, :);

            % 创建球体标记
            [x, y, z] = sphere(20);
            % 4. 动态计算标记大小（基于大脑尺寸）
            brainSize = max(max(vertices) - min(vertices));
            r = brainSize * 0.01; % 大脑尺寸的2%

            hold(app.handles1.axes, 'on');

            marker = surf(app.handles1.axes, ...
                x*r + vertexPos(1), ...
                y*r + vertexPos(2), ...
                z*r + vertexPos(3), ...
                'FaceColor', 'r', ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.9, ...
                'Tag', 'SelectionMarker', ...
                'HitTest', 'off', ...
                'PickableParts', 'none');       % 对 uiaxes 生效，彻底透明

            hold(app.handles1.axes, 'off');
            app.handles1.singlemarker = marker;
        end

        function plotAtlasRegions(app, checkednodes)

    vc = get(app.handles1.hp,'Vertices');
    tri = app.cortexData.Faces;   % 用完整 faces 做过滤基准，避免受 Left/Right 下 hp 已被过滤的影响

    % 清除之前的区域面片
    if isfield(app.handles1, 'regionPatches') && ~isempty(app.handles1.regionPatches)
        for i = 1:length(app.handles1.regionPatches)
            if isvalid(app.handles1.regionPatches(i))
                delete(app.handles1.regionPatches(i));
            end
        end
        app.handles1.regionPatches = [];
    end

    % 获取当前坐标轴
    ax = app.handles1.axes;

    % 为每个分区创建独立面片
    scouts = app.cortexData.Atlas(app.AtlasIdx).Scouts;
    regionPatches = gobjects(length(scouts), 1);

    % 计算偏移距离（基于大脑尺寸的动态比例）
    brainSize = max(max(vc) - min(vc));
    offsetDistance = brainSize * 0.01;
    offsetVertices = vc;

    for i = 1:length(checkednodes)
        nodeIndex = checkednodes(i).NodeData;
        regionVertices = scouts(nodeIndex).Vertices;

        % 验证顶点索引
        if isempty(regionVertices)
            fprintf('Scout %d has no vertices, skipping.\n', nodeIndex);
            continue;
        end

        % 移除无效索引
        validIdx = regionVertices > 0 & regionVertices <= size(vc, 1);
        if sum(validIdx) ~= numel(regionVertices)
            fprintf('Scout %d has %d invalid vertices, skipping.\n', ...
                nodeIndex, numel(regionVertices)-sum(validIdx));
            regionVertices = regionVertices(validIdx);
        end

        if isempty(regionVertices)
            fprintf('No valid vertices for scout %d, skipping.\n', nodeIndex);
            continue;
        end

        % 获取包含这些顶点的面
        [FacesInRegion, localFaces] = ismember(tri, regionVertices);
        FacesInRegion = any(FacesInRegion, 2);

        regionFaces = tri(FacesInRegion, :);
        localFaces = localFaces(FacesInRegion, :);
        if isempty(regionFaces)
            fprintf('No faces found for scout %d, skipping.\n', nodeIndex);
            continue;
        end
        InnerIdx = all(localFaces, 2);
        EdgeVertices = localFaces(~InnerIdx, :);
        EdgeVertices = unique(EdgeVertices(:));
        EdgeVertices(EdgeVertices==0) = [];

        try
            scoutNormals = app.cortexData.Normals(regionVertices, :);
        catch
            scoutNormals = 0;
        end
        offsetVertices(regionVertices, :) = vc(regionVertices, :) + ...
            scoutNormals * offsetDistance;

        if strcmp(app.hemisphere,'structrue')
            app.sepregionVertices{nodeIndex} = offsetVertices;
        else
            app.OffsetVertices{nodeIndex} = offsetVertices;
        end

        % === Left/Right 模式下：按 Faces 过滤到目标半球，彻底弃用 alpha 过滤 ===
        if strcmp(app.hemisphere,'left') || strcmp(app.hemisphere,'right')
            keepF = all(ismember(regionFaces, app.plothemIdx), 2);
            regionFaces = regionFaces(keepF, :);
            if isempty(regionFaces)
                continue;
            end
        end

        alphaData = ones(size(vc,1), 1) * app.alphasignal;
        alphaData(EdgeVertices) = app.edgealphasignal;

        app.regioncolorset = repmat(scouts(nodeIndex).Color, size(vc,1), 1);

        % 创建分区面片
        set(ax, 'SortMethod', 'depth');
        regionPatches(nodeIndex) = patch(ax, 'Vertices', offsetVertices, 'Faces', regionFaces, ...
            'FaceVertexCData', app.regioncolorset, 'FaceColor', 'flat', ...
            'EdgeColor', 'none', 'alphaDataMapping', 'none', 'edgeAlpha', 'interp', ...
            'FaceAlpha', 'flat', 'FaceVertexAlphaData', alphaData, ...
            'Tag', sprintf('RegionPatch_%d', nodeIndex), 'facelighting', 'gouraud', ...
            'specularstrength', 0.2, 'ambientstrength', 0.5, 'diffusestrength', 0.5, ...
            'BackfaceLighting', 'lit', 'AmbientStrength', 0.5, 'SpecularExponent', 1, ...
            'SpecularColorReflectance', 0.5, 'EdgeLighting', 'gouraud');

        set(regionPatches(nodeIndex), ...
            'HitTest', 'on', 'PickableParts', 'all', ...
            'ButtonDownFcn', @(s,e) startDrag(app.handles1.h, ax, s, e));
        material dull
    end

    % 保存区域面片句柄
    app.handles1.regionPatches = regionPatches;
end

        function Valpha = getAlphaVelue(app,side)
            vc=app.cortexData.Vertices;
            tIdx = strcmp('Structures',string(struct2table(app.cortexData.Atlas).Name));
            if strcmp(side,'right')
                app.plothemIdx = app.cortexData.Atlas(tIdx).Scouts(2).Vertices;
            else
                app.plothemIdx = app.cortexData.Atlas(tIdx).Scouts(1).Vertices;
            end
            Valpha = zeros(size(vc,1),1);
            Valpha(app.plothemIdx) = 0.4;
        end

        function resectPatch(app, V, F, cuts, side)
            xcut=cuts(1); ycut=cuts(2); zcut=cuts(3);
            axesN = [1 0 0; 0 1 0; 0 0 1];
            cuts  = [xcut ycut zcut];

            % 取现有颜色数据（若为按顶点色）
            C = [];
            try
                C = get(app.handles1.hp,'FaceVertexCData');
                if size(C,1)~=size(V,1)   % 不是按顶点色，则不插值
                    C = [];
                end
            catch
            end

            V2 = V; F2 = F; C2 = C;

            for a = 1:3
                d = cuts(a); sgn = side(a);
                if ~isnan(d) && sgn~=0
                    [V2,F2,C2] = app.clipPlane_once(V2,F2,C2,axesN(a,:),d,sgn);
                end
            end

            set(app.handles1.hp,'Vertices',V2,'Faces',F2);           % 更新几何
            if ~isempty(C2)                                           % 若有按顶点色，回写
                set(app.handles1.hp,'FaceVertexCData',C2,'FaceColor','interp');
            end
        end

        % 对单一平面做一次裁剪
        function [V2,F2,C2] = clipPlane_once(app,V,F,C,n,d,keepSign)
            % 保留 keepSign*(n·x - d) >= 0 一侧；C 为空则不处理颜色
            n = n(:).'/norm(n);
            s = keepSign * (V*n.' - d);     % 顶点到平面的有符号距离（统一成“>=0 保留”）

            V2 = V; F2 = zeros(0,3);
            if isempty(C), C2 = []; else, C2 = C; end

            for k = 1:size(F,1)
                idx = F(k,:); vv = V(idx,:); ss = s(idx);
                pos = ss >= 0; m = sum(pos);

                if m==3                      % 全保留
                    F2(end+1,:) = idx;
                elseif m==0                  % 全丢弃
                    % skip
                else                         % 跨平面：补交点
                    % 重排：把保留侧排前面，便于统一处理
                    [~,ord] = sort(pos,'descend'); idx = idx(ord); vv = vv(ord,:); ss = ss(ord); pos = pos(ord);
                    if m==2
                        % 两保留一丢弃：四边形 -> 两三角
                        A=vv(1,:); B=vv(2,:); Cx=vv(3,:);
                        sA=ss(1);  sB=ss(2);  sC=ss(3);

                        tAC = sA/(sA - sC);  E = A + tAC*(Cx - A);
                        tBC = sB/(sB - sC);  G = B + tBC*(Cx - B);

                        eIdx = size(V2,1)+1; gIdx = eIdx+1;
                        V2 = [V2; E; G];

                        if ~isempty(C2)
                            CE = (1-tAC)*C2(idx(1),:) + tAC*C2(idx(3),:);
                            CG = (1-tBC)*C2(idx(2),:) + tBC*C2(idx(3),:);
                            C2 = [C2; CE; CG];
                        end

                        F2 = [F2; idx(1) idx(2) gIdx; idx(1) gIdx eIdx];

                    elseif m==1
                        % 一保留两丢弃：一个三角
                        A=vv(1,:); B=vv(2,:); Cx=vv(3,:);
                        sA=ss(1);  sB=ss(2);  sC=ss(3);

                        tAB = sA/(sA - sB);  E = A + tAB*(B - A);
                        tAC = sA/(sA - sC);  G = A + tAC*(Cx - A);

                        eIdx = size(V2,1)+1; gIdx = eIdx+1;
                        V2 = [V2; E; G];

                        if ~isempty(C2)
                            CE = (1-tAB)*C2(idx(1),:) + tAB*C2(idx(2),:);
                            CG = (1-tAC)*C2(idx(1),:) + tAC*C2(idx(3),:);
                            C2 = [C2; CE; CG];
                        end

                        F2 = [F2; idx(1) eIdx gIdx];
                    end
                end
            end
        end

        function closeFigureCallback(app, ~, ~)
            % 防止递归调用
            if app.isClosing, return; end

            % 标记关闭状态
            app.isClosing = true;
            app.stopPositionTimer();
            app.deleteSourceFigure();
            app.deleteControlFigure();
        end

        function closeAppCallback(app, ~, ~)
            % 防止递归调用
            if app.isClosing, return; end

            % 标记关闭状态
            app.isClosing = true;
            app.stopPositionTimer();
            app.deleteSourceFigure();
            app.deleteControlFigure();
        end    
        
        function checkedNodes = popregionnode(app,checkedNodes)
            indice = find(checkedNodes == app.RegionsNode);
            if ~isempty(indice)
                checkedNodes(indice) = [];
            end
        end

        function patchButtonDown(app, src, event)
            % 1）只在“选点模式”下才响应
            if ~app.SelectionMode   % 下面会加这个属性
                return
            end

            % 2）只响应左键（Button=1），中键/右键啥也不做
            if isprop(event, 'Button') && event.Button ~= 1
                return;
            end

            % 3）真正的选点逻辑
            app.vertexClicked(src, event);
        end

        function concatenated_data = concatenate_trials(app, data_3d)
            % 输入：data_3d - 三维数组（通道数 × 数据点 × 试次）
            % 输出：concatenated_data - 二维拼接矩阵（通道数 × 总数据点）
            % 说明：试次索引信息将存储在 app.TrialIndices 属性中

            % 获取数组尺寸
            [n_channels, n_points_per_trial, n_trials] = size(data_3d);

            % 预分配拼接后的二维矩阵
            total_points = n_points_per_trial * n_trials;
            concatenated_data = zeros(n_channels, total_points);

%             % 初始化索引记录结构体
%             app.TrialIndices = struct(...
%                 'trial_number', num2cell(1:n_trials), ...
%                 'start_index', num2cell(zeros(1, n_trials)), ...
%                 'end_index', num2cell(zeros(1, n_trials)) ...
%                 );


            % 循环处理每个试次
            for trial = 1:n_trials
                % 计算当前试次的起始和结束索引
                start_idx = (trial - 1) * n_points_per_trial + 1;
                end_idx = trial * n_points_per_trial;

                % 提取当前试次数据
                current_trial_data = data_3d(:, :, trial);

                % 将数据放入拼接矩阵
                concatenated_data(:, start_idx:end_idx) = current_trial_data;
% 
%                 % 记录索引信息
%                 app.TrialIndices(trial).start_index = start_idx;
%                 app.TrialIndices(trial).end_index = end_idx;
            end

%             % 添加元数据到索引结构
%             app.TrialIndices(1).n_channels = n_channels;
%             app.TrialIndices(1).points_per_trial = n_points_per_trial;
%             app.TrialIndices(1).total_trials = n_trials;
%             app.TrialIndices(1).total_points = total_points;
        end
    
        function checkPositionChange(app)
            % timer 回调：检测哪个窗口移动了
            if app.isSyncing
                return;
            end

            try
                if ~isvalid(app) || ~isvalid(app.handles1.h)
                    stop(app.posTimer);
                    return;
                end
            catch
                return;
            end

            figPos = app.handles1.h.Position;
            appPos = app.UIFigure.Position;

            figMoved = ~isempty(app.lastFigPos) && any(figPos(1:2) ~= app.lastFigPos(1:2));
            appMoved = ~isempty(app.lastAppPos) && any(appPos(1:2) ~= app.lastAppPos(1:2));

            if figMoved && ~appMoved
                % 图像窗口动了 → App 跟随
                app.isSyncing = true;
                newAppPos = figPos(1:2) - app.figOffset;
                app.UIFigure.Position(1:2) = newAppPos;
                app.isSyncing = false;

            elseif appMoved && ~figMoved
                % App 动了 → 图像窗口跟随
                app.isSyncing = true;
                newFigPos = appPos(1:2) + app.figOffset;
                app.handles1.h.Position(1:2) = newFigPos;
                app.isSyncing = false;
            end

            % 更新记录
            app.lastFigPos = app.handles1.h.Position;
            app.lastAppPos = app.UIFigure.Position;
        end
    
        function [cuts, sides] = computeCutsAndSides(app)
            % 把 3 个滑块值映射为切割面坐标 cuts 和保留方向 sides
            % 滑块=0 → 不切；>0 → 从正端切入；<0 → 从负端切入
            xrange = [min(app.smoothedVC(:,1)) max(app.smoothedVC(:,1))];
            yrange = [min(app.smoothedVC(:,2)) max(app.smoothedVC(:,2))];
            zrange = [min(app.smoothedVC(:,3)) max(app.smoothedVC(:,3))];
            ranges = {xrange, yrange, zrange};
            vals   = [app.SliderX.Value, app.SliderY.Value, app.SliderZ.Value];

            cuts  = zeros(1,3);
            sides = zeros(1,3);
            for a = 1:3
                r = ranges{a}; v = vals(a);
                if v > 0
                    cuts(a)  = r(2) - v;   % 从正端往里推 v
                    sides(a) = -1;         % 保留 x <= cuts(a)
                elseif v < 0
                    cuts(a)  = r(1) - v;   % 等价于 r(1)+|v|，从负端往里推 |v|
                    sides(a) = +1;         % 保留 x >= cuts(a)
                else
                    cuts(a)  = 0;
                    sides(a) = 0;          % 不切
                end
            end
        end

        function showOutOfTrial(app, selectTime)
            % 在标题里显示 out-of-trial 提示;脑图保持上次的
            if nargin >= 2 && ~isempty(selectTime) && isnumeric(selectTime)
                msg = sprintf('%s  — out of trial (cursor t=%.3fs)', ...
                    app.BaseTitle, selectTime);
            else
                msg = sprintf('%s  — out of trial', app.BaseTitle);
            end
            if ~strcmp(msg, app.OutOfTrialMsg)   % 避免反复刷
                app.UIFigure.Name = msg;
                app.OutOfTrialMsg = msg;
            end
        end

        function showPartialMsg(app)
            msg = sprintf('%s  — partial range', app.BaseTitle);
            if ~strcmp(msg, app.OutOfTrialMsg)
                app.UIFigure.Name = msg;
                app.OutOfTrialMsg = msg;
            end
        end

        function restoreTitleIfNeeded(app)
            if ~isempty(app.OutOfTrialMsg)
                app.UIFigure.Name = app.BaseTitle;
                app.OutOfTrialMsg = '';
            end
        end

    end

    methods (Access = public)

        function TimePointValueChanged(app, selectTime)
            % 非 trial 分片结果:走原逻辑,不做范围判断
            if ~app.IsTrialSplit
                colIdx = app.timeToColumn(selectTime);
                app.redrawSourceAtColumn(colIdx);
                if ~isempty(app.TimeVector) && numel(app.TimeVector) >= colIdx
                    app.TimePointEditField.Value = app.TimeVector(colIdx);
                end
                app.restoreTitleIfNeeded();
                return;
            end

            % trial 分片结果:判断是否在本 trial 范围内
            sGlobal = round(selectTime * app.srate) + 1;   % 全局采样点,1-based
            s1 = app.TrialRange(1);
            s2 = app.TrialRange(2);

            if sGlobal < s1 || sGlobal > s2
                % 不在范围内:脑图不动,只更新标题提示
                app.showOutOfTrial(selectTime);
                return;
            end

            % 在范围内:换算成本 trial 内的局部列索引
            localIdx = sGlobal - s1 + 1;
            if localIdx < 1, localIdx = 1; end
            if localIdx > size(app.resultData,2), localIdx = size(app.resultData,2); end

            app.redrawSourceData(app.resultData(:, localIdx));

            app.restoreTitleIfNeeded();
        end

        function datachanged(app, PoriodLeft, PoriodRight)
            % 非 trial 分片:走原逻辑
            if ~app.IsTrialSplit
                L = max(1, min(PoriodLeft,  size(app.resultData,2)));
                R = max(1, min(PoriodRight, size(app.resultData,2)));
                if L > R, [L,R] = deal(R,L); end

                subMatrix = app.resultData(:, L:R);
                app.PoriodData = sum(subMatrix.^2, 2);
                app.redrawSourceData(app.PoriodData, [], false);
                app.restoreTitleIfNeeded();
                return;
            end

            % trial 分片:和本 trial 范围求交集
            s1 = app.TrialRange(1);
            s2 = app.TrialRange(2);

            a = min(PoriodLeft, PoriodRight);
            b = max(PoriodLeft, PoriodRight);

            interA = max(a, s1);
            interB = min(b, s2);

            if interA > interB
                % 无交集:脑图不动,提示
                app.showOutOfTrial([]);
                return;
            end

            % 有交集:换算成本 trial 内的局部列索引
            L = interA - s1 + 1;
            R = interB - s1 + 1;
            L = max(1, min(L, size(app.resultData,2)));
            R = max(1, min(R, size(app.resultData,2)));

            subMatrix = app.resultData(:, L:R);
            app.PoriodData = sum(subMatrix.^2, 2);
            app.redrawSourceData(app.PoriodData, [], false);

            % 如果交集不等于全选范围,加个提示但不影响显示
            if interA ~= a || interB ~= b
                app.showPartialMsg();
            else
                app.restoreTitleIfNeeded();
            end
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            addpath(genpath('CoreFunctions'));
            app.UIFigure.Name='Show Brain Model';
            app.CallingApp=mainapp;
            
            uiNode = app.CallingApp.FileTree.SelectedNodes;
            if numel(uiNode) ~= 1
                app.close();
                return
            end

            node = uiNode.NodeData;
            screenSize = get(0, 'ScreenSize');
            
            if isa(node, "CortexNode")
                app.cortexData=node.data;
                app.handles1=PlotSource([],app.cortexData,'smooth',app.smoothvalue, ...
                    'Position',[screenSize(3)*0.3,screenSize(4)*0.4 800, 600]);
                app.jump2signalButton.Enable="off";
                app.SourceDataThreshold = 0.5;
                app.syncSourceThresholdControls();
                app.setSourceThresholdControlsEnabled(false);
                app.smoothedVC=seal_getSmoothedVertices(app.cortexData.Vertices, app.cortexData.VertConn, app.smoothvalue);
                xrange = [min(app.smoothedVC(:,1)) max(app.smoothedVC(:,1))];
                yrange = [min(app.smoothedVC(:,2)) max(app.smoothedVC(:,2))];
                zrange = [min(app.smoothedVC(:,3)) max(app.smoothedVC(:,3))];

                xspan = diff(xrange);
                yspan = diff(yrange);
                zspan = diff(zrange);

                app.SliderX.Limits = [-xspan xspan];
                app.SliderY.Limits = [-yspan yspan];
                app.SliderZ.Limits = [-zspan zspan];
                app.SliderX.Value = 0; app.SliderY.Value = 0; app.SliderZ.Value = 0;
                app.SliderX.Value=0;app.SliderY.Value=0;app.SliderZ.Value=0;
            else
               
                data = node;
                session = data.parent;
                protocol = session.parent;
                cortex = protocol.cortexNode;
                app.srate=double(data.params.srate);
                app.SliderX.Enable="off";app.SliderY.Enable="off";app.SliderZ.Enable="off";
                app.TimePointEditField.ValueChangedFcn = createCallbackFcn(app, @TimePointValueChanged, true);
                if ndims(node.sourceActivity)==3
                    app.resultData = app.concatenate_trials(node.sourceActivity);
                else
                    app.resultData = node.sourceActivity;
                end
                app.cortexData = cortex.data;
                app.smoothedVC=seal_getSmoothedVertices(app.cortexData.Vertices, app.cortexData.VertConn, app.smoothvalue);
                % === 读取 trial 分片元数据 ===
                app.BaseTitle = app.UIFigure.Name;

                rInfo = [];
                if isprop(node, 'resultInfo') && ~isempty(node.resultInfo)
                    rInfo = node.resultInfo;
                end

                % 关键:trial 元数据实际存放在 resultInfo.params 下
                params = [];
                if ~isempty(rInfo) && isprop(rInfo,'params') && ~isempty(rInfo.params)
                    params = rInfo.params;
                end

                if ~isempty(params) && isstruct(params) && ...
                        isfield(params,'IsTrialSplit') && params.IsTrialSplit

                    app.IsTrialSplit = true;
                    app.TrialRange   = params.TrialRange;
                    app.TrialIndex   = params.TrialIndex;
                    if isfield(params,'TrialCount'), app.TrialCount = params.TrialCount; end
                    if isfield(params,'BatchID'),    app.BatchID    = params.BatchID;    end

                    if ~isempty(app.TrialCount)
                        app.BaseTitle = sprintf('%s  [trial %d/%d]', app.BaseTitle, ...
                            app.TrialIndex, app.TrialCount);
                    else
                        app.BaseTitle = sprintf('%s  [trial %d]', app.BaseTitle, app.TrialIndex);
                    end
                    app.UIFigure.Name = app.BaseTitle;
                else
                    app.IsTrialSplit = false;
                end

                app.configureSourceDisplay(data);
                app.syncSourceThresholdControls();
                app.setSourceThresholdControlsEnabled(true);
                app.TimeVector = app.resolveResultTimeVector(data, params, size(app.resultData, 2));
                if ~isempty(app.TimeVector) && any(app.TimeVector <= 0) && any(app.TimeVector >= 0)
                    [~, initCol] = min(abs(app.TimeVector));
                else
                    initCol = 1;
                end
                app.TimePointEditField.Value = app.TimeVector(initCol);

                initRawSource = app.resultData(:, initCol);
                [initSource, initMap, initCLim, initThresh] = app.sourceDisplayArgs(initRawSource);
                app.handles1=PlotSource(initSource,app.cortexData,'smooth',app.smoothvalue, ...
                    'displaymode', app.SourceDisplayMode, ...
                    'colormap', initMap, ...
                    'clim', initCLim, ...
                    'thresh', initThresh, ...
                    'Position',[screenSize(3)*0.3,screenSize(4)*0.4 800, 600]);
                app.CurrentSourceData = initRawSource;
                app.CurrentSourceRGB = [];
                app.CurrentUseGlobalLimit = true;
                app.updateSourceColorbarLabel();
            end
            seal_bind3DInteraction(app.handles1.h, app.handles1.axes,  @(src,evt) app.vertexClicked(src,evt));
            app.handles1.h.UserData = 0;%0表示false
            app.UIFigure.Position=[screenSize(3)*0.3+800,screenSize(4)*0.4+274,480 326];
         
            disp(['Cortexdata Successfully loaded ', num2str(length(fieldnames(app.cortexData))), ...
                ' variables from the MainApp']);
            app.SmoothSlider.Value=30;app.smoothedit.Value=30;%初始化smooth
            app.TranspSlider.Value=1;app.Transpedit.Value=1;%初始化transp
            app.transp = 1;
            app.plothemIdx = 1:size(app.cortexData.Vertices,1);
         
            app.Label.Text=sprintf('Vertices: %d', size(app.smoothedVC,1));
            app.Label_2.Text=sprintf('Faces: %d',size(app.cortexData.Faces,1));

            if ~isfield(app.handles1, 'hp') || ~isvalid(app.handles1.hp)
                error('Brain patch object not found or invalid');
            end


            %atlas
            if isfield(app.cortexData, 'Atlas')
                num_als = numel(app.cortexData.Atlas);
                if num_als > 0
                    scouts = app.cortexData.Atlas(1).Scouts;

                    for i = 1:length(scouts)
                        if isfield(scouts(i),'Region')
                            regionNames = scouts(i).Region;
                            node = uitreenode(app.RegionsNode);
                            node.Text=regionNames;
                            node.NodeData = i;
                        end
                    end
                    app.RegionTree.CheckedNodes = [];
                    expand(app.RegionsNode);
                    als_list = {}; als_data = [];
                    for ials = 1:num_als
                        if isfield(app.cortexData.Atlas(ials),'Name')
                            als_list{ials} = app.cortexData.Atlas(ials).Name;
                            als_data(ials) = ials;
                        end
                    end
                    app.AtlasDrop.Items = als_list;
                    app.AtlasDrop.ItemsData = als_data;

                end


            end

            app.figOffset = app.handles1.h.Position(1:2) - app.UIFigure.Position(1:2);
            app.lastFigPos = app.handles1.h.Position;
            app.lastAppPos = app.UIFigure.Position;
            % figure 没有原生的移动事件，用 timer 轮询
            app.posTimer = timer(...
                'ExecutionMode', 'fixedRate', ...
                'Period', 0.05, ...
                'TimerFcn', @(src, evt) seal_showCortexCheckPosition(src, evt));
            app.posTimer.UserData = struct( ...
                'sourceFigure', app.handles1.h, ...
                'controlFigure', app.UIFigure, ...
                'figOffset', app.figOffset, ...
                'lastFigPos', app.lastFigPos, ...
                'lastAppPos', app.lastAppPos, ...
                'isSyncing', false);
            setappdata(app.handles1.h, 'SEALPeerFigure', app.UIFigure);
            setappdata(app.UIFigure, 'SEALPeerFigure', app.handles1.h);
            setappdata(app.handles1.h, 'SEALPositionTimer', app.posTimer);
            setappdata(app.UIFigure, 'SEALPositionTimer', app.posTimer);

            % 设置双向关闭回调
            % 当任一窗口关闭时，关闭另一个窗口并清理 timer
            set(app.handles1.h, 'CloseRequestFcn', @(src, evt) seal_showCortexClose(src, evt));
            set(app.UIFigure, 'CloseRequestFcn', @(src, evt) seal_showCortexClose(src, evt));
            start(app.posTimer);

            app.UIFigure.Tag = 'SEAL_ResultWindow';
            
            app.UIFigure.UserData = app;
        end

        % Value changed function: SmoothSlider
        function SmoothSliderValueChanged(app, event)
            %传递新smoothvalue
            app.smoothvalue = app.SmoothSlider.Value;
            app.smoothedit.Value=app.SmoothSlider.Value;
            %计算新的顶点坐标
            app.smoothedVC=seal_getSmoothedVertices(app.cortexData.Vertices, app.cortexData.VertConn, app.smoothvalue);
            %判断顶点坐标计算是否正确
            if any(isnan(app.smoothedVC(:))) || any(isinf(app.smoothedVC(:)))
                error('The smoothed vertices contain NaN or Inf and cannot be rendered');
            end
            % 更新顶点，不新建窗口
            set(app.handles1.hp, 'Vertices', app.smoothedVC);


        end

        % Value changed function: smoothedit
        function smootheditValueChanged(app, event)
            %传递新smoothvalue
            app.smoothvalue=app.smoothedit.Value;
            app.SmoothSlider.Value=app.smoothedit.Value;
            %计算新的顶点坐标
            app.smoothedVC=seal_getSmoothedVertices(app.cortexData.Vertices, app.cortexData.VertConn, app.smoothvalue);
            %判断顶点坐标计算是否正确
            if any(isnan(app.smoothedVC(:))) || any(isinf(app.smoothedVC(:)))
                error('The smoothed vertices contain NaN or Inf and cannot be rendered');
            end
            % 更新顶点，不新建窗口
            set(app.handles1.hp, 'Vertices', app.smoothedVC);

        end

        % Value changed function: ThresholdSlider
        function ThresholdSliderValueChanged(app, event)
            app.setSourceThreshold(app.ThresholdSlider.Value, 'slider');
        end

        % Value changed function: ThresholdEdit
        function ThresholdEditValueChanged(app, event)
            app.setSourceThreshold(app.ThresholdEdit.Value, 'edit');
        end

        % Value changed function: TranspSlider
        function TranspSliderValueChanged(app, event)
            app.transp = app.TranspSlider.Value;
            app.Transpedit.Value = app.TranspSlider.Value;
            Valpha = zeros(size(app.smoothedVC, 1), 1);
            Valpha(app.plothemIdx) = app.transp;
            set(app.handles1.hp, 'FaceVertexAlphaData', Valpha, ...
                'facealpha', 'interp' ,'AlphaDataMapping','none');
            %             set(app.handles1.hp, 'facealpha', app.transp);
        end

        % Value changed function: Transpedit
        function TranspeditValueChanged(app, event)
            app.transp = app.Transpedit.Value;
            app.TranspSlider.Value   =  app.Transpedit.Value;
            Valpha = zeros(size(app.smoothedVC, 1), 1);
            Valpha(app.plothemIdx) = app.transp;
            set(app.handles1.hp, 'FaceVertexAlphaData', Valpha, ...
                'facealpha', 'interp' ,'AlphaDataMapping','none');

        end

        % Button pushed function: RightButton
        function RightButtonPushed(app, event)

            app.hemisphere='right';
            app.RightorLeft(app.hemisphere);
        end

        % Button pushed function: LeftButton
        function LeftButtonPushed(app, event)

            app.hemisphere='left';
            app.RightorLeft(app.hemisphere);
        end

        % Button pushed function: StructureButton
        function StructureButtonPushed(app, event)

            app.hemisphere='structrue';

            app.RightorLeft(app.hemisphere);

            app.alreadystructure =1;
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)

            app.hemisphere='reset';
            app.RightorLeft(app.hemisphere);
            app.SliderX.Value=0; app.SliderY.Value=0; app.SliderZ.Value=0;
            app.resectPatch(app.smoothedVC,app.cortexData.Faces, [0,0,0], [0,0,0]) ;
        end

        % Callback function: RegionTree
        function RegionTreeCheckedNodesChanged(app, event)
            checkedNodes = app.RegionTree.CheckedNodes;
            checkedNodes=app.popregionnode(checkedNodes);
            app.plotAtlasRegions(checkedNodes);
        end

        % Value changed function: AtlasDrop
        function AtlasDropValueChanged(app, event)
            value = app.AtlasDrop.Value;
            scouts = app.cortexData.Atlas(value).Scouts;
            app.RegionsNode.Children.delete;
            for i = 1:length(scouts)
                regionNames = scouts(i).Region;
                node = uitreenode(app.RegionsNode);
                node.Text=regionNames;
                node.NodeData = i;
            end
            expand(app.RegionsNode);
            app.AtlasIdx = value;
            RegionTreeCheckedNodesChanged(app);

        end

        % Button pushed function: ColorButton
        function ColorButtonPushed(app, event)
            import javax.swing.JColorChooser;
            import java.awt.Color;

            defaultColor = java.awt.Color(0.3, 0.4, 0.8);  % 默认颜色
            chosenColor = JColorChooser.showDialog([], '选择颜色', defaultColor);

            if ~isempty(chosenColor)
                rgb = [chosenColor.getRed, chosenColor.getGreen, chosenColor.getBlue] / 255;
                colIdx = app.timeToColumn(app.TimePointEditField.Value);
                app.redrawSourceAtColumn(colIdx, rgb);

            end

        end

        % Value changed function: EdgeButton
        function EdgeButtonValueChanged(app, event)
            value = app.EdgeButton.Value;
            if value==1
                set(app.handles1.hp,'edgecolor',[0.4,0.4,0.4]);
            else
                set(app.handles1.hp,'edgecolor','none');
            end
        end

        % Value changed function: SliderX
        function SliderXValueChanged(app, event)
            [cuts, sides] = app.computeCutsAndSides();
            app.resectPatch(app.smoothedVC, app.cortexData.Faces, cuts, sides);
        end

        % Value changed function: SliderY
        function SliderYValueChanged(app, event)
            [cuts, sides] = app.computeCutsAndSides();
            app.resectPatch(app.smoothedVC, app.cortexData.Faces, cuts, sides);
        end

        % Value changed function: SliderZ
        function SliderZValueChanged(app, event)
            [cuts, sides] = app.computeCutsAndSides();
            app.resectPatch(app.smoothedVC, app.cortexData.Faces, cuts, sides);
        end

        % Selection changed function: ButtonGroup
        function ButtonGroupSelectionChanged(app, event)
            selectedButton = app.ButtonGroup.SelectedObject;
            % 清除之前的区域面片
            if isfield(app.handles1, 'regionPatches') && ~isempty(app.handles1.regionPatches)
                if isequal(app.NoColorButton,selectedButton)
                    app.alphasignal=0;
                elseif isequal(app.FullColorButten,selectedButton)
                    app.alphasignal=1;
                elseif isequal(app.HalfColorButton,selectedButton)
                    app.alphasignal=0.4;
                end
            end
            checkedNodes = app.RegionTree.CheckedNodes;
            checkedNodes=app.popregionnode(checkedNodes);
            app.plotAtlasRegions(checkedNodes);
        end

        % Selection changed function: ButtonGroup_2
        function ButtonGroup_2SelectionChanged(app, event)
            selectedButton = app.ButtonGroup_2.SelectedObject;
            if isfield(app.handles1, 'regionPatches') && ~isempty(app.handles1.regionPatches)
                if isequal(app.ShowFrameButten,selectedButton)
                    app.edgealphasignal=1;
                elseif isequal(app.HideFrameButton,selectedButton)
                    app.edgealphasignal=0;
                end
            end
            checkedNodes = app.RegionTree.CheckedNodes;
            checkedNodes=app.popregionnode(checkedNodes);
            app.plotAtlasRegions(checkedNodes);
        end

        % Value changed function: SelectpointButton
        function SelectpointButtonValueChanged(app, event)
            value = app.SelectpointButton.Value;
            if value==1
               
                app.handles1.h.UserData = 1;
                set(app.handles1.h, 'Pointer', 'crosshair');
            else
                app.handles1.h.UserData = 0;%0表示false，选点模式
                
                set(app.handles1.h, 'Pointer', 'arrow');
            end
        end

        % Button pushed function: jump2signalButton
        function jump2signalButtonPushed(app, event)
             if isempty(app.selectedVertices)
                return
            end
            selcetdata = app.resultData(app.selectedVertices, :);
            % Build a minimal struct compatible with SEAL_showEEG startupFcn
            mockNode.data = selcetdata;
            mockNode.srate = app.srate;
            mockNode.name = 'Source Signal';
            mockNode.parent.parent.path = '';
            mockNode.dataInfo.chanlocs = [];
            mockNode.dataInfo.event = [];
            SEAL_showEEG(mockNode, app.UIFigure.Position);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 487 324];
            app.UIFigure.Name = 'MATLAB App';

            % Create SurfaceoprionsPanel
            app.SurfaceoprionsPanel = uipanel(app.UIFigure);
            app.SurfaceoprionsPanel.Title = 'Surface oprions';
            app.SurfaceoprionsPanel.Position = [155 157 331 167];

            % Create smoothedit
            app.smoothedit = uieditfield(app.SurfaceoprionsPanel, 'numeric');
            app.smoothedit.ValueChangedFcn = createCallbackFcn(app, @smootheditValueChanged, true);
            app.smoothedit.Position = [259 121 44 22];

            % Create Transpedit
            app.Transpedit = uieditfield(app.SurfaceoprionsPanel, 'numeric');
            app.Transpedit.ValueChangedFcn = createCallbackFcn(app, @TranspeditValueChanged, true);
            app.Transpedit.Position = [260 53 44 22];

            % Create ThresholdEdit
            app.ThresholdEdit = uieditfield(app.SurfaceoprionsPanel, 'numeric');
            app.ThresholdEdit.Limits = [0 1];
            app.ThresholdEdit.ValueChangedFcn = createCallbackFcn(app, @ThresholdEditValueChanged, true);
            app.ThresholdEdit.Position = [260 87 44 22];
            app.ThresholdEdit.Value = 0.5;

            % Create SmoothLabel
            app.SmoothLabel = uilabel(app.SurfaceoprionsPanel);
            app.SmoothLabel.HorizontalAlignment = 'right';
            app.SmoothLabel.Position = [20 111 66 22];
            app.SmoothLabel.Text = 'Smooth(%)';

            % Create SmoothSlider
            app.SmoothSlider = uislider(app.SurfaceoprionsPanel);
            app.SmoothSlider.Limits = [1 100];
            app.SmoothSlider.MajorTicks = [1 20 40 60 80 100];
            app.SmoothSlider.ValueChangedFcn = createCallbackFcn(app, @SmoothSliderValueChanged, true);
            app.SmoothSlider.Position = [99 130 150 3];
            app.SmoothSlider.Value = 1;

            % Create ThresholdLabel
            app.ThresholdLabel = uilabel(app.SurfaceoprionsPanel);
            app.ThresholdLabel.HorizontalAlignment = 'right';
            app.ThresholdLabel.Position = [13 87 73 22];
            app.ThresholdLabel.Text = 'Threshold';

            % Create ThresholdSlider
            app.ThresholdSlider = uislider(app.SurfaceoprionsPanel);
            app.ThresholdSlider.Limits = [0 1];
            app.ThresholdSlider.MajorTicks = [0 0.25 0.5 0.75 1];
            app.ThresholdSlider.ValueChangedFcn = createCallbackFcn(app, @ThresholdSliderValueChanged, true);
            app.ThresholdSlider.Position = [100 97 150 3];
            app.ThresholdSlider.Value = 0.5;

            % Create TranspLabel
            app.TranspLabel = uilabel(app.SurfaceoprionsPanel);
            app.TranspLabel.HorizontalAlignment = 'right';
            app.TranspLabel.Position = [25 53 61 22];
            app.TranspLabel.Text = 'Transp(%)';

            % Create TranspSlider
            app.TranspSlider = uislider(app.SurfaceoprionsPanel);
            app.TranspSlider.Limits = [0 1];
            app.TranspSlider.ValueChangedFcn = createCallbackFcn(app, @TranspSliderValueChanged, true);
            app.TranspSlider.Position = [100 63 150 3];

            % Create ColorButton
            app.ColorButton = uibutton(app.SurfaceoprionsPanel, 'push');
            app.ColorButton.ButtonPushedFcn = createCallbackFcn(app, @ColorButtonPushed, true);
            app.ColorButton.Position = [31 10 100 23];
            app.ColorButton.Text = 'Color';

            % Create EdgeButton
            app.EdgeButton = uibutton(app.SurfaceoprionsPanel, 'state');
            app.EdgeButton.ValueChangedFcn = createCallbackFcn(app, @EdgeButtonValueChanged, true);
            app.EdgeButton.Text = 'Edge';
            app.EdgeButton.Position = [183 9 100 23];

            % Create RegionsPanel
            app.RegionsPanel = uipanel(app.UIFigure);
            app.RegionsPanel.Title = 'Regions';
            app.RegionsPanel.Position = [1 -1 154 325];

            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.RegionsPanel);
            app.ButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.ButtonGroup.Position = [0 203 155 102];

            % Create FullColorButten
            app.FullColorButten = uitogglebutton(app.ButtonGroup);
            app.FullColorButten.Icon = 'circle3.png';
            app.FullColorButten.Text = '';
            app.FullColorButten.Position = [7 69 25 25];

            % Create HalfColorButton
            app.HalfColorButton = uitogglebutton(app.ButtonGroup);
            app.HalfColorButton.Icon = 'circle2.png';
            app.HalfColorButton.Text = '';
            app.HalfColorButton.Position = [7 38 25 25];
            app.HalfColorButton.Value = true;

            % Create NoColorButton
            app.NoColorButton = uitogglebutton(app.ButtonGroup);
            app.NoColorButton.Icon = 'circle1.png';
            app.NoColorButton.IconAlignment = 'right';
            app.NoColorButton.Text = '';
            app.NoColorButton.Position = [7 7 25 26];

            % Create AtlasDrop
            app.AtlasDrop = uidropdown(app.ButtonGroup);
            app.AtlasDrop.Items = {'Atlas'};
            app.AtlasDrop.ValueChangedFcn = createCallbackFcn(app, @AtlasDropValueChanged, true);
            app.AtlasDrop.Position = [38 80 116 21];
            app.AtlasDrop.Value = 'Atlas';

            % Create ButtonGroup_2
            app.ButtonGroup_2 = uibuttongroup(app.RegionsPanel);
            app.ButtonGroup_2.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_2SelectionChanged, true);
            app.ButtonGroup_2.Position = [0 128 155 76];

            % Create ShowFrameButten
            app.ShowFrameButten = uitogglebutton(app.ButtonGroup_2);
            app.ShowFrameButten.Icon = 'alphaframe.png';
            app.ShowFrameButten.Text = '';
            app.ShowFrameButten.Position = [7 9 25 25];
            app.ShowFrameButten.Value = true;

            % Create HideFrameButton
            app.HideFrameButton = uitogglebutton(app.ButtonGroup_2);
            app.HideFrameButton.Icon = 'solidframe.png';
            app.HideFrameButton.Text = '';
            app.HideFrameButton.Position = [7 41 25 25];

            % Create ButtonGroup_3
            app.ButtonGroup_3 = uibuttongroup(app.RegionsPanel);
            app.ButtonGroup_3.Position = [0 64 155 65];

            % Create SelectpointButton
            app.SelectpointButton = uibutton(app.ButtonGroup_3, 'state');
            app.SelectpointButton.ValueChangedFcn = createCallbackFcn(app, @SelectpointButtonValueChanged, true);
            app.SelectpointButton.Icon = 'point.png';
            app.SelectpointButton.Text = '';
            app.SelectpointButton.Position = [7 35 25 23];

            % Create jump2signalButton
            app.jump2signalButton = uibutton(app.ButtonGroup_3, 'push');
            app.jump2signalButton.ButtonPushedFcn = createCallbackFcn(app, @jump2signalButtonPushed, true);
            app.jump2signalButton.Icon = 'S.png';
            app.jump2signalButton.Position = [7 7 25 23];
            app.jump2signalButton.Text = '';

            % Create ResectXYZPanel
            app.ResectXYZPanel = uipanel(app.UIFigure);
            app.ResectXYZPanel.Title = 'Resect [X,Y,Z]';
            app.ResectXYZPanel.Position = [155 30 331 128];

            % Create LeftButton
            app.LeftButton = uibutton(app.ResectXYZPanel, 'push');
            app.LeftButton.ButtonPushedFcn = createCallbackFcn(app, @LeftButtonPushed, true);
            app.LeftButton.Position = [14 19 49 22];
            app.LeftButton.Text = 'Left';

            % Create RightButton
            app.RightButton = uibutton(app.ResectXYZPanel, 'push');
            app.RightButton.ButtonPushedFcn = createCallbackFcn(app, @RightButtonPushed, true);
            app.RightButton.Position = [80 19 48 22];
            app.RightButton.Text = 'Right';

            % Create StructureButton
            app.StructureButton = uibutton(app.ResectXYZPanel, 'push');
            app.StructureButton.ButtonPushedFcn = createCallbackFcn(app, @StructureButtonPushed, true);
            app.StructureButton.Position = [144 19 70 22];
            app.StructureButton.Text = 'Structure';

            % Create ResetButton
            app.ResetButton = uibutton(app.ResectXYZPanel, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [237 19 75 22];
            app.ResetButton.Text = 'Reset';

            % Create SliderX
            app.SliderX = uislider(app.ResectXYZPanel);
            app.SliderX.MajorTicks = [];
            app.SliderX.ValueChangedFcn = createCallbackFcn(app, @SliderXValueChanged, true);
            app.SliderX.MinorTicks = [];
            app.SliderX.Position = [15 73 83 3];

            % Create SliderY
            app.SliderY = uislider(app.ResectXYZPanel);
            app.SliderY.MajorTicks = [];
            app.SliderY.ValueChangedFcn = createCallbackFcn(app, @SliderYValueChanged, true);
            app.SliderY.MinorTicks = [];
            app.SliderY.Position = [124 73 83 3];

            % Create SliderZ
            app.SliderZ = uislider(app.ResectXYZPanel);
            app.SliderZ.MajorTicks = [];
            app.SliderZ.ValueChangedFcn = createCallbackFcn(app, @SliderZValueChanged, true);
            app.SliderZ.MinorTicks = [];
            app.SliderZ.Position = [230 73 83 3];

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.Position = [186 6 137 22];

            % Create Label_2
            app.Label_2 = uilabel(app.UIFigure);
            app.Label_2.Position = [335 6 124 22];

            % Create RegionTree
            app.RegionTree = uitree(app.UIFigure, 'checkbox');
            app.RegionTree.Position = [40 -1 116 285];

            % Create RegionsNode
            app.RegionsNode = uitreenode(app.RegionTree);
            app.RegionsNode.Text = 'Regions';

            % Assign Checked Nodes
            app.RegionTree.CheckedNodes = [app.RegionsNode];
            % Assign Checked Nodes
            app.RegionTree.CheckedNodesChangedFcn = createCallbackFcn(app, @RegionTreeCheckedNodesChanged, true);

            % Create EditFieldLabel
            app.EditFieldLabel = uilabel(app.UIFigure);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Position = [142 -37 55 22];
            app.EditFieldLabel.Text = 'Edit Field';

            % Create TimePointEditField
            app.TimePointEditField = uieditfield(app.UIFigure, 'numeric');
            app.TimePointEditField.Position = [212 -37 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SEAL_showCortex(varargin)

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

            app.stopPositionTimer();
            app.deleteSourceFigure();
            app.deleteControlFigure();
        end
    end
end

function seal_showCortexClose(src, ~)
try
    if isempty(src) || ~isvalid(src)
        return;
    end

    peer = [];
    posTimer = [];
    try
        peer = getappdata(src, 'SEALPeerFigure');
    catch
    end
    try
        posTimer = getappdata(src, 'SEALPositionTimer');
    catch
    end

    if ~isempty(posTimer) && isvalid(posTimer)
        stop(posTimer);
        delete(posTimer);
    end

    if isvalid(src)
        set(src, 'CloseRequestFcn', '');
        delete(src);
    end
    if ~isempty(peer) && isvalid(peer)
        set(peer, 'CloseRequestFcn', '');
        delete(peer);
    end
catch
end
end

function seal_showCortexCheckPosition(timerObj, ~)
try
    if isempty(timerObj) || ~isvalid(timerObj)
        return;
    end
    state = timerObj.UserData;
    if ~isstruct(state) || ~isfield(state, 'sourceFigure') || ~isfield(state, 'controlFigure')
        return;
    end
    if ~isfield(state, 'isSyncing')
        state.isSyncing = false;
    end
    if state.isSyncing
        return;
    end

    hFig = state.sourceFigure;
    hApp = state.controlFigure;
    if isempty(hFig) || isempty(hApp) || ~isvalid(hFig) || ~isvalid(hApp)
        stop(timerObj);
        return;
    end

    figPos = hFig.Position;
    appPos = hApp.Position;
    figMoved = ~isempty(state.lastFigPos) && any(figPos(1:2) ~= state.lastFigPos(1:2));
    appMoved = ~isempty(state.lastAppPos) && any(appPos(1:2) ~= state.lastAppPos(1:2));

    if figMoved && ~appMoved
        state.isSyncing = true;
        hApp.Position(1:2) = figPos(1:2) - state.figOffset;
        state.isSyncing = false;
    elseif appMoved && ~figMoved
        state.isSyncing = true;
        hFig.Position(1:2) = appPos(1:2) + state.figOffset;
        state.isSyncing = false;
    end

    state.lastFigPos = hFig.Position;
    state.lastAppPos = hApp.Position;
    timerObj.UserData = state;
catch
end
end
