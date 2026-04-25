classdef SEAL_showEEG < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        FiltersMenu                     matlab.ui.container.Menu
        FFTMenu                         matlab.ui.container.Menu
        LowPassMenu                     matlab.ui.container.Menu
        HighPassMenu                    matlab.ui.container.Menu
        BandPassMenu                    matlab.ui.container.Menu
        BandSuppressionMenu             matlab.ui.container.Menu
        FiltersMenu_2                   matlab.ui.container.Menu

        % --- 顶层网格 ---
        MainGrid                        matlab.ui.container.GridLayout

        % --- 主要面板 ---
        ChennalsPanel                   matlab.ui.container.Panel
        ChannelsGrid                    matlab.ui.container.GridLayout
        EventListPanel                  matlab.ui.container.Panel
        EventListGrid                   matlab.ui.container.GridLayout
        MarkneweventPanel               matlab.ui.container.Panel
        MarkEventGrid                   matlab.ui.container.GridLayout
        ControlPanel                    matlab.ui.container.Panel
        ControlGrid                     matlab.ui.container.GridLayout

        % --- 控制面板内的子网格 ---
        SelectionGrid                   matlab.ui.container.GridLayout
        WindowGrid                      matlab.ui.container.GridLayout
        NavGrid                         matlab.ui.container.GridLayout

        % --- 中央绘图区域容器 ---
        CenterPanel                     matlab.ui.container.Panel
        CenterGrid                      matlab.ui.container.GridLayout

        % --- 控件 ---
        RangeEditField                  matlab.ui.control.EditField
        RangeEditFieldLabel             matlab.ui.control.Label
        ScaleEditField                  matlab.ui.control.NumericEditField
        ScaleEditFieldLabel             matlab.ui.control.Label
        WinlenthEditField               matlab.ui.control.NumericEditField
        WinlenthEditFieldLabel          matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        rangeend                        matlab.ui.control.Label
        rangestart                      matlab.ui.control.Label
        biggerButton                    matlab.ui.control.Button
        smallerButton                   matlab.ui.control.Button
        LLButton                        matlab.ui.control.Button
        RRButton                        matlab.ui.control.Button
        RButton                         matlab.ui.control.Button
        LButton                         matlab.ui.control.Button
        SelectedTimeEditField           matlab.ui.control.NumericEditField
        SelectedTimeEditField_2Label_2  matlab.ui.control.Label
        StackingmodeButton              matlab.ui.control.StateButton
        EventLabelEditField             matlab.ui.control.EditField
        EventNameLabel                  matlab.ui.control.Label
        AddEventButton                  matlab.ui.control.Button
        EventListListBox                matlab.ui.control.ListBox
        SelectAllButton                 matlab.ui.control.Button
        ClearButton                     matlab.ui.control.Button
        ChannelTree                     matlab.ui.container.CheckBoxTree
        UIAxes                          matlab.ui.control.UIAxes
    end


    properties (Access = private)
        eegdata         % EEG 数据
        srate = 256     % 采样率
        spacing = 20    % 通道间距
        winlength = 1   % 窗口长度（秒）
        time =   0      % 当前时间(区间左)
        chans = 0       % 通道数
        frames = 0      % 时间点数，后续初始化
        dispchans = 0   % 显示通道数
        elecoffset = 0  % 垂直偏移
        DialogApp       %调用滤波器app
        submean= 'off'
        eloc_file       %通道文件
        YLabels         %生成通道序号
        timeStart       %区间左
        timeEnd         %区间右
        dataNode
        ChannelNode={}
        chanlocs
        CheckedData=[]          %选择显示通道的索引
        eventStrings            %标记列表元素名称   
        duidie=1                %堆叠标签
        SaveIndex = 1;          % 文件保存计数器
        SaveIndexrowdata = 1;   % 文件保存计数器
        SavePath = 'Data';      % 默认保存路径
        CallingAppResult             %联动结果显示
        TrialIndices
        datapoints
        SelectedChannel = 0;       % 当前选中的通道索引
        EventMarkers = struct('latency', {}, 'type', {}); % 事件标记      
        EventData = struct('latency', {}, 'type', {}); % 自定义事件标记     
        LineHandles = [];          % 通道线条句柄
        TextHandles = [];          % 通道标签句柄
        EventLineHandles = [];     % 事件线句柄
        EventTextHandles = [];     % 事件文本句柄
        TimeSelectionLine;         % 时间选择线
        ChannelSelectionText;      % 通道选择文本
        OriginalYLim2              %大y
        segmentInfo                 %信号分割信息
        StartChan (1,1) double = 1        % 当前起始通道（从1计）
        ChannelSlider matlab.ui.control.Slider  % 垂直滚动条引用

        % === 可调的幅值缩放参数 ===
        ScaleMode char = 'global';   % 选项: 'none'|'std'|'robust'|'minmax'|'zscore'|'percentile'
        ScaleFactor =0                    %缩放因子
        ScalePercentiles double = [2 98]; % 百分位缩放用到的上下界
        ScaleTarget double = [];          % 目标半幅值(纵向)，留空则按 spacing 自动确定
        SelDebounceTimer  % 一次性定时器，用于SelectionChanged去抖
        ClickHoldTimer            % 单次定时器：判断是否进入框选
        MouseDownTic uint64 = uint64(0);%记录按下时刻（tic）
        MouseDownPoint double = []% 记录按下时坐标 [x y]

        AutoScaleFactor = []   % 基于整段数据一次性算好的自动缩放系数
    end

    properties (Access = public)
        currentSaveIndex = 1; % 当前结构体保存索引
        historydatamanu
        filteredDataHistory = struct('data', {}, 'params', {}, 'index', {}); % 历史记录结构体数组
        SelectedTime = 0;               % 当前选中的时间点（秒）
        IsSelecting = false;            % 是否正在选择时间段
        SelectionStart = [];            % 选择开始的时间点
        SelectionEnd = [];              % 选择结束的时间点
        SelectionRectangle = [];        % 选择矩形的图形句柄
        SelectedTimeRange = [];         % 存储选择的时间范围
    end

    methods (Access = private)

        function redrawPlot(app)
            % 绑定到 UIAxes，清空现有内容
            cla(app.UIAxes, 'reset');
            hold(app.UIAxes, 'on');
            app.LineHandles = []; % 清空线条句柄
            app.TextHandles = []; % 清空文本句柄
            % 计算当前窗口数据范围（全局帧s→窗口内相对位置）
            [lowlim, highlim] = app.calcWindowRange();

            % 绘制 EEG 通道数据（上层，ZData=1）
            app.plotEEGChannels(lowlim, highlim);
            % 绘制自定义实践
            if ~isempty(app.EventData)
                for i = 1:length(app.EventData)
                    eventTime = app.EventData(i).latency;

                    if isempty(eventTime) || ~isnumeric(eventTime) || ~isscalar(eventTime) || ~isfinite(eventTime)
                        continue;
                    end

                    evType = app.EventData(i).type;
                    if isempty(evType)
                        evType = '';
                    elseif ~ischar(evType)
                        evType = char(string(evType));
                    end

                    if eventTime >= app.timeStart && eventTime <= app.timeEnd
                        lh = line(app.UIAxes, [eventTime, eventTime], app.UIAxes.YLim, ...
                            'Color', 'g', 'LineWidth', 1.5, 'LineStyle', '--');
                        app.EventLineHandles(end+1) = lh;

                        th2 = text(app.UIAxes, eventTime, app.UIAxes.YLim(2), evType, ...
                            'Color', 'g', 'FontSize', 8, 'VerticalAlignment', 'bottom','Clipping','off');
                        app.EventTextHandles(end+1) = th2;
                    end
                end
            end
        
        end

        function [lowlim, highlim] = calcWindowRange(app)   %计算点数的上下限
            lowlim = round(app.time*app.srate  + 1);
            highlim = round(min((app.time + app.winlength)*app.srate  + 1, app.frames));
        end

        function initializeChannelLabels(app)
            % 1. 确立绝对基准：总通道数永远等于"数据矩阵的行数"
            app.chans = size(app.eegdata, 1);

            % 2. 检查数据对象中是否带了 chanlocs 信息，并且包含 'labels' 字段
            if isempty(app.chanlocs) || ~isfield(app.chanlocs, 'labels')
                app.YLabels = arrayfun(@(x) sprintf('CH %d', x), 1:app.chans, 'UniformOutput', false);
            else
                extractedLabels = {app.chanlocs.labels};

                if length(extractedLabels) == app.chans
                    app.YLabels = extractedLabels;
                elseif length(extractedLabels) > app.chans
                    warning('通道信息(%d)多于数据通道(%d)，自动截取前 %d 个标签。', ...
                        length(extractedLabels), app.chans, app.chans);
                    app.YLabels = extractedLabels(1:app.chans);
                else
                    warning('通道信息数量不足！回退至默认编号。');
                    app.YLabels = arrayfun(@(x) sprintf('CH %d', x), 1:app.chans, 'UniformOutput', false);
                end
            end

            % 3. 构建左侧 UI 树节点 (ChannelTree)
            for i = 1:app.chans
                app.ChannelNode{i} = uitreenode(app.ChannelTree);
                app.ChannelNode{i}.Text = app.YLabels{i};
                app.ChannelTree.CheckedNodes = [app.ChannelTree.CheckedNodes; app.ChannelNode{i}];
                app.ChannelNode{i}.NodeData = i;
            end
        end

        function plotEEGChannels(app, lowlim, highlim)
            data = app.eegdata;
            dataRange = max(abs(data(:)));
            if dataRange < 1e-3
                data = data * 1e6;
            elseif dataRange < 1
                data = data * 1e3;
            end

            % ====== 决定可用通道索引 ======
            if isempty(app.CheckedData)
                availIdx = 1:app.chans;
                labelList = app.YLabels;
                badMask = false(1, app.chans);
                if ~isempty(app.eloc_file) && isfield(app.eloc_file,'badchan')
                    badMask = logical([app.eloc_file.badchan]);
                end
            else
                availIdx = app.CheckedData(:)';
                if isempty(app.eloc_file)
                    labelList = cellfun(@(x) num2str(x), num2cell(availIdx), 'UniformOutput', false);
                    badMask = false(1, numel(availIdx));
                else
                    labelList = app.YLabels(availIdx);
                    if isfield(app.eloc_file,'badchan')
                        badMask = logical([app.eloc_file(availIdx).badchan]);
                    else
                        badMask = false(1, numel(availIdx));
                    end
                end
                data = data(availIdx, :);
            end

            % ====== 基于 StartChan 选择"可见窗口"的通道 ======
            nAvail = numel(availIdx);
            app.dispchans = min(app.dispchans, nAvail);
            app.StartChan = min(max(1, app.StartChan), max(1, nAvail - app.dispchans + 1));
            visRange = app.StartChan : (app.StartChan + app.dispchans - 1);

            meandata = app.calcSubmean(data, lowlim, highlim);

            oldspacing = 580/app.dispchans;
            if app.duidie==1, app.spacing = oldspacing; else, app.spacing = 0; end

            app.timeStart = (lowlim - 1) / app.srate;
            app.timeEnd   = (highlim - 1) / app.srate;
            timePoints = linspace(app.timeStart, app.timeEnd, highlim - lowlim + 1);

            winDataVis = data(visRange, lowlim:highlim) - meandata(visRange);
            [winDataVisScaled, ~] = app.applyChannelScaling(winDataVis, app.ScaleFactor);
            [winDataVisRef,    ~] = app.applyChannelScaling(winDataVis, 0);

            app.LineHandles = gobjects(app.dispchans,1);
            app.TextHandles = gobjects(app.dispchans,1);

            for i = 1:app.dispchans
                row = app.dispchans - i + 1;
                yTrace = winDataVisScaled(row, :);
                if app.duidie==1
                    yOffset = i * app.spacing + (app.dispchans + 1)*(oldspacing - app.spacing)/2 ...
                        + app.elecoffset*(oldspacing - app.spacing);
                else
                    yOffset=0;
                end
                lh = plot(app.UIAxes, timePoints,  yTrace  + yOffset, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
                app.LineHandles(i) = lh;

                origIdx = availIdx( visRange(i) );
                lh.ButtonDownFcn = @(src, event) app.channelClicked(origIdx, event);

                %textX = app.UIAxes.XLim(1)-(0.02 * diff(app.UIAxes.XLim));
                textY = yOffset;
                textX = timePoints(1);
                th = text(app.UIAxes, textX, textY, [labelList{visRange(row)} ' '], ...
                    'HorizontalAlignment','right','Color','k','FontSize',8, ...
                    'Clipping','off');
                app.TextHandles(i) = th;
            end

            app.spacing = oldspacing;
            app.adjustUIAxes(lowlim, highlim, winDataVisRef);
            app.OriginalYLim2 = app.UIAxes.YLim;
            app.drawEventMarkers();

            if app.SelectedTime > 0
                app.drawTimeSelectionLine();
            end

            app.UIAxes.ButtonDownFcn = @(src, event) app.axesButtonDown(src, event);
            app.UIFigure.WindowButtonUpFcn = @(src, event) app.windowButtonUp(src, event);
            app.UIFigure.WindowButtonMotionFcn = @(src, event) app.windowButtonMotion(src, event);
            app.UIFigure.WindowScrollWheelFcn = @(src, event) app.windowScrollWheel(src, event);
        end
     
        function highlightSelectedChannel(app)
            if app.SelectedChannel > 0 && ~isempty(app.LineHandles)
                for i = 1:length(app.LineHandles)
                    if i == app.SelectedChannel
                        set(app.LineHandles(i), 'LineWidth', 2, 'Color', 'r');
                    else
                        set(app.LineHandles(i), 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
                    end
                end
            end
        end

        function channelClicked(app, channelIndex, event)
            cp = event.IntersectionPoint;
            app.SelectedTime = cp(1);

            app.drawTimeSelectionLine();
            x = cp(1); y = cp(2);

            if ~strcmp(get(app.UIFigure,'SelectionType'),'normal')
                return;
            end

            app.MouseDownTic   = tic;
            app.MouseDownPoint = [x y];

            if ~isempty(app.ClickHoldTimer) && isvalid(app.ClickHoldTimer)
                stop(app.ClickHoldTimer);
                delete(app.ClickHoldTimer);
            end
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end
            app.IsSelecting = false;
            app.SelectionStart = x;
            app.SelectionEnd   = x;

            app.ClickHoldTimer = timer( ...
                'ExecutionMode','singleShot', ...
                'StartDelay',0.3, ...
                'TimerFcn',@(~,~)beginDragSelect(app));
            start(app.ClickHoldTimer);
        end

        function recoverChannels(app)
            for i = 1:length(app.LineHandles)
                set(app.LineHandles(i), 'LineWidth', 1,'Color', [0.5 0.5 0.5]);
            end
        end

        function addTrialStartEvents(app)
            if isempty(app.TrialIndices)
                return;
            end

            nTrials = length(app.TrialIndices);
            baseUr  = numel(app.EventMarkers);

            for trialIdx = 1:nTrials
                startSample = app.TrialIndices(trialIdx).start_index;
                startTime   = (startSample - 1) / app.srate;

                newEvent = struct( ...
                    'type',     sprintf('TrialStart_%d', trialIdx), ...
                    'latency',  startTime, ...
                    'duration', 0, ...
                    'urevent',  baseUr + trialIdx, ...
                    'epoch',    trialIdx ...
                    );

                if ~isempty(app.EventMarkers)
                    latAll  = [app.EventMarkers.latency];
                    typeAll = {app.EventMarkers.type};
                    if any(abs(latAll - startTime) < 1e-9 & strcmp(typeAll, newEvent.type))
                        continue;
                    end
                end

                app.EventMarkers = app.appendEvent(app.EventMarkers, newEvent);
            end

            if ~isempty(app.EventMarkers)
                [~, idx] = sort([app.EventMarkers.latency]);
                app.EventMarkers = app.EventMarkers(idx);
            end

            app.updateEventList();
            app.drawEventMarkers();
        end

        function updateEventList(app)
            app.eventStrings = cell(1, length(app.EventMarkers));
            for i = 1:length(app.EventMarkers)
                app.eventStrings{i} = sprintf('%.3f s: %s', ...
                    app.EventMarkers(i).latency, ...
                    app.EventMarkers(i).type);
            end

            app.EventListListBox.Items = app.eventStrings;
        end

        function drawEventMarkers(app)
            delete(app.EventLineHandles(ishandle(app.EventLineHandles)));
            delete(app.EventTextHandles(ishandle(app.EventTextHandles)));
            app.EventLineHandles = [];
            app.EventTextHandles = [];

            if isempty(app.EventMarkers)
                return;
            end

            app.EventLineHandles = gobjects(length(app.EventMarkers),1);
            app.EventTextHandles = gobjects(length(app.EventMarkers),1);

            for i = 1:length(app.EventMarkers)
                eventTime = app.EventMarkers(i).latency;

                if isempty(eventTime) || ~isnumeric(eventTime) || ~isscalar(eventTime) || ~isfinite(eventTime)
                    continue;
                end

                evType = app.EventMarkers(i).type;
                if isempty(evType)
                    evType = '';
                elseif ~ischar(evType)
                    evType = char(string(evType));
                end

                if eventTime >= app.timeStart && eventTime <= app.timeEnd
                    lh = line(app.UIAxes, [eventTime, eventTime], app.UIAxes.YLim, ...
                        'Color', 'b', 'LineWidth', 1.5, 'LineStyle', '--');
                    app.EventLineHandles(i) = lh;

                    th2 = text(app.UIAxes, eventTime, app.UIAxes.YLim(2), evType, ...
                        'Color', 'b', 'FontSize', 8, 'VerticalAlignment', 'bottom','Clipping','off');
                    app.EventTextHandles(i) = th2;
                end
            end
        end

        function adjustUIAxes(app, ~,~, winDataVisRef)
            app.UIAxes.XLim = [app.timeStart, app.timeEnd];
            xticks = app.timeStart:0.2:app.timeEnd;
            if xticks(end) ~= app.timeEnd
                xticks = [xticks, app.timeEnd];
            end
            app.UIAxes.XTick = unique(xticks);
            app.UIAxes.YTick = [];

            if app.dispchans > 0
                if app.duidie == 1
                    totalHeight = app.spacing * (app.dispchans + 1);
                    app.UIAxes.YLim = [0, totalHeight];
                else
                    lo = prctile(winDataVisRef(:), 2);
                    hi = prctile(winDataVisRef(:), 98);
                    pad = 0.1 * (hi - lo + eps);
                    if hi > lo
                        app.UIAxes.YLim = [lo - pad, hi + pad];
                    end
                end
            end

            % 注意:UIAxes 现在被放进了 GridLayout,不再设 normalized 位置
            app.UIAxes.Box = 'off';
        end

        function drawTimeSelectionLine(app)
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end
            app.TimeSelectionLine = line(app.UIAxes, [app.SelectedTime, app.SelectedTime], app.OriginalYLim2, ...
                'Color', 'r', 'LineWidth', 1.5);
            app.SelectedTimeEditField.Value=app.SelectedTime;

            timepoint=round(app.SelectedTime*app.srate) + 1;

            if timepoint>app.datapoints
                return
            end

            targetFigs = findall(0, 'Type', 'figure', 'Tag', 'SEAL_ResultWindow');
            
            for i = 1:length(targetFigs)
                try
                    app.CallingAppResult = targetFigs(i).UserData;
                    
                    if isvalid(app.CallingAppResult)
                        app.CallingAppResult.TimePointValueChanged(app.SelectedTime);
                    end
                catch
                    disp('Datachange failed.Please select the data in first trial and try again.');
                end
            end
        end

        function axesButtonDown(app, src, ~)        
            cp = src.CurrentPoint;
            x = cp(1,1); y = cp(1,2);
            if ~strcmp(get(app.UIFigure,'SelectionType'),'normal')
                return;
            end

            app.MouseDownTic   = tic;
            app.MouseDownPoint = [x y];

            if ~isempty(app.ClickHoldTimer) && isvalid(app.ClickHoldTimer)
                stop(app.ClickHoldTimer);
                delete(app.ClickHoldTimer);
            end
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end
            app.IsSelecting = false;
            app.SelectionStart = x;
            app.SelectionEnd   = x;

            app.ClickHoldTimer = timer( ...
                'ExecutionMode','singleShot', ...
                'StartDelay',0.3, ...
                'TimerFcn',@(~,~)beginDragSelect(app));
            start(app.ClickHoldTimer);
        end

        function windowButtonUp(app,~, ~)
            cp_fig = app.UIFigure.CurrentPoint;
            axPos  = getpixelposition(app.UIAxes, true);

            isInAxes = ...
                cp_fig(1) >= axPos(1) && ...
                cp_fig(1) <= axPos(1)+axPos(3) && ...
                cp_fig(2) >= axPos(2) && ...
                cp_fig(2) <= axPos(2)+axPos(4);

            if ~isInAxes
                return
            end
            if ~isempty(app.ClickHoldTimer) && isvalid(app.ClickHoldTimer)
                stop(app.ClickHoldTimer);
                delete(app.ClickHoldTimer);
            end
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end
            dt = 0;
            if app.MouseDownTic ~= 0
                dt = toc(app.MouseDownTic);
            end
            app.MouseDownTic = uint64(0);

            if ~app.IsSelecting && dt < 0.3
                cp = app.UIAxes.CurrentPoint;
                x = cp(1,1);
                app.SelectedTime = x;
                if ishandle(app.TimeSelectionLine)
                    delete(app.TimeSelectionLine);
                end
                app.TimeSelectionLine = line(app.UIAxes, [x x], app.UIAxes.YLim, ...
                    'Color','r','LineWidth',1.5);
                app.SelectedTimeEditField.Value = app.SelectedTime;
            elseif app.IsSelecting
                startTime = min(app.SelectionStart, app.SelectionEnd);
                endTime   = max(app.SelectionStart, app.SelectionEnd);
                app.SelectedTimeRange = [startTime endTime];
                app.highlightSelectedTimeRange();
            end
            app.IsSelecting      = false;
            app.MouseDownTic     = uint64(0);
            app.MouseDownPoint   = [];
            app.SelectionStart   = [];
            app.SelectionEnd     = [];
            drawnow;
        end

        function highlightSelectedTimeRange(app)
            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end

            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.ChannelSelectionText)
                delete(app.ChannelSelectionText);
            end

            app.recoverChannels();
            if ~isempty(app.SelectedTimeRange)
                startTime = app.SelectedTimeRange(1);
                endTime = app.SelectedTimeRange(2);
                yLimits = app.UIAxes.YLim;

                x = [startTime, endTime, endTime, startTime];
                y = [yLimits(1), yLimits(1), yLimits(2), yLimits(2)];
                app.SelectionRectangle = patch(app.UIAxes, x, y, [0.9 0.9 0.9], ...
                    'EdgeColor', 'r', ...
                    'FaceAlpha', 0.3);
                app.rangestart.Text=sprintf('%.4f s', startTime);
                app.rangeend.Text=sprintf('%.4f s', endTime);
                app.SelectedTimeEditField.Value = 0;
            end
            starttimepoint=round(startTime*app.srate)+1;
            endtimepoint=round(endTime*app.srate)+1;
            if starttimepoint>app.datapoints||endtimepoint>app.datapoints
                return
            end
            
            targetFigs = findall(0, 'Type', 'figure', 'Tag', 'SEAL_ResultWindow');
            
            for i = 1:length(targetFigs)
                try
                    app.CallingAppResult = targetFigs(i).UserData;
                    
                    if isvalid(app.CallingAppResult)
                        app.CallingAppResult.datachanged(starttimepoint,endtimepoint);
                    end
                catch
                    disp('Datachange failed. Please select the data in first trial and try again.');
                end
            end
        end

        function windowButtonMotion(app, ~, ~)
            if ~app.IsSelecting, return; end

            cp = app.UIAxes.CurrentPoint;
            x = cp(1,1);
            yLimits = app.UIAxes.YLim;

            app.SelectionEnd = x;
            if ishandle(app.SelectionRectangle)
                xs = sort([app.SelectionStart, app.SelectionEnd]);
                set(app.SelectionRectangle, 'XData', [xs(1) xs(2) xs(2) xs(1)], ...
                    'YData', [yLimits(1) yLimits(1) yLimits(2) yLimits(2)]);
            end
        end

        function windowScrollWheel(app, ~, event)
            mousePos = get(app.UIFigure, 'CurrentPoint');
            axPos = getpixelposition(app.UIAxes, true);

            isMouseInAxes = mousePos(1) >= axPos(1) && ...
                mousePos(1) <= axPos(1) + axPos(3) && ...
                mousePos(2) >= axPos(2) && ...
                mousePos(2) <= axPos(2) + axPos(4);

            if ~isMouseInAxes
                return;
            end

            cp = app.UIAxes.CurrentPoint(1,1:2);
            if any(isnan(cp)) || isempty(cp), return; end

            baseFactor = 1.1;
            currentXRange = diff(app.UIAxes.XLim);
            if currentXRange < 1
                baseFactor = 1.05;
            end

            zoomFactor = baseFactor;
            if event.VerticalScrollCount < 0
                zoomFactor = 1 / baseFactor;
            end

            xRange = diff(app.UIAxes.XLim);
            yRange = diff(app.UIAxes.YLim);
            xRatio = (cp(1) - app.UIAxes.XLim(1)) / xRange;
            yRatio = (cp(2) - app.UIAxes.YLim(1)) / yRange;

            newXLim = cp(1) + [-xRatio, (1-xRatio)] * xRange * zoomFactor;
            newYLim = cp(2) + [-yRatio, (1-yRatio)] * yRange * zoomFactor;
            newXLim = app.constrainRange(newXLim, [0, app.timeEnd]);
            newYLim = app.constrainRange(newYLim, app.OriginalYLim2);

            app.UIAxes.XLim = newXLim;
            app.UIAxes.YLim = newYLim;

            xLeft = app.UIAxes.XLim(1);
            %xOffset = 0.02 * diff(app.UIAxes.XLim);
            for i = 1:numel(app.TextHandles)
                if isvalid(app.TextHandles(i))
                    pos = app.TextHandles(i).Position;
                    pos(1) = xLeft;
                    app.TextHandles(i).Position = pos;
                end
            end

            for i = 1:numel(app.EventTextHandles)
                if isa(app.EventTextHandles(i),'matlab.graphics.primitive.Text')
                    if isvalid(app.EventTextHandles(i))
                        pos = app.EventTextHandles(i).Position;
                        pos(2) = app.UIAxes.YLim(2);
                        app.EventTextHandles(i).Position = pos;
                    end
                end
            end
        end

        function lim = constrainRange(app,lim, validRange)
            minRange = 0.1 * diff(validRange);
            lim(1) = max(validRange(1), min(lim(1), validRange(2)-minRange));
            lim(2) = min(validRange(2), max(lim(2), validRange(1)+minRange));
        end
       
        function meandata = calcSubmean(app, data, lowlim, highlim)
            if strcmpi(app.submean, 'on')
                minData = data(:, lowlim:highlim);
                meandata = mean(minData, 2);
            else
                meandata = zeros(app.chans, 1);
            end
        end

        function historyDataMenuSelected(app, src,~)
            filename=src.Text;
            data = load(src.UserData);
            disp(['Selected file: ', filename]);
            app.eegdata=data.Data;
            app.redrawPlot();
        end

        function concatenated_data = concatenate_trials(app, data_3d)
            [n_channels, n_points_per_trial, n_trials] = size(data_3d);

            total_points = n_points_per_trial * n_trials;
            concatenated_data = zeros(n_channels, total_points);

            app.TrialIndices = struct(...
                'trial_number', num2cell(1:n_trials), ...
                'start_index', num2cell(zeros(1, n_trials)), ...
                'end_index', num2cell(zeros(1, n_trials)) ...
                );

            for trial = 1:n_trials
                start_idx = (trial - 1) * n_points_per_trial + 1;
                end_idx = trial * n_points_per_trial;

                current_trial_data = data_3d(:, :, trial);

                concatenated_data(:, start_idx:end_idx) = current_trial_data;

                app.TrialIndices(trial).start_index = start_idx;
                app.TrialIndices(trial).end_index = end_idx;
            end

            app.TrialIndices(1).n_channels = n_channels;
            app.TrialIndices(1).points_per_trial = n_points_per_trial;
            app.TrialIndices(1).total_trials = n_trials;
            app.TrialIndices(1).total_points = total_points;
        end

        function addSeamEventsFromSegmentInfo(app)
            if ~isprop(app, 'segmentInfo') || isempty(app.segmentInfo)
                return;
            end
            if size(app.segmentInfo,1) ~= 2 || size(app.segmentInfo,2) < 2
                return;
            end
            if isempty(app.srate) || app.srate <= 0
                error('app.srate 未设置或非法。');
            end

            seamSamples = app.segmentInfo(2, 2:end);
            seamTimes   = (seamSamples - 1) ./ app.srate;
            nSeams      = numel(seamTimes);

            baseUr = numel(app.EventMarkers);

            for k = 1:nSeams
                if ~isempty(app.EventMarkers)
                    latAll = [app.EventMarkers.latency];
                    if any(abs(latAll - seamTimes(k)) < 1e-9)
                        continue;
                    end
                end

                newEvent = struct( ...
                    'type',     sprintf('Seam_%d (Trial %d end)', k, k), ...
                    'latency',  seamTimes(k), ...
                    'duration', 0, ...
                    'urevent',  baseUr + k, ...
                    'epoch',    [] ...
                    );

                app.EventMarkers = app.appendEvent(app.EventMarkers, newEvent);
            end

            app.updateEventList();
            app.drawEventMarkers();
        end

        function onChannelSliderChanging(app, val)
            maxStart = app.chans-60+1;

            startChan = maxStart - val + 1 ;
            app.StartChan = seal_clamp(round(startChan), 1, maxStart);

            app.redrawPlot();
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end

            if ishandle(app.ChannelSelectionText)
                delete(app.ChannelSelectionText);
            end
            app.recoverChannels();
        end

        function onChannelSliderChanged(app, val)
            maxStart =app.chans-60+1;

            startChan = maxStart -val + 1;
            app.StartChan = seal_clamp(round(startChan), 1, maxStart);

            app.redrawPlot();
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end

            if ishandle(app.ChannelSelectionText)
                delete(app.ChannelSelectionText);
            end
            app.recoverChannels();
        end

        function target = getTargetAmplitude(app)
            if isempty(app.ScaleTarget)
                s = max(app.spacing, 580/app.dispchans);
                target = 0.35 * s;
            else
                target = app.ScaleTarget;
            end
        end

        function [scaledWin, scaleInfo] = applyChannelScaling(app, winData, f)
            target = app.getTargetAmplitude();
            nCh = size(winData, 1);

            if f == 0
                rel = 1;
            else
                rel = f;
            end

            switch lower(app.ScaleMode)
                case 'none'
                    scaledWin = winData * rel;
                    scaleInfo = rel * ones(nCh, 1);

                case 'std'
                    s = std(winData, 0, 2);  s(s == 0) = 1;
                    fVec = seal_clampScale((target ./ s) * rel);
                    scaledWin = winData .* fVec;
                    scaleInfo = fVec;

                case 'minmax'
                    mn = min(winData, [], 2); mx = max(winData, [], 2);
                    rng = (mx - mn) / 2;     rng(rng == 0) = 1;
                    centered = winData - (mn + mx) / 2;
                    fVec = seal_clampScale((target ./ rng) * rel);
                    scaledWin = centered .* fVec;
                    scaleInfo = fVec;

                case 'zscore'
                    mu = mean(winData, 2);
                    s  = std(winData, 0, 2);  s(s == 0) = 1;
                    fVec = seal_clampScale((target ./ s) * rel);
                    scaledWin = (winData - mu) .* fVec;
                    scaleInfo = fVec;

                case 'percentile'
                    lo = prctile(winData, app.ScalePercentiles(1), 2);
                    hi = prctile(winData, app.ScalePercentiles(2), 2);
                    c  = (hi + lo) / 2;
                    h  = (hi - lo) / 2;  h(h == 0) = 1;
                    fVec = seal_clampScale((target ./ h) * rel);
                    scaledWin = (winData - c) .* fVec;
                    scaleInfo = fVec;

                case 'global'
                    if ~isempty(app.AutoScaleFactor)
                        base = app.AutoScaleFactor;
                    else
                        allData = winData(:);
                        lo = prctile(allData, app.ScalePercentiles(1));
                        hi = prctile(allData, app.ScalePercentiles(2));
                        h  = (hi - lo) / 2;
                        if h == 0, h = 1; end
                        base = target / h;
                    end
                    fScalar = base * rel;
                    scaledWin = winData * fScalar;
                    scaleInfo = fScalar * ones(nCh, 1);

                otherwise
                    scaledWin = winData * rel;
                    scaleInfo = rel * ones(nCh, 1);
            end
        end

        function applySelectionHighlight(app, origIdx)
            checked = app.ChannelTree.CheckedNodes;
            if isempty(checked), return; end
            checkedIdx = cellfun(@(n) n.NodeData, num2cell(checked));
            if ~ismember(origIdx, checkedIdx), return; end

            if isempty(app.CheckedData)
                availIdx = 1:app.chans;
            else
                availIdx = app.CheckedData(:)';
            end

            posInAvail = find(availIdx == origIdx, 1);
            if isempty(posInAvail), return; end

            startPos = app.StartChan;
            endPos   = min(app.StartChan + app.dispchans - 1, numel(availIdx));
            if posInAvail < startPos || posInAvail > endPos
                return;
            end

            localIdx  = posInAvail - startPos + 1;
            screenIdx = app.dispchans - localIdx + 1;

            app.SelectedChannel = screenIdx;
            app.highlightSelectedChannel();
        end
   
    end

    methods (Access = public)

        function [freq, powerSpectrum] = computeChannelFFT(app, EEGdata,channelIdx)
            data = EEGdata(channelIdx, :);
            data = data - mean(data);

            if mod(length(data), 2) ~= 0
                data = data(1:end-1);
            end
            N = length(data);

            window = hann(N)';
            windowedData = data .* window;

            fftResult = fft(windowedData);
            P2 = abs(fftResult/N);

            halfN = floor(N/2);
            P1 = P2(1:halfN + 1);
            P1(2:end-1) = 2*P1(2:end-1);

            freq = app.srate*(0:halfN)/N;

            powerSpectrum = (P1.^2) / (sum(window.^2) * app.srate);

            validRange = freq >= 0.01 & freq <= 300;
            freq = freq(validRange);
            powerSpectrum = powerSpectrum(validRange);
        end

        function plotMultiChannelSpectrum(app, EEGdata,channelIndices)
            fig = figure('Name', 'Multi-Channel Power Spectrum', ...
                'Position', [100 100 800 600]);
            ax = axes(fig);
            hold(ax, 'on');

            colors = lines(length(channelIndices));

            for i = 1:length(channelIndices)
                chanIdx = channelIndices(i);

                [freq, spectrum] = app.computeChannelFFT(EEGdata,chanIdx);

                plot(ax, freq, 10*log10(spectrum), ...
                    'Color', colors(i,:), ...
                    'LineWidth', 1.5, ...
                    'DisplayName', sprintf('Ch%d: %s', chanIdx, app.YLabels{chanIdx}));
            end

            grid(ax, 'on');
            xlabel(ax, 'Frequency (Hz)');
            ylabel(ax, 'Power (dB)');
            title(ax, 'Multi-Channel EEG Power Spectrum');
            legend(ax, 'show', 'Location', 'northeastoutside');

            bandColors = {'r', 'g', 'm'};
            bands = [1 4; 4 8; 8 13];
            for b = 1:size(bands,1)
                x = [bands(b,1), bands(b,1), bands(b,2), bands(b,2)];
                y = [ax.YLim(1), ax.YLim(2), ax.YLim(2), ax.YLim(1)];
                patch(x, y, bandColors{b}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            end

            text(2.5, ax.YLim(2)-3, 'Delta', 'Color', 'r');
            text(6, ax.YLim(2)-3, 'Theta', 'Color', 'g');
            text(10.5, ax.YLim(2)-3, 'Alpha', 'Color', 'm');
        end

        function moveTimeSelection(app, delta)
            if isempty(app.SelectedTime)
                app.SelectedTime = app.time;
            end

            newTime = app.SelectedTime + delta;

            if newTime < 0
                newTime = 0;
            elseif newTime > app.frames/app.srate
                newTime = app.frames/app.srate;
            end

            app.SelectedTime = newTime;

            app.drawTimeSelectionLine();
        end

        function beginDragSelect(app)
            if isempty(app.MouseDownPoint) || (exist('app.MouseDownTic','var') && app.MouseDownTic==0)
                return;
            end
            app.IsSelecting = true;

            if ishandle(app.SelectionRectangle), delete(app.SelectionRectangle); end
            if ishandle(app.TimeSelectionLine),  delete(app.TimeSelectionLine);  end
            if ishandle(app.ChannelSelectionText), delete(app.ChannelSelectionText); end

            yLimits = app.UIAxes.YLim;
            x0 = app.SelectionStart;
            app.SelectionRectangle = patch(app.UIAxes, ...
                [x0 x0 x0 x0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
                'r', 'EdgeColor','r', 'LineWidth',1, 'LineStyle','--', 'FaceColor','none');
        end


        function computeAutoScaleFactor(app)
            data = app.eegdata;
            dataRange = max(abs(data(:)));
            if dataRange < 1e-3
                data = data * 1e6;
            elseif dataRange < 1
                data = data * 1e3;
            end

            target = app.getTargetAmplitude();
            switch lower(app.ScaleMode)
                case 'none'
                    f = 1;
                case 'std'
                    s = std(data, 0, 2); s(s==0) = 1;
                    f = target ./ s;
                case 'minmax'
                    rng = (max(data,[],2) - min(data,[],2))/2;  rng(rng==0)=1;
                    f = target ./ rng;
                case 'percentile'
                    lo = prctile(data, app.ScalePercentiles(1), 2);
                    hi = prctile(data, app.ScalePercentiles(2), 2);
                    h = (hi-lo)/2; h(h==0)=1;
                    f = target ./ h;
                case 'zscore'
                    s = std(data,0,2); s(s==0)=1;
                    f = target ./ s;
                case 'global'
                    lo = prctile(data(:), app.ScalePercentiles(1));
                    hi = prctile(data(:), app.ScalePercentiles(2));
                    h = (hi-lo)/2; if h==0, h=1; end
                    f = target / h;
                otherwise
                    f = 1;
            end
            app.AutoScaleFactor = f;
        end


        function importEEGLABEvents(app, eventStruct)
            if isempty(eventStruct) || ~isfield(eventStruct,'latency') || ~isfield(eventStruct,'type')
                return;
            end

            hasEpochField = isfield(eventStruct, 'epoch');
            is3D = ~isempty(app.TrialIndices) && ...
                isfield(app.TrialIndices, 'start_index') && ...
                ~isempty([app.TrialIndices.start_index]);

            baseUr = numel(app.EventMarkers);

            for i = 1:numel(eventStruct)
                lat = eventStruct(i).latency;
                if isempty(lat) || ~isnumeric(lat) || ~isscalar(lat) || ~isfinite(lat)
                    continue;
                end

                epNum = [];
                if hasEpochField && ~isempty(eventStruct(i).epoch)
                    epNum = double(eventStruct(i).epoch);
                end

                if is3D && ~isempty(epNum) && epNum >= 1 && epNum <= numel(app.TrialIndices)
                    nPts = app.TrialIndices(1).points_per_trial;
                    if double(lat) <= nPts
                        lat = app.TrialIndices(epNum).start_index + double(lat) - 1;
                    end
                end

                latSec = (double(lat) - 1) / app.srate;

                t = eventStruct(i).type;
                if isnumeric(t)
                    typeStr = num2str(t);
                elseif ischar(t) || isstring(t)
                    typeStr = char(t);
                else
                    typeStr = 'Event';
                end

                durVal = 0;
                if isfield(eventStruct, 'duration') && ~isempty(eventStruct(i).duration)
                    d = eventStruct(i).duration;
                    if isnumeric(d) && isscalar(d) && isfinite(d)
                        durVal = double(d) / app.srate;
                    end
                end

                urVal = baseUr + i;
                if isfield(eventStruct, 'urevent') && ~isempty(eventStruct(i).urevent)
                    urVal = double(eventStruct(i).urevent);
                end

                newEvent = struct( ...
                    'type',     typeStr, ...
                    'latency',  latSec, ...
                    'duration', durVal, ...
                    'urevent',  urVal, ...
                    'epoch',    epNum ...
                    );

                if ~isempty(app.EventMarkers)
                    latAll  = [app.EventMarkers.latency];
                    typeAll = {app.EventMarkers.type};
                    if any(abs(latAll - latSec) < 1e-9 & strcmp(typeAll, typeStr))
                        continue;
                    end
                end

                app.EventMarkers = app.appendEvent(app.EventMarkers, newEvent);
            end
        end

        function jumpToTime(app, targetTime, rangeEnd)
            if nargin < 3, rangeEnd = []; end

            maxTime = app.frames / app.srate;
            targetTime = max(0, min(targetTime, maxTime));

            if ~isempty(rangeEnd) && rangeEnd > targetTime
                rangeEnd  = max(0, min(rangeEnd, maxTime));
                app.winlength = rangeEnd - targetTime;
                app.WinlenthEditField.Value = app.winlength;
                app.time = targetTime;
                app.SelectedTime = 0;
                app.redrawPlot();
            else
                newStart = targetTime - app.winlength/2;
                newStart = max(0, min(newStart, max(0, maxTime - app.winlength)));
                app.time = newStart;
                app.SelectedTime = targetTime;
                app.redrawPlot();
                app.drawTimeSelectionLine();
            end
        end

        function vals = askFilterParams(app, dlgTitle, labels, defaults)
            n    = numel(labels);
            rowH = 32;
            figW = 320;
            figH = 50 + n*rowH + 50;

            sz = get(0,'ScreenSize');
            d  = uifigure('Name', dlgTitle, ...
                'Position', [sz(3)/2-figW/2, sz(4)/2-figH/2, figW, figH], ...
                'WindowStyle','modal','Resize','off');

            eds = gobjects(n,1);
            for i = 1:n
                y = figH - 30 - i*rowH;
                uilabel(d, 'Text', labels{i}, ...
                    'HorizontalAlignment','right', ...
                    'Position', [15, y, 150, 22]);
                eds(i) = uieditfield(d, 'numeric', ...
                    'Value', defaults(i), ...
                    'Position', [175, y, 125, 22]);
            end

            vals = [];
            uibutton(d,'Text','OK',    'Position',[60,15,80,28], ...
                'ButtonPushedFcn',@(~,~) onOK());
            uibutton(d,'Text','Cancel','Position',[180,15,80,28], ...
                'ButtonPushedFcn',@(~,~) onCancel());
            d.CloseRequestFcn = @(~,~) onCancel();

            uiwait(d);
            if isvalid(d), delete(d); end

            function onOK()
                vals = arrayfun(@(e) e.Value, eds);
                uiresume(d);
            end
            function onCancel()
                vals = [];
                uiresume(d);
            end
        end

        function applyFilterAndUpdate(app, firCoeff, paramStr, cutoffs, N, Rs)
            checkedNodes = app.ChannelTree.CheckedNodes;
            if isempty(checkedNodes)
                chIdx = 1:app.chans;
            else
                chIdx = arrayfun(@(n) n.NodeData, checkedNodes);
            end

            filtered = app.eegdata;
            for k = 1:numel(chIdx)
                ch = chIdx(k);
                filtered(ch,:)= filtfilt(firCoeff, 1, double(app.eegdata(ch,:)));
            end

            app.computeAutoScaleFactor();
            app.redrawPlot();

            app.visualizeFilterResponse(firCoeff, app.srate, paramStr, cutoffs);
            app.plotMultiChannelSpectrum(filtered, chIdx);
        end


        function visualizeFilterResponse(~, b, Fs, titleStr, cutoffs)
            nfft = 1024;
            [H, W] = freqz(b, 1, nfft, Fs);

            fig = figure('Name', titleStr, 'Color', 'w');

            ax1 = subplot(2,1,1, 'Parent', fig);
            magdB = 20*log10(abs(H) + eps);
            plot(ax1, W, magdB, 'LineWidth', 1.2);
            grid(ax1, 'on'); hold(ax1, 'on');
            xlabel(ax1, 'Frequency (Hz)');
            ylabel(ax1, 'Magnitude (dB)');
            title(ax1, titleStr, 'Interpreter', 'none');
            xlim(ax1, [0 Fs/2]);

            yl = [max(min(magdB), -120), max(magdB) + 5];
            ylim(ax1, yl);

            for k = 1:numel(cutoffs)
                plot(ax1, [cutoffs(k) cutoffs(k)], yl, 'r--');
                text(ax1, cutoffs(k), yl(1) + 0.1*diff(yl), ...
                    sprintf(' %.2f Hz', cutoffs(k)), 'Color', 'r');
            end

            ax2 = subplot(2,1,2, 'Parent', fig);
            plot(ax2, W, unwrap(angle(H)) * 180/pi, 'LineWidth', 1.2);
            grid(ax2, 'on');
            xlabel(ax2, 'Frequency (Hz)');
            ylabel(ax2, 'Phase (deg)');
            xlim(ax2, [0 Fs/2]);
        end

  
        function evArr = appendEvent(~, evArr, newEvent)
            if isempty(evArr)
                evArr = newEvent;
                return;
            end

            oldFields = fieldnames(evArr);
            newFields = fieldnames(newEvent);
            allFields = union(oldFields, newFields, 'stable');

            for i = 1:numel(allFields)
                f = allFields{i};
                if ~isfield(evArr, f)
                    [evArr.(f)] = deal([]);
                end
                if ~isfield(newEvent, f)
                    newEvent.(f) = [];
                end
            end

            evArr    = orderfields(evArr, allFields);
            newEvent = orderfields(newEvent, allFields);
            evArr(end+1) = newEvent;
        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, dataNode, position)
            myPos = app.UIFigure.Position;
            w = myPos(3); h = myPos(4);

            if nargin >= 3 && ~isempty(position) && numel(position) == 4
                mainX = position(1);
                mainY = position(2);
                mainW = position(3);
                mainH = position(4);

                newX = mainX + (mainW - w) / 2;
                newY = mainY + mainH - 10;

                screenSize = get(0,'ScreenSize');
                screenW = screenSize(3); screenH = screenSize(4);
                if newY + h > screenH
                    newY = max(1, screenH - h - 40);
                end
                newX = max(1, min(newX, screenW - w));

                app.UIFigure.Position = [newX, newY, w, h];
            else
                screenSize = get(0,'ScreenSize');
                app.UIFigure.Position(1) = screenSize(3) * 0.4;
                app.UIFigure.Position(2) = screenSize(3) * 0.3;
            end

            app.dataNode = dataNode;
            app.SavePath=dataNode.parent.parent.path;
            app.UIFigure.Name=dataNode.name;
            app.srate=double(dataNode.srate);
            tempdata = dataNode.data; 
            

            if ndims(tempdata)==3
                app.eegdata= app.concatenate_trials(tempdata);
            else
                app.eegdata= tempdata;
            end

            app.chans     = size(app.eegdata, 1);
            app.frames    = size(app.eegdata, 2);
            app.datapoints = app.frames;
            app.time      = 0;
            app.winlength = min(app.frames / app.srate, 10);
            app.timeStart = 0;
            app.timeEnd   = app.winlength;


            if ~isempty(dataNode.dataInfo.chanlocs)
                app.chanlocs = dataNode.dataInfo.chanlocs;
            end
            app.chans=size(app.eegdata,1);
            app.dispchans=min(app.chans,60);
            app.frames=size(app.eegdata,2);
            app.datapoints=size(app.eegdata,2);

           
            evt = dataNode.dataInfo.event;

            if ~isempty(evt)
                app.importEEGLABEvents(evt);
                
            end


            if isfield(dataNode.data,'segmentInfo')
                app.segmentInfo=dataNode.data.segmentInfo;
                app.addSeamEventsFromSegmentInfo();
            end
            if ndims(tempdata)==3
                app.addTrialStartEvents();
            end

            app.initializeChannelLabels();

            % 创建滑动条 - 放进 CenterGrid 的右侧第二列
            if app.chans>60
                app.ChannelSlider = uislider(app.CenterGrid, ...
                    'Orientation','vertical', ...
                    'Limits',[0, max(0, app.chans-60)+1], ...
                    'Value',max(0, app.chans-60+1), ...
                    'MinorTicks',[],'MajorTicks',[],...
                    'ValueChangingFcn',@(src,evt) onChannelSliderChanging(app, evt.Value), ...
                    'ValueChangedFcn', @(src,evt) onChannelSliderChanged(app, src.Value));
                app.ChannelSlider.Layout.Row = 1;
                app.ChannelSlider.Layout.Column = 2;
            end

            if ~isempty(app.EventMarkers)
                [~, idx] = sort([app.EventMarkers.latency]);
                app.EventMarkers = app.EventMarkers(idx);
            end

            app.updateEventList();
            app.redrawPlot();
            

            app.WinlenthEditField.Value=app.UIAxes.XLim(2)-app.UIAxes.XLim(1);
            app.winlength= app.WinlenthEditField.Value;

            app.computeAutoScaleFactor();
        end

        % Button pushed function: LLButton
        function LLButtonPushed(app, event)
            step = app.winlength ;
            app.time =max(0, app.time - step);
            app.redrawPlot;
        end

        % Button pushed function: smallerButton
        function smallerButtonPushed(app, event)
            if app.winlength==app.frames/app.srate
                app.winlength=fix(app.winlength*10)/10;
                app.WinlenthEditField.Value = app.winlength;
                app.redrawPlot;
                return;
            end
            if app.winlength>0.1
                app.winlength=app.winlength-0.1;
                app.WinlenthEditField.Value = app.winlength;
                app.redrawPlot;
            else
                return;
            end
        end

        % Button pushed function: RRButton
        function RRButtonPushed(app, event)
            step = app.winlength ;
            rightbountry = app.time + step;
            rightbountry = min(app.frames/app.srate,rightbountry+step);
            app.time = rightbountry - step;

            app.redrawPlot;
        end

        % Button pushed function: biggerButton
        function biggerButtonPushed(app, event)
            app.winlength=app.winlength+0.1;
            if app.winlength<app.frames/app.srate
                app.WinlenthEditField.Value = app.winlength;
                app.redrawPlot;
            else
                app.winlength=app.frames/app.srate;
                app.WinlenthEditField.Value = app.winlength;
                app.redrawPlot;
            end
        end

        % Button pushed function: RButton
        function RButtonPushed(app, event)
            step = app.winlength/2 ;
            rightbountry = app.time + step;
            rightbountry = min(app.frames/app.srate,rightbountry+step);
            app.time = rightbountry - step;

            app.redrawPlot;
        end

        % Button pushed function: LButton
        function LButtonPushed(app, event)
            step = app.winlength/2 ;
            app.time = max(0, app.time - step);

            app.redrawPlot;
        end

        % Button pushed function: AddEventButton
        function AddEventButtonPushed(app, event)
            eventTimeSec = app.SelectedTimeEditField.Value;
            eventLabel   = strtrim(app.EventLabelEditField.Value);

            maxT = (app.frames - 1) / app.srate;
            if eventTimeSec < 0 || eventTimeSec > maxT
                uialert(app.UIFigure, ...
                    sprintf('Latency %.3fs 超出数据范围 [0, %.3f]s', eventTimeSec, maxT), ...
                    'Invalid latency');
                return;
            end

            if isempty(eventLabel)
                uialert(app.UIFigure, 'Event type 不能为空', 'Invalid type');
                return;
            end

            epNum = [];
            if ~isempty(app.TrialIndices) && isfield(app.TrialIndices,'points_per_trial')
                sampleIdx = round(eventTimeSec * app.srate) + 1;
                nPts = app.TrialIndices(1).points_per_trial;
                epCandidate = ceil(sampleIdx / nPts);
                if epCandidate >= 1 && epCandidate <= numel(app.TrialIndices)
                    epNum = epCandidate;
                end
            end

            newEvent = struct( ...
                'type',     eventLabel, ...
                'latency',  eventTimeSec, ...
                'duration', 0, ...
                'urevent',  numel(app.EventMarkers) + 1, ...
                'epoch',    epNum ...
                );

            app.EventMarkers = app.appendEvent(app.EventMarkers, newEvent);
            app.EventData    = app.appendEvent(app.EventData,    newEvent);

            if ~isempty(app.EventMarkers)
                [~, idx]  = sort([app.EventMarkers.latency]);
                app.EventMarkers = app.EventMarkers(idx);
            end
            if ~isempty(app.EventData)
                [~, idx2] = sort([app.EventData.latency]);
                app.EventData = app.EventData(idx2);
            end
            app.EventLabelEditField.Value = '';

            app.updateEventList();
            app.redrawPlot();
        end

        % Value changed function: StackingmodeButton
        function StackingmodeButtonValueChanged(app, event)
            value = app.StackingmodeButton.Value;
            if value==1
                app.duidie=0;
                app.redrawPlot;
                for i = 1:length(app.TextHandles)
                    set(app.TextHandles(i), 'Visible', 'off');
                end
            else
                app.duidie=1;
                app.redrawPlot();
            end
        end

        % Callback function: ChannelTree
        function ChannelTreeCheckedNodesChanged(app, event)
            checkedNodes = app.ChannelTree.CheckedNodes;
            app.CheckedData=[];
            app.CheckedData=[app.CheckedData,checkedNodes.NodeData];
            if ~ismember(app.ChannelTree.SelectedNodes.NodeData,app.CheckedData)
                app.SelectedChannel=0;
            end
          
            app.dispchans=min(60,length(checkedNodes));
            app.spacing= 580/app.dispchans;
            app.redrawPlot();
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end

            if ishandle(app.ChannelSelectionText)
                delete(app.ChannelSelectionText);
            end
        end

        % Menu selected function: FFTMenu
        function FFTMenuSelected(app, event)
            checkedNodes = app.ChannelTree.CheckedNodes;
            if isempty(checkedNodes)
                chIdx = 1:app.chans;
            else
                chIdx = arrayfun(@(n) n.NodeData, checkedNodes);
            end
            app.CheckedData = chIdx;

            if isempty(chIdx)
                uialert(app.UIFigure, '请先勾选至少一个通道。', 'No channel selected');
                return;
            end

            app.plotMultiChannelSpectrum(app.eegdata, chIdx);
        end

        % Menu selected function: LowPassMenu
        function LowPassMenuSelected(app, event)
            p = app.askFilterParams('Low Pass Filter', ...
                {'Fc (Hz)','Rs (dB)','Order N (even)'}, ...
                [30, 40, 100]);
            if isempty(p), return; end
            Fc = p(1); Rs = p(2); N = round(p(3));
            if mod(N,2) ~= 0, N = N + 1; end
            try
                b = designFIRLowpass(app.srate, Fc, Rs, N);
            catch ME
                uialert(app.UIFigure, ME.message, 'Low pass design failed'); return;
            end
            app.applyFilterAndUpdate(b, sprintf('Low Pass  Fc=%.2f Hz', Fc), Fc, N, Rs);
        end

        % Menu selected function: HighPassMenu
        function HighPassMenuSelected(app, event)
            p = app.askFilterParams('High Pass Filter', ...
                {'Fc (Hz)','Rs (dB)','Order N (even)'}, ...
                [1, 40, 100]);
            if isempty(p), return; end
            Fc = p(1); Rs = p(2); N = round(p(3));
            if mod(N,2) ~= 0, N = N + 1; end
            try
                b = designFIRHighpass(app.srate, Fc, Rs, N);
            catch ME
                uialert(app.UIFigure, ME.message, 'High pass design failed'); return;
            end
            app.applyFilterAndUpdate(b, sprintf('High Pass  Fc=%.2f Hz', Fc), Fc, N, Rs);
        end

        % Menu selected function: BandPassMenu
        function BandPassMenuSelected(app, event)
            p = app.askFilterParams('Band Pass Filter', ...
                {'Flow (Hz)','Fhigh (Hz)','Rs (dB)','Order N (even)'}, ...
                [1, 40, 40, 100]);
            if isempty(p), return; end
            Flow = p(1); Fhigh = p(2); Rs = p(3); N = round(p(4));
            if mod(N,2) ~= 0, N = N + 1; end
            try
                b = designFIRBandpass(app.srate, Flow, Fhigh, Rs, N);
            catch ME
                uialert(app.UIFigure, ME.message, 'Band pass design failed'); return;
            end
            app.applyFilterAndUpdate(b, ...
                sprintf('Band Pass  [%.2f-%.2f] Hz', Flow, Fhigh), ...
                [Flow Fhigh], N, Rs);
        end

        % Menu selected function: BandSuppressionMenu
        function BandSuppressionMenuSelected(app, event)
            p = app.askFilterParams('Band Suppression (Notch)', ...
                {'Flow (Hz)','Fhigh (Hz)','Rs (dB)','Order N (even)'}, ...
                [49, 51, 40, 200]);
            if isempty(p), return; end
            Flow = p(1); Fhigh = p(2); Rs = p(3); N = round(p(4));
            if mod(N,2) ~= 0, N = N + 1; end
            try
                b = designFIRBandstop(app.srate, Flow, Fhigh, Rs, N);
            catch ME
                uialert(app.UIFigure, ME.message, 'Band stop design failed'); return;
            end
            app.applyFilterAndUpdate(b, ...
                sprintf('Band Stop  [%.2f-%.2f] Hz', Flow, Fhigh), ...
                [Flow Fhigh], N, Rs);
        end

        % Menu selected function: FiltersMenu_2
        function FiltersMenu_2Selected(app, event)
            data=app.eegdata;
            if ~isempty(app.CheckedData)
                data = data(app.CheckedData, :);
            end
            checkedNodes = app.ChannelTree.CheckedNodes;
            app.CheckedData=[];
            app.CheckedData=[app.CheckedData,checkedNodes.NodeData];
            app.DialogApp=Filters(app,app.srate,data,app.CheckedData);

            waitfor(app.DialogApp.UIFigure);

            filename = sprintf('FilteredEEG_%03d',app.currentSaveIndex);
            app.historydatamanu=uimenu(app.HistoryMenu);
            app.historydatamanu.Text=filename;
            app.historydatamanu.UserData=app.filteredDataHistory(end);
            app.historydatamanu.MenuSelectedFcn = @(src,event) historyDataMenuSelected(app, src, event);
            app.currentSaveIndex = app.currentSaveIndex + 1;
        end

        % Key press function: UIFigure
        function UIFigureKeyPress(app, event)
            isAltPressed = any(strcmp(event.Modifier, 'alt'));

            if isAltPressed
                switch event.Key
                    case 'leftarrow'
                        app.moveTimeSelection(-0.001);

                    case 'rightarrow'
                        app.moveTimeSelection(0.001);
                end
            else
                switch event.Key
                    case 'leftarrow'
                        app.moveTimeSelection(-0.01);

                    case 'rightarrow'
                        app.moveTimeSelection(0.01);
                end
            end
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            app.ChannelTree.CheckedNodes=[];

            app.CheckedData=[];

            app.dispchans=0;
            if ishandle(app.TimeSelectionLine)
                delete(app.TimeSelectionLine);
            end

            if ishandle(app.SelectionRectangle)
                delete(app.SelectionRectangle);
            end
            if ishandle(app.ChannelSelectionText)
                delete(app.ChannelSelectionText);
            end
            app.redrawPlot();
        end

        % Button pushed function: SelectAllButton
        function SelectAllButtonPushed(app, event)
            app.ChannelTree.CheckedNodes=[];
            app.CheckedData=[];
            for i=1:length(app.ChannelNode)
                app.ChannelTree.CheckedNodes=[app.ChannelTree.CheckedNodes;app.ChannelNode{i}];
                app.CheckedData=[app.CheckedData,app.ChannelTree.CheckedNodes(i).NodeData];
            end
            app.dispchans=min(length(app.ChannelTree.CheckedNodes),60);
            disp(app.dispchans);
            app.spacing= 580/app.dispchans;
            app.redrawPlot();
        end

        % Selection changed function: ChannelTree
        function ChannelTreeSelectionChanged(app, event)
            selNodes = app.ChannelTree.SelectedNodes;
            if isempty(selNodes), return; end
            selNode  = selNodes(1);
            origIdx  = selNode.NodeData;

            if ~isempty(app.SelDebounceTimer) && isvalid(app.SelDebounceTimer)
                stop(app.SelDebounceTimer);
                delete(app.SelDebounceTimer);
            end

            app.SelDebounceTimer = timer( ...
                'ExecutionMode','singleShot', ...
                'StartDelay', 0.2, ...
                'TimerFcn', @(~,~) app.applySelectionHighlight(origIdx));

            start(app.SelDebounceTimer);
        end

        % Value changed function: WinlenthEditField
        function WinlenthEditFieldValueChanged(app, event)
            app.winlength = app.WinlenthEditField.Value;
            if app.winlength<app.frames/app.srate
                app.redrawPlot;
            else
                app.winlength=app.frames/app.srate;
                app.WinlenthEditField.Value = app.winlength;
                app.redrawPlot;
            end
        end

        % Value changed function: ScaleEditField
        function ScaleEditFieldValueChanged(app, event)
            value = app.ScaleEditField.Value;
            app.ScaleFactor = value;
            if app.ScaleFactor == 0
                app.computeAutoScaleFactor();
            end
            app.redrawPlot();
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            if isobject(app.dataNode)
                evtToSave = app.EventMarkers;
                for i = 1:numel(evtToSave)
                    evtToSave(i).latency  = evtToSave(i).latency  * app.srate + 1;
                    if isfield(evtToSave, 'duration') && ~isempty(evtToSave(i).duration)
                        evtToSave(i).duration = evtToSave(i).duration * app.srate;
                    end
                end
                app.dataNode.dataInfo.event = evtToSave;
                app.dataNode.save();
            end
            delete(app)
        end

        % Value changed function: RangeEditField
        function RangeEditFieldValueChanged(app, event)
            raw = strtrim(app.RangeEditField.Value);
            if isempty(raw), return; end

            raw = regexprep(raw, '[,\-~:]', ' ');
            nums = sscanf(raw, '%f');

            if isempty(nums)
                uialert(app.UIFigure, ...
                    'Invalid format. Use "t" for a point or "t1-t2" for a range.', ...
                    'Invalid Input');
                return;
            end

            if isscalar(nums)
                app.jumpToTime(nums(1));
            else
                a = min(nums(1), nums(2));
                b = max(nums(1), nums(2));
                app.jumpToTime(a, b);
            end
        end

        % Double-clicked callback: EventListListBox
        function EventListListBoxDoubleClicked(app, event)
            selectedItem = app.EventListListBox.Value;
            if isempty(selectedItem), return; end

            idx = find(strcmp(app.EventListListBox.Items, selectedItem), 1);
            if isempty(idx) || idx > numel(app.EventMarkers), return; end

            targetTime = app.EventMarkers(idx).latency;
            app.jumpToTime(targetTime);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % ============ 1. 主窗口 ============
            app.UIFigure = uifigure('Visible', 'off');
            % 注意: 不再设置 AutoResizeChildren='off',允许子元素跟随缩放
            app.UIFigure.Position = [100 100 1158 694];
            app.UIFigure.Name = 'MATLAB App';
            
        
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.UIFigure.KeyPressFcn = createCallbackFcn(app, @UIFigureKeyPress, true);

            % ============ 2. 菜单 ============
            app.FiltersMenu = uimenu(app.UIFigure);
            app.FiltersMenu.Text = 'Filters';

            app.FFTMenu = uimenu(app.FiltersMenu);
            app.FFTMenu.MenuSelectedFcn = createCallbackFcn(app, @FFTMenuSelected, true);
            app.FFTMenu.Text = 'FFT';

            app.LowPassMenu = uimenu(app.FiltersMenu);
            app.LowPassMenu.MenuSelectedFcn = createCallbackFcn(app, @LowPassMenuSelected, true);
            app.LowPassMenu.Text = 'Low Pass';

            app.HighPassMenu = uimenu(app.FiltersMenu);
            app.HighPassMenu.MenuSelectedFcn = createCallbackFcn(app, @HighPassMenuSelected, true);
            app.HighPassMenu.Text = 'High Pass';

            app.BandPassMenu = uimenu(app.FiltersMenu);
            app.BandPassMenu.MenuSelectedFcn = createCallbackFcn(app, @BandPassMenuSelected, true);
            app.BandPassMenu.Text = 'Band Pass';

            app.BandSuppressionMenu = uimenu(app.FiltersMenu);
            app.BandSuppressionMenu.MenuSelectedFcn = createCallbackFcn(app, @BandSuppressionMenuSelected, true);
            app.BandSuppressionMenu.Text = 'Band Suppression ';

            app.FiltersMenu_2 = uimenu(app.FiltersMenu);
            app.FiltersMenu_2.MenuSelectedFcn = createCallbackFcn(app, @FiltersMenu_2Selected, true);
            app.FiltersMenu_2.Text = 'Filters';

            % ============ 3. 顶层网格布局 ============
            % 列: [左侧通道树 | 中央绘图 | 右侧事件区]
            % 行: [主显示区  | 控制面板]
            app.MainGrid = uigridlayout(app.UIFigure);
            app.MainGrid.ColumnWidth   = {140, '1x', 165};   % 左右固定,中间弹性
            app.MainGrid.RowHeight     = {'1x', 140};        % 上方弹性,下方控制面板固定
            app.MainGrid.ColumnSpacing = 4;
            app.MainGrid.RowSpacing    = 4;
            app.MainGrid.Padding       = [4 4 4 4];

            % ============ 4. 左侧:通道面板 ============
            app.ChennalsPanel = uipanel(app.MainGrid);
            app.ChennalsPanel.Title = 'Chennals';
            app.ChennalsPanel.Layout.Row    = 1;
            app.ChennalsPanel.Layout.Column = 1;

            app.ChannelsGrid = uigridlayout(app.ChennalsPanel);
            app.ChannelsGrid.ColumnWidth = {'1x','1x'};
            app.ChannelsGrid.RowHeight   = {'1x', 26};
            app.ChannelsGrid.ColumnSpacing = 2;
            app.ChannelsGrid.RowSpacing    = 2;
            app.ChannelsGrid.Padding       = [2 2 2 2];

            app.ChannelTree = uitree(app.ChannelsGrid, 'checkbox');
            app.ChannelTree.SelectionChangedFcn   = createCallbackFcn(app, @ChannelTreeSelectionChanged, true);
            app.ChannelTree.CheckedNodesChangedFcn = createCallbackFcn(app, @ChannelTreeCheckedNodesChanged, true);
            app.ChannelTree.Layout.Row    = 1;
            app.ChannelTree.Layout.Column = [1 2];

            app.SelectAllButton = uibutton(app.ChannelsGrid, 'push');
            app.SelectAllButton.ButtonPushedFcn = createCallbackFcn(app, @SelectAllButtonPushed, true);
            app.SelectAllButton.Text = 'Select All';
            app.SelectAllButton.Layout.Row    = 2;
            app.SelectAllButton.Layout.Column = 1;

            app.ClearButton = uibutton(app.ChannelsGrid, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.Text = 'Clear';
            app.ClearButton.Layout.Row    = 2;
            app.ClearButton.Layout.Column = 2;

            % ============ 5. 中央:绘图区 + (可选)滚动条 ============
            app.CenterPanel = uipanel(app.MainGrid);
            app.CenterPanel.BorderType = 'none';
            app.CenterPanel.Layout.Row    = 1;
            app.CenterPanel.Layout.Column = 2;

            app.CenterGrid = uigridlayout(app.CenterPanel);
            app.CenterGrid.ColumnWidth = {'1x', 18};   % 主图弹性,右侧滚动条 18px
            app.CenterGrid.RowHeight   = {'1x'};
            app.CenterGrid.ColumnSpacing = 2;
            app.CenterGrid.RowSpacing    = 0;
            app.CenterGrid.Padding       = [20 2 2 20];

            app.UIAxes = uiaxes(app.CenterGrid);
            title(app.UIAxes, '')
            xlabel(app.UIAxes, '')
            ylabel(app.UIAxes, '')
            app.UIAxes.Layout.Row    = 1;
            app.UIAxes.Layout.Column = 1;

            % 滚动条在 startupFcn 里按需创建,占第 2 列

            % ============ 6. 右侧:事件列表 ============
            app.EventListPanel = uipanel(app.MainGrid);
            app.EventListPanel.Title = 'Event List';
            app.EventListPanel.Layout.Row    = 1;
            app.EventListPanel.Layout.Column = 3;

            app.EventListGrid = uigridlayout(app.EventListPanel);
            app.EventListGrid.ColumnWidth = {'1x'};
            app.EventListGrid.RowHeight   = {'1x'};
            app.EventListGrid.Padding     = [2 2 2 2];

            app.EventListListBox = uilistbox(app.EventListGrid);
            app.EventListListBox.Items = {};
            app.EventListListBox.DoubleClickedFcn = createCallbackFcn(app, @EventListListBoxDoubleClicked, true);
            app.EventListListBox.Value = {};
            app.EventListListBox.Layout.Row    = 1;
            app.EventListListBox.Layout.Column = 1;

            % ============ 7. 底部:控制面板 (跨满 3 列) ============
            app.ControlPanel = uipanel(app.MainGrid);
            app.ControlPanel.Title = 'Control';
            app.ControlPanel.Layout.Row    = 2;
            app.ControlPanel.Layout.Column = [1 2];   % 跨左+中

            % 控制面板内部网格:把控件按"列段"组织,横向自适应
            % 列布局:[Stacking | Selected/Range | WinLen/Scale | Nav按钮 | Range跳转]
            app.ControlGrid = uigridlayout(app.ControlPanel);
            app.ControlGrid.ColumnWidth = {110, '1x', '1x', 220, '1x'};
            app.ControlGrid.RowHeight   = {'1x'};
            app.ControlGrid.ColumnSpacing = 6;
            app.ControlGrid.Padding       = [6 6 6 6];

            % --- 7.1 Stacking 按钮 ---
            app.StackingmodeButton = uibutton(app.ControlGrid, 'state');
            app.StackingmodeButton.ValueChangedFcn = createCallbackFcn(app, @StackingmodeButtonValueChanged, true);
            app.StackingmodeButton.Text = 'Stacking mode';
            app.StackingmodeButton.Layout.Row    = 1;
            app.StackingmodeButton.Layout.Column = 1;

            % --- 7.2 选中时间 + 范围显示 ---
            app.SelectionGrid = uigridlayout(app.ControlGrid);
            app.SelectionGrid.ColumnWidth = {90, '1x'};
            app.SelectionGrid.RowHeight   = {22, 22, 22};
            app.SelectionGrid.ColumnSpacing = 4;
            app.SelectionGrid.RowSpacing    = 2;
            app.SelectionGrid.Padding       = [0 0 0 0];
            app.SelectionGrid.Layout.Row    = 1;
            app.SelectionGrid.Layout.Column = 2;

            app.SelectedTimeEditField_2Label_2 = uilabel(app.SelectionGrid);
            app.SelectedTimeEditField_2Label_2.Text = 'Selected Time';
            app.SelectedTimeEditField_2Label_2.HorizontalAlignment = 'right';
            app.SelectedTimeEditField_2Label_2.Layout.Row    = 1;
            app.SelectedTimeEditField_2Label_2.Layout.Column = 1;

            app.SelectedTimeEditField = uieditfield(app.SelectionGrid, 'numeric');
            app.SelectedTimeEditField.Layout.Row    = 1;
            app.SelectedTimeEditField.Layout.Column = 2;

            % "范围:start - end" 这一行用一个内嵌网格放3个label
            rngRow = uigridlayout(app.SelectionGrid);
            rngRow.ColumnWidth = {'1x', 12, '1x'};
            rngRow.RowHeight   = {'1x'};
            rngRow.ColumnSpacing = 2;
            rngRow.Padding     = [0 0 0 0];
            rngRow.Layout.Row  = 2;
            rngRow.Layout.Column = [1 2];

            app.rangestart = uilabel(rngRow);
            app.rangestart.Text = '';
            app.rangestart.HorizontalAlignment = 'right';
            app.rangestart.Layout.Row    = 1;
            app.rangestart.Layout.Column = 1;

            app.Label = uilabel(rngRow);
            app.Label.Text = '-';
            app.Label.HorizontalAlignment = 'center';
            app.Label.Layout.Row    = 1;
            app.Label.Layout.Column = 2;

            app.rangeend = uilabel(rngRow);
            app.rangeend.Text = '';
            app.rangeend.HorizontalAlignment = 'left';
            app.rangeend.Layout.Row    = 1;
            app.rangeend.Layout.Column = 3;

            % --- 7.3 窗长 / 缩放 ---
            app.WindowGrid = uigridlayout(app.ControlGrid);
            app.WindowGrid.ColumnWidth = {70, '1x'};
            app.WindowGrid.RowHeight   = {22, 22};
            app.WindowGrid.ColumnSpacing = 4;
            app.WindowGrid.RowSpacing    = 4;
            app.WindowGrid.Padding       = [0 0 0 0];
            app.WindowGrid.Layout.Row    = 1;
            app.WindowGrid.Layout.Column = 3;

            app.WinlenthEditFieldLabel = uilabel(app.WindowGrid);
            app.WinlenthEditFieldLabel.Text = 'Winlenth';
            app.WinlenthEditFieldLabel.HorizontalAlignment = 'right';
            app.WinlenthEditFieldLabel.Layout.Row    = 1;
            app.WinlenthEditFieldLabel.Layout.Column = 1;

            app.WinlenthEditField = uieditfield(app.WindowGrid, 'numeric');
            app.WinlenthEditField.ValueChangedFcn = createCallbackFcn(app, @WinlenthEditFieldValueChanged, true);
            app.WinlenthEditField.Layout.Row    = 1;
            app.WinlenthEditField.Layout.Column = 2;

            app.ScaleEditFieldLabel = uilabel(app.WindowGrid);
            app.ScaleEditFieldLabel.Text = 'Scale';
            app.ScaleEditFieldLabel.HorizontalAlignment = 'right';
            app.ScaleEditFieldLabel.Layout.Row    = 2;
            app.ScaleEditFieldLabel.Layout.Column = 1;

            app.ScaleEditField = uieditfield(app.WindowGrid, 'numeric');
            app.ScaleEditField.ValueChangedFcn = createCallbackFcn(app, @ScaleEditFieldValueChanged, true);
            app.ScaleEditField.Layout.Row    = 2;
            app.ScaleEditField.Layout.Column = 2;

            % --- 7.4 导航按钮区 (<<  <  >  >>  ;  -  +) ---
            app.NavGrid = uigridlayout(app.ControlGrid);
            app.NavGrid.ColumnWidth = {'1x','1x','1x','1x'};
            app.NavGrid.RowHeight   = {22, 22};
            app.NavGrid.ColumnSpacing = 4;
            app.NavGrid.RowSpacing    = 4;
            app.NavGrid.Padding       = [0 0 0 0];
            app.NavGrid.Layout.Row    = 1;
            app.NavGrid.Layout.Column = 4;

            app.LLButton = uibutton(app.NavGrid, 'push');
            app.LLButton.ButtonPushedFcn = createCallbackFcn(app, @LLButtonPushed, true);
            app.LLButton.Text = '<<';
            app.LLButton.Layout.Row    = 1;
            app.LLButton.Layout.Column = 1;

            app.LButton = uibutton(app.NavGrid, 'push');
            app.LButton.ButtonPushedFcn = createCallbackFcn(app, @LButtonPushed, true);
            app.LButton.Text = '<';
            app.LButton.Layout.Row    = 1;
            app.LButton.Layout.Column = 2;

            app.RButton = uibutton(app.NavGrid, 'push');
            app.RButton.ButtonPushedFcn = createCallbackFcn(app, @RButtonPushed, true);
            app.RButton.Text = '>';
            app.RButton.Layout.Row    = 1;
            app.RButton.Layout.Column = 3;

            app.RRButton = uibutton(app.NavGrid, 'push');
            app.RRButton.ButtonPushedFcn = createCallbackFcn(app, @RRButtonPushed, true);
            app.RRButton.Text = '>>';
            app.RRButton.Layout.Row    = 1;
            app.RRButton.Layout.Column = 4;

            app.biggerButton = uibutton(app.NavGrid, 'push');
            app.biggerButton.ButtonPushedFcn = createCallbackFcn(app, @biggerButtonPushed, true);
            app.biggerButton.Text = '-';
            app.biggerButton.Layout.Row    = 2;
            app.biggerButton.Layout.Column = [1 2];

            app.smallerButton = uibutton(app.NavGrid, 'push');
            app.smallerButton.ButtonPushedFcn = createCallbackFcn(app, @smallerButtonPushed, true);
            app.smallerButton.Text = '+';
            app.smallerButton.Layout.Row    = 2;
            app.smallerButton.Layout.Column = [3 4];

            % --- 7.5 Range 跳转输入框 ---
            rangeJumpGrid = uigridlayout(app.ControlGrid);
            rangeJumpGrid.ColumnWidth = {50, '1x'};
            rangeJumpGrid.RowHeight   = {22};
            rangeJumpGrid.ColumnSpacing = 4;
            rangeJumpGrid.Padding       = [0 0 0 0];
            rangeJumpGrid.Layout.Row    = 1;
            rangeJumpGrid.Layout.Column = 5;

            app.RangeEditFieldLabel = uilabel(rangeJumpGrid);
            app.RangeEditFieldLabel.Text = 'Range';
            app.RangeEditFieldLabel.HorizontalAlignment = 'right';
            app.RangeEditFieldLabel.Layout.Row    = 1;
            app.RangeEditFieldLabel.Layout.Column = 1;

            app.RangeEditField = uieditfield(rangeJumpGrid, 'text');
            app.RangeEditField.ValueChangedFcn = createCallbackFcn(app, @RangeEditFieldValueChanged, true);
            app.RangeEditField.Layout.Row    = 1;
            app.RangeEditField.Layout.Column = 2;

            % ============ 8. 右下角:Mark new event 面板 ============
            app.MarkneweventPanel = uipanel(app.MainGrid);
            app.MarkneweventPanel.Title = 'Mark new event';
            app.MarkneweventPanel.Layout.Row    = 2;
            app.MarkneweventPanel.Layout.Column = 3;

            app.MarkEventGrid = uigridlayout(app.MarkneweventPanel);
            app.MarkEventGrid.ColumnWidth = {'1x'};
            app.MarkEventGrid.RowHeight   = {22, 24, 26};
            app.MarkEventGrid.RowSpacing  = 4;
            app.MarkEventGrid.Padding     = [6 6 6 6];

            app.EventNameLabel = uilabel(app.MarkEventGrid);
            app.EventNameLabel.Text = 'Event Name';
            app.EventNameLabel.HorizontalAlignment = 'center';
            app.EventNameLabel.Layout.Row    = 1;
            app.EventNameLabel.Layout.Column = 1;

            app.EventLabelEditField = uieditfield(app.MarkEventGrid, 'text');
            app.EventLabelEditField.Value = 'Event';
            app.EventLabelEditField.HorizontalAlignment = 'center';
            app.EventLabelEditField.Layout.Row    = 2;
            app.EventLabelEditField.Layout.Column = 1;

            app.AddEventButton = uibutton(app.MarkEventGrid, 'push');
            app.AddEventButton.ButtonPushedFcn = createCallbackFcn(app, @AddEventButtonPushed, true);
            app.AddEventButton.Text = 'AddEvent';
            app.AddEventButton.Layout.Row    = 3;
            app.AddEventButton.Layout.Column = 1;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
   
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SEAL_showEEG(varargin)

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
