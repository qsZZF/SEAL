classdef SEAL_Preprocess < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        SaveMenu                        matlab.ui.container.Menu
        ResetMenu                       matlab.ui.container.Menu
        LaplacianPanel                  matlab.ui.container.Panel
        LapMEditField                   matlab.ui.control.NumericEditField
        SphericalsplineordermEditField_2Label  matlab.ui.control.Label
        LapLambdaEditField              matlab.ui.control.NumericEditField
        LambdaEditFieldLabel            matlab.ui.control.Label
        LapUseCSDCheckBox               matlab.ui.control.CheckBox
        LapEdgeDropdown                 matlab.ui.control.DropDown
        EdgehandlingLabel               matlab.ui.control.Label
        LapOKButten                     matlab.ui.control.Button
        InterpolatePanel                matlab.ui.container.Panel
        InterpRadiusEditField           matlab.ui.control.NumericEditField
        InterpRadiusEditFieldLabel      matlab.ui.control.Label
        InterpKEditField                matlab.ui.control.NumericEditField
        InterpKEditFieldLabel           matlab.ui.control.Label
        UsetemplatemontageCheckBox      matlab.ui.control.CheckBox
        InterpolateMethodDropDown       matlab.ui.control.DropDown
        MethodDropDown_2Label           matlab.ui.control.Label
        interpolateOKButton             matlab.ui.control.Button
        BadchannelsEditField            matlab.ui.control.EditField
        BadchannelsEditFieldLabel       matlab.ui.control.Label
        RereferencingPanel              matlab.ui.container.Panel
        RemovethereferencechannelCheckBox  matlab.ui.control.CheckBox
        inputexamplesPanel              matlab.ui.container.Panel
        AverageofspecificchannelgroupsLabel  matlab.ui.control.Label
        AveragereferenceLabel           matlab.ui.control.Label
        SingleelectrodereferenceLabel   matlab.ui.control.Label
        DualelectrodereferenceLabel     matlab.ui.control.Label
        orA1A2Label                     matlab.ui.control.Label
        AverageLabel                    matlab.ui.control.Label
        A1A2Label                       matlab.ui.control.Label
        A1Label                         matlab.ui.control.Label
        rereferenceOKButton             matlab.ui.control.Button
        MethodDropDown                  matlab.ui.control.DropDown
        MethodDropDownLabel             matlab.ui.control.Label
        ChooseparamsEditField           matlab.ui.control.EditField
        ChooseparamsEditFieldLabel      matlab.ui.control.Label
        TimeFrequencyDecompositionPanel  matlab.ui.container.Panel
        STFTOKButton                    matlab.ui.control.Button
        WTOKButton                      matlab.ui.control.Button
        DownsamplingFacorEditField      matlab.ui.control.NumericEditField
        DownsamplingFacorEditFieldLabel  matlab.ui.control.Label
        ResamplingLabel                 matlab.ui.control.Label
        DSOKButton                      matlab.ui.control.Button
        VoicesperOctaveEditField        matlab.ui.control.NumericEditField
        VoicesperOctaveEditFieldLabel   matlab.ui.control.Label
        HzLabel_4                       matlab.ui.control.Label
        toLabel                         matlab.ui.control.Label
        ScalesEditField_2               matlab.ui.control.NumericEditField
        ScalesEditField                 matlab.ui.control.NumericEditField
        ScalesEditFieldLabel            matlab.ui.control.Label
        WaveletFamliyDropDown           matlab.ui.control.DropDown
        WaveletFamliyDropDownLabel      matlab.ui.control.Label
        NFFTWindowLengthEditField       matlab.ui.control.NumericEditField
        NFFTWindowLengthEditFieldLabel  matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        OverlapRatioEditField           matlab.ui.control.NumericEditField
        OverlapRatioEditFieldLabel      matlab.ui.control.Label
        WindowLengthEditField           matlab.ui.control.NumericEditField
        WindowLengthEditFieldLabel      matlab.ui.control.Label
        WindowDropDown                  matlab.ui.control.DropDown
        WindowDropDownLabel             matlab.ui.control.Label
        WaveletTransformLabel           matlab.ui.control.Label
        STFTLabel                       matlab.ui.control.Label
        ICAPanel                        matlab.ui.container.Panel
        ComponentsButton                matlab.ui.control.Button
        RemoveComponengtsIndexes        matlab.ui.control.EditField
        RemoveComponentsButton          matlab.ui.control.Button
        ICAOKButton                     matlab.ui.control.Button
        ChanneltypesorchannelindicesEditField  matlab.ui.control.EditField
        ChanneltypesorchannelindicesEditFieldLabel  matlab.ui.control.Label
        CommandlineoptiionsEditField    matlab.ui.control.EditField
        CommandlineoptiionsEditFieldLabel  matlab.ui.control.Label
        ICAalgorithmDropDown            matlab.ui.control.DropDown
        ICAalgorithmDropDownLabel       matlab.ui.control.Label
        ReordercomponentsbyvarianceCheckBox  matlab.ui.control.CheckBox
        SpectrumPanel                   matlab.ui.container.Panel
        UIAxes                          matlab.ui.control.UIAxes
        SpatialFilteringPanel           matlab.ui.container.Panel
        ERPrelatedPanel                 matlab.ui.container.Panel
        ExtractTrialsLabel              matlab.ui.control.Label
        EpochingButton                  matlab.ui.control.Button
        AverageTrialsLabel              matlab.ui.control.Label
        BaselineCorrectLabel            matlab.ui.control.Label
        BadChannelDetectionLabel        matlab.ui.control.Label
        AveragetrialsSetButton_3        matlab.ui.control.Button
        Baselinesetbutten               matlab.ui.control.Button
        ExtracttraialsSetButton         matlab.ui.control.Button
        MethodsButtonGroup              matlab.ui.container.ButtonGroup
        LaplacianButton                 matlab.ui.control.RadioButton
        InterpolateButton               matlab.ui.control.RadioButton
        RereferenceButton               matlab.ui.control.RadioButton
        ICAButton                       matlab.ui.control.RadioButton
        FrequencydomainFilteringPanel   matlab.ui.container.Panel
        FrequencydomainOKButton         matlab.ui.control.Button
        frequencydominSetbutten         matlab.ui.control.Button
        NotchEditField                  matlab.ui.control.NumericEditField
        NotchEditFieldLabel             matlab.ui.control.Label
        LowpassEditField                matlab.ui.control.NumericEditField
        LowpassEditFieldLabel           matlab.ui.control.Label
        HighpassEditField               matlab.ui.control.NumericEditField
        HighpassEditFieldLabel          matlab.ui.control.Label
        HzLabel_3                       matlab.ui.control.Label
        HzLabel_2                       matlab.ui.control.Label
        HzLabel                         matlab.ui.control.Label
        NotchCheckBox                   matlab.ui.control.CheckBox
        LowpassCheckBox                 matlab.ui.control.CheckBox
        HighpassCheckBox                matlab.ui.control.CheckBox
    end


    properties (Access = private)
        dataNode
        protocolNode
        sessionNode
        defaultpath
        rawdata
        ModifiedData
        srate
        Type
        channelCount
        SaveIndex=1  %保存数据索引
        Callingapp
        icaresult
        chanlocs
        chanLabels
        change
        TF_Result = []     
        NoiseCovMat = [] %噪声协方差矩阵
    end

    methods (Access = private)

        function plotMultiChannelSpectrum(app, EEGdata)
            % 绘制多通道功率谱
            ax = app.UIAxes;
            cla(ax); hold(ax, 'on');
            title(ax, 'Multi-Channel Power Spectrum');
            xlabel(ax, 'Frequency (Hz)');
            ylabel(ax, 'Power (dB \muV^2/Hz)'); % 单位：分贝(微伏平方/赫兹)
            grid(ax, 'on');

            if app.isEpoched(EEGdata)
                dataToPlot = double(EEGdata(:, :, 1));  % 取第一个试次预览
            else
                dataToPlot = double(EEGdata);
            end

            [~, nT] = size(dataToPlot);

            % 2. 预处理：去直流并确保偶数点
            dataToPlot = dataToPlot - mean(dataToPlot, 2);
            if mod(nT, 2) ~= 0
                dataToPlot = dataToPlot(:, 1:end-1);
                nT = nT - 1;
            end

            % 3. 向量化加窗 (利用隐式扩展，瞬间完成所有通道)
            win = hann(nT)'; % 1 x nT 的行向量
            windowedData = dataToPlot .* win; % nCh x nT

            % 4. 并行 FFT 计算 (沿第二个维度 '2'，即时间维度)
            fftResult = fft(windowedData, [], 2);

            % 5. 提取单边谱并应用正确的 PSD 缩放公式
            halfN = floor(nT/2);
            xdft = fftResult(:, 1:halfN+1);

            % 核心数学公式: (1 / (Fs * sum(w^2))) * |X(f)|^2
            psdMatrix = (1 / (app.srate * sum(win.^2))) * abs(xdft).^2;

            % 只有非直流和非奈奎斯特频率的能量需要乘以 2
            psdMatrix(:, 2:end-1) = 2 * psdMatrix(:, 2:end-1);

            % 生成频率轴
            freq = app.srate * (0:halfN) / nT;

            % 6. 限制频率范围 (0.01 到 300 Hz)
            validRange = freq >= 0.01 & freq <= 300;
            freq = freq(validRange);
            psdMatrix = psdMatrix(:, validRange);
            unit = app.dataNode.dataInfo.unit;
            switch unit
                case "V"
                    scale = 1e12;    % V² → μV²
                case "mV"
                    scale = 1e6;     % mV² → μV²
                case "μV"
                    scale = 1;       % 已经是 μV²
                case "nV"
                    scale = 1e-6;    % nV² → μV²
                otherwise
                    warning('SEAL:UnknownUnit', '未知单位 %s，按 μV 处理', unit);
                    scale = 1;
            end
            % 7. 单位转换：假定原数据为伏特(V)，乘 1e12 转为微伏(μV^2)
            % 加上 eps 彻底杜绝 log10(0) 报错
            psd_dB = 10 * log10(psdMatrix * scale + eps);

            % 8. 一键绘制所有通道 (需转置矩阵，让每一列代表一个通道画一条线)
            plot(ax, freq, psd_dB', 'LineWidth', 1.5);

            axis(ax, 'tight');
            hold(ax, 'off');
        end

        function closeall(app)
            app.ICAPanel.Visible='off';
            app.TimeFrequencyDecompositionPanel.Visible='off';
            app.RereferencingPanel.Visible='off';
            app.InterpolatePanel.Visible='off';
            app.LaplacianPanel.Visible='off';
        end

        function ica = runICA_FromGUI(app, data3d, alg, cmdopt_str, doReorder, chan_text)%需要EEGLAB中的runica和bunica函数
            % data3d: 通道 × 时间 × 试次

            disp('[STEP 1/6] 通道选择');
            nCh = size(data3d,1);
            chanind = app.parseChannelText(chan_text, nCh);
            disp(['  -> 已选择通道数：', num2str(numel(chanind))]);

            disp('[STEP 2/6] 展开为二维矩阵');
            X = reshape(data3d(chanind,:,:), numel(chanind), []);

            disp('[STEP 3/6] 估计数值秩');
            rankX = app.computeNumericRank(X, []);
            disp(['  -> rank = ', num2str(rankX)]);

            disp('[STEP 4/6] 解析命令行参数');
            extraArgs = app.parseNameValueString(cmdopt_str);
            disp(['  -> 附加参数个数：', num2str(numel(extraArgs))]);

            disp('[STEP 5/6] 运行 ICA');
            tic;
            ica = seal_runICA_OnMatrix(X, alg, rankX, extraArgs);
            ica.alg = alg;
            ica.chanind = chanind;
            ica.rank = rankX;
            toc;

            if doReorder
                disp('[STEP 6/6] 按方差重排');
                ica = seal_reorderIC_by_variance_matrix(ica);
                disp('  -> 已重排');
            else
                disp('[STEP 6/6] 跳过重排');
            end

            % 重塑到三维形式
            nT = size(data3d,2);
            nTr = size(data3d,3);
            nIC = size(ica.S,1);
            ica.S3d = reshape(ica.S, nIC, nT, nTr);

            disp('=== ICA 完成 ===');
        end

        function chanind = parseChannelText(app, chan_text, nMax)
            % 空字符串 => 全部通道
            if isempty(chan_text)
                chanind = 1:nMax; return;
            end
            chan_text = string(chan_text);
            parts = strtrim(split(chan_text, ','));
            idx = [];
            for k=1:numel(parts)
                tk = strtrim(parts{k});
                if contains(tk,'-')
                    ab = str2double(split(tk,'-'));
                    if numel(ab)==2 && all(isfinite(ab))
                        idx = [idx, ab(1):ab(2)]; 
                    end
                else
                    v = str2double(tk);
                    if isfinite(v), idx = [idx, v]; end 
                end
            end
            idx = unique(idx);
            chanind = idx(idx>=1 & idx<=nMax);
        end

        function idx = parseIndexString(app, s, nMax)
            parts = strtrim(split(s,','));
            idx = [];
            for k=1:numel(parts)
                tk = parts{k};
                if contains(tk,'-')
                    ab = str2double(split(tk,'-'));
                    if numel(ab)==2 && all(isfinite(ab))
                        idx = [idx, ab(1):ab(2)]; 
                    end
                else
                    v = str2double(tk);
                    if isfinite(v), idx = [idx, v]; end 
                end
            end
            idx = unique(idx);
            idx = idx(idx>=1 & idx<=nMax);
        end

        function r = computeNumericRank(app, X, tol)
            if nargin<2 || isempty(tol), tol = 1e-7; end
            [~,S,~] = svd(X, 'econ');
            s = diag(S);
            r = sum(s > tol * max(s));
        end

        function C = parseNameValueString(app, s)
            C = {};
            s = strtrim(s);
            if isempty(s), return; end
            parts = app.splitCommasRespectingQuotes(s);
            for k=1:numel(parts)
                tk = strtrim(parts{k});
                if (startsWith(tk,'''') && endsWith(tk,''''))
                    C{end+1} = tk(2:end-1); 
                elseif (startsWith(tk,'"') && endsWith(tk,'"'))
                    C{end+1} = tk(2:end-1); 
                else
                    val = str2double(tk);
                    if ~isnan(val)
                        C{end+1} = val; 
                    else
                        C{end+1} = tk; 
                    end
                end
            end
        end

        function parts = splitCommasRespectingQuotes(app, s)
            parts = {};
            buf = '';
            inq = false; qch = '';
            for i=1:numel(s)
                c = s(i);
                if ~inq && (c=='''' || c=='"')
                    inq = true; qch = c; buf = [buf c]; 
                elseif inq && c==qch
                    inq = false; buf = [buf c]; 
                elseif ~inq && c==','
                    parts{end+1} = buf; 
                    buf = '';
                else
                    buf = [buf c]; 
                end
            end
            if ~isempty(buf), parts{end+1} = buf; end
        end


        function wname = mapWaveletName(app,famUI)
            fam = lower(string(famUI));
            switch fam
                case {'morl','morlet'}   % GUI 里的 'morl' 映射为 MATLAB 的 'amor'
                    wname = 'amor';      % analytic Morlet
                case 'amor'
                    wname = 'amor';
                case 'morse'
                    wname = 'morse';
                case 'bump'
                    wname = 'bump';
                otherwise
                    warning('未知 Wavelet Family: %s，改用 amor', fam);
                    wname = 'amor';
            end
        end

        function verifyICACloseReq(app, fig, clean_data3d, chanind, bad_comps)
            % 当用户尝试关闭预览窗口时触发此函数
            
            % 弹出 EEGLAB 风格的确认框
            selection = uiconfirm(fig, ...
                sprintf('您确定要剔除这 %d 个成分吗？\n(这将会覆盖当前工作区的数据)', length(bad_comps)), ...
                '确认应用 ICA 清洗', ...
                'Options', {'Yes, apply changes', 'No, cancel'}, ...
                'DefaultOption', 1, 'CancelOption', 2);
                
            if strcmp(selection, 'Yes, apply changes')
                % === 用户选择“是”：正式更新数据 ===
                
                app.ModifiedData(chanind, :, :) = clean_data3d;
                app.change = 1;
                
                % 刷新主界面的频谱图 (利用我们之前优化好的方法)
                app.plotMultiChannelSpectrum(app.ModifiedData);
                
                msg = sprintf('成功剔除了 %d 个成分。', length(bad_comps));
                disp(['=== ', msg, ' ===']);
                uialert(app.UIFigure, msg, '清理完成');
            else
                % === 用户选择“否”：什么也不做 ===
                disp('=== [ICA] 用户取消了成分剔除操作，数据未更改。 ===');
            end
            
            % 最后，销毁预览对比图窗
            delete(fig);
        end
   
        function showSTFTPreview(app, result)
            %SHOWSTFTPREVIEW 弹出独立窗口显示 STFT 时频图（第一通道第一 trial）

            fig = figure('Name', 'STFT Preview', 'Color', 'w', ...
                'Position', [200, 200, 800, 500]);
            ax = axes(fig);

            % 推荐改法:
            S_psd = abs(result.S(:, :, 1, 1)).^2;
            S_psd(2:end-1, :) = 2 * S_psd(2:end-1, :);   % 单边谱补偿

            % 单位换算(假设数据是 V,想显示 μV²/Hz)
            unit = app.dataNode.dataInfo.unit;
            switch unit
                case "V",  scale = 1e12;
                case "mV", scale = 1e6;
                case "μV", scale = 1;
                case "nV", scale = 1e-6;
                otherwise, scale = 1;
            end
            S_dB = 10 * log10(S_psd * scale + eps);

            imagesc(ax, result.T, result.F, S_dB);
            axis(ax, 'xy');
            colormap(ax, 'jet');
            hBar = colorbar(ax);
            ylabel(hBar, 'PSD (dB μV²/Hz)');

            xlabel(ax, 'Time (s)');
            ylabel(ax, 'Frequency (Hz)');
            title(ax, sprintf('STFT Spectrogram | Channel 1, Trial 1 | Window: %s, Length: %d', ...
                result.params.window, result.params.winLen));
        end

        function showCWTPreview(app, result)
            %SHOWCWTPREVIEW 显示 CWT 时频图(频段归一化后的相对功率)

            fig = figure('Name', 'CWT Preview', 'Color', 'w', ...
                'Position', [200, 200, 900, 500]);
            ax = axes(fig);

            % ===== 时间轴 =====
            nTpts = size(result.WT, 2);
            if isfield(result, 'T') && numel(result.T) == nTpts
                tSec = result.T;
            else
                tSec = (0:nTpts - 1) / result.params.srate;
            end
            if isduration(tSec), tSec = seconds(tSec); end
            tSec = tSec(:).';

            F = result.F(:);

            % ===== 第一通道第一 trial =====
            S_power = abs(result.WT(:, :, 1, 1)).^2;

            % ===== 频段归一化(每个频率减去自己的时间平均) =====
            S_baseline = mean(S_power, 2);
            S_relative = S_power ./ (S_baseline + eps);
            S_dB = 10 * log10(S_relative + eps);

            % ===== 绘图 =====
            imagesc(ax, tSec, F, S_dB);
            set(ax, 'YDir', 'normal');
            set(ax, 'YScale', 'log');

            colormap(ax, 'turbo');
            hBar = colorbar(ax);
            ylabel(hBar, 'Power (dB, vs frequency mean)');

            clim(ax, [-10, 10]);   % 对称范围,绿色为平均水平

            % ===== COI 遮罩(若有) =====
            if isfield(result, 'COI') && ~isempty(result.COI)
                coi = result.COI;
                if isduration(coi), coi = seconds(coi); end
                coi = coi(:);

                if numel(coi) == numel(tSec)
                    hold(ax, 'on');
                    fLow = max(min(F), 1e-10);
                    xPatch = [tSec(:);  flipud(tSec(:))];
                    yPatch = [coi(:);   fLow * ones(numel(coi), 1)];
                    valid = isfinite(yPatch);
                    patch(ax, xPatch(valid), yPatch(valid), [0.3 0.3 0.3], ...
                        'FaceAlpha', 0.45, 'EdgeColor', 'none');
                    plot(ax, tSec, coi, 'w--', 'LineWidth', 1.2);
                    hold(ax, 'off');
                end
            end

            xlabel(ax, 'Time (s)');
            ylabel(ax, 'Frequency (Hz)');
            title(ax, sprintf(['CWT Scalogram (frequency-normalized) | Ch1, Trial1\n' ...
                'Wavelet: %s, VPO: %d'], ...
                result.params.wavelet, result.params.voicesPerOctave));

            ylim(ax, [min(F), max(F)]);
            xlim(ax, [tSec(1), tSec(end)]);
            grid(ax, 'on');
            ax.Layer = 'top';
        end

        function tf = isEpoched(~, data)
            % 判断数据是否为分段(epoched)格式
            % 2D 数据或第三维为1的数据都视为连续数据
            tf = size(data, 3) > 1;
        end

        function data2d = ensureContinuous(~, data)
            % 如果是 3D 分段数据，沿第三维拼接为 2D
            % 如果已经是 2D，原样返回
            if size(data, 3) > 1
                [nCh, nT, nTr] = size(data);
                data2d = reshape(data, nCh, nT * nTr);
            else
                data2d = data;
            end
        end
        function saveTrialsSeparately(app, data3d, parentPath, baseName)
            % 按 trial 把 3D 数据拆成多个 2D 单文件
            % data3d:    channels × time × trials
            % parentPath: 保存目录
            % baseName:  主文件基础名(不含扩展名)

            [nCh, nT, nTr] = size(data3d);

            % 创建子文件夹
            trialDir = fullfile(parentPath, [baseName, '_trials']);
            if ~isfolder(trialDir)
                mkdir(trialDir);
            end

            % 预取公共元信息(避免循环里反复拷贝)
            sratee    = app.srate;
            chanlocss = [];
            unitStr  = '';
            fullEvent = [];

            if ~isempty(app.dataNode.dataInfo)
                if isfield(app.dataNode.dataInfo, 'chanlocs')
                    chanlocss = app.dataNode.dataInfo.chanlocs;
                end
                if isfield(app.dataNode.dataInfo, 'unit')
                    unitStr = app.dataNode.dataInfo.unit;
                end
                if isfield(app.dataNode.dataInfo, 'event')
                    fullEvent = app.dataNode.dataInfo.event;
                end
            end

            % 进度条
            hWait = uiprogressdlg(app.UIFigure, ...
                'Title', '保存单 trial 文件', ...
                'Message', sprintf('共 %d 个 trial...', nTr), ...
                'Cancelable', 'on');

            % 数位数用于文件名补零(100 个 trial → 3 位,1000 个 → 4 位)
            nDigits = max(3, numel(num2str(nTr)));
            fmt = sprintf('trial_%%0%dd.mat', nDigits);

            try
                for tr = 1:nTr
                    if hWait.CancelRequested
                        break;
                    end

                    trialStruct = struct();
                    trialStruct.data       = data3d(:, :, tr);   % 2D: ch × t
                    trialStruct.srate      = sratee;
                    trialStruct.trialIndex = tr;
                    trialStruct.nChannels  = nCh;
                    trialStruct.nSamples   = nT;

                    if ~isempty(chanlocss)
                        trialStruct.chanlocs = chanlocss;
                    end
                    if ~isempty(unitStr)
                        trialStruct.unit = unitStr;
                    end

                    % 筛出属于这个 trial 的事件(如果 event 有 epoch 字段)
                    if ~isempty(fullEvent) && isstruct(fullEvent) && isfield(fullEvent, 'epoch')
                        mask = arrayfun(@(e) isequal(e.epoch, tr), fullEvent);
                        trialStruct.event = fullEvent(mask);
                    end

                    trialFile = fullfile(trialDir, sprintf(fmt, tr));
                    saveData(trialFile, trialStruct);
            
            hWait.Value   = tr / nTr;
            hWait.Message = sprintf('已保存 %d / %d', tr, nTr);
        end
    catch ME
        close(hWait);
        rethrow(ME);
    end
    
    close(hWait);
    
    fprintf('[Save] 已将 %d 个 trial 分别保存到:\n  %s\n', nTr, trialDir);
        end

    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.Callingapp = mainapp;
            app.UIFigure.Name='Preprocess';
            app.UIFigure.Resize = 'off';   
                %% 调整app窗口位置
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
            app.TimeFrequencyDecompositionPanel.Visible='off';
            app.HighpassEditField.Enable='off';
            app.LowpassEditField.Enable='off';
            app.NotchEditField.Enable='off';
            app.dataNode = app.Callingapp.FileTree.SelectedNodes.NodeData;
            app.sessionNode = app.dataNode.parent;
            app.protocolNode = app.dataNode.parent.parent;
            app.defaultpath = app.protocolNode.path;
            rawDataTemp = app.dataNode.data;
            
            if ismatrix(rawDataTemp) || ndims(rawDataTemp) == 3
                % 2D 连续数据保持 2D，3D 分段数据保持 3D，不做强制变形
                app.rawdata = rawDataTemp;
            else
                uialert(app.UIFigure, '不支持的数据维度，仅支持 2D 连续或 3D 分段数据。', '数据错误');
                return;
            end
            app.ModifiedData = app.rawdata;
            app.Type = app.dataNode.type;
            app.srate = double(app.dataNode.srate);
            app.channelCount = size(app.rawdata,1);
            app.plotMultiChannelSpectrum(app.rawdata)
            app.change = 0;

            if ~isempty(app.dataNode.dataInfo.chanlocs)
                app.chanlocs = app.dataNode.dataInfo.chanlocs;
                if isfield(app.chanlocs, 'labels')
                    app.chanLabels = {app.chanlocs.labels};  % 打包成 cell 数组
                end
            end
            app.ComponentsButton.Enable = 'off';
            app.RemoveComponengtsIndexes.Enable = 'off';
            app.RemoveComponentsButton.Enable = 'off';
            %默认参数设置
            defaultWinSec = 0.25;  % 默认 0.5 秒窗
            defaultWinLen = round(app.srate * defaultWinSec);
            app.WindowLengthEditField.Value = defaultWinLen;
            app.NFFTWindowLengthEditField.Value = 2^nextpow2(defaultWinLen);


        end

        % Button pushed function: frequencydominSetbutten
        function frequencydominSetbuttenPushed(app, event)
            app.closeall();
            app.TimeFrequencyDecompositionPanel.Visible='on';
            app.TimeFrequencyDecompositionPanel.Position=[628 1 439 531];
        end

        % Selection changed function: MethodsButtonGroup
        function MethodsButtonGroupSelectionChanged(app, event)
            selectedButton = app.MethodsButtonGroup.SelectedObject;
            app.closeall();
            if selectedButton==app.ICAButton
                app.ICAPanel.Visible='on';
                app.ICAPanel.Position=[628 312 438 221];
            end
            if  selectedButton==app.RereferenceButton
                app.RereferencingPanel.Visible='on';
                app.RereferencingPanel.Position=[628 1 439 531];
            end
            if  selectedButton==app.InterpolateButton
                app.InterpolatePanel.Visible='on';
                app.InterpolatePanel.Position=[628 1 439 531];
            end
            if  selectedButton==app.LaplacianButton
                app.LaplacianPanel.Visible='on';
                app.LaplacianPanel.Position=[628 1 439 531];
            end
        end

        % Value changed function: HighpassCheckBox
        function HighpassCheckBoxValueChanged(app, event)
            value = app.HighpassCheckBox.Value;
            if value==1
                app.HighpassEditField.Enable='on';
            else
                app.HighpassEditField.Enable='off';
            end
        end

        % Value changed function: LowpassCheckBox
        function LowpassCheckBoxValueChanged(app, event)
            value = app.LowpassCheckBox.Value;
            if value==1
                app.LowpassEditField.Enable='on';
            else
                app.LowpassEditField.Enable='off';
            end

        end

        % Value changed function: NotchCheckBox
        function NotchCheckBoxValueChanged(app, event)
            value = app.NotchCheckBox.Value;
            if value==1
                app.NotchEditField.Enable='on';
            else
                app.NotchEditField.Enable='off';
            end
        end

        % Button pushed function: FrequencydomainOKButton
        function FrequencydomainOKButtonPushed(app, event)
            try
                % 1. 初始化参数结构体
                filterParams = struct();
                
                % 2. 读取并校验 GUI 参数
                if app.HighpassCheckBox.Value == 1
                    hp_val = app.HighpassEditField.Value;
                    if hp_val <= 0
                        uialert(app.UIFigure, 'The cut-off frequency must be greater than 0.', 'Invalid Parameter');
                        return;
                    end
                    filterParams.highpass = hp_val;
                end
                
                if app.LowpassCheckBox.Value == 1
                    lp_val = app.LowpassEditField.Value;
                    if lp_val <= 0
                        uialert(app.UIFigure, 'The cut-off frequency must be greater than 0.', 'Invalid Parameter');
                        return;
                    end
                    filterParams.lowpass = lp_val;
                end
                
                if app.NotchCheckBox.Value == 1
                    filterParams.notch = app.NotchEditField.Value;
                end
                
                % 3. 交叉校验 (高通必须小于低通)
                if isfield(filterParams, 'highpass') && isfield(filterParams, 'lowpass')
                    if filterParams.highpass >= filterParams.lowpass
                        uialert(app.UIFigure, 'Highpass frequency must be strictly less than lowpass frequency.', 'Logic Error');
                        return;
                    end
                end
                
                % 确保至少选了一个
                if isempty(fieldnames(filterParams))
                    uialert(app.UIFigure, 'You did not choose any filter method.', 'No Method Selected');
                    return;
                end
                
                % 4. 调用外部核心算法进行极速滤波
                % 注意：这里输入的是 ModifiedData，保证预处理的流水线是连续的
                app.UIFigure.Pointer = 'watch'; drawnow; % 鼠标变沙漏，防假死
                
                if app.isEpoched(app.ModifiedData)
                    % 逐试次滤波，避免拼接点的边界效应
                    nTr = size(app.ModifiedData, 3);
                    for tr = 1:nTr
                        app.ModifiedData(:,:,tr) = seal_filter_data( ...
                            app.ModifiedData(:,:,tr), app.srate, filterParams);
                    end
                else
                    app.ModifiedData = seal_filter_data( ...
                        app.ModifiedData, app.srate, filterParams);
                end

                app.change = 1;
                app.plotMultiChannelSpectrum(app.ModifiedData);
                app.UIFigure.Pointer = 'arrow';

            catch ME
                app.UIFigure.Pointer = 'arrow';
                if strcmp(ME.identifier, 'SEAL:NoFilter')
                    uialert(app.UIFigure, ME.message, '参数错误');
                else
                    disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                    uialert(app.UIFigure, '滤波过程中发生意外错误，请查看命令行。', '运行出错');
                end
            end
            
            

        end

        % Button pushed function: rereferenceOKButton
        function rereferenceOKButtonPushed(app, event)
            try
                disp('=== [Re-reference] START ===');

                % ✅ 从流水线的当前状态读取
                tempdata = app.ModifiedData;

                % 读取 GUI
                refText = string(app.ChooseparamsEditField.Value);
                method  = string(app.MethodDropDown.Value);
                removeRef = app.RemovethereferencechannelCheckBox.Value;   % ← 新增

                if isempty(tempdata) || ~isnumeric(tempdata)
                    uialert(app.UIFigure, '数据为空或类型不合法', '错误');
                    return;
                end
                % 通道标签兜底
                [nCh, ~, ~] = size(tempdata);
                if isempty(app.chanLabels)
                    app.chanLabels = arrayfun(@(x) sprintf('Ch%d', x), 1:nCh, ...
                        'UniformOutput', false);
                    disp('[Info] chanLabels 为空，已自动生成默认通道名 Ch1...ChN');
                end

                % 解析参考参数
                refSpec = seal_parseRefSpec(refText, nCh, app.chanLabels);

                opts = struct();
                opts.removeRefChans = removeRef;

                % 执行重参考
                app.UIFigure.Pointer = 'watch'; drawnow;
                tic;
                [data_new, removedIdx, info] = seal_applyRereference( ...
                    tempdata, refSpec, method, opts);
                elapsed = toc;

                app.ModifiedData = data_new;
                app.change = 1;

                if ~isempty(removedIdx)
                    keepMask = true(nCh, 1);
                    keepMask(removedIdx) = false;

                    % 更新通道标签
                    if ~isempty(app.chanLabels)
                        app.chanLabels = app.chanLabels(keepMask);
                    end

                    % 更新通道位置
                    if ~isempty(app.chanlocs) && length(app.chanlocs) == nCh
                        app.chanlocs = app.chanlocs(keepMask);
                    end

                    fprintf('[Re-ref] 数据从 %d 通道 → %d 通道(剔除了参考通道)\n', ...
                        info.nChIn, info.nChOut);
                end


                app.plotMultiChannelSpectrum(app.ModifiedData);

                app.UIFigure.Pointer = 'arrow';

                if isempty(removedIdx)
                    msg = sprintf('重参考完成:Method=%s | 类型=%s | 通道数保持 %d | 用时 %.2fs', ...
                        method, refSpec.type, info.nChOut, elapsed);
                else
                    msg = sprintf(['重参考完成:Method=%s | 类型=%s\n' ...
                        '剔除 %d 个参考通道,剩余 %d 通道 | 用时 %.2fs'], ...
                        method, refSpec.type, numel(removedIdx), info.nChOut, elapsed);
                end
                disp(msg);
                uialert(app.UIFigure, msg, '完成');

                disp('=== [Re-reference] DONE ===');

            catch ME
                app.UIFigure.Pointer = 'arrow';
                disp('[ERR] Re-reference 出错：');
                disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), '运行出错');
            end
        end

        % Button pushed function: ICAOKButton
        function ICAOKButtonPushed(app, event)
            try
                disp('=== [ICA] 启动 ===');
                % 读取 GUI 值
                alg        = app.ICAalgorithmDropDown.Value;      % 'runica'|'binica'|'fastica'
                cmdopt_str = strtrim(app.CommandlineoptiionsEditField.Value);   % 例：'extended',1
                doReorder  = app.ReordercomponentsbyvarianceCheckBox.Value;         % true/false
                if doReorder==0  %for dispaly
                    doReorder='false';
                else
                    doReorder='true';
                end
                chan_text  = strtrim(app.ChanneltypesorchannelindicesEditField.Value);     % 'EEG' or '1-64'
                disp(['[GUI] Algorithm: ', alg]);
                disp(['[GUI] Cmd options: ', cmdopt_str]);
                disp(['[GUI] Reorder by variance: ', doReorder]);
                disp(['[GUI] Channel text: ', chan_text]);

                % 打印基本信息
                [nCh, nT, nTr] = size(app.ModifiedData);
                fprintf('[DATA] chans=%d, time=%d, trials=%d\n', nCh, nT, nTr);
                tic;
                app.icaresult = app.runICA_FromGUI(app.ModifiedData, alg, cmdopt_str, doReorder, chan_text);
                t = toc;

                % 提示完成
                msg = sprintf('ICA 完成：算法=%s | IC=%d | 用时 %.2fs', ...
                    app.icaresult.alg, size(app.icaresult.S,1), t);
                disp(msg);
                uialert(app.UIFigure, msg, '完成');
                app.change = 1;

                app.ComponentsButton.Enable = 'on';
                app.RemoveComponengtsIndexes.Enable = 'on';
                app.RemoveComponentsButton.Enable = 'on';

            catch ME
                % 友好错误提示
                disp(getReport(ME,'extended','hyperlinks','off'));

                uialert(app.UIFigure, getReport(ME,'basic','hyperlinks','off'), '运行出错');
            end


        end

        % Button pushed function: STFTOKButton
        function STFTOKButtonPushed(app, event)
            try
                disp('=== [STFT] 开始读取参数 ===');

                % ✅ 基于流水线当前数据
                if isempty(app.ModifiedData)
                    uialert(app.UIFigure, '没有可用数据。', '错误');
                    return;
                end

                % 收集 UI 参数
                winName   = string(app.WindowDropDown.Value);
                winLenInp = str2double(string(app.WindowLengthEditField.Value));
                ovlpPct   = str2double(string(app.OverlapRatioEditField.Value));
                nfftInp   = str2double(string(app.NFFTWindowLengthEditField.Value));

                % 默认值
                if isempty(winLenInp) || isnan(winLenInp), winLenInp = 128; end
                if isempty(ovlpPct)   || isnan(ovlpPct),   ovlpPct   = 50;  end
                if isempty(nfftInp)   || isnan(nfftInp),   nfftInp   = 256; end

                % 参数合理性校验
                if nfftInp < winLenInp
                    uialert(app.UIFigure, ...
                        sprintf('NFFT (%d) 必须不小于窗长 (%d)', nfftInp, winLenInp), ...
                        '参数错误');
                    return;
                end
                if ovlpPct < 0 || ovlpPct >= 100
                    uialert(app.UIFigure, ...
                        '重叠比例必须在 [0, 100) 范围内', '参数错误');
                    return;
                end

                disp('=== [STFT] 调用核心计算引擎 ===');
                app.UIFigure.Pointer = 'watch'; drawnow;

                % ✅ 使用 ModifiedData
                tic;
                [S, F, T] = seal_compute_stft(app.ModifiedData, app.srate, ...
                    winName, winLenInp, ovlpPct, nfftInp);
                elapsed = toc;

                app.UIFigure.Pointer = 'arrow';

                % ✅ 存储结果，带参数记录便于追溯
                stftResult = struct();
                stftResult.S = S;
                stftResult.F = F;
                stftResult.T = T;
                stftResult.params = struct(...
                    'window', winName, ...
                    'winLen', winLenInp, ...
                    'overlap', ovlpPct, ...
                    'nfft', nfftInp, ...
                    'srate', app.srate);
                stftResult.timestamp = datetime('now');

                app.TF_Result = stftResult;


                % ✅ 可视化：展示第一个通道第一个 trial 的时频图
                app.showSTFTPreview(stftResult);

                msg = sprintf('STFT 计算完成！耗时: %.2fs\n维度: Freq=%d × Time=%d × Ch=%d × Trials=%d', ...
                    elapsed, size(S,1), size(S,2), size(S,3), size(S,4));
                disp(msg);
                uialert(app.UIFigure, msg, '成功');

            catch ME
                app.UIFigure.Pointer = 'arrow';
                disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                uialert(app.UIFigure, ME.message, 'STFT 参数或计算错误');
            end
        end

        % Button pushed function: WTOKButton
        function WTOKButtonPushed(app, event)
            if ~license('test', 'Wavelet_Toolbox') || ~exist('cwt', 'file')
                uialert(app.UIFigure, ...
                    sprintf(['此功能需要 MATLAB Wavelet Toolbox。\n\n', ...
                    '当前环境中未检测到该工具箱。\n', ...
                    '请安装 Wavelet Toolbox 后重试，或使用 STFT 作为替代。']), ...
                    '缺少工具箱', 'Icon', 'warning');
                return;
            end

            try
                disp('=== [CWT] START ===');

                % ✅ 用流水线当前数据
                if isempty(app.ModifiedData)
                    uialert(app.UIFigure, '没有可用数据。', '错误');
                    return;
                end
                tempdata = app.ModifiedData;

                % ---- 读取 GUI ----
                famUI = string(app.WaveletFamliyDropDown.Value);
                fLoIn = str2double(string(app.ScalesEditField.Value));
                fHiIn = str2double(string(app.ScalesEditField_2.Value));
                vpoIn = str2double(string(app.VoicesperOctaveEditField.Value));

                % ---- 基本检查 ----
                if ismatrix(tempdata)
                    [nCh, nT] = size(tempdata);
                    tempdata = reshape(tempdata, [nCh, nT, 1]);
                elseif ndims(tempdata) ~= 3
                    uialert(app.UIFigure, '数据维度不支持', '错误');
                    return;
                end

                if isempty(app.srate) || ~isscalar(app.srate)
                    uialert(app.UIFigure, '缺少采样率信息', '错误');
                    return;
                end

                fs = double(app.srate);
                [nCh, nT, nTr] = size(tempdata);
                fprintf('[DATA] ch=%d, time=%d, trials=%d, fs=%.3f Hz\n', nCh, nT, nTr, fs);

                % ---- 规范参数 ----
                wname = app.mapWaveletName(famUI);
                if isempty(vpoIn) || ~isfinite(vpoIn) || vpoIn <= 0
                    vpo = 32;
                else
                    vpo = round(vpoIn);
                end

                % 自动频率范围
                [fminAuto, fmaxAuto] = cwtfreqbounds(nT, fs, 'Wavelet', wname, 'VoicesPerOctave', vpo);
                if isempty(fLoIn) || ~isfinite(fLoIn) || fLoIn <= 0
                    fLo = fminAuto;
                else
                    fLo = max(fLoIn, fminAuto);   % 不能低于理论最小值
                end
                if isempty(fHiIn) || ~isfinite(fHiIn) || fHiIn <= 0
                    fHi = fmaxAuto;
                else
                    fHi = min(fHiIn, fmaxAuto);   % 不能高于理论最大值
                end
                if fLo >= fHi
                    warning('频率区间无效，改用自动区间。');
                    fLo = fminAuto;
                    fHi = fmaxAuto;
                end
                fprintf('[CWT] wavelet=%s, VPO=%d, F=[%.3f %.3f] Hz\n', wname, vpo, fLo, fHi);

                % ---- 预跑获取尺寸 ----
                app.UIFigure.Pointer = 'watch'; drawnow;

                x0 = double(tempdata(1, :, 1));
                [WT0, F, COI] = cwt(x0, fs, wname, ...     % ← 第三个返回值改名为 COI
                    'VoicesPerOctave', vpo, ...
                    'FrequencyLimits', [fLo fHi]);
                nF = numel(F);
                nSeg = nT;                                   % ← CWT 输出的时间维 = 输入信号长度

                % ✅ 自己构造时间向量
                T = (0:nT-1) / fs;                           % 单位秒,nT 个点

                % ---- 预分配 ----
                WT = complex(zeros(nF, nSeg, nCh, nTr, 'like', WT0));

                % ---- 进度条 ----
                hWait = uiprogressdlg(app.UIFigure, ...
                    'Title', 'CWT 计算中...', ...
                    'Message', sprintf('处理 %d 通道 × %d 试次', nCh, nTr), ...
                    'Cancelable', 'on');
                totalIter = nCh * nTr;
                iter = 0;

                % ---- 主循环 ----
                tic;
                for ch = 1:nCh
                    if hWait.CancelRequested
                        close(hWait);
                        app.UIFigure.Pointer = 'arrow';
                        disp('=== [CWT] 用户取消 ===');
                        return;
                    end

                    for tr = 1:nTr
                        x = double(tempdata(ch, :, tr));
                        WTi = cwt(x, fs, wname, ...
                            'VoicesPerOctave', vpo, ...
                            'FrequencyLimits', [fLo fHi]);
                        WT(:, :, ch, tr) = WTi;

                        iter = iter + 1;
                        hWait.Value = iter / totalIter;
                    end
                    hWait.Message = sprintf('已完成 %d/%d 通道', ch, nCh);
                end
                elapsed = toc;

                close(hWait);
                app.UIFigure.Pointer = 'arrow';

                % ✅ 结果封装并保存
                cwtResult = struct();
                cwtResult.WT = WT;
                cwtResult.F = F;
                cwtResult.T = T;
                cwtResult.COI = COI;
                cwtResult.params = struct(...
                    'wavelet', wname, ...
                    'voicesPerOctave', vpo, ...
                    'freqRange', [fLo, fHi], ...
                    'srate', fs);
                cwtResult.timestamp = datetime('now');

                app.TF_Result = cwtResult;

                % ✅ 可视化
                app.showCWTPreview(cwtResult);

                msg = sprintf('CWT 完成 | F=%d × T=%d × Ch=%d × Tr=%d | 用时 %.2fs', ...
                    nF, nSeg, nCh, nTr, elapsed);
                disp(msg);
                uialert(app.UIFigure, msg, '完成');

            catch ME
                if exist('hWait', 'var') && isvalid(hWait)
                    close(hWait);
                end
                app.UIFigure.Pointer = 'arrow';
                disp('[ERR] CWT 出错：');
                disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), '运行出错');
            end
        end

        % Button pushed function: DSOKButton
        function DSOKButtonPushed(app, event)
             try
                disp('=== [Resampling] START ===');

                % ✅ 用流水线数据
                if isempty(app.ModifiedData)
                    uialert(app.UIFigure, '没有可用数据。', '错误');
                    return;
                end
                tempdata = app.ModifiedData;

                % 读取 GUI
                D = str2double(string(app.DownsamplingFacorEditField.Value));
                if isempty(D) || ~isfinite(D) || D <= 1
                    uialert(app.UIFigure, ...
                        'Downsampling Factor 必须是 > 1 的数（例如 2、4、8）', '参数错误');
                    return;
                end

                if ndims(tempdata) ~= 3
                    [nCh, nT] = size(tempdata);
                    tempdata = reshape(tempdata, [nCh, nT, 1]);
                end

                % 新采样率预检查：防止低于奈奎斯特限制
                newSrate = double(app.srate) / D;

                % 在 GUI 层
                nyquist_old = app.srate / 2;
                nyquist_new = newSrate / 2;
                selection = uiconfirm(app.UIFigure, ...
                    sprintf(['降采样将把最高保留频率从 %.1f Hz 降为 %.1f Hz。\n' ...
                    '高于 %.1f Hz 的信号成分将被滤除。\n' ...
                    '是否继续?'], ...
                    nyquist_old, nyquist_new, nyquist_new), ...
                    '降采样确认', 'Options', {'继续', '取消'});
                if strcmp(selection, '取消'), return; end
                % 执行降采样
                app.UIFigure.Pointer = 'watch'; drawnow;
                tic;
                [data_ds, fs_ds] = seal_applyDownsampleResample(tempdata, app.srate, D);
                elapsed = toc;
                app.UIFigure.Pointer = 'arrow';

                % ✅ 写入流水线
                app.ModifiedData = data_ds;
                app.srate = fs_ds;
                app.change = 1;

                % ✅ 同步 DataNode 元数据
                app.dataNode.dataInfo.srate = fs_ds;
                app.dataNode.dataInfo.size = size(data_ds);

                % ✅ 刷新频谱
                app.plotMultiChannelSpectrum(app.ModifiedData);

                % 反馈
                msg = sprintf('Downsample 完成：D=%g | 新采样率=%.6g Hz | 新长度=%d | 用时 %.2fs', ...
                    D, fs_ds, size(data_ds, 2), elapsed);
                disp(msg);
                uialert(app.UIFigure, msg, '完成');
                disp('=== [Resampling] DONE ===');

            catch ME
                app.UIFigure.Pointer = 'arrow';
                disp('[ERR] Downsample 出错：');
                disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), '运行出错');
            end
        end

        % Button pushed function: interpolateOKButton
        function interpolateOKButtonPushed(app, event)
            try
                % ✅ 用流水线数据
                if isempty(app.ModifiedData)
                    uialert(app.UIFigure, '没有可用数据。', '错误');
                    return;
                end
                tempdata = app.ModifiedData;

                % ✅ 放宽:支持 2D 连续数据和 3D 分段数据
                assert(~isempty(tempdata) && ismember(ndims(tempdata), [2,3]), ...
                    'app.data 需为 2D (ch×t) 或 3D (ch×t×tr)');
                nCh = size(tempdata, 1);

                % 若没标签,自动生成
                if ~isprop(app,'chanLabels') || isempty(app.chanLabels)
                    app.chanLabels = arrayfun(@(k)sprintf('Ch%d',k),1:nCh,'uni',0);
                end

                % 若没 3D 坐标,尝试从极坐标补
                assert(isprop(app,'chanlocs') && ~isempty(app.chanlocs), ...
                    '缺少 chanlocs(X/Y/Z 或 theta/phi)');

                % -------- Interpolate --------
                badText = string(app.BadchannelsEditField.Value);
                method  = lower(string(app.InterpolateMethodDropDown.Value));
                k       = max(1, round(str2double(string(app.InterpKEditField.Value))));
                radius  = str2double(string(app.InterpRadiusEditField.Value));
                useTpl  = app.UsetemplatemontageCheckBox.Value;

                badIdx  = seal_parseChanList(badText, nCh, app.chanLabels);
                if isempty(badIdx)
                    uialert(app.UIFigure, '未检测到需要修复的坏道,无需插值。', '提示');
                    return;
                end

                app.ModifiedData = seal_runInterpolation(tempdata, app.chanlocs, ...
                    badIdx, method, k, radius, useTpl);

                uialert(app.UIFigure, ...
                    sprintf('Interpolate 完成(%s),修复坏道数=%d', method, numel(badIdx)), ...
                    '完成');
                app.change = 1;

            catch ME
                disp(' 出错:');
                disp(getReport(ME,'extended','hyperlinks','off'));
                uialert(app.UIFigure, getReport(ME,'basic','hyperlinks','off'), '运行出错');
            end
        end

        % Button pushed function: LapOKButten
        function LapOKButtenButtonPushed(app, event)
            try
                % ===== 数据可用性检查 =====
                if isempty(app.ModifiedData)
                    uialert(app.UIFigure, '没有可用数据。', '错误');
                    return;
                end
                tempdata = app.ModifiedData;

                % ✅ 维度放宽:支持 2D (ch×t) 和 3D (ch×t×tr)
                assert(~isempty(tempdata) && ismember(ndims(tempdata), [2,3]), ...
                    'app.data 需为 2D (ch×t) 或 3D (ch×t×tr)');
                nCh = size(tempdata, 1);

                % 若没标签,自动生成
                if ~isprop(app, 'chanLabels') || isempty(app.chanLabels)
                    app.chanLabels = arrayfun(@(k) sprintf('Ch%d', k), 1:nCh, 'uni', 0);
                end

                % 必须有通道坐标(Laplacian 需要空间几何信息)
                assert(isprop(app, 'chanlocs') && ~isempty(app.chanlocs), ...
                    '缺少 app.chanlocs(X/Y/Z 或 theta/phi)');
                if length(app.chanlocs) ~= nCh
                    uialert(app.UIFigure, ...
                        sprintf('通道坐标数(%d)与数据通道数(%d)不一致', ...
                        length(app.chanlocs), nCh), ...
                        '通道不匹配');
                    return;
                end

                % ===== 读取 GUI 参数(带范围保护) =====
                m = str2double(string(app.LapMEditField.Value));
                if ~isfinite(m) || m < 2 || m > 10
                    warning('SEAL:Laplacian:BadM', 'm 非法或超出常用范围 [2,10],改用默认值 4');
                    m = 4;
                end

                lambda = str2double(string(app.LapLambdaEditField.Value));
                if ~isfinite(lambda) || lambda < 0
                    warning('SEAL:Laplacian:BadLambda', 'λ 非法,改用默认值 1e-5');
                    lambda = 1e-5;
                end

                useCSD = app.LapUseCSDCheckBox.Value;
                edge   = string(app.LapEdgeDropdown.Value);

                % ===== 执行 =====
                app.UIFigure.Pointer = 'watch'; drawnow;
                tic;
                app.ModifiedData = seal_runLaplacian(tempdata, app.chanlocs, ...
                    m, lambda, useCSD, edge);
                elapsed = toc;
                app.UIFigure.Pointer = 'arrow';

                app.change = 1;

                % ✅ 刷新频谱图(Laplacian 改变了信号,频谱也变了)
                app.plotMultiChannelSpectrum(app.ModifiedData);

                % ===== 反馈 =====
                msg = sprintf(['Laplacian 完成\n' ...
                    '──────────\n' ...
                    'm = %g, λ = %g\n' ...
                    '用 CSD 工具箱: %s\n' ...
                    '边缘处理: %s\n' ...
                    '用时: %.2f s\n\n' ...
                    '注意: 信号单位已变为 μV/cm²,后续不要与未 Laplacian 的数据混用'], ...
                    m, lambda, mat2str(useCSD), edge, elapsed);
                uialert(app.UIFigure, msg, '完成');

            catch ME
                app.UIFigure.Pointer = 'arrow';   % ✅ 异常时也要复位指针
                disp('[ERR] Laplacian 出错:');
                disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    '运行出错');
            end

        end

        % Menu selected function: SaveMenu
        function SaveMenuSelected(app, event)
            SavePath = fullfile(app.sessionNode.path, "Seal_outputData");
            if ~isfolder(SavePath)
                mkdir(SavePath);
            end
            Path = uigetdir(SavePath, '选择结果保存文件夹');
            app.Callingapp.bringToFront();
            if isequal(Path, 0)
                return;  % 用户取消
            end
            timestamp = datestr(now, 'yyyymmdd_HHMMss');
            fileName = sprintf('%s_preprocessed_%s.mat', app.dataNode.name, timestamp);
            fullPath = fullfile(Path, fileName);

            savestruct.data = app.ModifiedData;
            savestruct.srate = app.srate;

            saveData(fullPath, savestruct);
            datanode = app.sessionNode.openDataFromData(fullPath);
            if ~isempty(app.dataNode.dataInfo) && isprop(app.dataNode.dataInfo, 'event')
                datanode.dataInfo.event = app.dataNode.dataInfo.event;
            end
            app.Callingapp.flushFileTree();
            if ~isempty(app.NoiseCovMat)
                datanode.dataInfo.NoiseCovMat = app.NoiseCovMat;
            end
            % ===== 2. 额外按 trial 拆分保存(仅当数据是分段的)=====
            if app.isEpoched(app.ModifiedData)
                app.saveTrialsSeparately(app.ModifiedData, Path, baseName);
            end

            app.sessionNode.save();
            app.change = 0;
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            
            if app.change == 1
                selection = uiconfirm(app.UIFigure, ...
                    '是否保存？', ...
                    '确认', ...
                    'Options', {'Yes','No'}, ...
                    'DefaultOption', 1);
                if strcmp(selection,'Yes')
                    SavePath = fullfile(app.sessionNode.path, "Seal_outputData");
                    if ~isfolder(SavePath)
                        mkdir(SavePath);
                    end
                    Path = uigetdir(SavePath, '选择结果保存文件夹');
                    app.Callingapp.bringToFront();
                    if isequal(Path, 0)
                        return;  % 用户取消
                    end
                    timestamp = datestr(now, 'yyyymmdd_HHMMss');
                    fileName = sprintf('%s_preprocessed_%s.mat', app.dataNode.name, timestamp);
                    fullPath = fullfile(Path, fileName);

                    savestruct.data = app.ModifiedData;
                    savestruct.srate = app.srate;

                    saveData(fullPath, savestruct);
                    datanode = app.sessionNode.openDataFromData(fullPath);
                    if ~isempty(app.dataNode.dataInfo) && isprop(app.dataNode.dataInfo, 'event')
                        datanode.dataInfo.event = app.dataNode.dataInfo.event;
                    end

                    if ~isempty(app.NoiseCovMat)
                        datanode.dataInfo.NoiseCovMat = app.NoiseCovMat;
                    end

                    if app.isEpoched(app.ModifiedData)
                        app.saveTrialsSeparately(app.ModifiedData, Path, ...
                            sprintf('%s_preprocessed_%s', app.dataNode.name, timestamp));
                    end

                    app.Callingapp.flushFileTree();
                    app.sessionNode.save();
                end
            end
            delete(app)
        end

        % Button pushed function: RemoveComponentsButton
        function RemoveComponentsButtonPushed(app, event)
            try
                disp('=== [ICA] 开始计算并生成预览对比图 ===');

                % 1. 检查前提条件
                if isempty(app.icaresult) || ~isfield(app.icaresult, 'A') || ~isfield(app.icaresult, 'S')
                    uialert(app.UIFigure, '尚未找到完整的 ICA 结果 (缺少矩阵 A 或 S)，请先运行 ICA！', '错误');
                    return;
                end

                % 2. 获取用户输入
                bad_comps_str = strtrim(app.RemoveComponengtsIndexes.Value);
                if isempty(bad_comps_str)
                    uialert(app.UIFigure, '请输入需要剔除的成分序号！', '提示');
                    return;
                end

                bad_comps = unique(round(str2num(bad_comps_str))); %#ok<ST2NM>
                if isempty(bad_comps)
                    uialert(app.UIFigure, '成分输入格式无法解析', '格式错误');
                    return;
                end

                % 3. 维度检查与矩阵计算
                [~, nT, nTr] = size(app.ModifiedData);
                ica = app.icaresult;
                total_comps = size(ica.S, 1);

                if any(bad_comps < 1) || any(bad_comps > total_comps)
                    uialert(app.UIFigure, sprintf('成分序号越界！有效范围为 1 至 %d。', total_comps), '越界错误');
                    return;
                end

                good_comps = setdiff(1:total_comps, bad_comps);
                if isempty(good_comps)
                    uialert(app.UIFigure, '不能剔除所有成分，信号将变为空白！', '操作被拒绝');
                    return;
                end

                % === 核心计算，但不覆盖 app.ModifiedData ===
                A_good = ica.A(:, good_comps);
                S_good = ica.S(good_comps, :);
                clean_X = A_good * S_good;
                num_ica_chans = length(ica.chanind);
                clean_data3d = reshape(clean_X, [num_ica_chans, nT, nTr]);

                % ==========================================
                % 生成 EEGLAB 风格的波形对比预览图
                % ==========================================
                % 创建临时 UIFigure
                cmpFig = uifigure('Name', 'Verify Component Removal', 'Position', [150, 150, 900, 600]);
                ax = uiaxes(cmpFig, 'Position', [50, 50, 800, 500]);
                hold(ax, 'on');

                %% 为了防止画面太乱，最多只画出前 5 个通道的第 1 个试次（Trial）
                num_plot = min(5, num_ica_chans);
                time_pts = 1:nT;

                % 动态计算偏移量，防止通道重叠 (取这几个通道标准差的最大值的 6 倍作为间距)
                offset_base = max(std(app.ModifiedData(ica.chanind(1:num_plot), :, 1), 0, 2)) * 6;

                % 逐通道绘制
                for c = 1:num_plot
                    % 计算纵向偏移量（通道序号越小，画得越高）
                    offset = (num_plot - c) * offset_base;

                    % 原始信号 (深蓝色)
                    orig_sig = app.ModifiedData(ica.chanind(c), :, 1) + offset;
                    % 清洗后的信号 (红色)
                    clean_sig = clean_data3d(c, :, 1) + offset;

                    % 绘制对比（保留句柄用于图例）
                    p1 = plot(ax, time_pts, orig_sig, 'Color', [0 0 0.5], 'LineWidth', 1.2); % 深蓝色
                    p2 = plot(ax, time_pts, clean_sig, 'Color', [0.8 0 0], 'LineWidth', 1);  % 红色
                end

                % 美化图表
                title(ax, sprintf('ICA Cleaning Preview (Showing Trial 1, up to %d channels)\nClose this window to Confirm or Cancel.', num_plot));
                xlabel(ax, 'Time Points');
                ylabel(ax, 'Amplitude');
                ax.YTick = []; % 隐藏毫无意义的绝对数值刻度
                legend(ax, [p1, p2], {'Original (Dark Blue)', 'Cleaned (Red)'}, 'Location', 'northeast');
                hold(ax, 'off');
                %%
                cmpFig.CloseRequestFcn = @(src, ev) app.verifyICACloseReq(src, clean_data3d, ica.chanind, bad_comps);

            catch ME
                disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), '预览生成出错');
            end
        end

        % Menu selected function: ResetMenu
        function ResetMenuSelected(app, event)
            app.ModifiedData = app.rawdata;
            app.Callingapp.appendOutput('The data has been reseted.')
        end

        % Button pushed function: ComponentsButton
        function ComponentsButtonPushed(app, event)
             try
        disp('=== [ICA] 准备绘制独立成分空间特征 ===');

        % ===== 前置检查 =====
        if isempty(app.icaresult) || ~isfield(app.icaresult, 'A')
            uialert(app.UIFigure, '尚未找到完整的 ICA 结果,请先运行 ICA!', '错误');
            return;
        end

        ica = app.icaresult;
        A = ica.A;                     % 混合矩阵 (通道数 × 成分数)
        num_ics = size(A, 2);

        % 限制一次最多绘制的数量
        plot_ics = min(num_ics, 36);
        if num_ics > 36
            disp('[提示] 成分较多,本次仅绘制前 36 个成分。');
        end

        cols = ceil(sqrt(plot_ics));
        rows = ceil(plot_ics / cols);

        % ===== 检测可用的 topoplot =====
        % 优先使用 EEGLAB 的 topoplot(效果更好)
        useTopoplot = exist('topoplot', 'file') == 2 && ...
                      ~isempty(app.chanlocs) && ...
                      isfield(app.chanlocs, 'theta') && ...
                      isfield(app.chanlocs, 'radius');

        if ~useTopoplot
            if isempty(app.chanlocs) || ~isfield(app.chanlocs, 'theta')
                disp('[警告] 缺少有效的通道坐标信息(chanlocs),将降级绘制柱状图');
            else
                disp('[警告] 未检测到 topoplot 函数,将降级绘制柱状图');
            end
        end

        % ===== 创建图窗 =====
        figTopo = figure('Name', 'ICA Component Spatial Weights', ...
                         'NumberTitle', 'off', ...
                         'Color', 'w', ...
                         'Position', [150, 150, 900, 800]);

        app.UIFigure.Pointer = 'watch'; drawnow;

        % ===== 循环绘制每个成分 =====
        for i = 1:plot_ics
            ax = subplot(rows, cols, i);
            ic_topo = A(:, i);

            draw_bar = false;

            if useTopoplot
                locs = app.chanlocs(ica.chanind);
                locs = ensureFullChanlocs(locs); 
                if length(ic_topo) == length(locs)
                    try
                        % ✅ 切到目标 axes
                        axes(ax);
                        % ✅ 调用 EEGLAB 的 topoplot
                        topoplot(ic_topo, locs, ...
                            'electrodes', 'off', ...
                            'maplimits', 'absmax', ...     % 对称色标
                            'style',     'both', ...        % 填色 + 等高线
                            'shading',   'interp', ...      % 平滑
                            'whitebk',   'on', ...          % 白底
                            'verbose',   'off');            % 关闭刷屏
                    catch ME
                        % topoplot 偶发失败:降级
                        fprintf('[IC %d] topoplot 失败: %s\n', i, ME.message);
                        cla(ax);
                        draw_bar = true;
                    end
                else
                    draw_bar = true;
                end
            else
                draw_bar = true;
            end

            % ===== 降级:绘制柱状图权重 =====
            if draw_bar
                bar(ax, ic_topo, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
                axis(ax, 'tight');
                y_lim = max(abs(ic_topo)) * 1.1;
                if y_lim > 0, ylim(ax, [-y_lim, y_lim]); end
                hold(ax, 'on');
                plot(ax, xlim(ax), [0 0], 'k-', 'LineWidth', 0.5);
                hold(ax, 'off');
                set(ax, 'XTick', [], 'YTick', []);
            end

            title(ax, sprintf('IC %d', i), 'FontSize', 10, 'FontWeight', 'bold');
        end

        % ===== 统一颜色映射 =====
        colormap(figTopo, jet);

        app.UIFigure.Pointer = 'arrow';
        disp('=== [ICA] 空间特征绘制完成 ===');

    catch ME
        app.UIFigure.Pointer = 'arrow';
        disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
        uialert(app.UIFigure, '绘制过程发生错误,请查看控制台日志。', '绘图出错');
    end
        end

        % Button pushed function: ExtracttraialsSetButton
        function ExtracttraialsSetButtonPushed(app, event)
            try
                tempdata = app.ModifiedData;
                assert(~isempty(tempdata) && ndims(tempdata)<=3, '数据维度不合法');

                % 确保标签存在
                [nCh,~,~] = size(tempdata);
                if ~isprop(app,'chanLabels') || isempty(app.chanLabels)
                    app.chanLabels = arrayfun(@(k)sprintf('Ch%d',k),1:nCh,'uni',0);
                end

                % 读取用户在界面上设置的参数（如果没有控件，可以先用默认值 4 和 0.4）
                prompt = {'请输入方差 Z-score 阈值 (推荐 3.0 - 5.0，越小越严格):', ...
                    '请输入最小相关性阈值 (推荐 0.4，取值 0-1):'};
                dlgtitle = '自动坏导检测参数设置';
                dims = [1 55]; % 输入框的行数和字符宽度
                definput = {'4.0', '0.4'}; % 默认参数

                % 阻塞式弹窗，等待用户输入
                answer = inputdlg(prompt, dlgtitle, dims, definput);

                % 检查用户是否点击了“取消”或关闭了窗口
                if isempty(answer)
                    return; % 悄悄退出，不执行后续操作
                end

                % 将用户输入的字符串转换为数字
                zThresh = str2double(answer{1});
                corrThresh = str2double(answer{2});

                % 格式校验：如果输入的不是合法数字，拦截报错
                if isnan(zThresh) || isnan(corrThresh) || zThresh <= 0 || corrThresh < 0 || corrThresh >= 1
                    uialert(app.UIFigure, '输入的参数无效，请输入正确的数字！', '参数错误');
                    return;
                end

                % 执行自动剔除
                [app.ModifiedData, app.chanlocs, app.chanLabels, badIdx, badNames] = ...
                    seal_runAutoRejectChannels(tempdata, app.chanlocs, app.chanLabels, zThresh, corrThresh);

                % 弹窗报告结果
                if isempty(badIdx)
                    uialert(app.UIFigure, '检测完成。当前数据质量很好，未发现坏导。', '完成');
                else
                    msg = sprintf('检测并剔除了 %d 个坏导。\n剔除的通道：\n%s', ...
                        numel(badIdx), strjoin(string(badNames), ', '));
                    uialert(app.UIFigure, msg, '完成');
                end

                app.change = 1;

            catch ME
                disp(' 出错：');
                disp(getReport(ME,'extended','hyperlinks','off'));
                uialert(app.UIFigure, getReport(ME,'basic','hyperlinks','off'), '运行出错');
            end
        end

        % Button pushed function: Baselinesetbutten
        function BaselinesetbuttenButtonPushed(app, event)
            try
                tempdata = app.ModifiedData;
                % 兼容 3D (分段数据) 和 2D (连续数据)
                assert(~isempty(tempdata) && ndims(tempdata)<=3, '缺少有效数据，需为 ch×t×tr 或 ch×t 格式');

                [~, nTime, ~] = size(tempdata);

                % ==========================================
                % 弹出参数输入框
                % ==========================================
                prompt = {sprintf('请输入基线起点 (点数，范围 1 - %d):', nTime), ...
                    sprintf('请输入基线终点 (点数，范围 1 - %d):', nTime)};
                dlgtitle = '基线校准参数设置';
                dims = [1 55];

                % 默认值：推荐前 200 个点，如果总数据不到 200 点，则取 1/10
                defaultEnd = num2str(min(200, floor(nTime/10)));
                definput = {'1', defaultEnd};

                % 阻塞式弹窗，等待用户输入
                answer = inputdlg(prompt, dlgtitle, dims, definput);

                if isempty(answer)
                    return; % 用户点击取消，悄悄退出
                end

                % 解析输入的数字（四舍五入防手滑输入小数）
                baseStart = round(str2double(answer{1}));
                baseEnd = round(str2double(answer{2}));

                % 格式和范围校验
                if isnan(baseStart) || isnan(baseEnd) || baseStart < 1 || baseEnd > nTime || baseStart > baseEnd
                    uialert(app.UIFigure, '输入的基线区间无效，请检查数值范围！', '参数错误');
                    return;
                end
                % ==========================================

                % 执行底层基线校准算法
                app.ModifiedData = seal_runBaselineCorrection(tempdata, baseStart, baseEnd);

                % 弹窗报告结果
                msg = sprintf('基线校准完成！\n选取的基线区间为：第 %d 点 到 第 %d 点。', baseStart, baseEnd); 
                uialert(app.UIFigure, msg, '完成');

                app.change = 1;

            catch ME
                disp(' 出错：');
                disp(getReport(ME,'extended','hyperlinks','off'));
                uialert(app.UIFigure, getReport(ME,'basic','hyperlinks','off'), '运行出错');
            end
        end

        % Button pushed function: AveragetrialsSetButton_3
        function AveragetrialsSetButton_3Pushed(app, event)
            try
                tempdata = app.ModifiedData;

                % 维度检查：必须是分段后的 3D 数据才能叠加
                if isempty(tempdata) || ndims(tempdata) < 3 || size(tempdata, 3) == 1
                    uialert(app.UIFigure, '当前数据不满足叠加条件。请确保数据已分段 (Epoching) 且包含多个 trials。', '操作无效');
                    return;
                end

                nTrials = size(tempdata, 3);

                % ==========================================
                % 弹出选择框，询问平均策略
                % ==========================================
                qString = sprintf('当前数据包含 %d 个 trials。\n请选择叠加平均的策略：\n\n[传统均值]：标准做法\n[稳健修剪]：自动剔除极端伪迹后求均值(推荐)\n[中位数]：绝对抗噪，但波形不平滑', nTrials);
                choice = questdlg(qString, '叠加平均设置', ...
                    '传统均值', '稳健修剪', '中位数', '传统均值');

                % 用户点击关闭或取消
                if isempty(choice)
                    return;
                end

                % 将中文选项映射为底层方法字符串
                switch choice
                    case '传统均值', method = 'mean';
                    case '稳健修剪', method = 'trim';
                    case '中位数',   method = 'median';
                end
                % ==========================================

                % 执行叠加平均
                app.ModifiedData = seal_runEpochAveraging(tempdata, method);

                % 弹窗报告结果
                uialert(app.UIFigure, sprintf('叠加平均完成！\n使用策略：%s\n数据已从 3D (单次 trial) 降维为 2D (ERP 波形)。', choice), '完成');

                app.change = 1;

            catch ME
                disp(' 出错：');
                disp(getReport(ME,'extended','hyperlinks','off'));
                uialert(app.UIFigure, getReport(ME,'basic','hyperlinks','off'), '运行出错');
            end
        end

        % Button pushed function: EpochingButton
        function EpochingButtonPushed(app, event)
            try
                % ===== Step 1:校验数据状态 =====
                tempdata = app.ModifiedData;
                if isempty(tempdata)
                    uialert(app.UIFigure, '没有可用数据。', '错误');
                    return;
                end
                if ~ismatrix(tempdata)
                    uialert(app.UIFigure, ...
                        sprintf(['分段功能需要 2D 连续数据,当前数据维度 %d,' ...
                        '可能已经是分段数据。\n' ...
                        '若需重新分段,请先 Reset 恢复原始数据。'], ndims(tempdata)), ...
                        '数据维度错误');
                    return;
                end

                [nCh, nT] = size(tempdata);
                totalSec = nT / app.srate;

                % ===== Step 2:选择事件源 =====
                sourceChoice = questdlg(...
                    sprintf(['当前数据:%d 通道 × %d 样本(%.1f 秒 @ %g Hz)\n\n' ...
                    '请选择事件来源:'], nCh, nT, totalSec, app.srate), ...
                    '分段 Step 1/2:事件来源', ...
                    '事件通道(从数据通道检测脉冲)', ...
                    '固定间隔(静息态)', ...
                    'events 结构体(高级)', ...
                    '事件通道(从数据通道检测脉冲)');

                if isempty(sourceChoice), return; end

                eventSpec = struct();

                switch sourceChoice
                    case '事件通道(从数据通道检测脉冲)'
                        eventSpec.source = 'event_channel';
                        prompt = {...
                            sprintf('事件通道索引(1 - %d):', nCh), ...
                            '脉冲检测阈值(留空为自动 = 最大值 × 0.5):'};
                        definput = {num2str(nCh), ''};  % 最后一个通道常被用作触发
                        answer = inputdlg(prompt, '事件通道参数', [1 50], definput);
                        if isempty(answer), return; end

                        eventSpec.data = round(str2double(answer{1}));
                        thresh = str2double(answer{2});
                        if isnan(thresh) || isempty(answer{2})
                            eventSpec.threshold = NaN;  % 让底层自动决定
                        else
                            eventSpec.threshold = thresh;
                        end

                        if ~isfinite(eventSpec.data) || eventSpec.data < 1 || eventSpec.data > nCh
                            uialert(app.UIFigure, '事件通道索引无效', '参数错误'); return;
                        end

                    case '固定间隔(静息态)'
                        eventSpec.source = 'fixed_interval';
                        answer = inputdlg({'间隔秒数(如 2 表示每 2 秒切一段):'}, ...
                            '固定间隔参数', [1 50], {'2'});
                        if isempty(answer), return; end
                        eventSpec.data = str2double(answer{1});
                        if ~isfinite(eventSpec.data) || eventSpec.data <= 0
                            uialert(app.UIFigure, '间隔秒数必须为正', '参数错误'); return;
                        end

                    case 'events 结构体(高级)'

                        ev = [];
                        if isprop(app.dataNode.dataInfo, 'event')
                            ev = app.dataNode.dataInfo.event;
                        end
                        if isempty(ev)
                            uialert(app.UIFigure, ...
                                '当前数据节点未附带 events 结构体。\n请使用其他事件来源方式。', ...
                                '无 events 数据');
                            return;
                        end
                        eventSpec.source = 'events_struct';
                        eventSpec.data = ev;

                        % ===== 格式化事件类型显示 =====
                        allTypesRaw = {ev.type};

                        % 先把 type 字段统一到 cell of string
                        allTypesStr = cell(size(allTypesRaw));
                        for i = 1:numel(allTypesRaw)
                            t = allTypesRaw{i};
                            if isnumeric(t)
                                allTypesStr{i} = num2str(t);
                            elseif ischar(t) || isstring(t)
                                allTypesStr{i} = char(t);
                            else
                                allTypesStr{i} = '<unknown>';
                            end
                        end

                        % 去重 + 计数
                        [uniqueTypes, ~, idx] = unique(allTypesStr);
                        counts = accumarray(idx, 1);

                        % 组装显示字符串:每行一个类型 + 计数
                        lines = cell(numel(uniqueTypes), 1);
                        for i = 1:numel(uniqueTypes)
                            lines{i} = sprintf('  %s  (× %d)', uniqueTypes{i}, counts(i));
                        end
                        typeStr = strjoin(lines, newline);

                        promptMsg = sprintf(['检测到 %d 种事件类型:\n\n%s\n\n' ...
                            '输入要提取的类型(逗号分隔,空 = 全部):'], ...
                            numel(uniqueTypes), typeStr);

                        answer = inputdlg({promptMsg}, '事件类型筛选', [10 60], {''});
                        if isempty(answer), return; end

                        if ~isempty(strtrim(answer{1}))
                            parts = strtrim(strsplit(answer{1}, ','));
                            parts = parts(~cellfun(@isempty, parts));

                            % 探测 ev.type 的实际类型(以第一个非空事件为准)
                            refType = 'string';   % 默认当字符串处理
                            for i = 1:numel(ev)
                                t = ev(i).type;
                                if ~isempty(t)
                                    if isnumeric(t)
                                        refType = 'numeric';
                                    else
                                        refType = 'string';
                                    end
                                    break;
                                end
                            end

                            % 按 ev.type 的实际类型来转换用户输入
                            switch refType
                                case 'numeric'
                                    % event.type 是数字 → 把用户输入全部转成 double 数组
                                    nums = str2double(parts);
                                    if any(isnan(nums))
                                        badIdx = find(isnan(nums), 1);
                                        uialert(app.UIFigure, sprintf(...
                                            ['事件类型为数字型,但输入的 "%s" 不是有效数字。\n' ...
                                            '请输入数字(如 1, 2, 3)。'], parts{badIdx}), ...
                                            '类型不匹配');
                                        return;
                                    end
                                    eventSpec.typesWanted = nums;

                                case 'string'
                                    % event.type 是字符串 → 全部当成字符串 cell
                                    eventSpec.typesWanted = parts;   % 已经是 cell of char
                            end
                        end
                end

                % ===== Step 3:询问 epoch 窗口和质量筛查参数 =====

                % 根据事件源和数据长度推荐默认窗口
                switch eventSpec.source
                    case 'fixed_interval'
                        % 固定间隔:窗口 = [0, 间隔秒数]
                        defPre  = '0';
                        defPost = num2str(eventSpec.data);
                        hintStr = sprintf('(建议 [0, %g],对应固定间隔长度)', eventSpec.data);

                    otherwise
                        % 事件驱动:根据总时长给出建议
                        if totalSec < 2
                            % 极短数据,给一个很窄的窗口
                            defPre  = '-0.05';
                            defPost = num2str(min(0.2, totalSec/2 - 0.05), '%.2f');
                            hintStr = '(数据较短,建议窄窗口)';
                        elseif totalSec < 10
                            defPre  = '-0.1';
                            defPost = '0.5';
                            hintStr = '(短数据默认窗口)';
                        else
                            defPre  = '-0.2';
                            defPost = '0.8';
                            hintStr = '(ERP 常用默认窗口)';
                        end
                end

                prompt = {...
                    sprintf('Epoch 起点(秒,刺激前为负):\n数据总长 %.2f s %s', totalSec, hintStr), ...
                    'Epoch 终点(秒,刺激后为正):', ...
                    '振幅阈值(μV,超过则剔除该 trial,留空=不筛查):', ...
                    '边界处理(reject / zeropad / mirror):'};
                dlgtitle = sprintf('分段 Step 2/2:窗口与筛查 [数据 %.2f s]', totalSec);
                definput = {defPre, defPost, '100', 'reject'};
                answer = inputdlg(prompt, dlgtitle, [1 60], definput);
                if isempty(answer), return; end

                tPre  = str2double(answer{1});
                tPost = str2double(answer{2});
                rejAmp = str2double(answer{3});
                if isnan(rejAmp) || isempty(strtrim(answer{3}))
                    rejAmp = Inf;
                end
                boundaryMode = lower(strtrim(answer{4}));

                % ===== 窗口合法性检查(增强版)=====
                if ~isfinite(tPre) || ~isfinite(tPost) || tPre >= tPost
                    uialert(app.UIFigure, 'Epoch 窗口无效(需 t_pre < t_post)', '参数错误');
                    return;
                end

                winDur = tPost - tPre;
                if winDur > totalSec
                    uialert(app.UIFigure, sprintf(...
                        ['Epoch 窗口长度(%.2f s)超过数据总长(%.2f s),无法分段。\n' ...
                        '请缩小窗口范围。'], winDur, totalSec), '窗口过大');
                    return;
                end

                if winDur > totalSec / 2
                    choice = questdlg(sprintf(...
                        ['Epoch 窗口长度(%.2f s)超过数据总长的一半(%.2f s),\n' ...
                        '可能只能提取很少的 trial 或全部被边界剔除。\n\n是否继续?'], ...
                        winDur, totalSec), ...
                        '窗口偏大', '继续', '取消', '取消');
                    if ~strcmp(choice, '继续'), return; end
                end

                if ~ismember(boundaryMode, {'reject','zeropad','mirror'})
                    uialert(app.UIFigure, '边界处理方式无效', '参数错误');
                    return;
                end
                % ===== Step 4:组装 opts 并调用底层 =====
                opts = struct();
                opts.rejectAmp    = rejAmp;
                opts.boundaryMode = boundaryMode;
                opts.ampUnit      = 'uV';
                if isfield(app.dataNode.dataInfo, 'unit')
                    opts.ampUnit = char(app.dataNode.dataInfo.unit);
                end
                opts.verbose = true;

                % ===== Step 5:执行分段 =====
                app.UIFigure.Pointer = 'watch'; drawnow;
                tic;
                [dataEpoched, info] = seal_runEpoching( ...
                    tempdata, app.srate, eventSpec, [tPre, tPost], opts);
                elapsed = toc;
                app.UIFigure.Pointer = 'arrow';

                % ===== Step 6:写回 app 状态 =====
                app.ModifiedData = dataEpoched;
                app.change = 1;

      
                app.plotMultiChannelSpectrum(app.ModifiedData);

                % ===== Step 7:弹窗反馈 =====
                msg = sprintf([...
                    '分段完成!\n' ...
                    '──────────────\n' ...
                    '总事件数:     %d\n' ...
                    '成功提取:     %d\n' ...
                    '边界剔除:     %d\n' ...
                    '振幅剔除:     %d\n' ...
                    '──────────────\n' ...
                    '输出:         %d 通道 × %d 样本 × %d trials\n' ...
                    'Epoch 窗口:   [%.3f, %.3f] s\n' ...
                    '用时:         %.2f s'], ...
                    info.nTrialsFound, info.nTrialsExtracted, ...
                    info.nRejectedBound, info.nRejectedAmp, ...
                    size(dataEpoched, 1), size(dataEpoched, 2), size(dataEpoched, 3), ...
                    info.epochWindow(1), info.epochWindow(2), elapsed);
                uialert(app.UIFigure, msg, '完成');

            catch ME
                app.UIFigure.Pointer = 'arrow';
                disp('[ERR] Extract Trials 出错:');
                disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), '运行出错');
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1066 532];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create SaveMenu
            app.SaveMenu = uimenu(app.FileMenu);
            app.SaveMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveMenuSelected, true);
            app.SaveMenu.Text = 'Save';

            % Create ResetMenu
            app.ResetMenu = uimenu(app.FileMenu);
            app.ResetMenu.MenuSelectedFcn = createCallbackFcn(app, @ResetMenuSelected, true);
            app.ResetMenu.Text = 'Reset';

            % Create FrequencydomainFilteringPanel
            app.FrequencydomainFilteringPanel = uipanel(app.UIFigure);
            app.FrequencydomainFilteringPanel.Title = 'Frequency-domain Filtering';
            app.FrequencydomainFilteringPanel.FontWeight = 'bold';
            app.FrequencydomainFilteringPanel.Position = [1 243 260 290];

            % Create HighpassCheckBox
            app.HighpassCheckBox = uicheckbox(app.FrequencydomainFilteringPanel);
            app.HighpassCheckBox.ValueChangedFcn = createCallbackFcn(app, @HighpassCheckBoxValueChanged, true);
            app.HighpassCheckBox.Text = '';
            app.HighpassCheckBox.Position = [13 225 25 22];

            % Create LowpassCheckBox
            app.LowpassCheckBox = uicheckbox(app.FrequencydomainFilteringPanel);
            app.LowpassCheckBox.ValueChangedFcn = createCallbackFcn(app, @LowpassCheckBoxValueChanged, true);
            app.LowpassCheckBox.Text = '';
            app.LowpassCheckBox.Position = [13 185 25 22];

            % Create NotchCheckBox
            app.NotchCheckBox = uicheckbox(app.FrequencydomainFilteringPanel);
            app.NotchCheckBox.ValueChangedFcn = createCallbackFcn(app, @NotchCheckBoxValueChanged, true);
            app.NotchCheckBox.Text = '';
            app.NotchCheckBox.Position = [13 144 25 22];

            % Create HzLabel
            app.HzLabel = uilabel(app.FrequencydomainFilteringPanel);
            app.HzLabel.Position = [229 225 25 22];
            app.HzLabel.Text = 'Hz';

            % Create HzLabel_2
            app.HzLabel_2 = uilabel(app.FrequencydomainFilteringPanel);
            app.HzLabel_2.Position = [229 185 25 22];
            app.HzLabel_2.Text = 'Hz';

            % Create HzLabel_3
            app.HzLabel_3 = uilabel(app.FrequencydomainFilteringPanel);
            app.HzLabel_3.Position = [229 144 25 22];
            app.HzLabel_3.Text = 'Hz';

            % Create HighpassEditFieldLabel
            app.HighpassEditFieldLabel = uilabel(app.FrequencydomainFilteringPanel);
            app.HighpassEditFieldLabel.HorizontalAlignment = 'right';
            app.HighpassEditFieldLabel.Position = [39 225 60 22];
            app.HighpassEditFieldLabel.Text = 'High-pass';

            % Create HighpassEditField
            app.HighpassEditField = uieditfield(app.FrequencydomainFilteringPanel, 'numeric');
            app.HighpassEditField.Position = [114 225 100 22];

            % Create LowpassEditFieldLabel
            app.LowpassEditFieldLabel = uilabel(app.FrequencydomainFilteringPanel);
            app.LowpassEditFieldLabel.HorizontalAlignment = 'right';
            app.LowpassEditFieldLabel.Position = [42 185 57 22];
            app.LowpassEditFieldLabel.Text = 'Low-pass';

            % Create LowpassEditField
            app.LowpassEditField = uieditfield(app.FrequencydomainFilteringPanel, 'numeric');
            app.LowpassEditField.Position = [114 185 100 22];

            % Create NotchEditFieldLabel
            app.NotchEditFieldLabel = uilabel(app.FrequencydomainFilteringPanel);
            app.NotchEditFieldLabel.HorizontalAlignment = 'right';
            app.NotchEditFieldLabel.Position = [52 144 37 22];
            app.NotchEditFieldLabel.Text = 'Notch';

            % Create NotchEditField
            app.NotchEditField = uieditfield(app.FrequencydomainFilteringPanel, 'numeric');
            app.NotchEditField.Position = [114 144 100 22];

            % Create frequencydominSetbutten
            app.frequencydominSetbutten = uibutton(app.FrequencydomainFilteringPanel, 'push');
            app.frequencydominSetbutten.ButtonPushedFcn = createCallbackFcn(app, @frequencydominSetbuttenPushed, true);
            app.frequencydominSetbutten.Position = [30 41 203 22];
            app.frequencydominSetbutten.Text = 'Time-Frequency Decomposition >>';

            % Create FrequencydomainOKButton
            app.FrequencydomainOKButton = uibutton(app.FrequencydomainFilteringPanel, 'push');
            app.FrequencydomainOKButton.ButtonPushedFcn = createCallbackFcn(app, @FrequencydomainOKButtonPushed, true);
            app.FrequencydomainOKButton.Position = [171 85 69 23];
            app.FrequencydomainOKButton.Text = 'Apply';

            % Create SpatialFilteringPanel
            app.SpatialFilteringPanel = uipanel(app.UIFigure);
            app.SpatialFilteringPanel.Title = 'Spatial Filtering';
            app.SpatialFilteringPanel.FontWeight = 'bold';
            app.SpatialFilteringPanel.Position = [260 242 369 291];

            % Create MethodsButtonGroup
            app.MethodsButtonGroup = uibuttongroup(app.SpatialFilteringPanel);
            app.MethodsButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @MethodsButtonGroupSelectionChanged, true);
            app.MethodsButtonGroup.Title = 'Methods';
            app.MethodsButtonGroup.Position = [1 155 370 117];

            % Create ICAButton
            app.ICAButton = uiradiobutton(app.MethodsButtonGroup);
            app.ICAButton.Text = '     ICA';
            app.ICAButton.FontSize = 16;
            app.ICAButton.Position = [95 59 112 23];
            app.ICAButton.Value = true;

            % Create RereferenceButton
            app.RereferenceButton = uiradiobutton(app.MethodsButtonGroup);
            app.RereferenceButton.Text = '     Re-reference';
            app.RereferenceButton.FontSize = 16;
            app.RereferenceButton.Position = [95 40 167 22];

            % Create InterpolateButton
            app.InterpolateButton = uiradiobutton(app.MethodsButtonGroup);
            app.InterpolateButton.Text = '     Interpolate';
            app.InterpolateButton.FontSize = 16;
            app.InterpolateButton.Position = [95 20 120 22];

            % Create LaplacianButton
            app.LaplacianButton = uiradiobutton(app.MethodsButtonGroup);
            app.LaplacianButton.Text = '     Laplacian';
            app.LaplacianButton.FontSize = 16;
            app.LaplacianButton.Position = [95 0 113 22];

            % Create ERPrelatedPanel
            app.ERPrelatedPanel = uipanel(app.SpatialFilteringPanel);
            app.ERPrelatedPanel.Title = 'ERP-related';
            app.ERPrelatedPanel.FontWeight = 'bold';
            app.ERPrelatedPanel.Position = [0 0 369 154];

            % Create ExtracttraialsSetButton
            app.ExtracttraialsSetButton = uibutton(app.ERPrelatedPanel, 'push');
            app.ExtracttraialsSetButton.ButtonPushedFcn = createCallbackFcn(app, @ExtracttraialsSetButtonPushed, true);
            app.ExtracttraialsSetButton.Position = [223 102 75 22];
            app.ExtracttraialsSetButton.Text = 'Set >>';

            % Create Baselinesetbutten
            app.Baselinesetbutten = uibutton(app.ERPrelatedPanel, 'push');
            app.Baselinesetbutten.ButtonPushedFcn = createCallbackFcn(app, @BaselinesetbuttenButtonPushed, true);
            app.Baselinesetbutten.Position = [223 71 75 22];
            app.Baselinesetbutten.Text = 'Set >>';

            % Create AveragetrialsSetButton_3
            app.AveragetrialsSetButton_3 = uibutton(app.ERPrelatedPanel, 'push');
            app.AveragetrialsSetButton_3.ButtonPushedFcn = createCallbackFcn(app, @AveragetrialsSetButton_3Pushed, true);
            app.AveragetrialsSetButton_3.Position = [223 41 75 22];
            app.AveragetrialsSetButton_3.Text = 'Set >>';

            % Create BadChannelDetectionLabel
            app.BadChannelDetectionLabel = uilabel(app.ERPrelatedPanel);
            app.BadChannelDetectionLabel.Position = [56 102 128 22];
            app.BadChannelDetectionLabel.Text = 'Bad Channel Detection';

            % Create BaselineCorrectLabel
            app.BaselineCorrectLabel = uilabel(app.ERPrelatedPanel);
            app.BaselineCorrectLabel.Position = [56 71 95 22];
            app.BaselineCorrectLabel.Text = 'Baseline-Correct';

            % Create AverageTrialsLabel
            app.AverageTrialsLabel = uilabel(app.ERPrelatedPanel);
            app.AverageTrialsLabel.Position = [56 41 82 22];
            app.AverageTrialsLabel.Text = 'Average Trials';

            % Create EpochingButton
            app.EpochingButton = uibutton(app.ERPrelatedPanel, 'push');
            app.EpochingButton.ButtonPushedFcn = createCallbackFcn(app, @EpochingButtonPushed, true);
            app.EpochingButton.Position = [223 8 75 23];
            app.EpochingButton.Text = 'Epoching';

            % Create ExtractTrialsLabel
            app.ExtractTrialsLabel = uilabel(app.ERPrelatedPanel);
            app.ExtractTrialsLabel.Position = [56 12 74 22];
            app.ExtractTrialsLabel.Text = 'Extract Trials';

            % Create SpectrumPanel
            app.SpectrumPanel = uipanel(app.UIFigure);
            app.SpectrumPanel.Title = 'Spectrum';
            app.SpectrumPanel.FontWeight = 'bold';
            app.SpectrumPanel.Position = [1 1 628 242];

            % Create UIAxes
            app.UIAxes = uiaxes(app.SpectrumPanel);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [3 0 625 220];

            % Create ICAPanel
            app.ICAPanel = uipanel(app.UIFigure);
            app.ICAPanel.Title = 'ICA';
            app.ICAPanel.FontWeight = 'bold';
            app.ICAPanel.Position = [628 284 438 249];

            % Create ReordercomponentsbyvarianceCheckBox
            app.ReordercomponentsbyvarianceCheckBox = uicheckbox(app.ICAPanel);
            app.ReordercomponentsbyvarianceCheckBox.Text = 'Reorder components by variance';
            app.ReordercomponentsbyvarianceCheckBox.Position = [64 112 200 22];

            % Create ICAalgorithmDropDownLabel
            app.ICAalgorithmDropDownLabel = uilabel(app.ICAPanel);
            app.ICAalgorithmDropDownLabel.HorizontalAlignment = 'right';
            app.ICAalgorithmDropDownLabel.Position = [60 179 78 22];
            app.ICAalgorithmDropDownLabel.Text = 'ICA algorithm';

            % Create ICAalgorithmDropDown
            app.ICAalgorithmDropDown = uidropdown(app.ICAPanel);
            app.ICAalgorithmDropDown.Items = {'runica', 'binica'};
            app.ICAalgorithmDropDown.Position = [197 179 180 22];
            app.ICAalgorithmDropDown.Value = 'runica';

            % Create CommandlineoptiionsEditFieldLabel
            app.CommandlineoptiionsEditFieldLabel = uilabel(app.ICAPanel);
            app.CommandlineoptiionsEditFieldLabel.HorizontalAlignment = 'right';
            app.CommandlineoptiionsEditFieldLabel.Position = [60 146 124 22];
            app.CommandlineoptiionsEditFieldLabel.Text = 'Commandline optiions';

            % Create CommandlineoptiionsEditField
            app.CommandlineoptiionsEditField = uieditfield(app.ICAPanel, 'text');
            app.CommandlineoptiionsEditField.HorizontalAlignment = 'center';
            app.CommandlineoptiionsEditField.Position = [197 146 180 22];
            app.CommandlineoptiionsEditField.Value = '''extended'',1';

            % Create ChanneltypesorchannelindicesEditFieldLabel
            app.ChanneltypesorchannelindicesEditFieldLabel = uilabel(app.ICAPanel);
            app.ChanneltypesorchannelindicesEditFieldLabel.HorizontalAlignment = 'right';
            app.ChanneltypesorchannelindicesEditFieldLabel.Position = [60 82 191 22];
            app.ChanneltypesorchannelindicesEditFieldLabel.Text = 'Channel type(s) or channel indices';

            % Create ChanneltypesorchannelindicesEditField
            app.ChanneltypesorchannelindicesEditField = uieditfield(app.ICAPanel, 'text');
            app.ChanneltypesorchannelindicesEditField.Position = [277 82 100 22];

            % Create ICAOKButton
            app.ICAOKButton = uibutton(app.ICAPanel, 'push');
            app.ICAOKButton.ButtonPushedFcn = createCallbackFcn(app, @ICAOKButtonPushed, true);
            app.ICAOKButton.Position = [277 40 100 22];
            app.ICAOKButton.Text = 'OK';

            % Create RemoveComponentsButton
            app.RemoveComponentsButton = uibutton(app.ICAPanel, 'push');
            app.RemoveComponentsButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveComponentsButtonPushed, true);
            app.RemoveComponentsButton.Position = [221 6 156 23];
            app.RemoveComponentsButton.Text = 'Remove Components';

            % Create RemoveComponengtsIndexes
            app.RemoveComponengtsIndexes = uieditfield(app.ICAPanel, 'text');
            app.RemoveComponengtsIndexes.Position = [83 6 100 22];

            % Create ComponentsButton
            app.ComponentsButton = uibutton(app.ICAPanel, 'push');
            app.ComponentsButton.ButtonPushedFcn = createCallbackFcn(app, @ComponentsButtonPushed, true);
            app.ComponentsButton.Position = [137 39 100 23];
            app.ComponentsButton.Text = 'Components';

            % Create TimeFrequencyDecompositionPanel
            app.TimeFrequencyDecompositionPanel = uipanel(app.UIFigure);
            app.TimeFrequencyDecompositionPanel.Title = 'Time-Frequency Decomposition';
            app.TimeFrequencyDecompositionPanel.FontWeight = 'bold';
            app.TimeFrequencyDecompositionPanel.Position = [1 -572 439 531];

            % Create STFTLabel
            app.STFTLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.STFTLabel.Position = [20 478 35 22];
            app.STFTLabel.Text = 'STFT';

            % Create WaveletTransformLabel
            app.WaveletTransformLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.WaveletTransformLabel.Position = [20 290 194 22];
            app.WaveletTransformLabel.Text = 'Wavelet Transform';

            % Create WindowDropDownLabel
            app.WindowDropDownLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.WindowDropDownLabel.HorizontalAlignment = 'right';
            app.WindowDropDownLabel.Position = [108 447 48 22];
            app.WindowDropDownLabel.Text = 'Window';

            % Create WindowDropDown
            app.WindowDropDown = uidropdown(app.TimeFrequencyDecompositionPanel);
            app.WindowDropDown.Items = {'hann', 'hamming', 'gaussian', 'blackman'};
            app.WindowDropDown.Position = [291 447 100 22];
            app.WindowDropDown.Value = 'hann';

            % Create WindowLengthEditFieldLabel
            app.WindowLengthEditFieldLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.WindowLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.WindowLengthEditFieldLabel.Position = [95 418 88 22];
            app.WindowLengthEditFieldLabel.Text = 'Window Length';

            % Create WindowLengthEditField
            app.WindowLengthEditField = uieditfield(app.TimeFrequencyDecompositionPanel, 'numeric');
            app.WindowLengthEditField.Position = [291 418 100 22];

            % Create OverlapRatioEditFieldLabel
            app.OverlapRatioEditFieldLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.OverlapRatioEditFieldLabel.HorizontalAlignment = 'right';
            app.OverlapRatioEditFieldLabel.Position = [93 386 79 22];
            app.OverlapRatioEditFieldLabel.Text = 'Overlap Ratio';

            % Create OverlapRatioEditField
            app.OverlapRatioEditField = uieditfield(app.TimeFrequencyDecompositionPanel, 'numeric');
            app.OverlapRatioEditField.Position = [291 386 100 22];
            app.OverlapRatioEditField.Value = 50;

            % Create Label
            app.Label = uilabel(app.TimeFrequencyDecompositionPanel);
            app.Label.Position = [397 386 14 22];
            app.Label.Text = '%';

            % Create NFFTWindowLengthEditFieldLabel
            app.NFFTWindowLengthEditFieldLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.NFFTWindowLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.NFFTWindowLengthEditFieldLabel.Position = [73 353 134 22];
            app.NFFTWindowLengthEditFieldLabel.Text = 'NFFT(>Window Length)';

            % Create NFFTWindowLengthEditField
            app.NFFTWindowLengthEditField = uieditfield(app.TimeFrequencyDecompositionPanel, 'numeric');
            app.NFFTWindowLengthEditField.Position = [291 353 100 22];

            % Create WaveletFamliyDropDownLabel
            app.WaveletFamliyDropDownLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.WaveletFamliyDropDownLabel.HorizontalAlignment = 'right';
            app.WaveletFamliyDropDownLabel.Position = [93 239 87 22];
            app.WaveletFamliyDropDownLabel.Text = 'Wavelet Famliy';

            % Create WaveletFamliyDropDown
            app.WaveletFamliyDropDown = uidropdown(app.TimeFrequencyDecompositionPanel);
            app.WaveletFamliyDropDown.Items = {'morl', 'mesh', 'cmor'};
            app.WaveletFamliyDropDown.Position = [289 239 100 22];
            app.WaveletFamliyDropDown.Value = 'morl';

            % Create ScalesEditFieldLabel
            app.ScalesEditFieldLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.ScalesEditFieldLabel.HorizontalAlignment = 'right';
            app.ScalesEditFieldLabel.Position = [93 200 42 22];
            app.ScalesEditFieldLabel.Text = 'Scales';

            % Create ScalesEditField
            app.ScalesEditField = uieditfield(app.TimeFrequencyDecompositionPanel, 'numeric');
            app.ScalesEditField.Position = [164 200 100 22];
            app.ScalesEditField.Value = 1;

            % Create ScalesEditField_2
            app.ScalesEditField_2 = uieditfield(app.TimeFrequencyDecompositionPanel, 'numeric');
            app.ScalesEditField_2.Position = [289 200 100 22];
            app.ScalesEditField_2.Value = 80;

            % Create toLabel
            app.toLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.toLabel.Position = [275 200 25 22];
            app.toLabel.Text = 'to';

            % Create HzLabel_4
            app.HzLabel_4 = uilabel(app.TimeFrequencyDecompositionPanel);
            app.HzLabel_4.Position = [397 200 25 22];
            app.HzLabel_4.Text = 'Hz';

            % Create VoicesperOctaveEditFieldLabel
            app.VoicesperOctaveEditFieldLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.VoicesperOctaveEditFieldLabel.HorizontalAlignment = 'right';
            app.VoicesperOctaveEditFieldLabel.Position = [86 157 103 22];
            app.VoicesperOctaveEditFieldLabel.Text = 'Voices per Octave';

            % Create VoicesperOctaveEditField
            app.VoicesperOctaveEditField = uieditfield(app.TimeFrequencyDecompositionPanel, 'numeric');
            app.VoicesperOctaveEditField.Position = [289 157 100 22];
            app.VoicesperOctaveEditField.Value = 32;

            % Create DSOKButton
            app.DSOKButton = uibutton(app.TimeFrequencyDecompositionPanel, 'push');
            app.DSOKButton.ButtonPushedFcn = createCallbackFcn(app, @DSOKButtonPushed, true);
            app.DSOKButton.Position = [289 27 100 22];
            app.DSOKButton.Text = 'OK';

            % Create ResamplingLabel
            app.ResamplingLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.ResamplingLabel.Position = [23 93 69 22];
            app.ResamplingLabel.Text = 'Resampling';

            % Create DownsamplingFacorEditFieldLabel
            app.DownsamplingFacorEditFieldLabel = uilabel(app.TimeFrequencyDecompositionPanel);
            app.DownsamplingFacorEditFieldLabel.HorizontalAlignment = 'right';
            app.DownsamplingFacorEditFieldLabel.Position = [76 58 122 22];
            app.DownsamplingFacorEditFieldLabel.Text = ' Downsampling Facor';

            % Create DownsamplingFacorEditField
            app.DownsamplingFacorEditField = uieditfield(app.TimeFrequencyDecompositionPanel, 'numeric');
            app.DownsamplingFacorEditField.Position = [289 58 100 22];

            % Create WTOKButton
            app.WTOKButton = uibutton(app.TimeFrequencyDecompositionPanel, 'push');
            app.WTOKButton.ButtonPushedFcn = createCallbackFcn(app, @WTOKButtonPushed, true);
            app.WTOKButton.Position = [289 123 100 22];
            app.WTOKButton.Text = 'OK';

            % Create STFTOKButton
            app.STFTOKButton = uibutton(app.TimeFrequencyDecompositionPanel, 'push');
            app.STFTOKButton.ButtonPushedFcn = createCallbackFcn(app, @STFTOKButtonPushed, true);
            app.STFTOKButton.Position = [291 322 100 22];
            app.STFTOKButton.Text = 'OK';

            % Create RereferencingPanel
            app.RereferencingPanel = uipanel(app.UIFigure);
            app.RereferencingPanel.Title = 'Re-referencing';
            app.RereferencingPanel.FontWeight = 'bold';
            app.RereferencingPanel.Position = [468 -572 439 531];

            % Create ChooseparamsEditFieldLabel
            app.ChooseparamsEditFieldLabel = uilabel(app.RereferencingPanel);
            app.ChooseparamsEditFieldLabel.HorizontalAlignment = 'right';
            app.ChooseparamsEditFieldLabel.Position = [49 439 94 22];
            app.ChooseparamsEditFieldLabel.Text = 'Choose params ';

            % Create ChooseparamsEditField
            app.ChooseparamsEditField = uieditfield(app.RereferencingPanel, 'text');
            app.ChooseparamsEditField.Position = [179 439 212 22];

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.RereferencingPanel);
            app.MethodDropDownLabel.HorizontalAlignment = 'right';
            app.MethodDropDownLabel.Position = [106 224 46 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create MethodDropDown
            app.MethodDropDown = uidropdown(app.RereferencingPanel);
            app.MethodDropDown.Items = {'direct', 'median', 'trimmed'};
            app.MethodDropDown.Position = [257 224 100 22];
            app.MethodDropDown.Value = 'direct';

            % Create rereferenceOKButton
            app.rereferenceOKButton = uibutton(app.RereferencingPanel, 'push');
            app.rereferenceOKButton.ButtonPushedFcn = createCallbackFcn(app, @rereferenceOKButtonPushed, true);
            app.rereferenceOKButton.Position = [257 185 100 22];
            app.rereferenceOKButton.Text = 'OK';

            % Create inputexamplesPanel
            app.inputexamplesPanel = uipanel(app.RereferencingPanel);
            app.inputexamplesPanel.Title = 'input examples';
            app.inputexamplesPanel.Position = [19 290 381 126];

            % Create A1Label
            app.A1Label = uilabel(app.inputexamplesPanel);
            app.A1Label.FontName = 'Segoe UI Semilight';
            app.A1Label.Position = [226 75 112 22];
            app.A1Label.Text = ''' A1 ''';

            % Create A1A2Label
            app.A1A2Label = uilabel(app.inputexamplesPanel);
            app.A1A2Label.FontName = 'Segoe UI Semilight';
            app.A1A2Label.Position = [226 51 130 22];
            app.A1A2Label.Text = '{ '' A1 '' , '' A2 '' }';

            % Create AverageLabel
            app.AverageLabel = uilabel(app.inputexamplesPanel);
            app.AverageLabel.FontName = 'Segoe UI Semilight';
            app.AverageLabel.Position = [226 29 95 22];
            app.AverageLabel.Text = ''' Average ''';

            % Create orA1A2Label
            app.orA1A2Label = uilabel(app.inputexamplesPanel);
            app.orA1A2Label.FontName = 'Segoe UI Semilight';
            app.orA1A2Label.Position = [226 8 146 22];
            app.orA1A2Label.Text = '[ 1 , 3 , 5 ] or { '' A1 '' , '' A2 '' }';

            % Create DualelectrodereferenceLabel
            app.DualelectrodereferenceLabel = uilabel(app.inputexamplesPanel);
            app.DualelectrodereferenceLabel.FontName = 'Segoe UI Semilight';
            app.DualelectrodereferenceLabel.Position = [18 52 138 22];
            app.DualelectrodereferenceLabel.Text = 'Dual-electrode reference';

            % Create SingleelectrodereferenceLabel
            app.SingleelectrodereferenceLabel = uilabel(app.inputexamplesPanel);
            app.SingleelectrodereferenceLabel.FontName = 'Segoe UI Semilight';
            app.SingleelectrodereferenceLabel.Position = [18 75 146 22];
            app.SingleelectrodereferenceLabel.Text = 'Single-electrode reference';

            % Create AveragereferenceLabel
            app.AveragereferenceLabel = uilabel(app.inputexamplesPanel);
            app.AveragereferenceLabel.FontName = 'Segoe UI Semilight';
            app.AveragereferenceLabel.Position = [18 30 104 22];
            app.AveragereferenceLabel.Text = 'Average reference';

            % Create AverageofspecificchannelgroupsLabel
            app.AverageofspecificchannelgroupsLabel = uilabel(app.inputexamplesPanel);
            app.AverageofspecificchannelgroupsLabel.FontName = 'Segoe UI Semilight';
            app.AverageofspecificchannelgroupsLabel.Position = [18 8 193 22];
            app.AverageofspecificchannelgroupsLabel.Text = 'Average of specific channel groups';

            % Create RemovethereferencechannelCheckBox
            app.RemovethereferencechannelCheckBox = uicheckbox(app.RereferencingPanel);
            app.RemovethereferencechannelCheckBox.Text = 'Remove the reference channel';
            app.RemovethereferencechannelCheckBox.Position = [205 254 186 22];

            % Create InterpolatePanel
            app.InterpolatePanel = uipanel(app.UIFigure);
            app.InterpolatePanel.Title = 'Interpolate';
            app.InterpolatePanel.FontWeight = 'bold';
            app.InterpolatePanel.Position = [952 -572 439 531];

            % Create BadchannelsEditFieldLabel
            app.BadchannelsEditFieldLabel = uilabel(app.InterpolatePanel);
            app.BadchannelsEditFieldLabel.HorizontalAlignment = 'right';
            app.BadchannelsEditFieldLabel.Position = [85 447 78 22];
            app.BadchannelsEditFieldLabel.Text = 'Bad channels';

            % Create BadchannelsEditField
            app.BadchannelsEditField = uieditfield(app.InterpolatePanel, 'text');
            app.BadchannelsEditField.Position = [269 447 100 22];

            % Create interpolateOKButton
            app.interpolateOKButton = uibutton(app.InterpolatePanel, 'push');
            app.interpolateOKButton.ButtonPushedFcn = createCallbackFcn(app, @interpolateOKButtonPushed, true);
            app.interpolateOKButton.Position = [269 239 100 22];
            app.interpolateOKButton.Text = 'OK';

            % Create MethodDropDown_2Label
            app.MethodDropDown_2Label = uilabel(app.InterpolatePanel);
            app.MethodDropDown_2Label.HorizontalAlignment = 'right';
            app.MethodDropDown_2Label.Position = [104 407 46 22];
            app.MethodDropDown_2Label.Text = 'Method';

            % Create InterpolateMethodDropDown
            app.InterpolateMethodDropDown = uidropdown(app.InterpolatePanel);
            app.InterpolateMethodDropDown.Items = {'spherical-spline', 'idw', 'nearest'};
            app.InterpolateMethodDropDown.Position = [244 407 149 22];
            app.InterpolateMethodDropDown.Value = 'spherical-spline';

            % Create UsetemplatemontageCheckBox
            app.UsetemplatemontageCheckBox = uicheckbox(app.InterpolatePanel);
            app.UsetemplatemontageCheckBox.Text = 'Use template montage';
            app.UsetemplatemontageCheckBox.Position = [91 290 143 22];

            % Create InterpKEditFieldLabel
            app.InterpKEditFieldLabel = uilabel(app.InterpolatePanel);
            app.InterpKEditFieldLabel.HorizontalAlignment = 'right';
            app.InterpKEditFieldLabel.Position = [106 370 44 22];
            app.InterpKEditFieldLabel.Text = 'InterpK';

            % Create InterpKEditField
            app.InterpKEditField = uieditfield(app.InterpolatePanel, 'numeric');
            app.InterpKEditField.Position = [269 370 100 22];
            app.InterpKEditField.Value = 4;

            % Create InterpRadiusEditFieldLabel
            app.InterpRadiusEditFieldLabel = uilabel(app.InterpolatePanel);
            app.InterpRadiusEditFieldLabel.HorizontalAlignment = 'right';
            app.InterpRadiusEditFieldLabel.Position = [91 332 74 22];
            app.InterpRadiusEditFieldLabel.Text = 'InterpRadius';

            % Create InterpRadiusEditField
            app.InterpRadiusEditField = uieditfield(app.InterpolatePanel, 'numeric');
            app.InterpRadiusEditField.Position = [269 332 100 22];
            app.InterpRadiusEditField.Value = 25;

            % Create LaplacianPanel
            app.LaplacianPanel = uipanel(app.UIFigure);
            app.LaplacianPanel.Title = 'Laplacian';
            app.LaplacianPanel.FontWeight = 'bold';
            app.LaplacianPanel.Position = [1158 1 439 531];

            % Create LapOKButten
            app.LapOKButten = uibutton(app.LaplacianPanel, 'push');
            app.LapOKButten.ButtonPushedFcn = createCallbackFcn(app, @LapOKButtenButtonPushed, true);
            app.LapOKButten.Position = [269 292 100 22];
            app.LapOKButten.Text = 'OK';

            % Create EdgehandlingLabel
            app.EdgehandlingLabel = uilabel(app.LaplacianPanel);
            app.EdgehandlingLabel.HorizontalAlignment = 'right';
            app.EdgehandlingLabel.Position = [89 375 82 22];
            app.EdgehandlingLabel.Text = 'Edge handling';

            % Create LapEdgeDropdown
            app.LapEdgeDropdown = uidropdown(app.LaplacianPanel);
            app.LapEdgeDropdown.Items = {'spherical-spline', 'idw', 'nearest'};
            app.LapEdgeDropdown.Position = [244 375 149 22];
            app.LapEdgeDropdown.Value = 'spherical-spline';

            % Create LapUseCSDCheckBox
            app.LapUseCSDCheckBox = uicheckbox(app.LaplacianPanel);
            app.LapUseCSDCheckBox.Text = 'Use CSD toolbox if available';
            app.LapUseCSDCheckBox.Position = [63 327 174 22];

            % Create LambdaEditFieldLabel
            app.LambdaEditFieldLabel = uilabel(app.LaplacianPanel);
            app.LambdaEditFieldLabel.HorizontalAlignment = 'right';
            app.LambdaEditFieldLabel.Position = [101 417 49 22];
            app.LambdaEditFieldLabel.Text = 'Lambda';

            % Create LapLambdaEditField
            app.LapLambdaEditField = uieditfield(app.LaplacianPanel, 'numeric');
            app.LapLambdaEditField.Position = [269 417 100 22];
            app.LapLambdaEditField.Value = 1e-05;

            % Create SphericalsplineordermEditField_2Label
            app.SphericalsplineordermEditField_2Label = uilabel(app.LaplacianPanel);
            app.SphericalsplineordermEditField_2Label.HorizontalAlignment = 'right';
            app.SphericalsplineordermEditField_2Label.Position = [61 457 135 22];
            app.SphericalsplineordermEditField_2Label.Text = 'Spherical spline order m';

            % Create LapMEditField
            app.LapMEditField = uieditfield(app.LaplacianPanel, 'numeric');
            app.LapMEditField.Position = [269 457 100 22];
            app.LapMEditField.Value = 4;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SEAL_Preprocess(varargin)

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