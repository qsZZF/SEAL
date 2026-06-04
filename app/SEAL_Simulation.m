classdef SEAL_Simulation < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        LoadRealisticNoiseButton     matlab.ui.control.Button
        StandardAnswerButton         matlab.ui.control.Button
        OKButton                     matlab.ui.control.Button
        SourceListPanel              matlab.ui.container.Panel
        SourceTable                  matlab.ui.control.Table
        SensorLevelSettingPanel      matlab.ui.container.Panel
        SetBackgroundActivityButton  matlab.ui.control.Button
        SimulateNoiseButton          matlab.ui.control.Button
        SourceLevelSettingPanel      matlab.ui.container.Panel
        AddSourceButton              matlab.ui.control.Button
        SourcePatchGenerationPanel   matlab.ui.container.Panel
        SeedVoxelsEditField          matlab.ui.control.NumericEditField
        SeedVoxelsEditField_2Label   matlab.ui.control.Label
        PowerDecayEditField          matlab.ui.control.EditField
        PowerDecayEditFieldLabel     matlab.ui.control.Label
        DecayModelDropDown           matlab.ui.control.DropDown
        DecayModelDropDownLabel      matlab.ui.control.Label
        PatchAreamm2EditField        matlab.ui.control.EditField
        PatchAreamm2EditFieldLabel   matlab.ui.control.Label
        SelectfromcortexButton       matlab.ui.control.Button
        SelectfromAtlasButton        matlab.ui.control.Button
        LoadButton                   matlab.ui.control.Button
        RandButton                   matlab.ui.control.Button
        SourceTimeSeriesPanel        matlab.ui.container.Panel
        SamplingRateEditField        matlab.ui.control.NumericEditField
        SamplingRateEditFieldLabel   matlab.ui.control.Label
        PreviewButton                matlab.ui.control.Button
        AddCheckBox                  matlab.ui.control.CheckBox
        Pink_noiseCheckBox           matlab.ui.control.CheckBox
        ARCheckBox                   matlab.ui.control.CheckBox
        Rest_2CheckBox               matlab.ui.control.CheckBox
        Rest_sinCheckBox             matlab.ui.control.CheckBox
        NMM_spikeCheckBox            matlab.ui.control.CheckBox
        Shap_spikeCheckBox           matlab.ui.control.CheckBox
        Gauss_PulseCheckBox          matlab.ui.control.CheckBox
        ERP_GaborCheckBox            matlab.ui.control.CheckBox
        ERSPCheckBox                 matlab.ui.control.CheckBox
        ERPCheckBox                  matlab.ui.control.CheckBox
        TimeRangeEditField           matlab.ui.control.EditField
        TimeRangeEditFieldLabel      matlab.ui.control.Label
        SourceListMenu               matlab.ui.container.ContextMenu
        deleteMenu                   matlab.ui.container.Menu
    end

    
    properties (Access = private)
        Callingapp
        ProtocolNode
        L_fixed             %leadfieldmatrix
        NumOrientations     %方向数1/3
        numChannels         %通道数
        numSources          %等于Vertices
        PreviewSignal = []       %生成信号预览
        TimeRange           %生成信号范围
        srate               %生成波形采样率   
        freq                %生成波形频率
        simulatedEEG        %生成的纯净的EEG信号
        Noise            %生成的带噪声的EEG信号
        CleanEEG_cache = []   % computeCleanEEG 的缓存结果 (不含 V_background)
        s_real                  %用于标准显示
        SensorSNR_dB     = []    % SimulateNoise 使用的 SNR
        BackgroundSNR_dB = []    % SetBackgroundActivity 使用的 SNR
        activeIndices     =[]      %实际激活点

        Params_Rest_sin
        Params_ERP
        Params_Gauss_Pulse
        Params_Shap_spike
        Params_NMM_spike
        Params_ERSP
        Params_ERP_Gabor
        Params_Rest_2
        Params_AR
        Params_Pink_noise

        V_background

        SourceData = struct('Vertices', {}, 'Waveform', {}, 'SourceName', {});%源活动列表
    end
    
    properties (Access = public)
        handles1            %Cortex图形句柄
        selectedVertices = []
        
        selectedVerticesNum
    end

    methods (Access = private)

        function enforceSingleSelection(app, activeBox)
            % 1. 把你所有的 10 个刺激复选框打包成一个数组
            % 注意：这里的名字必须和你 App Designer 里组件的名字一模一样
            allBoxes = [app.ERPCheckBox, app.ERSPCheckBox, app.Rest_sinCheckBox, ...
                app.Rest_2CheckBox, app.ARCheckBox, app.ERP_GaborCheckBox, ...
                app.Gauss_PulseCheckBox, app.Shap_spikeCheckBox, ...
                app.NMM_spikeCheckBox, app.Pink_noiseCheckBox];

            % 2. 遍历这个数组，除了当前激活的这个，其他的全部设为 false
            for i = 1:length(allBoxes)
                if allBoxes(i) ~= activeBox
                    allBoxes(i).Value = false;
                end
            end
        end

        function clearAllSelections(app)
            % 1. 把你所有的 10 个刺激复选框打包成一个数组
            allBoxes = [app.ERPCheckBox, app.ERSPCheckBox, app.Rest_sinCheckBox, ...
                app.Rest_2CheckBox, app.ARCheckBox, app.ERP_GaborCheckBox, ...
                app.Gauss_PulseCheckBox, app.Shap_spikeCheckBox, ...
                app.NMM_spikeCheckBox, app.Pink_noiseCheckBox];

            % 2. 遍历这个数组，无差别地把所有选项设为 false
            for i = 1:length(allBoxes)
                allBoxes(i).Value = false;
            end

            % 3. (极其重要) 清空底层的参数和预览信号
            app.PreviewSignal = [];
            % 如果有预览图，在这里清空坐标轴
            if isprop(app, 'UIAxes')
                cla(app.UIAxes);
                title(app.UIAxes, 'Waiting for signal generation...');
            end
        end

        function updatePreviewSignal(app)

            t = app.TimeRange(1) : 1/app.srate : app.TimeRange(2) - 1/app.srate;
            totalSignal = zeros(size(t));
            isAddMode = app.AddCheckBox.Value;

            paramList = buildWaveformParamsList(app);
            for k = 1:length(paramList)
                s = computeWaveformFromParams(app, paramList{k}, t) .* (t >= 0);
                if isAddMode
                    totalSignal = totalSignal + s;
                else
                    totalSignal = s;     % 非 Add 模式只会有一个被勾选
                end
            end

            % 不再做 *1e-9 缩放：amplitude 就是 seal_generate_waveform 的总增益
            app.PreviewSignal = totalSignal;

        end

         function paramList = buildWaveformParamsList(app)
            % 把当前所有勾选的复选框打包成 waveform_params 的 cell 数组,
            % 字段命名严格对齐 seal_generate_waveform 的要求。
            %
            % 【单位约定】
            %   UI 输入: 幅值 nA·m (1 nA·m = 1e-9 A·m)
            %   打包输出: p.amplitude 以 A·m 为单位 (SI 单位)
            %   这样 L(V/Am) × S(Am) = V(伏特), 物理单位链闭合
            NAM_TO_AM = 1e-9;    % nA·m -> A·m
 
            paramList = {};
 
            if app.Rest_sinCheckBox.Value
                p.waveformtype = 'rest_sin_custom';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_Rest_sin.Amp * NAM_TO_AM;
                p.Freq         = app.Params_Rest_sin.Freq;
                paramList{end+1} = p; clear p;
            end
 
            if app.ERPCheckBox.Value
                p.waveformtype = 'erp';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_ERP.Amp * NAM_TO_AM;
                p.components   = struct( ...
                    'type','gaussian', ...
                    'latency_s', app.Params_ERP.Lat, ...
                    'amplitude_rel', 1.0, ...
                    'width_s', app.Params_ERP.Wid);
                paramList{end+1} = p; clear p;
            end
 
            if app.Gauss_PulseCheckBox.Value
                p.waveformtype        = 'gaussian_pulse';
                p.samplingRate        = app.srate;
                p.amplitude           = app.Params_Gauss_Pulse.Amp * NAM_TO_AM;
                p.pulse_center_time_s = app.Params_Gauss_Pulse.Lat;
                p.pulse_width_time_s  = app.Params_Gauss_Pulse.Wid;
                paramList{end+1} = p; clear p;
            end
 
            if app.Shap_spikeCheckBox.Value
                p.waveformtype = 'spike';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_Shap_spike.Amp * NAM_TO_AM;
                p.latency_s    = app.Params_Shap_spike.Lat;
                p.shape        = 'exponential_decay';
                p.duration_s   = 0.03;
                p.decay_tau_s  = 0.005;
                paramList{end+1} = p; clear p;
            end
 
            if app.NMM_spikeCheckBox.Value
                p.waveformtype = 'nmm_spike_custom';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_NMM_spike.Amp * NAM_TO_AM;
                p.Lat          = app.Params_NMM_spike.Lat;
                p.Tau          = app.Params_NMM_spike.Tau;
                paramList{end+1} = p; clear p;
            end
 
            if app.ERSPCheckBox.Value
                p.waveformtype = 'erp';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_ERSP.Amp * NAM_TO_AM;
                p.components   = struct( ...
                    'type','gabor', ...
                    'latency_s', app.Params_ERSP.Lat, ...
                    'amplitude_rel', 1.0, ...
                    'frequency_Hz', app.Params_ERSP.Freq, ...
                    'std_dev_s', app.Params_ERSP.Wid, ...
                    'phase_rad', 0);
                paramList{end+1} = p; clear p;
            end
 
            if app.ERP_GaborCheckBox.Value
                p.waveformtype = 'erp';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_ERP_Gabor.Amp * NAM_TO_AM;
                p.components   = struct( ...
                    'type','gabor', ...
                    'latency_s', app.Params_ERP_Gabor.Lat, ...
                    'amplitude_rel', 1.0, ...
                    'frequency_Hz', app.Params_ERP_Gabor.Freq, ...
                    'std_dev_s', app.Params_ERP_Gabor.Wid, ...
                    'phase_rad', 0);
                paramList{end+1} = p; clear p;
            end
 
            if app.Rest_2CheckBox.Value
                p.waveformtype = 'rest_2_custom';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_Rest_2.Amp * NAM_TO_AM;
                p.Freq1        = app.Params_Rest_2.Freq1;
                p.Freq2        = app.Params_Rest_2.Freq2;
                paramList{end+1} = p; clear p;
            end
 
            if app.ARCheckBox.Value
                p.waveformtype        = 'ar';
                p.samplingRate        = app.srate;
                p.amplitude           = app.Params_AR.Amp * NAM_TO_AM;
                p.ar_coeffs           = app.Params_AR.Coeff;
                p.innovation_variance = 1.0;
                paramList{end+1} = p; clear p;
            end
 
            if app.Pink_noiseCheckBox.Value
                p.waveformtype = 'resting_pinknoise';
                p.samplingRate = app.srate;
                p.amplitude    = app.Params_Pink_noise.Amp * NAM_TO_AM;
                paramList{end+1} = p; clear p;
            end
        end

        function w = computeWaveformFromParams(~, p, t)
            % seal_generate_waveform 不原生支持的三种波形在这里手算，
            % 其余统统转交给 seal_generate_waveform。amplitude 字段已经是总增益。
            switch p.waveformtype
                case 'rest_sin_custom'
                    w = p.amplitude * sin(2*pi*p.Freq*t);
                case 'nmm_spike_custom'
                    ts = t - p.Lat;
                    w  = zeros(size(t));
                    idx = ts >= 0;
                    w(idx) = p.amplitude * (ts(idx)/p.Tau) .* exp(1 - ts(idx)/p.Tau);
                case 'rest_2_custom'
                    w = (p.amplitude/2) * sin(2*pi*p.Freq1*t) + ...
                        (p.amplitude/2) * sin(2*pi*p.Freq2*t);
                otherwise
                    w = seal_generate_waveform(t, p.samplingRate, p.waveformtype, p);
            end
        end

        function VC = computeVertConnFromFaces(~, Faces, nVerts)
            % 当 cortex 没有 VertConn 时的兜底（seal_generate_source_activity 内部需要）
            rows = [Faces(:,1); Faces(:,2); Faces(:,2); Faces(:,3); Faces(:,3); Faces(:,1)];
            cols = [Faces(:,2); Faces(:,1); Faces(:,3); Faces(:,2); Faces(:,1); Faces(:,3)];
            VC = sparse(double(rows), double(cols), true, nVerts, nVerts);
        end

        function sourceWindowCloseCallback(app, src, ~)
            % src 是即将关闭的 figure 句柄

            % 1. 检查是否有选点，可以在这里做最后的确认或记录
            if ~isempty(app.selectedVertices)
                fprintf('窗口关闭，最终选中了 %d 个顶点。\n', length(app.selectedVertices));
                % 如果你需要把索引同步到主界面的某个特定地方，在这里写
                % app.SomeOtherField.Value = ...
            else
                disp('窗口关闭，未选中任何顶点。');
            end

            % 2. 恢复主界面的鼠标指针（如果之前改成了 crosshair）
            set(app.UIFigure, 'Pointer', 'arrow');

            % 3. 清理 app.handles1 中的引用，防止“无效句柄”报错
            % 这样主程序就知道 3D 窗口已经不在了
            app.handles1 = [];

            % 4. 必须手动执行 delete(src)，否则窗口关不掉
            delete(src);

            disp('3D 交互窗口已安全关闭。');
        end

      function clean_EEG = computeCleanEEG(app)
            if isempty(app.SourceData)
                clean_EEG = [];
                return;
            end
 
            % --- 1. 命中缓存就直接用 ---
            if ~isempty(app.CleanEEG_cache)
                clean_EEG =  app.CleanEEG_cache;
                return
            else
                % --- 2. 未命中: 按配方算一遍 ---
                cortex = app.ProtocolNode.cortexNode.data;
                if ~isfield(cortex, 'VertConn') || isempty(cortex.VertConn)
                    cortex.VertConn = computeVertConnFromFaces(app, cortex.Faces, ...
                        size(cortex.Vertices,1));
                end
 
                t = app.TimeRange(1) : 1/app.srate : app.TimeRange(2) - 1/app.srate;
                n_sources = size(cortex.Vertices, 1);
                S_total   = zeros(n_sources, length(t));
                allActiveIdx = {};   % 累积所有 ROI 的激活索引 (顺带修 #3 bug)
 
                for srcIdx = 1:length(app.SourceData)
                    src = app.SourceData(srcIdx);
 
                    patchArea_m2   = src.PatchArea_mm2 * 1e-6;
                    patchRadius_mm = sqrt(max(src.PatchArea_mm2, eps) / pi);
                    decay_m        = (patchRadius_mm * max(src.PowerDecay_pct, eps) / 100) * 1e-3;
 
                    % ✅ 统一使用 DecayModel 字段,不再用 useGaussian 逻辑开关
                    decay_model = lower(src.DecayModel);   % 'gaussian' | 'linear' | 'flat' | 'exponential'
 
                    for k = 1:length(src.ParamList)
                        wp = src.ParamList{k};
 
                        if app.isCustomFallback(wp.waveformtype)
                            [S_part, idx_cell] = runFallbackROIBatch(app, cortex, wp, t, ...
                                src.SeedVertices, patchArea_m2, decay_m, decay_model);
                        else
                            ROIs = [];
                            for s = 1:length(src.SeedVertices)
                                roi = struct();
                                roi.seedvox         = src.SeedVertices(s);
                                roi.extent          = patchArea_m2;
                                roi.waveform_params = wp;
                                % ✅ 透传衰减模型和特征距离到引擎
                                roi.waveform_params.decay_model = decay_model;
                                if ~strcmp(decay_model, 'flat')
                                    roi.waveform_params.decay = decay_m;
                                end
                                ROIs = [ROIs, roi]; %#ok<AGROW>
                            end
                            [S_part, idx_cell] = seal_generate_source_activity(cortex, ROIs, t);
                        end
 
                        S_total = S_total + S_part;
                        allActiveIdx = [allActiveIdx, idx_cell]; %#ok<AGROW>
                    end
                end
 
                app.s_real = S_total;
                app.activeIndices = allActiveIdx;
                app.CleanEEG_cache = app.L_fixed * S_total;
                clean_EEG = app.CleanEEG_cache;
            end         
        end

        function tf = isCustomFallback(app, typeStr)
            tf = ismember(typeStr, {'rest_sin_custom','rest_2_custom','nmm_spike_custom'});
        end

      function [S_part, idx_cell] = runFallbackROIBatch(app, cortex, wp, t, ...
                 seedVerts, patchArea_m2, decay_m, decay_model)
            % seal_generate_waveform 不支持的三种类型走这里。
            % 复用 seal_generate_source_activity 的 PatchGenerate (通过零幅值 dummy ROI),
            % 然后用我们手算的波形 + 与引擎严格一致的空间权重公式
            n_sources = size(cortex.Vertices, 1);
            S_part    = zeros(n_sources, length(t));
            idx_cell  = cell(1, length(seedVerts));
 
            waveform  = computeWaveformFromParams(app, wp, t);   % 1×T
 
            dummyWP = struct('waveformtype','gaussian_pulse', ...
                'samplingRate', wp.samplingRate, ...
                'amplitude', 0);
 
            for s = 1:length(seedVerts)
                seed = seedVerts(s);
                dummyROI = struct('seedvox', seed, ...
                    'extent',  patchArea_m2, ...
                    'waveform_params', dummyWP);
                [~, idx_one] = seal_generate_source_activity(cortex, dummyROI, t);
                patch_idx = idx_one{1};
                idx_cell{s} = patch_idx;
 
                % ✅ 与引擎层严格一致的空间权重计算
                cc = cortex.Vertices(seed, :);
                pc = cortex.Vertices(patch_idx, :);
                d  = sqrt(sum((pc - cc).^2, 2));
 
                switch lower(decay_model)
                    case 'flat'
                        w_spatial = ones(length(patch_idx), 1);
                    case 'gaussian'
                        sigma2 = decay_m^2 / log(2);
                        w_spatial = exp(-d.^2 / sigma2);
                    case 'linear'
                        w_spatial = max(0, 1 - d / decay_m);
                    case 'exponential'
                        tau = decay_m / log(2);
                        w_spatial = exp(-d / tau);
                    otherwise
                        sigma2 = decay_m^2 / log(2);
                        w_spatial = exp(-d.^2 / sigma2);
                end
 
                S_part(patch_idx, :) = S_part(patch_idx, :) + w_spatial * waveform;
            end
        end

        function sanityCheckUnits(app, final_EEG)
            % 打印源与 EEG 的数值量级, 用于判断单位链是否正确
            %
            % 生理合理范围:
            %   源偶极子矩峰值:  1e-10  ~  1e-6   A·m   (即 0.1 nA·m ~ 1000 nA·m)
            %   头皮 EEG 峰值:   1e-6   ~  1e-3   V     (即 1 μV ~ 1 mV)
 
            fprintf('\n╔══════════════════════════════════════════════════════╗\n');
            fprintf('║            单 位 链 自 检 (Unit Sanity Check)        ║\n');
            fprintf('╚══════════════════════════════════════════════════════╝\n');
 
            % ---- 源强度 ----
            if ~isempty(app.s_real)
                src_peak = max(abs(app.s_real(:)));
                fprintf(' 源偶极子峰值: %.3e A·m  (=%.2f nA·m)\n', src_peak, src_peak*1e9);
                if src_peak >= 1e-10 && src_peak <= 1e-6
                    fprintf('   ✅ 在生理合理范围 [0.1, 1000] nA·m 内\n');
                elseif src_peak == 0
                    fprintf('   ⚠️ 源为零 (可能未配置波形)\n');
                elseif src_peak > 1e-6
                    fprintf('   ❌ 偏大 (> 1000 nA·m), 可能忘了 nA·m→A·m 换算\n');
                else
                    fprintf('   ⚠️ 偏小 (< 0.1 nA·m), 信号可能被噪声淹没\n');
                end
            else
                fprintf(' 源偶极子峰值: (无数据)\n');
            end
 
            % ---- EEG ----
            eeg_peak = max(abs(final_EEG(:)));
            fprintf(' EEG 电位峰值: %.3e V    (=%.2f μV)\n', eeg_peak, eeg_peak*1e6);
            if eeg_peak >= 1e-6 && eeg_peak <= 1e-3
                fprintf('   ✅ 在生理合理范围 [1, 1000] μV 内\n');
            elseif eeg_peak == 0
                fprintf('   ⚠️ EEG 为零\n');
            elseif eeg_peak > 1e-3
                fprintf('   ❌ 偏大 (> 1 mV), 检查源幅值或导联场单位\n');
            else
                fprintf('   ⚠️ 偏小 (< 1 μV), 实际采集时会被噪声淹没\n');
            end
 
            % ---- 比值 ----
            if ~isempty(app.s_real) && src_peak > 0
                gain = eeg_peak / src_peak;
                fprintf(' 有效增益 |V|_peak / |S|_peak = %.3e V/(A·m)\n', gain);
                fprintf('   (典型导联场 Gain 元素量级应为 1e-3 ~ 1e2 V/(A·m))\n');
            end
 
            % ---- SNR (如果有) ----
            if ~isempty(app.SensorSNR_dB)
                fprintf(' 传感器噪声 SNR: %.1f dB\n', app.SensorSNR_dB);
            end
            if ~isempty(app.BackgroundSNR_dB)
                fprintf(' 背景活动 SNR:  %.1f dB\n', app.BackgroundSNR_dB);
            end
            fprintf('──────────────────────────────────────────────────────\n\n');
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp, ProtocolNode)
            app.UIFigure.Name='Simulation';
            app.Callingapp = mainapp;
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
            app.ProtocolNode = ProtocolNode;
            nColumns_L=size(ProtocolNode.leadfieldNode.data,2);
            nPhysicalSources = size(ProtocolNode.cortexNode.data.Vertices, 1);
            nd = nColumns_L / nPhysicalSources;
            if nd == 1
                app.L_fixed = seal_reduceLeadfield(ProtocolNode.leadfieldNode.data, 1);
            elseif nd == 3
                app.L_fixed = seal_reduceLeadfield(ProtocolNode.leadfieldNode.data, 3,ProtocolNode.leadfieldNode.leadfieldInfo.normals);
            elseif nd ~= 1 && nd ~= 3
                error('SEAL:Algorithm:DimensionMismatch', ...
                    '导联场维度异常！\n导联场列数 (%d) 必须是皮层节点数 (%d) 的 1 倍 (固定朝向) 或 3 倍 (自由朝向)。\n请检查头模型与导联场是否匹配！', ...
                    nColumns_L, nPhysicalSources);
            end
            app.NumOrientations = nd;
            
            app.numChannels = size(app.L_fixed,1);
            app.numSources = nPhysicalSources;
            timeRangeStr = app.TimeRangeEditField.Value;
            timeRange = str2num(timeRangeStr);

            % 防止用户乱输入，比如输入了 '[A, B]' 或者只输入了一个数字 '[-1]'
            if isempty(timeRange) || length(timeRange) ~= 2
                % 如果解析失败或长度不对，弹出错误提示并终止运行
                uialert(app.UIFigure, '时间范围格式错误！请使用形如 [-1, 1] 的格式。', '输入错误', 'Icon', 'error');
                return; % 终止后续的回调代码执行
            end
            % 确保起始时间小于结束时间
            if timeRange(1) >= timeRange(2)
                uialert(app.UIFigure, '起始时间必须小于结束时间！', '逻辑错误', 'Icon', 'warning');
                return;
            end
            app.TimeRange = timeRange;
            app.srate = app.SamplingRateEditField.Value;
            app.OKButton.Enable = 'off';
            app.StandardAnswerButton.Enable = "off";
        end

        % Value changed function: Rest_sinCheckBox
        function Rest_sinCheckBoxValueChanged(app, event)
            if app.Rest_sinCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    % 传入 app 和当前这个复选框的句柄
                    enforceSingleSelection(app, app.Rest_sinCheckBox);
                end

                % 2. 弹窗获取参数并存入 app 属性
                answer = inputdlg({'频率(Hz):','幅值（nAm）:'}, '设置', [1 50], {'10','5'});
                if ~isempty(answer)
                    app.Params_Rest_sin.Freq = str2double(answer{1});
                    app.Params_Rest_sin.Amp = str2double(answer{2});
                else
                    app.Rest_sinCheckBox.Value = false;
                end
            end

            % 3. 无论勾选还是取消，都让“大脑”重新算一遍
            updatePreviewSignal(app);
        end

        % Button pushed function: PreviewButton
        function PreviewButtonPushed(app, event)

            % 看看生成的源信号是什么样子的
            figure('Name', 'Source Waveform Preview', 'Color', 'w');
            t = app.TimeRange(1) : 1/app.srate : app.TimeRange(2) - 1/app.srate;
            % PreviewSignal 单位是 A·m, 显示时 ×1e9 换算成 nA·m
            plot(t, app.PreviewSignal * 1e9, 'LineWidth', 1.2);
            xlabel('Time (s)');
            ylabel('Amplitude (nA·m)');
            title('Preview: Source Waveform');
            grid on;

        end

        % Button pushed function: RandButton
        function RandButtonPushed(app, event)
            SeedNum = app.SeedVoxelsEditField.Value;
            seedVoxelIdx = randperm(app.numSources, SeedNum);
            app.selectedVertices = seedVoxelIdx;
            disp(['选中的随机顶点 ID: ', num2str(seedVoxelIdx(:)')]);     
        end

        % Button pushed function: SimulateNoiseButton
        function SimulateNoiseButtonPushed(app, event)
            clean_EEG = computeCleanEEG(app);
            if isempty(clean_EEG)
                uialert(app.UIFigure, '请先添加至少一个源！', '缺少参考信号');
                return;
            end

            disp('正在根据纯净信号总功率生成传感器噪声...');

   
            answer = inputdlg('请输入目标信噪比 (dB):', '噪声参数设置', [1 35], {'10'});
            if isempty(answer)
                return;  % 用户取消
            end

            snr_db = str2double(answer{1});
            if isnan(snr_db)
                uialert(app.UIFigure, '请输入有效的数值！', '输入错误');
                return;
            end

            [~, numPoints] = size(clean_EEG);
       

            % 3. 计算功率与缩放
            signal_power = mean(clean_EEG(:).^2);
            noise_power = signal_power / (10^(snr_db / 10));
            noise_std = sqrt(noise_power);


            app.Noise = noise_std * randn(app.numChannels, numPoints);
            app.SensorSNR_dB = snr_db;

            disp(['✅ 噪声生成完毕！目标 SNR: ', num2str(snr_db), ' dB']);

        end

        % Button pushed function: SetBackgroundActivityButton
        function SetBackgroundActivityButtonPushed(app, event)
            clean_EEG = computeCleanEEG(app);
            if isempty(clean_EEG)
                uialert(app.UIFigure, '请先添加至少一个源！', '缺少参考信号');
                return;
            end

            % 弹窗获取背景噪声 SNR
            answer = inputdlg(...
                {'背景脑活动信噪比 (dB):'}, ...
                '背景活动参数设置', [1 40], {'10'});
            if isempty(answer)
                return;
            end

            snr_db = str2double(answer{1});
            
            if isnan(snr_db) 
                uialert(app.UIFigure, '请输入有效的数值！', '输入错误');
                return;
            end

            % --- 1. 排除靶源顶点 ---
            allTargetVertices = [];
            if isprop(app, 'SourceData') && ~isempty(app.SourceData)
                for i = 1:length(app.SourceData)
                    allTargetVertices = [allTargetVertices, app.SourceData(i).SeedVertices(:)'];
                end
            end
            if isprop(app, 'selectedVertices') && ~isempty(app.selectedVertices)
                allTargetVertices = [allTargetVertices, app.selectedVertices(:)'];
            end
            allTargetVertices = unique(allTargetVertices);

            % --- 2. 构建候选池 ---
            availablePool = setdiff(1:app.numSources, allTargetVertices);
            numBgSources = max(round(0.2 * app.numSources), 100);
            numBgSources = min(numBgSources, length(availablePool));

            % --- 3. 随机采样背景源 ---
            randomIndex = randperm(length(availablePool), numBgSources);
            bgIndices = availablePool(randomIndex);
            L_bg = app.L_fixed(:, bgIndices);

            t = app.TimeRange(1) : 1/app.srate : app.TimeRange(2) - 1/app.srate;
            numTimePoints = length(t);
            
            pink_noise = zeros(numBgSources, numTimePoints);
            for k = 1:numBgSources
                pink_noise(k, :) = seal_generate_pink_noise_channel(numTimePoints);
            end
            pink_noise = pink_noise ./ (std(pink_noise, 0, 2) + eps);
            bg_all = L_bg * pink_noise;


            % --- 5. 根据 SNR 计算缩放系数 ---
            signal_power = mean(clean_EEG(:).^2);
            target_noise_power = signal_power / (10^(snr_db / 10));
            current_noise_power = mean(bg_all(:).^2);

            scale = sqrt(target_noise_power / current_noise_power);
            app.V_background = bg_all * scale;
            app.BackgroundSNR_dB = snr_db;  % ← 新增
            disp(['✅ 背景噪声生成完毕！目标 SNR: ', num2str(snr_db), ' dB']);
        end

        % Value changed function: ERPCheckBox
        function ERPCheckBoxValueChanged(app, event)
            value = app.ERPCheckBox.Value;
            if value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.ERPCheckBox);
                end
                answer = inputdlg({'潜伏期 Latency (s):', '波宽 Width (s):', '幅值 (nAm):'}, 'ERP 设置', [1 50], {'0.3', '0.05', '50'});
                if ~isempty(answer)
                    app.Params_ERP.Lat = str2double(answer{1});
                    app.Params_ERP.Wid = str2double(answer{2});
                    app.Params_ERP.Amp = str2double(answer{3});
                else
                    app.ERPCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);
        end

        % Value changed function: Gauss_PulseCheckBox
        function Gauss_PulseCheckBoxValueChanged(app, event)
            if app.Gauss_PulseCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.Gauss_PulseCheckBox);
                end
                answer = inputdlg({'潜伏期 (s):', '波宽 (s):', '幅值 (nAm):'}, '高斯脉冲设置', [1 50], {'0.3', '1000', '5'});
                if ~isempty(answer)
                    app.Params_Gauss_Pulse.Lat = str2double(answer{1});
                    app.Params_Gauss_Pulse.Wid = str2double(answer{2});
                    app.Params_Gauss_Pulse.Amp = str2double(answer{3});
                else
                    app.Gauss_PulseCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);
        end

        % Value changed function: Shap_spikeCheckBox
        function Shap_spikeCheckBoxValueChanged(app, event)
            if app.Shap_spikeCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.Shap_spikeCheckBox);
                end
                answer = inputdlg({'潜伏期 Latency (s):', '幅值 (nAm):'}, '尖峰放电设置', [1 50], {'0.1', '500'});
                if ~isempty(answer)
                    app.Params_Shap_spike.Lat = str2double(answer{1});
                    app.Params_Shap_spike.Amp = str2double(answer{2});
                else
                    app.Shap_spikeCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);
        end

        % Value changed function: NMM_spikeCheckBox
        function NMM_spikeCheckBoxValueChanged(app, event)
            if app.NMM_spikeCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.NMM_spikeCheckBox);
                end
                answer = inputdlg({'潜伏期 (s):', '时间常数 Tau (s):', '幅值 (nAm):'}, 'NMM 脉冲设置', [1 50], {'0.1', '0.02', '100'});
                if ~isempty(answer)
                    app.Params_NMM_spike.Lat = str2double(answer{1});
                    app.Params_NMM_spike.Tau = str2double(answer{2});
                    app.Params_NMM_spike.Amp = str2double(answer{3});
                else
                    app.NMM_spikeCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);
        end

        % Value changed function: ERSPCheckBox
        function ERSPCheckBoxValueChanged(app, event)
            if app.ERSPCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.ERSPCheckBox);
                end
                answer = inputdlg({'中心频率 (Hz):', '爆发时间 (s):', '包络宽度 (s):', '幅值 (nAm):'}, 'ERSP 设置', [1 50], {'40', '0.2', '0.1', '30'});
                if ~isempty(answer)
                    app.Params_ERSP.Freq = str2double(answer{1});
                    app.Params_ERSP.Lat  = str2double(answer{2});
                    app.Params_ERSP.Wid  = str2double(answer{3});
                    app.Params_ERSP.Amp  = str2double(answer{4});
                else
                    app.ERSPCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);

        end

        % Value changed function: ERP_GaborCheckBox
        function ERP_GaborCheckBoxValueChanged(app, event)
            if app.ERP_GaborCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.ERP_GaborCheckBox);
                end
                answer = inputdlg({'中心频率 (Hz):', '爆发时间 (s):', '高斯窗宽 (s):', '幅值 (nAm):'}, 'Gabor 斑块设置', [1 50], {'40', '0.2', '0.1', '30'});
                if ~isempty(answer)
                    app.Params_ERP_Gabor.Freq = str2double(answer{1});
                    app.Params_ERP_Gabor.Lat  = str2double(answer{2});
                    app.Params_ERP_Gabor.Wid  = str2double(answer{3});
                    app.Params_ERP_Gabor.Amp  = str2double(answer{4});
                else
                    app.ERP_GaborCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);

        end

        % Value changed function: Rest_2CheckBox
        function Rest_2CheckBoxValueChanged(app, event)
            if app.Rest_2CheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.Rest_2CheckBox);
                end
                answer = inputdlg({'频率 1 (Hz):', '频率 2 (Hz):', '总幅值 (nAm):'}, '双频静息态设置', [1 50], {'10', '20', '15'});
                if ~isempty(answer)
                    app.Params_Rest_2.Freq1 = str2double(answer{1});
                    app.Params_Rest_2.Freq2 = str2double(answer{2});
                    app.Params_Rest_2.Amp   = str2double(answer{3});
                else
                    app.Rest_2CheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);
        end

        % Value changed function: ARCheckBox
        function ARCheckBoxValueChanged(app, event)
            if app.ARCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.ARCheckBox);
                end
                answer = inputdlg({'AR 平滑系数 (0-1):', '幅值 (nAm):'}, 'AR 噪声设置', [1 50], {'0.9', '10'});
                if ~isempty(answer)
                    app.Params_AR.Coeff = str2double(answer{1});
                    app.Params_AR.Amp   = str2double(answer{2});
                else
                    app.ARCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);
        end

        % Value changed function: Pink_noiseCheckBox
        function Pink_noiseCheckBoxValueChanged(app, event)
            if app.Pink_noiseCheckBox.Value == true
                if app.AddCheckBox.Value == false
                    enforceSingleSelection(app, app.Pink_noiseCheckBox);
                end
                answer = inputdlg({'噪声幅值级别 (nAm):'}, '粉红噪声设置', [1 50], {'2'});
                if ~isempty(answer)
                    app.Params_Pink_noise.Amp = str2double(answer{1});
                else
                    app.Pink_noiseCheckBox.Value = false;
                end
            end
            updatePreviewSignal(app);
        end

        % Button pushed function: SelectfromcortexButton
        function SelectfromcortexButtonPushed(app, event)
            app.handles1=PlotSource([],app.ProtocolNode.cortexNode.data,'smooth',30, ...
                    'Position',[app.UIFigure.Position(1)+app.UIFigure.Position(3) app.UIFigure.Position(2)+app.UIFigure.Position(4)-600 800, 600]);

            seal_bind3DInteraction(app.handles1.h, app.handles1.axes,  @(src,evt) vertexClicked(app,src,evt));

            app.handles1.h.UserData = 1;
            app.handles1.h.CloseRequestFcn = @(src, event) sourceWindowCloseCallback(app, src, event);
            set(app.handles1.h, 'Pointer', 'crosshair');

        end

        % Button pushed function: AddSourceButton
        function AddSourceButtonPushed(app, event)
            % --- 1. 拦截 ---
            if isempty(app.selectedVertices)
                uialert(app.UIFigure, '请先在 3D 大脑上选择至少一个顶点！', '缺少空间信息');
                return;
            end
            paramList = buildWaveformParamsList(app);
            if isempty(paramList)
                uialert(app.UIFigure, '请先勾选并配置至少一种刺激波形！', '缺少时间信息');
                return;
            end

            % --- 2. 拼标签 ---
            selectedWaves = {};
            allBoxes = [app.ERPCheckBox, app.ERSPCheckBox, app.Rest_sinCheckBox, ...
                app.Rest_2CheckBox, app.ARCheckBox, app.ERP_GaborCheckBox, ...
                app.Gauss_PulseCheckBox, app.Shap_spikeCheckBox, ...
                app.NMM_spikeCheckBox, app.Pink_noiseCheckBox];
            for i = 1:length(allBoxes)
                if allBoxes(i).Value, selectedWaves{end+1} = allBoxes(i).Text; end %#ok<AGROW>
            end
            if isempty(selectedWaves), waveType = 'Custom';
            else, waveType = strjoin(selectedWaves, ' + '); end

            % --- 3. 打包配方 (不生成任何信号!) ---
            sourceCount = length(app.SourceData) + 1;
            newSource.No             = sourceCount;
            newSource.WaveType       = waveType;
            newSource.SeedVertices   = app.selectedVertices;
            newSource.PatchArea_mm2  = str2double(app.PatchAreamm2EditField.Value);
            newSource.DecayModel     = app.DecayModelDropDown.Value;
            newSource.PowerDecay_pct = str2double(app.PowerDecayEditField.Value);
            newSource.ParamList      = paramList;
            newSource.NumVertices    = length(app.selectedVertices);

            if isempty(app.SourceData), app.SourceData = newSource;
            else, app.SourceData(sourceCount) = newSource; end

            % --- 4. 刷新 UI (Vertices 列展示种子点数) ---
            displayData = cell(sourceCount, 3);
            for i = 1:sourceCount
                displayData{i,1} = num2str(app.SourceData(i).No);
                displayData{i,2} = num2str(length(app.SourceData(i).SeedVertices));
                displayData{i,3} = app.SourceData(i).WaveType;
            end
            app.SourceTable.Data = displayData;

            app.SamplingRateEditField.Enable = 'off';
            app.TimeRangeEditField.Enable    = 'off';
            app.OKButton.Enable              = 'on';
            app.handles1                     = [];

            app.selectedVertices          = [];
            app.CleanEEG_cache = [];   % 配方变了
            app.V_background =[];
            app.Noise =[];
            app.SeedVoxelsEditField.Value = 0;
            clearAllSelections(app);
        end

        % Menu selected function: deleteMenu
        function deleteMenuSelected(app, event)
            drawnow;
            % 1. 获取当前选中的是哪一行
            % 因为我们把表格设为了 'row' 选择模式，所以 Selection 存的就是行号
            selectedRow = app.SourceTable.Selection;

            % 防呆检查：如果没选中任何行，直接返回
            if isempty(selectedRow)
                return;
            end

            % 2. 友好的二次确认弹窗 (防止手滑)
            sourceName = app.SourceTable.Data{selectedRow, 3}; % 获取要删的刺激类型名字
            confirmMsg = sprintf('确定要删除第 %d 个源 (%s) 吗？', selectedRow, sourceName);
            selection = uiconfirm(app.UIFigure, confirmMsg, '确认删除', ...
                'Options', {'确定删除', '取消'}, 'DefaultOption', 2, 'Icon', 'warning');

            if strcmp(selection, '取消')
                return;
            end

            % --- 3. 核心：从后台数据仓库中精准击杀 ---
            app.SourceData(selectedRow) = []; % 直接把这一行赋值为空，MATLAB会自动缩减数组

            % --- 4. 刷新前端显示与序号重排 ---
            if isempty(app.SourceData)
                % 终极情况：如果所有的源都被删光了，彻底清空表格
                app.SourceTable.Data = cell(0, 3);

                app.SamplingRateEditField.Enable = 'on';
                app.TimeRangeEditField.Enable = 'on';
                app.OKButton.Enable = "off";
                disp('所有源已清空，全局参数已重新解锁。');
            else
                % 常规情况：删除了中间的某个源，需要把后面的源序号往前顶
                nSources = length(app.SourceData);
                displayData = cell(nSources, 3);

                for i = 1:nSources
                    % 1. 修正后台的序号记录
                    app.SourceData(i).No = i;
                    % 2. 重新构建前端展示用的字符串矩阵
                    displayData{i, 1} = num2str(i);
                    displayData{i, 2} = num2str(app.SourceData(i).NumVertices);
                    displayData{i, 3} = app.SourceData(i).WaveType;
                end

                % 推送到表格界面
                app.SourceTable.Data = displayData;
                disp(['已删除指定源，当前剩余 ', num2str(nSources), ' 个源。']);
            end
            app.SourceTable.Selection = [];
            app.CleanEEG_cache = [];            %source变了
            app.V_background =[];
            app.Noise =[];
        end

        % Button pushed function: OKButton
        function OKButtonPushed(app, event)

           % 1. 拿纯净靶源投影
            final_EEG = computeCleanEEG(app);
            if isempty(final_EEG)
                uialert(app.UIFigure, 'Source List 为空，请先配置并 Add Source！', '操作无效');
                return;
            end
 
            % 2. 叠加背景脑活动
            if isprop(app, 'V_background') && ~isempty(app.V_background)
                min_t = min(size(final_EEG, 2), size(app.V_background, 2));
                final_EEG(:, 1:min_t) = final_EEG(:, 1:min_t) + app.V_background(:, 1:min_t);
                disp('✅ 成功注入背景脑活动。');
            else
                disp('ℹ️ 未检测到背景活动。');
            end
 
            % 3. 叠加传感器噪声
            if isprop(app, 'Noise') && ~isempty(app.Noise)
                noise_trial1 = app.Noise(:,:,1);
                min_t = min(size(final_EEG, 2), size(noise_trial1, 2));
                final_EEG(:, 1:min_t) = final_EEG(:, 1:min_t) + noise_trial1(:, 1:min_t);
                disp('✅ 成功注入传感器噪声。');
            else
                disp('ℹ️ 未检测到传感器噪声。');
            end
 
            out = struct();
            out.data       = final_EEG;          % 单位: V (伏特)
            out.srate      = app.SamplingRateEditField.Value;
            out.time_range = app.TimeRange;
            out.nChannels  = app.numChannels;
            out.nSources   = app.numSources;
            out.unit       = 'V';                % 明确标注数据单位
 
            out.SourceData    = app.SourceData;
            out.activeIndices = app.activeIndices;
            out.s_real        = app.s_real;        % ground-truth source time series for temporal metrics

            out.noise.sensor_SNR_dB     = app.SensorSNR_dB;
            out.noise.background_SNR_dB = app.BackgroundSNR_dB;
            out.noise.has_sensor_noise  = ~isempty(app.Noise);
            out.noise.has_background    = ~isempty(app.V_background);
 
            app.simulatedEEG  = out;
 
            app.PreviewSignal = [];
 
            filterSpec = {'*.mat', 'MATLAB'; '*.*', 'All Files'};
            [filename, path] = uiputfile(filterSpec, "Select a location to save results.");
            if isequal(filename, 0)
                disp('用户取消保存。');
                return;
            end
            savePath = fullfile(path, filename);
            saveData(savePath, app.simulatedEEG);
            node = app.ProtocolNode.openDataFromData(savePath);
            node.dataInfo.isSimu = 1;
            node.save();
            app.Callingapp.flushFileTree();
 
            disp('🎉 脑电仿真全部完成！');
            app.Callingapp.bringToFront();
            % ---------- EEG 可视化 (V -> μV) ----------
            figure('Name', 'Final Simulated Scalp EEG', 'Color', 'w');
            t = linspace(app.TimeRange(1), app.TimeRange(2), size(final_EEG, 2));
            plot(t, final_EEG' * 1e6, 'LineWidth', 0.8);   % V -> μV
            set(gca, 'FontSize', 12);
            xlabel('Time (s)', 'FontWeight', 'bold');
            ylabel('Amplitude (μV)', 'FontWeight', 'bold');
            title('Final Simulated Scalp EEG', 'FontWeight', 'bold');
            grid on;
 
            % ============ 快速验证方法: 单位自检 ============
            app.sanityCheckUnits(final_EEG);
 
            app.StandardAnswerButton.Enable = "on";
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
         % 1. 定义当前头模型的最大顶点数 (作为越界检查的合法天花板)
            % 请确保 app.numSources 在导入头模型时已经正确赋值
            if ~isprop(app, 'numSources') || isempty(app.numSources)
                uialert(app.UIFigure, '请先加载头模型 (Cortex/Leadfield)！', '操作被拒绝', 'Icon', 'warning');
                return;
            end
            maxVertices = app.numSources;

            % 2. 弹窗要求用户输入
            prompt = {sprintf('请输入顶点索引 (1 - %d):\n支持空格/逗号分隔，或使用 1:10 语法', maxVertices)};
            dlgtitle = 'Load Vertices';
            dims = [1 50];
            definput = {'1'}; % 默认占位符
            
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            
            % 3. 防线一：用户点击了“取消”或关闭了弹窗
            if isempty(answer)
                return; % 静默退出，不打扰用户
            end
            
            inputStr = answer{1};
            
            % 4. 防线二：将字符串转换为数字数组
            % 极其推荐使用 str2num 而不是 str2double，因为 str2num 可以完美解析 '1:100' 或 '1,2,3'
            idxArray = str2num(inputStr); %#ok<ST2NM> 
            
            if isempty(idxArray)
                uialert(app.UIFigure, '输入格式不合法！\n请输入纯数字，例如: 15, 20, 25 或 10:50', ...
                    '解析失败', 'Icon', 'error');
                return;
            end
            
            % 5. 防线三：小数防御 (顶点 ID 必须是整数)
            if any(floor(idxArray) ~= idxArray)
                uialert(app.UIFigure, '顶点索引不能包含小数！', '输入非法', 'Icon', 'error');
                return;
            end
            
            % 6. 防线四：范围防御 (不能小于 1，不能大于总顶点数)
            if any(idxArray < 1) || any(idxArray > maxVertices)
                errorMsg = sprintf('索引越界！\n输入的顶点必须在 1 到 %d 之间。', maxVertices);
                uialert(app.UIFigure, errorMsg, '越界错误', 'Icon', 'error');
                return;
            end
            
            idxArray = unique(idxArray(:)'); 
            
            app.selectedVertices = idxArray; 
            
            % 可选：在命令窗口打印日志，方便科研调试
            disp(['✅ 成功加载 ', num2str(length(idxArray)), ' 个顶点。']);
            
    
        
        end

        % Button pushed function: SelectfromAtlasButton
        function SelectfromAtlasButtonPushed(app, event)
            if ~isprop(app, 'ProtocolNode') || isempty(app.ProtocolNode.cortexNode.data)
                uialert(app.UIFigure, '请先加载大脑皮层 (Cortex) 模型！', '未检测到模型', 'Icon', 'warning');
                return;
            end
            cortexData = app.ProtocolNode.cortexNode.data;

            % 2. 🛡️ 核心安检：检测是否存在 Atlas 字段且不为空
            if ~isfield(cortexData, 'Atlas') || isempty(cortexData.Atlas)
                uialert(app.UIFigure, '当前加载的皮层模型不包含 Atlas (脑图谱) 数据！\n请确保导入的是包含 Scout 划分的 Brainstorm Surface 文件。', ...
                    '无图谱数据', 'Icon', 'error');
                return;
            end

            atlasStruct = cortexData.Atlas;
            atlasNames = {atlasStruct.Name}; % 提取所有可用图谱的名称 (如 Desikan-Killiany, Destrieux)

            % 3. 图谱选择向导 (如果文件里打包了多个不同精度的图谱)
            if length(atlasNames) > 1
                [atlasIdx, tf] = listdlg('ListString', atlasNames, ...
                    'SelectionMode', 'single', ...
                    'PromptString', '请选择要使用的脑图谱划分标准:', ...
                    'Name', '选择图谱 (Atlas)', ...
                    'ListSize', [300, 150]);
                if ~tf
                    return; % 用户取消
                end
            else
                atlasIdx = 1; % 如果只有一个，直接默认使用
            end

            selectedAtlas = atlasStruct(atlasIdx);

            % 4. 🛡️ 次级安检：检测选中的图谱里是否有脑区 (Scouts)
            if ~isfield(selectedAtlas, 'Scouts') || isempty(selectedAtlas.Scouts)
                uialert(app.UIFigure, ['选中的图谱 "', selectedAtlas.Name, '" 中未划分具体的脑区 (Scouts)！'], ...
                    '图谱为空', 'Icon', 'error');
                return;
            end

            % 5. 脑区选择向导 (核心交互)
            scoutNames = {selectedAtlas.Scouts.Label}; % 提取所有脑区的英文/拼音名字
            [scoutIdx, tf2] = listdlg('ListString', scoutNames, ...
                'SelectionMode', 'multiple', ... % 🌟 允许按住 Ctrl 或 Shift 多选！
                'PromptString', {'请选择目标脑区 (ROI):', '(支持按住 Ctrl 多选合并)'}, ...
                'Name', ['选择脑区 - ' selectedAtlas.Name], ...
                'ListSize', [350, 400]);

            if ~tf2
                return; % 用户取消
            end

            % 6. 提取并合并所有选中脑区的顶点
            allSelectedVertices = [];
            for i = 1:length(scoutIdx)
                idx = scoutIdx(i);
                % Brainstorm 的 Scout.Vertices 默认就是该脑区包含的顶点 ID 数组
                allSelectedVertices = [allSelectedVertices, selectedAtlas.Scouts(idx).Vertices(:)'];
            end

            % 去重并转为行向量 (防止两个脑区在边界上有重叠顶点)
            allSelectedVertices = unique(allSelectedVertices);       
            app.selectedVertices = allSelectedVertices;
            
            % 在控制台打印详细日志，方便科研复查
            disp(['✅ 从 Atlas 成功提取脑区，共计 ', num2str(length(allSelectedVertices)), ' 个顶点。']);
        
        end

        % Value changed function: DecayModelDropDown
        function DecayModelDropDownValueChanged(app, event)
            value = app.DecayModelDropDown.Value;
            if   strcmp(value,'flat')
                app.PowerDecayEditField.Enable = 'off';
            else
                app.PowerDecayEditField.Enable = 'on';
            end
        end

        % Button pushed function: StandardAnswerButton
        function StandardAnswerButtonPushed(app, event)
            
           % 源能量分布图 (sum(s^2, 2) 单位: A²·m²·s, 只用于定位,不关心绝对量纲)
            PlotSource(sum(app.s_real.^2,2), app.ProtocolNode.cortexNode.data, ...
                'smooth', 10, 'thresh', 0.1, 'colormap', jet(256));
 
            t_start = app.TimeRange(1);
            t_end   = app.TimeRange(2) - 1/app.srate;
            t = app.TimeRange(1) : 1/app.srate : app.TimeRange(2) - 1/app.srate;
            
            figure('Name', 'Generated Time Courses at Seed Vertices', ...
                'NumberTitle', 'off', 'Color', 'w');
 
            totalPlots = 0;
            for i = 1:length(app.SourceData)
                totalPlots = totalPlots + length(app.SourceData(i).SeedVertices);
            end
 
            plotIdx = 0;
            for i = 1:length(app.SourceData)
                seeds = app.SourceData(i).SeedVertices;
                for j = 1:length(seeds)
                    plotIdx = plotIdx + 1;
                    subplot(totalPlots, 1, plotIdx);
                    % s_real 单位是 A·m, ×1e9 -> nA·m
                    plot(t, app.s_real(seeds(j), :) * 1e9, 'r', 'LineWidth', 1.5);
                    grid on;
                    xlabel('Time (s)');
                    ylabel('Amplitude (nA·m)');
                    xlim([t_start t_end]);
                    title(sprintf('Source %d — Seed vertex %d (%s)', ...
                        i, seeds(j), app.SourceData(i).WaveType));
                end
            end

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 398 490];
            app.UIFigure.Name = 'MATLAB App';

            % Create SourceLevelSettingPanel
            app.SourceLevelSettingPanel = uipanel(app.UIFigure);
            app.SourceLevelSettingPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.SourceLevelSettingPanel.TitlePosition = 'centertop';
            app.SourceLevelSettingPanel.Title = 'Source Level Setting';
            app.SourceLevelSettingPanel.FontName = 'Arial';
            app.SourceLevelSettingPanel.FontWeight = 'bold';
            app.SourceLevelSettingPanel.Position = [1 195 397 296];

            % Create SourceTimeSeriesPanel
            app.SourceTimeSeriesPanel = uipanel(app.SourceLevelSettingPanel);
            app.SourceTimeSeriesPanel.TitlePosition = 'centertop';
            app.SourceTimeSeriesPanel.Title = 'Source Time Series';
            app.SourceTimeSeriesPanel.FontName = 'Arial';
            app.SourceTimeSeriesPanel.FontWeight = 'bold';
            app.SourceTimeSeriesPanel.Position = [0 0 209 276];

            % Create TimeRangeEditFieldLabel
            app.TimeRangeEditFieldLabel = uilabel(app.SourceTimeSeriesPanel);
            app.TimeRangeEditFieldLabel.Position = [18 176 70 22];
            app.TimeRangeEditFieldLabel.Text = 'Time Range';

            % Create TimeRangeEditField
            app.TimeRangeEditField = uieditfield(app.SourceTimeSeriesPanel, 'text');
            app.TimeRangeEditField.HorizontalAlignment = 'right';
            app.TimeRangeEditField.Position = [115 176 76 22];
            app.TimeRangeEditField.Value = '[-1, 1]';

            % Create ERPCheckBox
            app.ERPCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.ERPCheckBox.ValueChangedFcn = createCallbackFcn(app, @ERPCheckBoxValueChanged, true);
            app.ERPCheckBox.Text = 'ERP';
            app.ERPCheckBox.Position = [22 137 47 22];

            % Create ERSPCheckBox
            app.ERSPCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.ERSPCheckBox.ValueChangedFcn = createCallbackFcn(app, @ERSPCheckBoxValueChanged, true);
            app.ERSPCheckBox.Text = 'ERSP';
            app.ERSPCheckBox.Position = [22 117 55 22];

            % Create ERP_GaborCheckBox
            app.ERP_GaborCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.ERP_GaborCheckBox.ValueChangedFcn = createCallbackFcn(app, @ERP_GaborCheckBoxValueChanged, true);
            app.ERP_GaborCheckBox.Text = 'ERP_Gabor';
            app.ERP_GaborCheckBox.Position = [96 137 87 22];

            % Create Gauss_PulseCheckBox
            app.Gauss_PulseCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.Gauss_PulseCheckBox.ValueChangedFcn = createCallbackFcn(app, @Gauss_PulseCheckBoxValueChanged, true);
            app.Gauss_PulseCheckBox.Text = 'Gauss_Pulse';
            app.Gauss_PulseCheckBox.Position = [96 117 93 22];

            % Create Shap_spikeCheckBox
            app.Shap_spikeCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.Shap_spikeCheckBox.ValueChangedFcn = createCallbackFcn(app, @Shap_spikeCheckBoxValueChanged, true);
            app.Shap_spikeCheckBox.Text = 'Shap_spike';
            app.Shap_spikeCheckBox.Position = [96 97 85 22];

            % Create NMM_spikeCheckBox
            app.NMM_spikeCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.NMM_spikeCheckBox.ValueChangedFcn = createCallbackFcn(app, @NMM_spikeCheckBoxValueChanged, true);
            app.NMM_spikeCheckBox.Text = 'NMM_spike';
            app.NMM_spikeCheckBox.Position = [96 77 85 22];

            % Create Rest_sinCheckBox
            app.Rest_sinCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.Rest_sinCheckBox.ValueChangedFcn = createCallbackFcn(app, @Rest_sinCheckBoxValueChanged, true);
            app.Rest_sinCheckBox.Text = 'Rest_sin';
            app.Rest_sinCheckBox.Position = [22 97 69 22];

            % Create Rest_2CheckBox
            app.Rest_2CheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.Rest_2CheckBox.ValueChangedFcn = createCallbackFcn(app, @Rest_2CheckBoxValueChanged, true);
            app.Rest_2CheckBox.Text = 'Rest_2';
            app.Rest_2CheckBox.Position = [22 77 60 22];

            % Create ARCheckBox
            app.ARCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.ARCheckBox.ValueChangedFcn = createCallbackFcn(app, @ARCheckBoxValueChanged, true);
            app.ARCheckBox.Text = 'AR';
            app.ARCheckBox.Position = [22 57 39 22];

            % Create Pink_noiseCheckBox
            app.Pink_noiseCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.Pink_noiseCheckBox.ValueChangedFcn = createCallbackFcn(app, @Pink_noiseCheckBoxValueChanged, true);
            app.Pink_noiseCheckBox.Text = 'Pink_noise';
            app.Pink_noiseCheckBox.Position = [96 57 81 22];

            % Create AddCheckBox
            app.AddCheckBox = uicheckbox(app.SourceTimeSeriesPanel);
            app.AddCheckBox.Text = 'Add';
            app.AddCheckBox.FontColor = [0.851 0.3255 0.098];
            app.AddCheckBox.Position = [68 35 71 22];

            % Create PreviewButton
            app.PreviewButton = uibutton(app.SourceTimeSeriesPanel, 'push');
            app.PreviewButton.ButtonPushedFcn = createCallbackFcn(app, @PreviewButtonPushed, true);
            app.PreviewButton.Position = [48 9 100 23];
            app.PreviewButton.Text = 'Preview';

            % Create SamplingRateEditFieldLabel
            app.SamplingRateEditFieldLabel = uilabel(app.SourceTimeSeriesPanel);
            app.SamplingRateEditFieldLabel.HorizontalAlignment = 'right';
            app.SamplingRateEditFieldLabel.Position = [12 211 84 22];
            app.SamplingRateEditFieldLabel.Text = 'Sampling Rate';

            % Create SamplingRateEditField
            app.SamplingRateEditField = uieditfield(app.SourceTimeSeriesPanel, 'numeric');
            app.SamplingRateEditField.Position = [115 211 76 22];
            app.SamplingRateEditField.Value = 1000;

            % Create SourcePatchGenerationPanel
            app.SourcePatchGenerationPanel = uipanel(app.SourceLevelSettingPanel);
            app.SourcePatchGenerationPanel.TitlePosition = 'centertop';
            app.SourcePatchGenerationPanel.Title = ' Source Patch Generation';
            app.SourcePatchGenerationPanel.FontName = 'Arial';
            app.SourcePatchGenerationPanel.FontWeight = 'bold';
            app.SourcePatchGenerationPanel.Position = [207 57 189 219];

            % Create RandButton
            app.RandButton = uibutton(app.SourcePatchGenerationPanel, 'push');
            app.RandButton.ButtonPushedFcn = createCallbackFcn(app, @RandButtonPushed, true);
            app.RandButton.Position = [13 141 48 22];
            app.RandButton.Text = 'Rand';

            % Create LoadButton
            app.LoadButton = uibutton(app.SourcePatchGenerationPanel, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Position = [13 115 48 22];
            app.LoadButton.Text = 'Load';

            % Create SelectfromAtlasButton
            app.SelectfromAtlasButton = uibutton(app.SourcePatchGenerationPanel, 'push');
            app.SelectfromAtlasButton.ButtonPushedFcn = createCallbackFcn(app, @SelectfromAtlasButtonPushed, true);
            app.SelectfromAtlasButton.Position = [66 141 106 22];
            app.SelectfromAtlasButton.Text = 'Select from Atlas';

            % Create SelectfromcortexButton
            app.SelectfromcortexButton = uibutton(app.SourcePatchGenerationPanel, 'push');
            app.SelectfromcortexButton.ButtonPushedFcn = createCallbackFcn(app, @SelectfromcortexButtonPushed, true);
            app.SelectfromcortexButton.Position = [66 115 106 22];
            app.SelectfromcortexButton.Text = 'Select from cortex';

            % Create PatchAreamm2EditFieldLabel
            app.PatchAreamm2EditFieldLabel = uilabel(app.SourcePatchGenerationPanel);
            app.PatchAreamm2EditFieldLabel.Position = [12 83 94 22];
            app.PatchAreamm2EditFieldLabel.Text = 'Patch Area/mm2';

            % Create PatchAreamm2EditField
            app.PatchAreamm2EditField = uieditfield(app.SourcePatchGenerationPanel, 'text');
            app.PatchAreamm2EditField.HorizontalAlignment = 'right';
            app.PatchAreamm2EditField.Position = [111 83 59 22];
            app.PatchAreamm2EditField.Value = '600';

            % Create DecayModelDropDownLabel
            app.DecayModelDropDownLabel = uilabel(app.SourcePatchGenerationPanel);
            app.DecayModelDropDownLabel.HorizontalAlignment = 'right';
            app.DecayModelDropDownLabel.Position = [7 53 75 22];
            app.DecayModelDropDownLabel.Text = 'Decay Model';

            % Create DecayModelDropDown
            app.DecayModelDropDown = uidropdown(app.SourcePatchGenerationPanel);
            app.DecayModelDropDown.Items = {'Gaussian', 'Linear', 'Flat', 'Exponential'};
            app.DecayModelDropDown.ValueChangedFcn = createCallbackFcn(app, @DecayModelDropDownValueChanged, true);
            app.DecayModelDropDown.Position = [93 53 77 22];
            app.DecayModelDropDown.Value = 'Gaussian';

            % Create PowerDecayEditFieldLabel
            app.PowerDecayEditFieldLabel = uilabel(app.SourcePatchGenerationPanel);
            app.PowerDecayEditFieldLabel.Position = [12 22 90 22];
            app.PowerDecayEditFieldLabel.Text = 'Power Decay/%';

            % Create PowerDecayEditField
            app.PowerDecayEditField = uieditfield(app.SourcePatchGenerationPanel, 'text');
            app.PowerDecayEditField.HorizontalAlignment = 'right';
            app.PowerDecayEditField.Position = [111 22 59 22];
            app.PowerDecayEditField.Value = '50';

            % Create SeedVoxelsEditField_2Label
            app.SeedVoxelsEditField_2Label = uilabel(app.SourcePatchGenerationPanel);
            app.SeedVoxelsEditField_2Label.HorizontalAlignment = 'right';
            app.SeedVoxelsEditField_2Label.Position = [8 172 72 22];
            app.SeedVoxelsEditField_2Label.Text = 'Seed Voxels';

            % Create SeedVoxelsEditField
            app.SeedVoxelsEditField = uieditfield(app.SourcePatchGenerationPanel, 'numeric');
            app.SeedVoxelsEditField.Position = [111 172 59 22];
            app.SeedVoxelsEditField.Value = 1;

            % Create AddSourceButton
            app.AddSourceButton = uibutton(app.SourceLevelSettingPanel, 'push');
            app.AddSourceButton.ButtonPushedFcn = createCallbackFcn(app, @AddSourceButtonPushed, true);
            app.AddSourceButton.FontWeight = 'bold';
            app.AddSourceButton.FontColor = [0.149 0.149 0.149];
            app.AddSourceButton.Position = [223 13 157 28];
            app.AddSourceButton.Text = 'Add Source';

            % Create SensorLevelSettingPanel
            app.SensorLevelSettingPanel = uipanel(app.UIFigure);
            app.SensorLevelSettingPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.SensorLevelSettingPanel.TitlePosition = 'centertop';
            app.SensorLevelSettingPanel.Title = 'Sensor Level Setting';
            app.SensorLevelSettingPanel.FontName = 'Arial';
            app.SensorLevelSettingPanel.FontWeight = 'bold';
            app.SensorLevelSettingPanel.Position = [1 136 397 60];

            % Create SimulateNoiseButton
            app.SimulateNoiseButton = uibutton(app.SensorLevelSettingPanel, 'push');
            app.SimulateNoiseButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateNoiseButtonPushed, true);
            app.SimulateNoiseButton.Position = [234 10 96 22];
            app.SimulateNoiseButton.Text = 'Simulate Noise';

            % Create SetBackgroundActivityButton
            app.SetBackgroundActivityButton = uibutton(app.SensorLevelSettingPanel, 'push');
            app.SetBackgroundActivityButton.ButtonPushedFcn = createCallbackFcn(app, @SetBackgroundActivityButtonPushed, true);
            app.SetBackgroundActivityButton.Position = [55 10 137 22];
            app.SetBackgroundActivityButton.Text = 'Set Background Activity';

            % Create SourceListPanel
            app.SourceListPanel = uipanel(app.UIFigure);
            app.SourceListPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.SourceListPanel.TitlePosition = 'centertop';
            app.SourceListPanel.Title = 'Source List';
            app.SourceListPanel.FontName = 'Arial';
            app.SourceListPanel.FontWeight = 'bold';
            app.SourceListPanel.Position = [1 28 397 109];

            % Create SourceTable
            app.SourceTable = uitable(app.SourceListPanel);
            app.SourceTable.ColumnName = {'No.'; 'Vertices'; 'Stimulus'};
            app.SourceTable.RowName = {};
            app.SourceTable.SelectionType = 'row';
            app.SourceTable.Position = [0 0 397 89];

            % Create OKButton
            app.OKButton = uibutton(app.UIFigure, 'push');
            app.OKButton.ButtonPushedFcn = createCallbackFcn(app, @OKButtonPushed, true);
            app.OKButton.Position = [279 4 100 23];
            app.OKButton.Text = 'OK';

            % Create StandardAnswerButton
            app.StandardAnswerButton = uibutton(app.UIFigure, 'push');
            app.StandardAnswerButton.ButtonPushedFcn = createCallbackFcn(app, @StandardAnswerButtonPushed, true);
            app.StandardAnswerButton.Position = [44 4 106 23];
            app.StandardAnswerButton.Text = 'Standard Answer';

            % Create LoadRealisticNoiseButton
            app.LoadRealisticNoiseButton = uibutton(app.UIFigure, 'push');
            app.LoadRealisticNoiseButton.Position = [499 136 125 22];
            app.LoadRealisticNoiseButton.Text = 'Load Realistic Noise';

            % Create SourceListMenu
            app.SourceListMenu = uicontextmenu(app.UIFigure);

            % Create deleteMenu
            app.deleteMenu = uimenu(app.SourceListMenu);
            app.deleteMenu.MenuSelectedFcn = createCallbackFcn(app, @deleteMenuSelected, true);
            app.deleteMenu.Text = 'delete';
            
            % Assign app.SourceListMenu
            app.SourceTable.ContextMenu = app.SourceListMenu;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SEAL_Simulation(varargin)

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
