function [dataOut, info] = seal_runEpoching(dataIn, srate, eventSpec, epochWin, opts)
%SEAL_RUNEPOCHING 将连续 EEG 数据切分为 trials
%
%   [dataOut, info] = seal_runEpoching(dataIn, srate, eventSpec, epochWin, opts)
%
%   Inputs:
%     dataIn     - ch×T 连续数据 (必须是 2D)
%     srate      - 采样率 (Hz)
%     eventSpec  - 事件规范结构体,含字段:
%                  .source   - 'event_channel' | 'events_struct' | 'fixed_interval'
%                  .data     - 根据 source 不同:
%                              * 'event_channel':  通道索引(标量,事件通道在 dataIn 中的行号)
%                              * 'events_struct':  结构体数组,含 .latency(样本索引)和 .type
%                              * 'fixed_interval': 标量,间隔秒数
%                  .typesWanted - (可选)要提取的事件类型,空表示全部
%                  .threshold   - (仅 event_channel 用)脉冲检测阈值
%     epochWin   - [t_pre, t_post] 单位秒,例如 [-0.2, 0.8]
%     opts       - 可选结构体:
%                  .boundaryMode - 'reject'(默认)| 'zeropad' | 'mirror'
%                  .rejectAmp    - 振幅阈值(μV),超过则剔除该 trial,默认 Inf(不剔除)
%                  .rejectGrad   - 梯度阈值(μV/样本),默认 Inf
%                  .rejectPtP    - 峰峰值阈值(μV),默认 Inf
%                  .ampUnit      - 数据振幅单位,'uV'(默认)| 'V' | 'mV' | 'nV'
%                  .verbose      - 是否打印详细信息,默认 true
%
%   Outputs:
%     dataOut    - ch × T_epoch × nTrials 分段后数据(3D)
%     info       - 元信息结构体:
%                  .nTrialsFound     总共找到的事件数
%                  .nTrialsExtracted 成功提取的 trials 数
%                  .nRejectedBound   因越界被剔除
%                  .nRejectedAmp     因振幅超限被剔除
%                  .trialLatencies   每个保留 trial 的 latency (样本索引)
%                  .trialTypes       每个保留 trial 的 type
%                  .epochTimes       epoch 内的时间向量 (秒)
%                  .rejectedInfo     被剔除 trial 的详情

    % ===== 默认 opts =====
    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'boundaryMode'), opts.boundaryMode = 'reject'; end
    if ~isfield(opts, 'rejectAmp'),    opts.rejectAmp    = Inf; end
    if ~isfield(opts, 'rejectGrad'),   opts.rejectGrad   = Inf; end
    if ~isfield(opts, 'rejectPtP'),    opts.rejectPtP    = Inf; end
    if ~isfield(opts, 'ampUnit'),      opts.ampUnit      = 'uV'; end
    if ~isfield(opts, 'verbose'),      opts.verbose      = true; end

    % ===== 输入校验 =====
    if isempty(dataIn) || ~isnumeric(dataIn)
        error('SEAL:Epoch:InvalidInput', '数据为空或类型不合法');
    end
    if ~ismatrix(dataIn)
        error('SEAL:Epoch:Need2D', ...
            '分段功能需要 2D 连续数据(ch×T),当前维度 %d,可能已分段', ndims(dataIn));
    end
    if ~isfinite(srate) || srate <= 0
        error('SEAL:Epoch:InvalidSrate', '采样率非法: %g', srate);
    end
    if numel(epochWin) ~= 2 || epochWin(1) >= epochWin(2)
        error('SEAL:Epoch:InvalidWindow', ...
            'epoch 窗口必须是 [t_pre, t_post] 且 t_pre < t_post,当前 %s', ...
            mat2str(epochWin));
    end

    [nCh, nT] = size(dataIn);

    % ===== 单位缩放(检测时统一到 μV) =====
    unitScale = getUnitScaleToMicroVolt(opts.ampUnit);

    % ===== Step 1:根据 eventSpec 得到事件 latency 和 type 列表 =====
    [evLatencies, evTypes] = resolveEvents(dataIn, srate, eventSpec, nT, opts.verbose);

    if isempty(evLatencies)
        error('SEAL:Epoch:NoEvents', '未检测到任何事件,请检查参数');
    end

    % ===== Step 2:按事件类型筛选 =====
    if isfield(eventSpec, 'typesWanted') && ~isempty(eventSpec.typesWanted)
        want = eventSpec.typesWanted;

        % 统一 evTypes 为字符串 cell
        evTypesStr = cell(size(evTypes));
        for i = 1:numel(evTypes)
            t = evTypes{i};
            if isnumeric(t)
                evTypesStr{i} = num2str(t);
            else
                evTypesStr{i} = char(t);
            end
        end

        % 统一 want 为字符串 cell
        if isnumeric(want)
            wantStr = arrayfun(@num2str, want, 'UniformOutput', false);
        elseif iscell(want)
            wantStr = cell(size(want));
            for i = 1:numel(want)
                if isnumeric(want{i})
                    wantStr{i} = num2str(want{i});
                else
                    wantStr{i} = char(want{i});
                end
            end
        else
            wantStr = {char(want)};
        end

        mask = ismember(evTypesStr, wantStr);
        evLatencies = evLatencies(mask);
        evTypes     = evTypes(mask);
        if opts.verbose
            fprintf('[Epoch] 按事件类型筛选后剩余 %d 个事件\n', numel(evLatencies));
        end
    end

    if isempty(evLatencies)
        error('SEAL:Epoch:NoMatchingEvents', '筛选后无匹配事件');
    end

    % ===== Step 3:计算 epoch 时间-样本映射 =====
    preSamples  = round(epochWin(1) * srate);   % 负数
    postSamples = round(epochWin(2) * srate) - 1;
    epochLen    = postSamples - preSamples + 1;

    if epochLen < 2
        error('SEAL:Epoch:WindowTooShort', ...
            'Epoch 窗口过短,只含 %d 个样本', epochLen);
    end

    epochTimes = (preSamples:postSamples) / srate;   % 秒

    if opts.verbose
        fprintf('[Epoch] epoch 窗口 = [%.3f, %.3f] s → 每段 %d 样本\n', ...
            epochWin(1), epochWin(2), epochLen);
    end

    % ===== Step 4:遍历事件,提取每个 trial =====
    nEvents = numel(evLatencies);
    validMask    = false(nEvents, 1);
    rejectReason = repmat("", nEvents, 1);

    % 预分配(按最大可能数)
    dataTmp = zeros(nCh, epochLen, nEvents, 'like', dataIn);

    for i = 1:nEvents
        eventSample = evLatencies(i);
        idxStart    = eventSample + preSamples;
        idxEnd      = eventSample + postSamples;

        % --- 边界处理 ---
        if idxStart < 1 || idxEnd > nT
            switch opts.boundaryMode
                case 'reject'
                    rejectReason(i) = "boundary";
                    continue;
                case 'zeropad'
                    trial = zeros(nCh, epochLen, 'like', dataIn);
                    src_start = max(1, idxStart);
                    src_end   = min(nT, idxEnd);
                    dst_start = src_start - idxStart + 1;
                    dst_end   = src_end   - idxStart + 1;
                    trial(:, dst_start:dst_end) = dataIn(:, src_start:src_end);
                case 'mirror'
                    trial = mirrorPadExtract(dataIn, idxStart, idxEnd);
                otherwise
                    error('未知 boundaryMode: %s', opts.boundaryMode);
            end
        else
            trial = dataIn(:, idxStart:idxEnd);
        end

        % --- 振幅/梯度/峰峰值筛查(仅当阈值有限时) ---
        if isfinite(opts.rejectAmp) || isfinite(opts.rejectGrad) || isfinite(opts.rejectPtP)
            trial_uV = double(trial) * unitScale;

            if isfinite(opts.rejectAmp) && max(abs(trial_uV), [], 'all') > opts.rejectAmp
                rejectReason(i) = "amplitude";
                continue;
            end
            if isfinite(opts.rejectGrad) && max(abs(diff(trial_uV, 1, 2)), [], 'all') > opts.rejectGrad
                rejectReason(i) = "gradient";
                continue;
            end
            if isfinite(opts.rejectPtP)
                ptp = max(trial_uV, [], 2) - min(trial_uV, [], 2);
                if max(ptp) > opts.rejectPtP
                    rejectReason(i) = "peak-to-peak";
                    continue;
                end
            end
        end

        dataTmp(:, :, i) = trial;
        validMask(i) = true;
    end

    % ===== Step 5:收集有效 trials =====
    dataOut = dataTmp(:, :, validMask);

    if isempty(dataOut)
        error('SEAL:Epoch:NoValidTrials', ...
            '所有 %d 个事件都被剔除,请放宽阈值或检查边界', nEvents);
    end

    % ===== 统计信息 =====
    info.nTrialsFound     = nEvents;
    info.nTrialsExtracted = sum(validMask);
    info.nRejectedBound   = sum(rejectReason == "boundary");
    info.nRejectedAmp     = sum(rejectReason == "amplitude");
    info.nRejectedGrad    = sum(rejectReason == "gradient");
    info.nRejectedPtP     = sum(rejectReason == "peak-to-peak");
    info.trialLatencies   = evLatencies(validMask);
    info.trialTypes       = evTypes(validMask);
    info.epochTimes       = epochTimes;
    info.epochWindow      = epochWin;
    info.srate            = srate;
    info.rejectedInfo     = struct( ...
        'indices', find(~validMask), ...
        'reasons', rejectReason(~validMask));

    if opts.verbose
        fprintf('\n===== [Epoch] 完成 =====\n');
        fprintf('  总事件数:   %d\n', info.nTrialsFound);
        fprintf('  成功提取:   %d\n', info.nTrialsExtracted);
        fprintf('  边界剔除:   %d\n', info.nRejectedBound);
        fprintf('  振幅剔除:   %d\n', info.nRejectedAmp);
        fprintf('  梯度剔除:   %d\n', info.nRejectedGrad);
        fprintf('  峰峰剔除:   %d\n', info.nRejectedPtP);
        fprintf('  输出维度:   %d × %d × %d (ch × T × trials)\n', ...
            size(dataOut, 1), size(dataOut, 2), size(dataOut, 3));
    end
end


% =========================================================================
% 辅助函数
% =========================================================================

function [latencies, types] = resolveEvents(dataIn, srate, eventSpec, nT, verbose)
%RESOLVEEVENTS 根据 eventSpec 提取事件 latency 和 type 列表

    switch lower(eventSpec.source)
        case 'event_channel'
            % 从专用事件通道中检测脉冲(触发)
            chIdx = eventSpec.data;
            if chIdx < 1 || chIdx > size(dataIn, 1)
                error('SEAL:Epoch:BadEventChannel', ...
                    '事件通道索引 %d 越界(有效 1..%d)', chIdx, size(dataIn,1));
            end

            thresh = eventSpec.threshold;
            if ~isfinite(thresh)
                % 自动阈值:信号最大值的 50%
                thresh = 0.5 * max(dataIn(chIdx, :));
            end

            % 检测上升沿(从低于阈值 → 高于阈值)
            sig = double(dataIn(chIdx, :));
            above = sig > thresh;
            latencies = find(diff([0, above]) == 1);

            % 事件通道的 type 通常就用阈值处的值(可以做多电平解码,这里简化)
            types = ones(1, numel(latencies));

            if verbose
                fprintf('[Epoch] 从通道 %d 检测到 %d 个事件(阈值 %.3g)\n', ...
                    chIdx, numel(latencies), thresh);
            end

        case 'events_struct'
            ev = eventSpec.data;
            if isempty(ev)
                latencies = [];
                types = [];
                return;
            end
            latencies = round([ev.latency]);

            % type 可能是数字、字符串,统一为 cell
            if isfield(ev, 'type')
                types = {ev.type};
                % 如果都是数字,转成数组更方便
                if all(cellfun(@(x) isnumeric(x) && isscalar(x), types))
                    types = cell2mat(types);
                end
            else
                types = ones(1, numel(latencies));
            end

            if verbose
                fprintf('[Epoch] 从 events 结构体读取 %d 个事件\n', numel(latencies));
            end

        case 'fixed_interval'
            intervalSec = eventSpec.data;
            if ~isfinite(intervalSec) || intervalSec <= 0
                error('SEAL:Epoch:BadInterval', '间隔秒数必须为正: %g', intervalSec);
            end
            intervalSamp = round(intervalSec * srate);
            % 从第一个间隔点开始放事件(留够 pre-epoch 空间)
            latencies = intervalSamp : intervalSamp : nT;
            types = ones(1, numel(latencies));

            if verbose
                fprintf('[Epoch] 固定间隔 %.3fs 生成 %d 个事件\n', ...
                    intervalSec, numel(latencies));
            end

        otherwise
            error('SEAL:Epoch:UnknownSource', ...
                '未知 eventSpec.source: %s', eventSpec.source);
    end

    latencies = latencies(:).';
    if iscell(types), types = types(:).'; else, types = types(:).'; end
end


function trial = mirrorPadExtract(dataIn, idxStart, idxEnd)
%MIRRORPADEXTRACT 越界部分用镜像延拓
    [nCh, nT] = size(dataIn);
    epochLen = idxEnd - idxStart + 1;
    trial = zeros(nCh, epochLen, 'like', dataIn);
    for k = 1:epochLen
        srcIdx = idxStart + k - 1;
        if srcIdx < 1
            srcIdx = 2 - srcIdx;          % 1 → 1, 0 → 2, -1 → 3 ...
        elseif srcIdx > nT
            srcIdx = 2 * nT - srcIdx;     % nT+1 → nT-1, nT+2 → nT-2 ...
        end
        srcIdx = max(1, min(nT, srcIdx)); % 最终夹紧(极端情况)
        trial(:, k) = dataIn(:, srcIdx);
    end
end


function scale = getUnitScaleToMicroVolt(unit)
%GETUNITSCALETOMICROVOLT 把数据振幅换算到 μV
    switch lower(string(unit))
        case "v",        scale = 1e6;
        case "mv",       scale = 1e3;
        case {"uv", "μv"}, scale = 1;
        case "nv",       scale = 1e-3;
        otherwise
            warning('SEAL:Epoch:UnknownUnit', '未知单位 %s,按 μV 处理', unit);
            scale = 1;
    end
end