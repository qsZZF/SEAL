function filteredData = seal_filter_data(data, srate, filterParams)
% SEAL_FILTER_DATA 频域滤波核心算法（无 GUI 依赖）
% 支持高通/低通(FIR)和陷波(IIR)，逐试次滤波避免边界伪迹
%
% 输入：
%   data         - [nCh × nT] 连续数据 或 [nCh × nT × nTr] 分段数据
%   srate        - 采样率 (Hz)
%   filterParams - 结构体，可含字段 .highpass .lowpass .notch
%
% 输出：
%   filteredData - 与输入同维度的滤波后数据

    % ========== 1. 构建滤波器 ==========
    filters = {};

    % --- 高通 (FIR) ---
    if isfield(filterParams, 'highpass') && filterParams.highpass > 0
        b_hp = designFIRHighpass(srate, filterParams.highpass, 40, 100);
        filters{end+1} = struct('b', b_hp, 'a', 1);
    end

    % --- 低通 (FIR) ---
    if isfield(filterParams, 'lowpass') && filterParams.lowpass > 0
        b_lp = designFIRLowpass(srate, filterParams.lowpass, 40, 10);
        filters{end+1} = struct('b', b_lp, 'a', 1);
    end

    % --- 陷波 (IIR Butterworth) ---
    if isfield(filterParams, 'notch') && filterParams.notch > 0
        width = 2;
        Nyq   = srate / 2;
        Flow  = max(filterParams.notch - width/2, 0.1);
        Fhigh = min(filterParams.notch + width/2, Nyq * 0.99);
        Wn    = [Flow, Fhigh] / Nyq;
        [b_notch, a_notch] = butter(2, Wn, 'stop');
        filters{end+1} = struct('b', b_notch, 'a', a_notch);
    end

    if isempty(filters)
        error('SEAL:NoFilter', '未指定任何有效的滤波器参数。');
    end

    % ========== 2. 计算最小数据长度要求 ==========
    minRequired = 0;
    for j = 1:length(filters)
        minRequired = max(minRequired, 3 * max(length(filters{j}.b), length(filters{j}.a)));
    end

    % ========== 3. 判断数据维度并逐试次滤波 ==========
    nCh = size(data, 1);
    nT  = size(data, 2);
    nTr = size(data, 3);   % 2D 时自动返回 1

    if nT <= minRequired
        error('SEAL:DataTooShort', ...
            '每个试次的时间点数(%d)太短，滤波要求至少 %d 个点。\n建议对未分段的连续数据滤波，或降低滤波器阶数。', ...
            nT, minRequired + 1);
    end

    filteredData = zeros(size(data), 'like', data);

    for tr = 1:nTr
        % 取出当前试次：[nCh × nT]
        trialData = data(:, :, tr);

        % 转置为 [nT × nCh]，filtfilt 沿第一维（时间）滤波
        trialData = trialData.';

        for j = 1:length(filters)
            trialData = filtfilt(filters{j}.b, filters{j}.a, trialData);
        end

        % 转置回 [nCh × nT] 写入结果
        filteredData(:, :, tr) = trialData.';
    end

    % 如果输入是 2D，输出也保持 2D
    if ismatrix(data)
        filteredData = squeeze(filteredData);
    end
end