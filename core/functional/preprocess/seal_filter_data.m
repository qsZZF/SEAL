function filteredData = seal_filter_data(data3d, srate, filterParams)
    % SEAL_FILTER_DATA 频域滤波核心算法 (无 GUI 依赖)
    % 支持高低通(FIR)和陷波(IIR)，自带向量化加速与数据长度保护

    % 1. 维度获取与重塑 (向量化核心)
    [nCh, nT, nTr] = size(data3d);
    data_perm = permute(data3d, [2, 1, 3]); 
    data2D = reshape(data_perm, nT, nCh * nTr); % 变为 Time x (Ch*Trials)
    
    % 2. 统一管理滤波器系数：将 b 和 a 存入结构体
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
    
    % --- 陷波 (改为 IIR 巴特沃斯：科学、极速、不卡长度) ---
    if isfield(filterParams, 'notch') && filterParams.notch > 0
        width = 2; % 陷波宽度 2Hz 非常足够
        Flow = max(filterParams.notch - width/2, 0.1);
        Fhigh = min(filterParams.notch + width/2, srate/2 * 0.99);
        
        % 计算归一化频率 (除以奈奎斯特频率)
        Nyq = srate / 2;
        Wn = [Flow, Fhigh] / Nyq;
        
        % MATLAB内置：设计二阶 IIR 带阻滤波器
        [b_notch, a_notch] = butter(2, Wn, 'stop');
        
        % IIR 滤波器有 b 也有 a
        filters{end+1} = struct('b', b_notch, 'a', a_notch);
    end
    
    if isempty(filters)
        error('SEAL:NoFilter', '未指定任何有效的滤波器参数。');
    end

    % 3. 极速向量化滤波执行
    for j = 1:length(filters)
        b = filters{j}.b;
        a = filters{j}.a;
        
        % ⚠️ 绝对防爆盾：检查数据长度是否满足 filtfilt 最小要求
        min_required_pts = 3 * max(length(b), length(a));
        if nT <= min_required_pts
            error('SEAL:DataTooShort', '数据点数(%d)太短！\n滤波要求长度大于 %d 点。\n请对未切段的连续数据滤波，或降低 FIR 阶数。', nT, min_required_pts);
        end
        
        % 执行零相位滤波 (完美兼容 FIR 的 a=1 和 IIR 的实际 a 值)
        data2D = filtfilt(b, a, data2D);
    end
    
    % 4. 还原三维矩阵
    data_perm_back = reshape(data2D, nT, nCh, nTr);
    filteredData = permute(data_perm_back, [2, 1, 3]);
end