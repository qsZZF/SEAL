<<<<<<< HEAD
<<<<<<< HEAD
function [data_out, fs_out] = applyDownsample_resample(data_in, fs_in, D)
% data_in: ch × time × trials
% D      : 降采样因子（>1，整数最佳）
% 输出   : data_out (ch × time_new × trials), fs_out = fs_in / D

    [nCh, nT, nTr] = size(data_in);

    % 若 D 不是整数，尝试转为有理数 p/q ≈ 1/D，以便 resample(p,q) 使用
    % 但本 GUI 语义是“降采样因子”，推荐整数；这里做个兜底：
    if abs(D - round(D)) < 1e-9
        p = 1; q = round(D);
    else
        % 用 rat 寻找 p/q ≈ 1/D（限制分母大小，避免过高阶滤波）
        [p,q] = rat(1/D, 1e-6);
        if q > 200
            % 分母太大就强制转为最近整数因子
            q = max(2, round(D)); p = 1;
            warning('非整数降采样因子已近似为 D=%d（p=1,q=%d）。', q, q);
        end
    end
    fs_out = fs_in * (p/q);   % 对于整数 D，fs_out = fs_in / D

    % 用第一条试验估算输出长度，以便一次性预分配
    x0 = double(squeeze(data_in(1,:,1)));
    y0 = resample(x0, p, q);          % FIR 抗混叠 + 延迟补偿
    nT_new = numel(y0);

    data_out = zeros(nCh, nT_new, nTr, class(y0));  % 预分配

    % 主循环
    for ch = 1:nCh
        for tr = 1:nTr
            x = double(squeeze(data_in(ch,:,tr)));  % 1 × nT
            y = resample(x, p, q);
            % 尺寸保护（极端参数下可能 ±1 点差异）
            nt = min(nT_new, numel(y));
            data_out(ch,1:nt,tr) = y(1:nt);
        end
        if mod(ch, max(1, round(nCh/10)))==0
            fprintf('  -> 已完成通道 %d/%d\n', ch, nCh);
        end
    end
end
=======
function [data_out, fs_out] = applyDownsample_resample(data_in, fs_in, D)
% data_in: ch × time × trials
% D      : 降采样因子（>1，整数最佳）
% 输出   : data_out (ch × time_new × trials), fs_out = fs_in / D

    [nCh, nT, nTr] = size(data_in);

    % 若 D 不是整数，尝试转为有理数 p/q ≈ 1/D，以便 resample(p,q) 使用
    % 但本 GUI 语义是“降采样因子”，推荐整数；这里做个兜底：
    if abs(D - round(D)) < 1e-9
        p = 1; q = round(D);
    else
        % 用 rat 寻找 p/q ≈ 1/D（限制分母大小，避免过高阶滤波）
        [p,q] = rat(1/D, 1e-6);
        if q > 200
            % 分母太大就强制转为最近整数因子
            q = max(2, round(D)); p = 1;
            warning('非整数降采样因子已近似为 D=%d（p=1,q=%d）。', q, q);
        end
    end
    fs_out = fs_in * (p/q);   % 对于整数 D，fs_out = fs_in / D

    % 用第一条试验估算输出长度，以便一次性预分配
    x0 = double(squeeze(data_in(1,:,1)));
    y0 = resample(x0, p, q);          % FIR 抗混叠 + 延迟补偿
    nT_new = numel(y0);

    data_out = zeros(nCh, nT_new, nTr, class(y0));  % 预分配

    % 主循环
    for ch = 1:nCh
        for tr = 1:nTr
            x = double(squeeze(data_in(ch,:,tr)));  % 1 × nT
            y = resample(x, p, q);
            % 尺寸保护（极端参数下可能 ±1 点差异）
            nt = min(nT_new, numel(y));
            data_out(ch,1:nt,tr) = y(1:nt);
        end
        if mod(ch, max(1, round(nCh/10)))==0
            fprintf('  -> 已完成通道 %d/%d\n', ch, nCh);
        end
    end
end
>>>>>>> origin
=======
function [data_out, fs_out] = applyDownsample_resample(data_in, fs_in, D)
% data_in: ch × time × trials
% D      : 降采样因子（>1，整数最佳）
% 输出   : data_out (ch × time_new × trials), fs_out = fs_in / D

    [nCh, nT, nTr] = size(data_in);

    % 若 D 不是整数，尝试转为有理数 p/q ≈ 1/D，以便 resample(p,q) 使用
    % 但本 GUI 语义是“降采样因子”，推荐整数；这里做个兜底：
    if abs(D - round(D)) < 1e-9
        p = 1; q = round(D);
    else
        % 用 rat 寻找 p/q ≈ 1/D（限制分母大小，避免过高阶滤波）
        [p,q] = rat(1/D, 1e-6);
        if q > 200
            % 分母太大就强制转为最近整数因子
            q = max(2, round(D)); p = 1;
            warning('非整数降采样因子已近似为 D=%d（p=1,q=%d）。', q, q);
        end
    end
    fs_out = fs_in * (p/q);   % 对于整数 D，fs_out = fs_in / D

    % 用第一条试验估算输出长度，以便一次性预分配
    x0 = double(squeeze(data_in(1,:,1)));
    y0 = resample(x0, p, q);          % FIR 抗混叠 + 延迟补偿
    nT_new = numel(y0);

    data_out = zeros(nCh, nT_new, nTr, class(y0));  % 预分配

    % 主循环
    for ch = 1:nCh
        for tr = 1:nTr
            x = double(squeeze(data_in(ch,:,tr)));  % 1 × nT
            y = resample(x, p, q);
            % 尺寸保护（极端参数下可能 ±1 点差异）
            nt = min(nT_new, numel(y));
            data_out(ch,1:nt,tr) = y(1:nt);
        end
        if mod(ch, max(1, round(nCh/10)))==0
            fprintf('  -> 已完成通道 %d/%d\n', ch, nCh);
        end
    end
end
>>>>>>> 3afa99beaec32453db712d4621fd05217f53fcfd
