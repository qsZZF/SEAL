function dataOut = seal_runBaselineCorrection(dataIn, baseStart, baseEnd)
% dataIn: ch x time x trials 的三维或二维脑电数据
% baseStart: 基线起始点索引
% baseEnd: 基线结束点索引
if ismatrix(dataIn) || size(dataIn, 3) == 1
    warning('SEAL:Baseline:ContinuousData', ...
        ['基线校准应用于连续数据(单 trial),将整段数据减去 [%d, %d] 区间的均值。\n' ...
        '若需要 ERP 风格的逐 trial 基线校准,请先执行 Extract Trials。'], ...
        baseStart, baseEnd);
end


[nCh, nTime, nTrials] = size(dataIn);

% 严格校验基线区间合法性
if baseStart < 1 || baseEnd > nTime || baseStart > baseEnd
    error('基线区间设置错误！当前数据时间点范围为 1 到 %d。', nTime);
end

nBasePoints = baseEnd - baseStart + 1;
if nBasePoints < 10
    warning('SEAL:Baseline:TooShort', ...
        '基线区间过短(%d 个样本),均值估计可能不稳定。建议至少 50ms 或 50 个样本。', ...
        nBasePoints);
end

% 初始化输出矩阵
dataOut = zeros(size(dataIn));

% 一次性提取所有 trial 的基线段:ch × baseLen × Tr
baseSegment = dataIn(:, baseStart:baseEnd, :);

% 沿时间维求均值:ch × 1 × Tr
baseMean = mean(baseSegment, 2);

% 隐式扩展减法:ch × T × Tr 减去 ch × 1 × Tr
dataOut = dataIn - baseMean;
end