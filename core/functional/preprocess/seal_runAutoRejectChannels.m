function [dataClean, locsClean, labelsClean, badIdx, badNames] = runAutoRejectChannels(dataIn, chanlocs, labels, zThresh, corrThresh)
    % 输入：
    % dataIn: ch x time x trials 的三维脑电数据
    % chanlocs: 结构体数组，包含通道坐标
    % labels: 字符串数组或元胞数组，通道名称
    % zThresh: 方差的稳健 Z-score 阈值（推荐 3 到 5）
    % corrThresh: 最小相关性阈值（推荐 0.4）
    
    [nCh, nTime, nTrials] = size(dataIn);
    
    % 1. 将数据展平为 2D (ch x total_time) 以便计算整体统计量
    data2D = double(reshape(dataIn, nCh, nTime * nTrials));
    
    % ==========================================
    % 指标 A：基于振幅/方差的异常检测 (Robust Z-score)
    % ==========================================
    % 计算每个通道的标准差
    chStd = std(data2D, 0, 2); 
    
    % 使用绝对中位差 (MAD) 计算稳健的 Z-score (防止极端坏道拉高整体均值)
    medStd = median(chStd);
    madStd = median(abs(chStd - medStd));
    
    % 稳健 Z-score 公式 (0.6745 是正态分布的缩放因子)
    zScore = 0.6745 * (chStd - medStd) ./ (madStd + eps); 
    
    % 找出极度吵闹 (Z > zThresh) 或 彻底死机 (Z < -zThresh) 的通道
    badVarIdx = find(abs(zScore) > zThresh);
    
    % ==========================================
    % 指标 B：基于跨通道相关性的异常检测
    % ==========================================
    % 为了计算速度，如果数据太长，可以随机抽取几万个点计算相关性
    maxPts = min(size(data2D, 2), 20000); 
    randIdx = randperm(size(data2D, 2), maxPts);
    
    % 计算通道间的相关系数矩阵
    corrMat = corrcoef(data2D(:, randIdx)'); 
    corrMat = abs(corrMat - eye(nCh)); % 取绝对值并去掉对角线的 1
    
    % 找出每个通道与其他通道的“最大相关性” (即它最好朋友的亲密程度)
    maxCorr = max(corrMat, [], 2);
    
    % 找出与所有通道都毫无默契的“孤岛”通道
    badCorrIdx = find(maxCorr < corrThresh);
    
    % ==========================================
    % 汇总坏道并执行剔除
    % ==========================================
    % 合并两类坏道并去重
    badIdx = unique([badVarIdx; badCorrIdx]);
    
    % 获取好道索引
    goodIdx = setdiff(1:nCh, badIdx);
    
    % 提取好道数据、坐标和标签
    dataClean   = dataIn(goodIdx, :, :);
    
    if ~isempty(chanlocs)
        locsClean = chanlocs(goodIdx);
    else
        locsClean = [];
    end
    
    if ~isempty(labels)
        labelsClean = labels(goodIdx);
        badNames    = labels(badIdx);
    else
        labelsClean = [];
        badNames    = arrayfun(@(x) sprintf('Ch%d', x), badIdx, 'UniformOutput', false);
    end
end