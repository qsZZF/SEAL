function C_reg = seal_runNoiseCovariance(dataIn, baseStart, baseEnd, lambda)
    % dataIn: ch x time x trials
    % baseStart, baseEnd: 基线期在时间维度上的起始和结束索引（点数）
    % lambda: 正则化系数 (0 到 1 之间)
    
    [nCh, nTime, nTrials] = size(dataIn);
    
    % 检查基线范围的合法性
    if ~isfinite(baseStart) || ~isfinite(baseEnd) || baseStart < 1 || baseEnd > nTime || baseStart >= baseEnd
        error('基线点范围设置不合法。当前数据总点数为 %d。', nTime);
    end
    
    % 截取基线期数据
    baselineData = dataIn(:, baseStart:baseEnd, :);
    nBasePoints = size(baselineData, 2);
    
    % 初始化协方差矩阵
    C_empirical = zeros(nCh, nCh);
    
    % 逐 trial 计算并累加协方差
    for tr = 1:nTrials
        % 取出当前 trial 的基线数据: ch x time
        tmp = double(baselineData(:, :, tr));
        
        % 去除每个通道基线期的直流偏移（Mean centering）
        tmp = tmp - mean(tmp, 2);
        
        % 计算协方差并累加 (除以 n-1 是无偏估计)
        C_empirical = C_empirical + (tmp * tmp') / (nBasePoints - 1);
    end
    
    % 取所有 trials 的平均
    C_empirical = C_empirical / nTrials;
    
    % ==========================================
    % 矩阵正则化 (Tikhonov / Shrinkage Regularization)
    % 解决插值或平均参考导致的矩阵不满秩问题
    % ==========================================
    if lambda > 0
        % 计算矩阵的迹的平均值（对角线元素的平均）
        meanTrace = trace(C_empirical) / nCh;
        % 将经验协方差矩阵向单位矩阵（按比例缩放后）收缩
        C_reg = (1 - lambda) * C_empirical + lambda * meanTrace * eye(nCh);
    else
        C_reg = C_empirical;
    end
end