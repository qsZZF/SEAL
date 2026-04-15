function dataOut = seal_runBaselineCorrection(dataIn, baseStart, baseEnd)
    % dataIn: ch x time x trials 的三维或二维脑电数据
    % baseStart: 基线起始点索引
    % baseEnd: 基线结束点索引
    
    [nCh, nTime, nTrials] = size(dataIn);
    
    % 严格校验基线区间合法性
    if baseStart < 1 || baseEnd > nTime || baseStart > baseEnd
        error('基线区间设置错误！当前数据时间点范围为 1 到 %d。', nTime);
    end
    
    % 初始化输出矩阵
    dataOut = zeros(size(dataIn));
    
    % 逐 trial 进行基线校准
    for tr = 1:nTrials
        % 取出当前 trial 数据 (ch x time)
        trialData = double(dataIn(:, :, tr));
        
        % 计算每个通道在基线期的平均值 (ch x 1)
        baseMean = mean(trialData(:, baseStart:baseEnd), 2);
        
        % 减去基线均值（利用 MATLAB 的隐式扩展，ch x 1 会自动与 ch x time 逐列相减）
        dataOut(:, :, tr) = trialData - baseMean;
    end
end