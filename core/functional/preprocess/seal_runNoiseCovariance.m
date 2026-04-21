function [C_reg, info] = seal_runNoiseCovariance(dataIn, baseStart, baseEnd, lambda)
%SEAL_RUNNOISECOVARIANCE 计算 EEG 噪声协方差矩阵
%   Inputs:
%     dataIn    - ch × time × trials (3D 必需)
%     baseStart - 基线起点(样本索引)
%     baseEnd   - 基线终点(样本索引)
%     lambda    - 收缩系数 [0, 1]
%   Outputs:
%     C_reg     - 正则化后协方差矩阵
%     info      - 元信息结构体(秩、样本数等)

    % ===== 输入校验 =====
    if ndims(dataIn) ~= 3 || size(dataIn, 3) < 2
        error('SEAL:NoiseCov:NeedsEpochs', ...
            '需要 3D 分段数据(ch×time×trials),且至少 2 个 trial');
    end
    
    if ~isfinite(lambda) || lambda < 0 || lambda > 1
        error('SEAL:NoiseCov:InvalidLambda', ...
            'lambda 必须在 [0, 1],当前 %g', lambda);
    end
    
    [nCh, nTime, nTrials] = size(dataIn);
    
    if ~isfinite(baseStart) || ~isfinite(baseEnd) || ...
       baseStart < 1 || baseEnd > nTime || baseStart >= baseEnd
        error('SEAL:NoiseCov:InvalidBaseline', ...
            '基线范围非法 [%g, %g],数据长度 %d', baseStart, baseEnd, nTime);
    end
    
    % ===== 提取并去均值(向量化) =====
    X = double(dataIn(:, baseStart:baseEnd, :));  % ch × T × Tr
    nBasePoints = size(X, 2);
    nSamples = nBasePoints * nTrials;
    
    % 逐 trial 去均值(隐式扩展)
    X = X - mean(X, 2);
    
    % 展平为 2D
    X2d = reshape(X, nCh, nSamples);
    
    % ===== Pooled 协方差(一次矩阵乘) =====
    C_empirical = (X2d * X2d') / (nSamples - 1);
    
    % ===== Tikhonov shrinkage =====
    if ~isfinite(lambda) || lambda < 0 || lambda > 1
        error('SEAL:NoiseCov:InvalidLambda', ...
            'lambda 必须在 [0, 1] 范围内,当前 %g', lambda);
    end

    if lambda > 0
        meanTrace = trace(C_empirical) / nCh;
        C_reg = (1 - lambda) * C_empirical + lambda * meanTrace * eye(nCh);
    else
        C_reg = C_empirical;
    end

    % ===== 元信息 =====
    if nargout > 1
        tol = 1e-8 * max(diag(C_empirical));
        info.rank_empirical = rank(C_empirical, tol);
        info.rank_reg = rank(C_reg);
        info.shrinkage = lambda;
        info.nSamples = nSamples;
        info.nChannels = nCh;
        info.nTrials = nTrials;
        info.baselineRange = [baseStart, baseEnd];
        
        if nSamples < 5 * nCh
            warning('SEAL:NoiseCov:Undersampled', ...
                '样本数(%d)/通道数(%d) 比值偏低,估计可能不稳定', ...
                nSamples, nCh);
        end
        
        if info.rank_empirical < nCh
            warning('SEAL:NoiseCov:RankDeficient', ...
                '经验协方差秩亏 (rank=%d / %d),通常源于平均参考、坏导插值或 ICA 剔除', ...
                info.rank_empirical, nCh);
        end
    end
end