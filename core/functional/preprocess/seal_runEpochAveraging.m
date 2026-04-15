function dataOut = seal_runEpochAveraging(dataIn, method)
    % dataIn: ch x time x trials 的三维脑电数据
    % method: 'mean' (传统平均), 'median' (中位数), 'trim' (修剪平均)
    
    if ndims(dataIn) < 3 || size(dataIn, 3) == 1
        warning('当前数据只有 1 个 trial 或已经是二维连续数据，无法进行叠加平均。');
        dataOut = dataIn;
        return;
    end
    
    switch lower(method)
        case 'mean'
            % 1. 传统均值 (速度最快，标准做法)
            dataOut = mean(dataIn, 3);
            
        case 'median'
            % 2. 中位数平均 (极度抗噪，完全无视极端伪迹，但波形可能略显起伏)
            dataOut = median(dataIn, 3);
            
        case 'trim'
            % 3. 修剪平均 / 截尾均值 (业界强烈推荐的前沿做法)
            % 剔除每个时间点上最高 10% 和最低 10% 的极值后，对剩下的 80% 求均值
            % 融合了 mean 的平滑性和 median 的稳健性
            dataOut = trimmean(dataIn, 20, 'round', 3);
            
        otherwise
            error('未知的叠加平均方法: %s', method);
    end
end