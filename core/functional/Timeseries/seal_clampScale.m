     function f = clampScale(f)
    %CLAMPSCALE 防止缩放系数过大或过小
    %
    % 策略：以所有通道的中位数缩放系数为基准，
    %        单个通道不允许偏离中位数太远
    
    medF = median(f);
    
    if medF == 0
        medF = 1;
    end
    
    % 允许的范围：中位数的 [1/maxRatio, maxRatio] 倍
    maxRatio = 5;  % 可以做成 app 属性让用户调
    
    loBound = medF / maxRatio;
    hiBound = medF * maxRatio;
    
    f = max(f, loBound);
    f = min(f, hiBound);
end
