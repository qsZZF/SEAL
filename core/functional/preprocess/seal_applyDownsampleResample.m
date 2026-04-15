function [data_out, fs_out] = seal_applyDownsampleResample(data_in, fs_in, D)
    % ✅ 类型防御
    fs_in = double(fs_in);
    data_in = double(data_in);
    
    [nCh, nT, nTr] = size(data_in);
    
    % 降采样因子 → 有理数 p/q
    if abs(D - round(D)) < 1e-9
        p = 1; q = round(D);
    else
        [p, q] = rat(1/D, 1e-6);
        if q > 200
            q = max(2, round(D)); p = 1;
            warning('非整数降采样因子已近似为 D=%d', q);
        end
    end
    
    fs_out = fs_in * (p/q);     % 现在 fs_in 是 double，结果也是 double
    
    % 预跑一次拿尺寸
    x0 = data_in(1, :, 1);
    y0 = resample(x0, p, q);
    nT_new = numel(y0);
    
    data_out = zeros(nCh, nT_new, nTr);
    
    % 主循环
    for ch = 1:nCh
        for tr = 1:nTr
            x = data_in(ch, :, tr);
            y = resample(x, p, q);
            nt = min(nT_new, numel(y));
            data_out(ch, 1:nt, tr) = y(1:nt);
        end
    end
end