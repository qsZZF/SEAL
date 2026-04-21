function [data_out, fs_out] = seal_applyDownsampleResample(data_in, fs_in, D)
    % ✅ 类型防御
    fs_in = double(fs_in);
    data_in = double(data_in);
    
    sz = size(data_in);
    nCh = sz(1);
    nT  = sz(2);
    nTr = size(data_in, 3);  % 显式取第 3 维,2D 时自动返回 1

    % 降采样因子 → 有理数 p/q
    if abs(D - round(D)) < 1e-9
        p = 1; q = round(D);
    else
        [p, q] = rat(1/D, 1e-6);
        if q > 200
            % 用 round(D) 作为更稳的近似
            q_approx = max(2, round(D));
            p_approx = 1;
            actual_D = q_approx / p_approx;
            warning('SEAL:Resample:InexactRatio', ...
                '原始 D=%.4f 的有理近似分母过大(q=%d),已退化为 D=%.4f', ...
                D, q, actual_D);
            p = p_approx;
            q = q_approx;
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