function [S, F, T] = seal_compute_stft(data3d, fs, winName, winLen, ovlpPct, nfft)
%SEAL_COMPUTE_STFT 手动实现的 STFT，不依赖 spectrogram
    
    % 类型防御
    fs = double(fs);
    data3d = double(data3d);
    winLen = double(winLen);
    nfft = double(nfft);
    
    if ndims(data3d) ~= 3
        [nCh, nT] = size(data3d);
        data3d = reshape(data3d, [nCh, nT, 1]);
    end
    
    [nCh, nT, nTr] = size(data3d);
    
    % 参数校验
    if winLen > nT
        error('STFT:WindowTooLong', '窗长 (%d) 不能大于信号长度 (%d)', winLen, nT);
    end
    if nfft < winLen
        error('STFT:NFFTTooSmall', 'NFFT (%d) 不能小于窗长 (%d)', nfft, winLen);
    end
    
    % 构造窗函数
    win = build_window_func(winName, winLen);
    win = win(:);   % 列向量
    
    % 计算步长和段数
    noverlap = floor((ovlpPct / 100) * winLen);
    hop = winLen - noverlap;
    nSeg = floor((nT - winLen) / hop) + 1;
    
    % 频率轴（单边谱，real 信号）
    nF = floor(nfft/2) + 1;
    F = (0:nF-1)' * (fs / nfft);
    
    % 时间轴：每段的中心时刻
    T = ((0:nSeg-1) * hop + winLen/2) / fs;
    T = T(:);
    
    % 预分配
    S = complex(zeros(nF, nSeg, nCh, nTr));
    
    % 主循环
    for ch = 1:nCh
        for tr = 1:nTr
            x = data3d(ch, :, tr);
            
            % 分段加窗 FFT
            for seg = 1:nSeg
                startIdx = (seg - 1) * hop + 1;
                endIdx = startIdx + winLen - 1;
                
                segment = x(startIdx:endIdx).' .* win;   % 列向量加窗
                
                % FFT（自动零填充到 nfft）
                X = fft(segment, nfft);
                
                % 只取单边谱
                S(:, seg, ch, tr) = X(1:nF);
            end
        end
    end
end

function w = build_window_func(name, N)
    switch lower(name)
        case {'hann', 'hanning'}
            w = hann(N, 'periodic');
        case 'hamming'
            w = hamming(N, 'periodic');
        case 'blackman'
            w = blackman(N, 'periodic');
        case 'gaussian'
            w = gausswin(N);
        case {'rect', 'rectwin', 'boxcar'}
            w = rectwin(N);
        otherwise
            warning('SEAL:STFT:UnknownWindow', '未知窗 %s，使用 hann', name);
            w = hann(N, 'periodic');
    end
end