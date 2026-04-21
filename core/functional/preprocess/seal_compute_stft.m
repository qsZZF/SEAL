function [S, F, T, info] = seal_compute_stft(data3d, fs, winName, winLen, ovlpPct, nfft, opts)
%SEAL_COMPUTE_STFT 手动实现的 STFT,不依赖 Signal Processing Toolbox 的 spectrogram
%
%   [S, F, T] = seal_compute_stft(data3d, fs, winName, winLen, ovlpPct, nfft)
%   [S, F, T, info] = seal_compute_stft(data3d, fs, winName, winLen, ovlpPct, nfft, opts)
%
%   Inputs:
%     data3d  - ch×T 或 ch×T×Tr 数据(2D 会自动升为 3D 处理)
%     fs      - 采样率 (Hz)
%     winName - 'hann' | 'hamming' | 'blackman' | 'gaussian' | 'rect'
%     winLen  - 窗长(样本数)
%     ovlpPct - 重叠比例 [0, 100)
%     nfft    - FFT 长度(≥ winLen)
%     opts    - 可选结构体:
%               .output       - 'complex'(默认,返回复数 STFT)| 'psd'(返回 PSD)
%               .precision    - 'double'(默认)| 'single'
%               .gaussAlpha   - gaussian 窗的 alpha 参数(默认 2.5)
%               .memWarnGB    - 内存预警阈值(GB,默认 4)
%
%   Outputs:
%     S    - 复数 STFT(归一化后): nF × nSeg × nCh × nTr
%            或 PSD(μV²/Hz 若输入单位为 μV): nF × nSeg × nCh × nTr
%     F    - 频率向量 (Hz): nF × 1
%     T    - 时间向量 (s): nSeg × 1, 每段中心时刻
%     info - 计算元信息(窗、归一化因子、丢弃样本数等)
%
%   Notes:
%     - 复数 STFT 已归一化,使得 |S|^2 后 ×2(除 DC/Nyquist)得到 PSD
%     - 具体公式: S_norm = fft(x.*w) / sqrt(fs * sum(w.^2))
%                 PSD = |S_norm|^2, 且 PSD(2:end-1) *= 2

    % ===== 可选参数处理 =====
    if nargin < 7, opts = struct(); end
    if ~isfield(opts, 'output'),     opts.output     = 'complex'; end
    if ~isfield(opts, 'precision'),  opts.precision  = 'double'; end
    if ~isfield(opts, 'gaussAlpha'), opts.gaussAlpha = 2.5; end
    if ~isfield(opts, 'memWarnGB'),  opts.memWarnGB  = 4; end

    % ===== 类型防御(只对参数,不强制转换数据) =====
    fs     = double(fs);
    winLen = double(winLen);
    nfft   = double(nfft);
    ovlpPct = double(ovlpPct);

    % ===== 维度归一化:2D → 3D =====
    if ndims(data3d) ~= 3
        [nCh, nT] = size(data3d);
        data3d = reshape(data3d, [nCh, nT, 1]);
    end
    [nCh, nT, nTr] = size(data3d);

    % ===== 参数校验 =====
    if ~isfinite(fs) || fs <= 0
        error('SEAL:STFT:InvalidFs', '采样率必须为正,当前 fs=%g', fs);
    end
    if winLen < 1 || winLen > nT
        error('SEAL:STFT:WindowTooLong', ...
            '窗长 (%d) 非法或超过信号长度 (%d)', winLen, nT);
    end
    if nfft < winLen
        error('SEAL:STFT:NFFTTooSmall', ...
            'NFFT (%d) 不能小于窗长 (%d)', nfft, winLen);
    end
    if ovlpPct < 0 || ovlpPct >= 100
        error('SEAL:STFT:InvalidOverlap', ...
            '重叠比例必须在 [0, 100),当前 %.2f%%', ovlpPct);
    end
    if ~ismember(opts.output, {'complex', 'psd'})
        error('SEAL:STFT:InvalidOutput', ...
            'opts.output 必须是 ''complex'' 或 ''psd''');
    end
    if ~ismember(opts.precision, {'double', 'single'})
        error('SEAL:STFT:InvalidPrecision', ...
            'opts.precision 必须是 ''double'' 或 ''single''');
    end

    % ===== 构造窗函数 =====
    win = build_window_func(winName, winLen, opts.gaussAlpha);
    win = win(:);   % 列向量
    winEnergy = sum(win.^2);

    % ===== 计算步长和段数 =====
    noverlap = floor((ovlpPct / 100) * winLen);
    hop      = winLen - noverlap;
    nSeg     = floor((nT - winLen) / hop) + 1;

    if nSeg < 1
        error('SEAL:STFT:NoSegment', ...
            '按当前参数(winLen=%d, hop=%d, nT=%d)无法生成任何完整段', ...
            winLen, hop, nT);
    end

    % ===== 末尾丢弃样本提醒 =====
    lastCovered = (nSeg - 1) * hop + winLen;
    nDropped    = nT - lastCovered;
    if nDropped > 0
        if nDropped >= 0.1 * nT
            warning('SEAL:STFT:SamplesDropped', ...
                '末尾 %d 个样本 (%.1f%% 的信号) 未被覆盖,可能丢失晚期成分。', ...
                nDropped, 100 * nDropped / nT);
        end
    end

    % ===== 频率轴(单边谱) =====
    nF = floor(nfft / 2) + 1;
    F  = (0:nF-1).' * (fs / nfft);

    % ===== 时间轴(每段中心) =====
    T = ((0:nSeg-1) * hop + winLen / 2).' / fs;

    % ===== 内存预估与预警 =====
    bytesPerElem = 16;  % complex double 默认
    if strcmp(opts.precision, 'single')
        bytesPerElem = 8;   % complex single
    end
    estBytes = nF * nSeg * nCh * nTr * bytesPerElem;
    if estBytes > opts.memWarnGB * 1e9
        warning('SEAL:STFT:LargeMemory', ...
            ['STFT 输出预计占用 %.2f GB 内存(超过 %.1f GB 预警阈值)\n' ...
             '建议:增大 hop(降低重叠)、减小 nfft、或对数据先降采样。\n' ...
             '或使用 opts.precision=''single'' 减半内存。'], ...
            estBytes / 1e9, opts.memWarnGB);
    end

    % ===== 归一化因子 =====
    % 复数 STFT 的归一化:除以 sqrt(fs * winEnergy)
    % 使得后续 |S|^2 直接等于 PSD(未×2单边谱补偿前)
    normFactor = 1 / sqrt(fs * winEnergy);

    % ===== 段索引矩阵(所有段的样本索引) =====
    % segIdx(i, k) = 第 k 段中第 i 个样本在原信号中的索引
    segIdx = (1:winLen).' + (0:nSeg-1) * hop;   % winLen × nSeg
    win3d  = cast(win, opts.precision);          % winLen × 1,与精度对齐

    % ===== 预分配输出 =====
    S = complex(zeros(nF, nSeg, nCh, nTr, opts.precision));

    % ===== 主循环(矩阵化:每次一个 ch×tr 的组合) =====
    for tr = 1:nTr
        for ch = 1:nCh
            % 提取该通道该 trial 的所有段:winLen × nSeg
            segments = data3d(ch, segIdx, tr);       % 1 × winLen × nSeg
            segments = reshape(segments, winLen, nSeg);

            % 精度对齐
            if ~strcmp(class(segments), opts.precision)
                segments = cast(segments, opts.precision);
            end

            % 加窗(隐式扩展)
            segments = segments .* win3d;            % winLen × nSeg

            % 一次 FFT 所有段
            X = fft(segments, nfft, 1);              % nfft × nSeg

            % 取单边谱并归一化
            S(:, :, ch, tr) = X(1:nF, :) * normFactor;
        end
    end

    % ===== 可选:输出 PSD 而不是复数 STFT =====
    if strcmp(opts.output, 'psd')
        S = abs(S).^2;
        % 单边谱补偿:除 DC 和 Nyquist 外的频点乘 2
        S(2:end-1, :, :, :) = 2 * S(2:end-1, :, :, :);
    end

    % ===== 元信息 =====
    if nargout > 3
        info.fs             = fs;
        info.winName        = winName;
        info.winLen         = winLen;
        info.win            = win;
        info.winEnergy      = winEnergy;
        info.nfft           = nfft;
        info.hop            = hop;
        info.overlapPct     = ovlpPct;
        info.nSegments      = nSeg;
        info.nDroppedSamples= nDropped;
        info.normFactor     = normFactor;
        info.outputType     = opts.output;
        info.precision      = opts.precision;
        info.freqResolution = fs / nfft;           % Hz
        info.timeResolution = winLen / fs;         % s(单段覆盖时长)
        info.hopDuration    = hop / fs;            % s(相邻段时间间隔)
    end
end


function w = build_window_func(name, N, gaussAlpha)
%BUILD_WINDOW_FUNC 构造 STFT 窗函数(periodic,适合分段重建)
    switch lower(name)
        case {'hann', 'hanning'}
            w = hann(N, 'periodic');
        case 'hamming'
            w = hamming(N, 'periodic');
        case 'blackman'
            w = blackman(N, 'periodic');
        case 'gaussian'
            % alpha 越大窗越窄(频率分辨率更好,时间定位更差)
            w = gausswin(N, gaussAlpha);
        case {'rect', 'rectwin', 'boxcar'}
            w = rectwin(N);
        otherwise
            warning('SEAL:STFT:UnknownWindow', '未知窗 ''%s'',使用 hann', name);
            w = hann(N, 'periodic');
    end
end