function peakMem = estimatePeakMemory(algName, nCh, nSrc, nT, nO)
    % 估算算法运行时的峰值内存（bytes）
    % 基本单位：一个 double = 8 bytes
    
    B = 8;              % bytes per element
    N = nCh;            % channels
    M = nO * nSrc;      % leadfield columns
    T = nT;             % time points
    
    % 基础：X (N×T), L (N×M), S (M×T)
    baseMem = (N*T + N*M + M*T) * B;
    
    switch algName
        % ========== L2 Tikhonov 类 ==========
        % 峰值：L*L' (N×N) + inv buffer (N×N) + L'*inv (M×N)
        case {'MNE','dSPM','sLORETA'}
            peakMem = baseMem + (2*N^2 + M*N) * B;
            
        % eLORETA 迭代，每步都维持 W (M×M 对角)、L*W*L' (N×N)
        case 'eLORETA'
            peakMem = baseMem + (M + 2*N^2) * B;
            
        % LORETA 需要 Laplacian L_lap (M×M)——这是最大的一项
        case {'LORETA','LAURA','LOGURA'}
            peakMem = baseMem + (M^2 + 2*N^2) * B;   % M×M 稀疏矩阵实际会小一些
            
        % ========== Beamformer ==========
        % 需要 inv(C_data) (N×N)，每个体素一个 filter
        case {'LCMV','MCMV'}
            peakMem = baseMem + (2*N^2 + M*N) * B;
            
        % ========== Sparse / Bayesian ==========
        % Champagne：维护每个源的 gamma (M×1)，sigma_y (N×N) 逆，源协方差
        case {'Champagne','TS_Champagne','SmoothChampagne','Smooth-Cham'}
            peakMem = baseMem + (M + 3*N^2 + M*N) * B;
            
        % BlockChampagne：块结构，内存和普通 Champagne 接近
        case 'BlockChampagne'
            peakMem = baseMem + (M + 3*N^2 + M*N) * B;
            
        % BESTIES：VB + MRF + TBF；主要缓存源空间 MxM 后验近似和若干 NxN 传感器矩阵
        case 'BESTIES'
            k = 4;   % 默认 NumTBFs
            peakMem = baseMem + (3*M^2 + (k+3)*N^2 + M*N + M*k) * B;
            
        % STARTS / sSTARTS / uSTAR：TBF 分解降维到 k 个基
        case {'STARTS','sSTARTS','uSTAR'}
            k = 4;   % 默认 NumTBFs，可从 prm 读更精确
            peakMem = baseMem + (M*k + k*T + 2*N^2) * B;
            
        % ========== L1 / L21 ==========
        % MxNE / irMxNE：迭代 proximal，需要维持 S_prev、gradient buffer
        case {'MxNE','irMxNE','debiased_GroupLasso'}
            peakMem = baseMem + (3*M*T + N*T) * B;

        % FOCUSS: source buffers plus a dense sensor-by-sensor weighted system
        case 'FOCUSS'
            peakMem = baseMem + (3*M*T + 2*N^2 + M*N) * B;
            
        case {'IRES','FASTIRES'}
            peakMem = baseMem + (4*M*T + M^2) * B;  % 含 Laplacian

        % VSSI-Lp: source, edge-domain primal/dual buffers, and sensor solve
        case {'VSSI_Lp','VSSI_GGD'}
            peakMem = baseMem + (12*M*T + 2*M*N + N^2) * B;

        case 'STOUT'
            peakMem = baseMem + (4*M*T + M*N + N^2) * B;
            
        case {'SISSY_L21','dMAP'}
            peakMem = baseMem + (3*M*T + M^2) * B;

        case 'dSPN'
            % Low-memory state-space solve: source trajectories plus Kalman gain buffers.
            peakMem = baseMem + (4*M*T + 3*M*N + 4*N^2 + M) * B;
            
        case 'TFMxNE'
            % 时频字典 K 倍膨胀，取 K=10 估
            K = 10;
            peakMem = baseMem + (K*M*T + 2*N*T) * B;
            
        case 'Lq_ADMM'
            peakMem = baseMem + (3*M*T + N^2) * B;
            
        otherwise
            % 未知算法，用保守估计：baseMem 的 3 倍
            peakMem = baseMem * 3;
    end
    
    % 加上 MATLAB 基础开销 (~500 MB) 和 20% 余量
    peakMem = peakMem * 1.2 + 500e6;
end
