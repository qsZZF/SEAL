function L_fixed = seal_reduceLeadfield(L, nd, normals)
    %REDUCELEADFIELD 将自由朝向导联场 [N×3M] 降维为固定朝向 [N×M]
    %
    % 输入:
    %   L       - 导联场矩阵
    %             若为 [N×3M]，执行降维
    %             若为 [N×M]，直接返回
    %   normals - (可选) 皮层法向量 [M×3]
    %             提供时用法线投影，否则用 SVD 最优朝向
    %
    % 输出:
    %   L_fixed - 固定朝向导联场 [N×M]
    %
    % 示例:
    %   L_fixed = reduceLeadfield(L);              % SVD 模式
    %   L_fixed = reduceLeadfield(L, normals);     % 法线投影模式

    [nChan, nCol] = size(L);

    % 固定朝向，无需处理
    if nd == 1
        L_fixed = L;
        return;
    end
    nSource = nCol/3;
    % ========== 降维核心 ==========
    L_fixed = zeros(nChan, nSource);
    useNormals = (nargin >= 3) && ~isempty(normals);

    if useNormals
        % ---------- 法线投影模式 ----------
        % 向量化: 把 L [N×3M] reshape 成 [N×3×M]，再与 normals 做批量点积
        L_3d = reshape(L, nChan, 3, nSource);          % [N × 3 × M]
        normals_t = permute(normals, [3 2 1]);          % [1 × 3 × M]
        normals_t = normals_t ./ vecnorm(normals_t, 2, 2); % 归一化
        L_fixed = squeeze(sum(L_3d .* normals_t, 2));   % [N × M]
        disp(['Using the academic default setting: Assumes the dipole orientation is perpendicular to the cortical surface.', newline, newline, ...
      'Note: This is a widely used physiological approximation that well represents the macroscopic discharge direction of pyramidal cells.']);
    else
        % ---------- SVD 最优朝向模式 ----------
        for k = 1:nSource
            cols = (k-1)*3 + (1:3);
            L_vertex = L(:, cols);                      % [N × 3]
            [~, ~, V] = svd(L_vertex, 'econ');
            L_fixed(:, k) = L_vertex * V(:, 1);        % 投影到第一主成分
        end
        disp(['Anatomical normal vectors not detected/used. Automatically switched to SVD mathematical optimal orientation.', newline, newline, ...
                   'Note: This mode extracts the principal component direction that generates the maximum scalp voltage. It maximizes signal energy but lacks real physiological structural constraints.']);
       
end