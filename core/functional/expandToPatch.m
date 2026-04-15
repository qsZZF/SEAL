function [patchVertices, patchWeights] = expandToPatch(seedVertices, cortexCoords, patchArea, decayModel, powerDecay)
    %EXPANDTOPATCH 将点源扩展为斑块源
    %
    % 输入:
    %   seedVertices  - 种子顶点 ID 数组 [1×K]
    %   cortexCoords  - 皮层网格所有顶点坐标 [V×3]，单位 mm
    %   patchArea     - 斑块面积 mm²
    %   decayModel    - 衰减模型: 'Gaussian' | 'Linear' | 'Flat'
    %   powerDecay    - 边缘功率衰减百分比 (如 10 表示边缘功率为中心的 10%)
    %
    % 输出:
    %   patchVertices - 斑块包含的所有顶点 ID [1×P]
    %   patchWeights  - 对应的幅度权重 [1×P]，中心=1，边缘按衰减模型递减


    coordRange = max(abs(cortexCoords(:)));

    if coordRange < 1
        % 单位是米，乘 1000 转 mm
        scaleFactor = 1000;
    elseif coordRange < 30
        % 单位是厘米，乘 10 转 mm
        scaleFactor = 10;
    else
        % 已经是毫米，不缩放
        scaleFactor = 1;
    end

    coords_mm = cortexCoords * scaleFactor;
    % 面积 → 半径 (假设近似圆形斑块)
    radius = sqrt(patchArea / pi);

    % 功率衰减 → 边缘幅度比
    % 功率 ∝ 幅度²，所以幅度比 = sqrt(功率比)
    edgeAmp = sqrt(powerDecay / 100);

    % ========== 搜索邻域顶点 ==========
    allPatchIdx = [];
    allPatchDist = [];

    for i = 1:length(seedVertices)
        seedCoord = coords_mm(seedVertices(i), :);

        % 欧氏距离（对小斑块足够精确，大斑块建议换测地距离）
        dists = vecnorm(coords_mm - seedCoord, 2, 2);

        % 筛选半径内的顶点
        inPatch = find(dists <= radius);
        allPatchIdx = [allPatchIdx; inPatch(:)];
        allPatchDist = [allPatchDist; dists(inPatch(:))];
    end

    % 去重：多个种子的斑块可能重叠，取距离最小值
    [patchVertices, ~, ic] = unique(allPatchIdx);
    patchDist = accumarray(ic, allPatchDist, [], @min);

    % ========== 计算衰减权重 ==========
    normDist = patchDist / radius;  % 归一化到 [0, 1]
    normDist = min(normDist, 1);    % 钳位

    switch lower(decayModel)
        case 'gaussian'
            % 高斯衰减：中心=1，边缘=edgeAmp
            % w = exp(-d²/(2σ²))，令 w(1)=edgeAmp → σ² = -1/(2·ln(edgeAmp))
            if edgeAmp > 0 && edgeAmp < 1
                sigma2 = -1 / (2 * log(edgeAmp));
                patchWeights = exp(-normDist.^2 / (2 * sigma2));
            else
                patchWeights = ones(size(normDist));
            end

        case 'linear'
            % 线性衰减：中心=1，边缘=edgeAmp
            patchWeights = 1 - (1 - edgeAmp) * normDist;

        case 'flat'
            % 无衰减：斑块内全部等权
            patchWeights = ones(size(normDist));

        otherwise
            error('SEAL:expandToPatch:UnknownModel', ...
                '未知的衰减模型: %s', decayModel);
    end

    % 转行向量
    patchVertices = patchVertices(:)';
    patchWeights = patchWeights(:)';
end