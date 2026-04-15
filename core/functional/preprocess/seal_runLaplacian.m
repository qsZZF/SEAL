function dataOut = seal_runLaplacian(dataIn, chanLocs, m, lambda, useCSD, edgeMode)
dataOut = dataIn;
% 统一获取 3D 笛卡尔坐标
P = ensureXYZ(chanLocs);
% 若可用 CSD/EEGLAB：优先使用（Perrin 球面样条，效果最佳）
if useCSD
    try


        % ========================================================
        % 分支 1：如果检测到 Kayser & Tenke CSD Toolbox
        % ========================================================
        if exist('CSD','file')==2 && exist('GetGH','file')==2
            % 将笛卡尔坐标转回极坐标系 (cart2sph 返回弧度)
            [th_rad, ph_rad, ~] = cart2sph(P(:,1), P(:,2), P(:,3));

            % 动态构建 CSD Toolbox 专属的 Montage 结构体 (M)
            % CSD 严格要求字段名为 theta 和 phi，且单位必须是【角度 (Degrees)】
            M = struct();
            for i = 1:numel(chanLocs)
                M(i).theta = rad2deg(th_rad(i));
                M(i).phi   = rad2deg(ph_rad(i));
                if isfield(chanLocs, 'labels') && ~isempty(chanLocs(i).labels)
                    M(i).label = chanLocs(i).labels;
                else
                    M(i).label = sprintf('E%d', i);
                end
            end

            % 计算 G 和 H 矩阵 (注意：GetGH 内部默认 m=4，也可手动传入)
            [G, H] = GetGH(M, m);

            % 逐 trial 运算
            for tr = 1:size(dataIn,3)
                V = double(squeeze(dataIn(:,:,tr)));  % CSD 通常要求 double 精度
                dataOut(:,:,tr) = CSD(V, G, H, lambda); % 【修复】：补上了缺失的 H 矩阵
            end
            return; % 成功则退出，跳过兜底算法

            % ========================================================
            % 分支 2：如果没有 CSD 库，但有 EEGLAB 自带的 laplacian_perrin
            % ========================================================
        elseif exist('laplacian_perrin','file')==2
            X = P(:,1); Y = P(:,2); Z = P(:,3);
            for tr = 1:size(dataIn,3)
                V = double(squeeze(dataIn(:,:,tr)));
                % EEGLAB 自带函数的接口不需要 G/H 矩阵，直接吃 X,Y,Z
                dataOut(:,:,tr) = laplacian_perrin(V, X, Y, Z, m, lambda);
            end
            return; % 成功则退出
        end

    catch ME
        % 若高级算法因为某些奇葩坐标或内存不足报错，拦截错误并交给下方的兜底算法
        warning('高级 CSD/Laplacian 库调用失败，将降级使用近邻差分兜底。错误信息: %s', ME.message);
    end
end

tri = convhull(P(:,1), P(:,2), P(:,3));
nCh = size(P,1);

% 2. 构建顶点邻接矩阵 (Adjacency Matrix)
adj = false(nCh, nCh);
for i = 1:size(tri, 1)
    v1 = tri(i,1); v2 = tri(i,2); v3 = tri(i,3);
    adj(v1, v2) = true; adj(v2, v1) = true;
    adj(v1, v3) = true; adj(v3, v1) = true;
    adj(v2, v3) = true; adj(v3, v2) = true;
end

% 3. 找出边缘电极 (邻居数少于3的视为边缘)
neighborCounts = sum(adj, 2);
edgeList = find(neighborCounts < 3);

% 4. 计算近邻差分拉普拉斯
for tr = 1:size(dataIn,3)
    V = double(squeeze(dataIn(:,:,tr)));  % ch × time
    L = zeros(size(V));

    for i = 1:nCh
        % 获取当前通道的所有邻居索引
        neiIdx = find(adj(i,:));
        if isempty(neiIdx)
            L(i,:) = 0;
            continue;
        end
        % 局部均值差（Laplacian 近似）
        mu = mean(V(neiIdx, :), 1);
        L(i,:) = V(i,:) - mu;
    end

    % 处理边缘电极
    if strcmpi(edgeMode,'shrink') && ~isempty(edgeList)
        L(edgeList,:) = 0;
    end
    dataOut(:,:,tr) = L;
end
end


function P = ensureXYZ(chanLocs)
n = numel(chanLocs);
P = nan(n,3);
hasXYZ = isfield(chanLocs,'X') && isfield(chanLocs,'Y') && isfield(chanLocs,'Z') ...
    && all(~cellfun(@isempty,{chanLocs.X}));
if hasXYZ
    for i=1:n, P(i,:)=[chanLocs(i).X, chanLocs(i).Y, chanLocs(i).Z]; end
else
    % 若只有极坐标 theta/phi（度）：转 XYZ
    for i=1:n
        th = deg2rad(chanLocs(i).theta);  % 水平角（EEGLAB 定义）
        ph = deg2rad(chanLocs(i).phi);    % 垂直角
        r  = 1;

        [x,y,z] = sph2cart(th, ph, r);
        P(i,:) = [x,y,z];
    end
end
% 单位化到单位球
P = P ./ vecnorm(P,2,2);
end


