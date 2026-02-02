function dataOut = runLaplacian(dataIn, chanLocs, m, lambda, useCSD, edgeMode)
    dataOut = dataIn;

    % 若可用 CSD/EEGLAB：优先使用（Perrin 球面样条，效果最佳）
    if useCSD && (exist('CSD','file')==2 || exist('laplacian_perrin','file')==2 || exist('eeg_laplacian','file')==2)
        try
            P = ensureXYZ(chanLocs);
            G = GetGH(P);               % 若你装了 CSD toolbox：GetGH.m / CSD.m
            % 下面伪代码：根据你装的实现调整调用
            for tr=1:size(dataIn,3)
                V = squeeze(dataIn(:,:,tr));   % ch×time
                Vcsd = CSD(V, G, m, lambda);   % 或 laplacian_perrin(V,P,m,lambda)
                dataOut(:,:,tr) = Vcsd;
            end
            return;
        catch
            warning('CSD/Laplacian 库调用失败，改用近邻差分兜底。');
        end
    end

%     % ---- 兜底：基于最近邻的离散拉普拉斯（球面三角剖分）----
%     P = ensureXYZ(chanLocs);
%     tri = delaunayTriangulation(P);      % 三角化
%     N  = neighbors(tri, (1:size(P,1))'); % cell：每点的邻居索引
% 
%     for tr=1:size(dataIn,3)
%         V = double(squeeze(dataIn(:,:,tr)));  % ch×time
%         L = zeros(size(V));
%         for i=1:size(P,1)
%             nei = N(i,:);
%             if isempty(nei), L(i,:) = 0; continue; end
%             mu  = mean(V(nei,:),1);
%             L(i,:) = V(i,:) - mu;             % 简单的局部均值差（Laplacian 近似）
%         end
%         if strcmpi(edgeMode,'shrink')
%             edgeList = find(cellfun(@numel,N)<3);
%             L(edgeList,:) = 0;
%         end
%         dataOut(:,:,tr) = L;
%     end
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
            [x,y,z] = sph2cart(th, deg2rad(90-chanLocs(i).sph_theta_besa), r); %#ok<*NASGU>
            % 不同数据源角度定义不一致，实际使用时请按你的坐标系修正
            % 这里保底：用简单极角->球坐标替换
            [x,y,z] = sph2cart(th, ph, r);
            P(i,:) = [x,y,z];
        end
    end
    % 单位化到单位球
    P = P ./ vecnorm(P,2,2);
end


