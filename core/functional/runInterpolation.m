function dataOut = runInterpolation(dataIn, chanLocs, badIdx, method, k, radius, useTemplate)
    dataOut = dataIn;
    goodIdx = setdiff(1:numel(chanLocs), badIdx);

    % 拿到 3D 单位球坐标
    P = ensureXYZ(chanLocs);    % ch×3，单位化；见下
    Pg = P(goodIdx,:); Pb = P(badIdx,:);

    switch method
        case 'spherical-spline'
            % 如果有 EEGLAB 的 eeg_interp，用它（最好效果）
            if exist('eeg_interp','file')==2
                % 包个最小 EEG 结构在内存中转一下
                EEGtmp.data = dataIn; EEGtmp.nbchan=size(dataIn,1); EEGtmp.trials=size(dataIn,3);
                EEGtmp.pnts = size(dataIn,2); EEGtmp.srate = 1;
                EEGtmp.chanlocs = chanLocs; EEGtmp.icachansind = 1:EEGtmp.nbchan;
                dataOut = eeg_interp(EEGtmp, badIdx, 'spherical').data;
            else
                % 纯 MATLAB 兜底：基于球面 RBF（thin-plate）近似
                dataOut = interp_rbf_sphere(dataIn, Pg, Pb, goodIdx, badIdx);
            end

        case 'idw'      % inverse distance weighting on sphere
            dataOut = interp_idw_sphere(dataIn, P, goodIdx, badIdx, k, radius);

        case 'nearest'
            dataOut = interp_nearest_sphere(dataIn, P, goodIdx, badIdx);

        otherwise
            error('未知插值方法：%s', method);
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

function dataOut = interp_nearest_sphere(dataIn, P, goodIdx, badIdx)
    dataOut = dataIn;
    % 球面距离≈夹角
    for b = badIdx
        ang = real(acos( max(-1,min(1, P(goodIdx,:)*P(b,:)')) ));
        [~,iMin] = min(ang);
        nn = goodIdx(iMin);
        dataOut(b,:,:) = dataIn(nn,:,:);
    end
end

function dataOut = interp_idw_sphere(dataIn, P, goodIdx, badIdx, k, radiusDeg)
    dataOut = dataIn;
    if ~isfinite(radiusDeg) || radiusDeg<=0, radiusDeg = 180; end
    R = deg2rad(radiusDeg);

    for b = badIdx
        ang = real(acos( max(-1,min(1, P(goodIdx,:)*P(b,:)')) ));
        inR = find(ang <= R);
        if isempty(inR), [~,inR] = mink(ang, k); end
        [~,ord] = sort(ang(inR));
        pick = inR(ord(1:min(k,numel(inR))));

        d = ang(pick) + eps;      % 权重 ~ 1/d
        w = (1./d).'/sum(1./d);
        x = double(reshape(dataIn(goodIdx(pick),:,:), numel(pick), [])); % k × (T*Tr)
        y = w * x;                % 1 × (T*Tr)
        dataOut(b,:,:) = reshape(y, 1, size(dataIn,2), size(dataIn,3));
    end
end

function dataOut = interp_rbf_sphere(dataIn, Pg, Pb, goodIdx, badIdx)
    % 简化的 RBF（thin-plate）在球面上的近似：phi(r)=r^2 log r
    dataOut = dataIn;
    G = acos( max(-1,min(1, Pg*Pg.')) );   % 夹角矩阵
    G(G==0) = eps;
    Phi = (G.^2) .* log(G);

    % 正则化避免病态
    lam = 1e-6;
    A = Phi + lam*eye(size(Phi));

    % 目标点到已知点的核
    B = acos( max(-1,min(1, Pg*Pb.')) );  B(B==0) = eps;
    K = (B.^2) .* log(B);                 % |good| × |bad|

    X = double(reshape(dataIn(goodIdx,:,:), numel(goodIdx), [])); % |good| × (T*Tr)
    W = A \ X;                           % 核权重
    Y = (K.' * W);                       % |bad| × (T*Tr)
    for i=1:numel(badIdx)
        dataOut(badIdx(i),:,:) = reshape(Y(i,:), 1, size(dataIn,2), size(dataIn,3));
    end
end
