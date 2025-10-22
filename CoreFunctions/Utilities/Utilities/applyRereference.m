<<<<<<< HEAD
<<<<<<< HEAD
function dataOut = applyRereference(dataIn, refSpec, method)
% data: ch x time x trials
    [nCh,~,~] = size(dataIn);
    dataOut = dataIn; % 先复制

    switch lower(method)
        case 'direct'
            fMean = @(X,dim) mean(X, dim);
        case 'median'
            fMean = @(X,dim) median(X, dim);
        case 'trimmed'
            % 20% 截尾均值
            fMean = @(X,dim) trimmean(X, 40, dim);  % 40=两侧各20%
        otherwise
            error('未知 Method：%s（可选：direct/median/trimmed）', method);
    end

    switch refSpec.type
        case 'average-all'
            % 1×time×trials 的参考信号
            ref = fMean(dataIn, 1);
            dataOut = dataIn - ref;   % 隐式扩展，所有通道都减同一参考

        case 'channels'
            idx = refSpec.idx;
            checkIndexRange(idx, nCh);
            ref = fMean(dataIn(idx,:,:), 1);
            dataOut = dataIn - ref;

        case 'groups'
            % 对每一组：组内通道减该组均值；不在任何组的通道保持不变
            covered = false(nCh,1);
            for g = 1:numel(refSpec.groups)
                gg = refSpec.groups{g};
                checkIndexRange(gg, nCh);
                covered(gg) = true;
                ref_g = fMean(dataIn(gg,:,:), 1);  % 1×T×Tr
                dataOut(gg,:,:) = dataIn(gg,:,:) - ref_g;
            end
            if ~all(covered)
                warnCh = find(~covered);
                fprintf('[Re-ref] 以下通道未在任何组里：'); fprintf('%d ', warnCh); fprintf('\n');
                % 如需对未分组的通道也做全局平均，可在此添加逻辑
            end

        otherwise
            error('未知 refSpec.type：%s', refSpec.type);
    end
end


function checkIndexRange(idx, nCh)
    assert(all(idx>=1 & idx<=nCh), '索引越界（1..%d）。', nCh);
end

=======
function dataOut = applyRereference(dataIn, refSpec, method)
% data: ch x time x trials
    [nCh,~,~] = size(dataIn);
    dataOut = dataIn; % 先复制

    switch lower(method)
        case 'direct'
            fMean = @(X,dim) mean(X, dim);
        case 'median'
            fMean = @(X,dim) median(X, dim);
        case 'trimmed'
            % 20% 截尾均值
            fMean = @(X,dim) trimmean(X, 40, dim);  % 40=两侧各20%
        otherwise
            error('未知 Method：%s（可选：direct/median/trimmed）', method);
    end

    switch refSpec.type
        case 'average-all'
            % 1×time×trials 的参考信号
            ref = fMean(dataIn, 1);
            dataOut = dataIn - ref;   % 隐式扩展，所有通道都减同一参考

        case 'channels'
            idx = refSpec.idx;
            checkIndexRange(idx, nCh);
            ref = fMean(dataIn(idx,:,:), 1);
            dataOut = dataIn - ref;

        case 'groups'
            % 对每一组：组内通道减该组均值；不在任何组的通道保持不变
            covered = false(nCh,1);
            for g = 1:numel(refSpec.groups)
                gg = refSpec.groups{g};
                checkIndexRange(gg, nCh);
                covered(gg) = true;
                ref_g = fMean(dataIn(gg,:,:), 1);  % 1×T×Tr
                dataOut(gg,:,:) = dataIn(gg,:,:) - ref_g;
            end
            if ~all(covered)
                warnCh = find(~covered);
                fprintf('[Re-ref] 以下通道未在任何组里：'); fprintf('%d ', warnCh); fprintf('\n');
                % 如需对未分组的通道也做全局平均，可在此添加逻辑
            end

        otherwise
            error('未知 refSpec.type：%s', refSpec.type);
    end
end


function checkIndexRange(idx, nCh)
    assert(all(idx>=1 & idx<=nCh), '索引越界（1..%d）。', nCh);
end

>>>>>>> origin
=======
function dataOut = applyRereference(dataIn, refSpec, method)
% data: ch x time x trials
    [nCh,~,~] = size(dataIn);
    dataOut = dataIn; % 先复制

    switch lower(method)
        case 'direct'
            fMean = @(X,dim) mean(X, dim);
        case 'median'
            fMean = @(X,dim) median(X, dim);
        case 'trimmed'
            % 20% 截尾均值
            fMean = @(X,dim) trimmean(X, 40, dim);  % 40=两侧各20%
        otherwise
            error('未知 Method：%s（可选：direct/median/trimmed）', method);
    end

    switch refSpec.type
        case 'average-all'
            % 1×time×trials 的参考信号
            ref = fMean(dataIn, 1);
            dataOut = dataIn - ref;   % 隐式扩展，所有通道都减同一参考

        case 'channels'
            idx = refSpec.idx;
            checkIndexRange(idx, nCh);
            ref = fMean(dataIn(idx,:,:), 1);
            dataOut = dataIn - ref;

        case 'groups'
            % 对每一组：组内通道减该组均值；不在任何组的通道保持不变
            covered = false(nCh,1);
            for g = 1:numel(refSpec.groups)
                gg = refSpec.groups{g};
                checkIndexRange(gg, nCh);
                covered(gg) = true;
                ref_g = fMean(dataIn(gg,:,:), 1);  % 1×T×Tr
                dataOut(gg,:,:) = dataIn(gg,:,:) - ref_g;
            end
            if ~all(covered)
                warnCh = find(~covered);
                fprintf('[Re-ref] 以下通道未在任何组里：'); fprintf('%d ', warnCh); fprintf('\n');
                % 如需对未分组的通道也做全局平均，可在此添加逻辑
            end

        otherwise
            error('未知 refSpec.type：%s', refSpec.type);
    end
end


function checkIndexRange(idx, nCh)
    assert(all(idx>=1 & idx<=nCh), '索引越界（1..%d）。', nCh);
end

>>>>>>> 3afa99beaec32453db712d4621fd05217f53fcfd
