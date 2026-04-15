function [dataOut,removedIdx] = Seal_applyRereference(dataIn, refSpec, method)
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
            
            % 计算参考信号
            if numel(idx) == 1
                ref = dataIn(idx, :, :);
            else
                ref = fMean(dataIn(idx,:,:), 1);
            end
            
            % 所有通道减参考
            dataOut = dataIn - ref;
            
            % 剔除参考通道本身（它们已经变成零）
            keepMask = true(nCh, 1);
            keepMask(idx) = false;
            dataOut = dataOut(keepMask, :, :);
            removedIdx = idx(:)';   % 返回被剔除的索引
            
            fprintf('[Re-ref] 已剔除 %d 个参考通道: %s\n', ...
                numel(idx), mat2str(idx));
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
                uncovered = find(~covered);
                % 明确策略：未分组通道用全局平均作为兜底
                warning('SEAL:Rereference:UncoveredChannels', ...
                    '%d 个通道未在任何组内，将使用全局平均作为参考: %s', ...
                    length(uncovered), mat2str(uncovered));
                global_ref = fMean(dataIn, 1);
                dataOut(uncovered,:,:) = dataIn(uncovered,:,:) - global_ref;
            end

        otherwise
            error('未知 refSpec.type：%s', refSpec.type);
    end
end


function checkIndexRange(idx, nCh)
    assert(all(idx>=1 & idx<=nCh), '索引越界（1..%d）。', nCh);
end

