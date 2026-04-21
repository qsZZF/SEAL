function [dataOut, removedIdx, info] = seal_applyRereference(dataIn, refSpec, method, opts)
%SEAL_APPLYREREFERENCE 对 EEG 数据应用重参考
%
%   [dataOut, removedIdx] = seal_applyRereference(dataIn, refSpec, method)
%   [dataOut, removedIdx, info] = seal_applyRereference(dataIn, refSpec, method, opts)
%
%   Inputs:
%     dataIn  - ch×T 或 ch×T×Tr 数据
%     refSpec - 参考规范结构体,含字段:
%               .type  - 'average-all' | 'channels' | 'groups'
%               .idx   - (当 type='channels')参考通道索引向量
%               .groups- (当 type='groups')cell 数组,每个 cell 是一个通道组
%     method  - 'direct'(均值)| 'median' | 'trimmed'
%     opts    - 可选结构体:
%               .removeRefChans - 是否剔除参考通道本身(默认 true,仅对 'channels' 生效)
%               .trimPct        - trimmed 方法的修剪百分比(默认 20,表示上下各去 10%)
%
%   Outputs:
%     dataOut    - 重参考后的数据(类型与输入相同)
%     removedIdx - 被剔除的通道索引(原始通道空间)
%     info       - 元信息结构体

    % ===== 默认 opts =====
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'removeRefChans'), opts.removeRefChans = true; end
    if ~isfield(opts, 'trimPct'),        opts.trimPct        = 20; end

    % ===== 输入防御 =====
    if isempty(dataIn) || ~isnumeric(dataIn)
        error('SEAL:Reref:InvalidInput', '数据为空或类型不合法');
    end

    if opts.trimPct < 0 || opts.trimPct >= 100
        error('SEAL:Reref:InvalidTrimPct', ...
            'trimPct 必须在 [0, 100),当前 %g', opts.trimPct);
    end

    % ===== 类型防御:非浮点升为 double(EEG 罕见 int,但保险) =====
    origClass = class(dataIn);
    if ~isfloat(dataIn)
        dataIn = double(dataIn);
    end

    if ismatrix(dataIn)
        dataIn = reshape(dataIn, size(dataIn,1), size(dataIn,2), 1);
    end
    [nCh, nT, nTr] = size(dataIn);

    dataOut    = dataIn;
    removedIdx = [];

    % ===== 选择聚合函数 =====
    switch lower(method)
        case 'direct'
            fAgg = @(X, dim) mean(X, dim);
        case 'median'
            fAgg = @(X, dim) median(X, dim);
        case 'trimmed'
            pct = opts.trimPct;
            fAgg = @(X, dim) trimmean(X, pct, dim);
        otherwise
            error('SEAL:Reref:UnknownMethod', ...
                '未知 Method:%s(可选:direct/median/trimmed)', method);
    end

    % ===== 按参考类型分派 =====
    switch refSpec.type
        case 'average-all'
            ref = fAgg(dataIn, 1);             % 1×T×Tr
            dataOut = dataIn - ref;

        case 'channels'
            idx = refSpec.idx(:).';
            checkIndexRange(idx, nCh);

            if numel(idx) == 1
                ref = dataIn(idx, :, :);
            else
                ref = fAgg(dataIn(idx, :, :), 1);
            end
            dataOut = dataIn - ref;

            % 是否剔除参考通道本身(由用户的 checkbox 决定)
            if opts.removeRefChans
                keepMask = true(nCh, 1);
                keepMask(idx) = false;
                dataOut = dataOut(keepMask, :, :);
                removedIdx = idx;
                fprintf('[Re-ref] 已剔除 %d 个参考通道:%s\n', ...
                    numel(idx), mat2str(idx));
            else
                fprintf('[Re-ref] 保留 %d 个参考通道(其值应近似为 0):%s\n', ...
                    numel(idx), mat2str(idx));
            end

        case 'groups'
            covered = false(nCh, 1);
            for g = 1:numel(refSpec.groups)
                gg = refSpec.groups{g}(:).';
                checkIndexRange(gg, nCh);

                % 组内是否重叠?
                if any(covered(gg))
                    warning('SEAL:Reref:OverlappingGroups', ...
                        '第 %d 组和之前的组有重叠通道:%s', ...
                        g, mat2str(intersect(find(covered), gg)));
                end
                covered(gg) = true;

                ref_g = fAgg(dataIn(gg, :, :), 1);
                dataOut(gg, :, :) = dataIn(gg, :, :) - ref_g;
            end

            % 严格模式:未覆盖通道直接报错,不做"全局平均补救"
            if ~all(covered)
                uncovered = find(~covered);
                error('SEAL:Reref:UncoveredChannels', ...
                    ['以下 %d 个通道未被任何组覆盖,无法确定其参考方案:\n%s\n' ...
                     '请补全分组(推荐把它们作为独立一组或添加到现有组),' ...
                     '或改用 ''average-all'' / ''channels'' 方案。'], ...
                    numel(uncovered), mat2str(uncovered));
            end

        otherwise
            error('SEAL:Reref:UnknownType', ...
                '未知 refSpec.type:%s', refSpec.type);
    end

    % ===== 还原维度(3D → 2D 若输入是 2D) =====
    if ismatrix(dataIn)
        dataOut = dataOut(:, :, 1);
    end

    % ===== 还原原始类型 =====
    if ~strcmp(class(dataOut), origClass) && ~isfloat(dataIn)
        % 只有在原本是整型时才转回(浮点的 single/double 保持计算精度)
        dataOut = cast(dataOut, origClass);
    end

    % ===== 元信息 =====
    if nargout > 2
        info.type        = refSpec.type;
        info.method      = method;
        info.trimPct     = opts.trimPct;
        info.removedIdx  = removedIdx;
        info.nChIn       = nCh;
        info.nChOut      = size(dataOut, 1);
        info.wasMatrix   = ismatrix(dataIn);
    end
end


function checkIndexRange(idx, nCh)
    if isempty(idx)
        error('SEAL:Reref:EmptyIndex', '参考通道索引不能为空');
    end
    if ~all(isfinite(idx)) || any(idx ~= round(idx))
        error('SEAL:Reref:NonIntegerIndex', '参考通道索引必须是正整数');
    end
    if any(idx < 1) || any(idx > nCh)
        error('SEAL:Reref:IndexOutOfRange', ...
            '索引越界(有效范围 1..%d),当前:%s', nCh, mat2str(idx));
    end
    if numel(unique(idx)) ~= numel(idx)
        error('SEAL:Reref:DuplicateIndex', ...
            '参考通道索引有重复:%s', mat2str(idx));
    end
end