function startDrag(fig, ax, src, event)
% fig:  figure 句柄
% ax :  axes / UIAxes 句柄
% src: 被点击的对象（可能是 patch / axes / 其他）
% event: 事件结构体，patch 点击时包含 IntersectionPoint

    % ===== ① 先处理“选点模式” =====

    isSelMode = fig.UserData;
  
 
    if isSelMode
        st = get(fig,'SelectionType');   % 'normal','extend','alt','open'
        % 只在左键点击 patch 时触发选点
        if (strcmp(st,'normal') || strcmp(st,'alt'))&& strcmp(get(src,'Type'),'patch')
            selHandler = [];
            if isappdata(fig,'SelectionHandler')
                selHandler = getappdata(fig,'SelectionHandler');
            end

            if ~isempty(selHandler)
                % 调用你在 App 里传进来的 vertexClicked
                selHandler(src, event);
            end
            % ⚠ 不再挂 WindowButtonMotionFcn，不发生旋转/平移
            return;
        end
    end

    % ===== ② 否则，走原来的旋转/平移逻辑 =====

    % 1) 记录起点（屏幕像素坐标）
    md.pt  = get(0,'PointerLocation');

    % 2) 判定模式：中键（'extend'）或 Shift+左键 => pan；否则 rotate
    st = get(fig,'SelectionType');           % 'normal','extend','alt','open'
    cm = get(fig,'CurrentModifier');         % cellstr，如 {'shift','control'} or {}
    isShift = ismember('shift', cm);

    if strcmp(st,'extend') || any(isShift)
        md.mode = 'pan';
    else
        md.mode = 'rotate';
    end

    % 3) 存状态，并挂上拖动回调
    setappdata(fig,'md', md);
    set(fig,'WindowButtonMotionFcn', @(s,e) doDrag(fig, ax));
end
