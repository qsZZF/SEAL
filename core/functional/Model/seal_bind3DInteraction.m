function bind3DInteraction(fig, ax, selectFcn)

% 关闭内置交互，避免抢事件
try, disableDefaultInteractivity(ax); end
zoom(fig,'off'); pan(fig,'off'); rotate3d(fig,'off');

% 命中设置：点在轴或子对象上都能触发
set(ax,'HitTest','on','PickableParts','all');
kids = allchild(ax);
set(kids, 'HitTest','on','PickableParts','all', ...
    'ButtonDownFcn', @(s,e) startDrag(fig, ax, s, e));

% 在轴上按下也能触发
set(ax,'ButtonDownFcn', @(s,e) startDrag(fig, ax, s, e));

% 全局抬起/滚轮
set(fig,'WindowButtonUpFcn',    @(s,e) stopDrag(fig), ...
    'WindowScrollWheelFcn', @(s,e) onScroll(ax,e));

% 可选：快捷键 r 复位、l 补光
set(fig,'KeyPressFcn',@(s,e) onKey(ax,e));

% ===== 在 fig 的 appdata 里保存“选点功能”相关信息 =====
if nargin >= 3 && ~isempty(selectFcn)
    setappdata(fig,'SelectionHandler', selectFcn);
else
    setappdata(fig,'SelectionHandler', []);
end

% 默认关闭“选点模式”
if ~isappdata(fig,'SelectionMode')
    setappdata(fig,'SelectionMode', false);
end
   
end