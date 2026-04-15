function vertexClicked(app, src, event)
% 获取点击位置
clickedPoint = event.IntersectionPoint;

%  获取所有顶点坐标
vertices = get(src, 'Vertices');

%  找到最近的顶点
[vertexIndex, ~] = findClosestVertex(vertices, clickedPoint);

isCtrlPressed = ismember('control', get(app.handles1.h, 'CurrentModifier'));
isAltPressed = ismember('alt', get(app.handles1.h, 'CurrentModifier'));

if isAltPressed
    % 如果当前这个点在已选集合中，就移除
    if ismember(vertexIndex, app.selectedVertices)
        app.selectedVertices = setdiff(app.selectedVertices, vertexIndex);
    end

    if isempty(app.selectedVertices)
        % 没有任何点被选中了，清空所有 marker
        if isfield(app.handles1, 'Multimarker') && ~isempty(app.handles1.Multimarker)
            for i = 1:length(app.handles1.Multimarker)
                if isvalid(app.handles1.Multimarker(i))
                    delete(app.handles1.Multimarker(i));
                end
            end
            app.handles1.Multimarker = [];
        end

        if isfield(app.handles1, 'singlemarker')
            delete(app.handles1.singlemarker);
            app.handles1.singlemarker = [];
        end
    else
        % 还有其他点，重新高亮剩余的多选点
        app.selectedVerticesNum = 0;
        highlightMultipleVertices(app,app.selectedVertices);
    end
    return;  % Alt 模式处理完就返回，不走后面的 Ctrl/单选逻辑
end

if isCtrlPressed
    % Ctrl+点击：多选模式

    % 检查是否已经选中了该顶点
    if ~ismember(vertexIndex, app.selectedVertices)
        app.selectedVertices = union(app.selectedVertices, vertexIndex);
    end


    % 高亮所有选中的顶点
    app.selectedVerticesNum = 0;
    highlightMultipleVertices(app,app.selectedVertices);
else
    % 清除现有的多选标记
    if isfield(app.handles1, 'Multimarker') && ~isempty(app.handles1.Multimarker)
        for i = 1:length(app.handles1.Multimarker)
            if isvalid(app.handles1.Multimarker(i))
                delete(app.handles1.Multimarker(i));
            end
        end
        app.handles1.Multimarker = [];
    end
    % 清除旧标记
    if isfield(app.handles1, 'singlemarker')
        delete(app.handles1.singlemarker);
    end

    app.selectedVertices=vertexIndex;
    addVisualMarker(app,vertexIndex);
    disp('vertice selected');
end
%  更新信息显示
% updateSelectionInfo(app, vertexIndex, clickedPoint);
end
function highlightMultipleVertices(app,vertexIndex)
% 获取顶点位置
vertices = get(app.handles1.hp, 'Vertices');
% 清除现有的多选标记
if isfield(app.handles1, 'Multimarker') && ~isempty(app.handles1.Multimarker)
    for i = 1:length(app.handles1.Multimarker)
        if isvalid(app.handles1.Multimarker(i))
            delete(app.handles1.Multimarker(i));
        end
    end
    app.handles1.Multimarker = [];
end

if isfield(app.handles1, 'singlemarker') && ~isempty(app.handles1.singlemarker)
        if isvalid(app.handles1.singlemarker)
            delete(app.handles1.singlemarker);
        end
        app.handles1.singlemarker = [];
end

% 创建球体标记
[x, y, z] = sphere(20);
% 4. 动态计算标记大小（基于大脑尺寸）
brainSize = max(max(vertices) - min(vertices));
r = brainSize * 0.01; % 大脑尺寸的2%

hold(app.handles1.axes, 'on');
for i=1:length(vertexIndex)
    vertexPos = vertices(vertexIndex(i), :);
    app.selectedVerticesNum = app.selectedVerticesNum + 1 ;
    marker(app.selectedVerticesNum) = surf(app.handles1.axes, ...
        x*r + vertexPos(1), ...
        y*r + vertexPos(2), ...
        z*r + vertexPos(3), ...
        'FaceColor', 'r', ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.9, ...
        'Tag', 'SelectionMarker', ...
        'HitTest', 'off', ...
        'PickableParts', 'none');       % 对 uiaxes 生效，彻底透明
end
hold(app.handles1.axes, 'off');
app.handles1.Multimarker = marker;
end

function addVisualMarker(app, vertexIndex)
% 获取顶点位置
vertices = get(app.handles1.hp, 'Vertices');
vertexPos = vertices(vertexIndex, :);

% 创建球体标记
[x, y, z] = sphere(20);
% 4. 动态计算标记大小（基于大脑尺寸）
brainSize = max(max(vertices) - min(vertices));
r = brainSize * 0.01; % 大脑尺寸的2%

hold(app.handles1.axes, 'on');

marker = surf(app.handles1.axes, ...
    x*r + vertexPos(1), ...
    y*r + vertexPos(2), ...
    z*r + vertexPos(3), ...
    'FaceColor', 'r', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.9, ...
    'Tag', 'SelectionMarker', ...
    'HitTest', 'off', ...
    'PickableParts', 'none');       % 对 uiaxes 生效，彻底透明

hold(app.handles1.axes, 'off');
app.handles1.singlemarker = marker;
end