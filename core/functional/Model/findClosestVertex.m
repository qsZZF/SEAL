function [vertexIndex, minDist] = findClosestVertex(vertices, point)
% 计算所有顶点到点击点的距离
dists = sum((vertices - point).^2, 2);

% 找到最近的顶点
[minDist, vertexIndex] = min(dists);
end