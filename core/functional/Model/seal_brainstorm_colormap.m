function cmap = seal_brainstorm_colormap(nColors)
%SEAL_BRAINSTORM_COLORMAP  Sequential source colormap similar to Brainstorm.
%
% Brainstorm displays source maps with an absolute-valued, sequential
% blue-purple-to-yellow colormap by default. This compact version uses a few
% control colors rather than copying Brainstorm's full table.

if nargin < 1 || isempty(nColors)
    nColors = 256;
end

anchors = [ ...
    0.04 0.06 0.61; ...
    0.18 0.12 0.55; ...
    0.36 0.16 0.51; ...
    0.56 0.15 0.47; ...
    0.72 0.16 0.34; ...
    0.81 0.28 0.12; ...
    0.83 0.62 0.05];

x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, nColors);
cmap = zeros(nColors, 3);
for k = 1:3
    cmap(:, k) = interp1(x, anchors(:, k), xi, 'pchip')';
end
cmap = min(max(cmap, 0), 1);

end
