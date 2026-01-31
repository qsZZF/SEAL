function onKey(ax,e)
switch lower(e.Key)
    case 'r'
        axis(ax,'vis3d'); view(ax,[-90 90]); camlookat(ax);
    case 'l'
        camlight(ax,'headlight','infinite');
end
end