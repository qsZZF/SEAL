function doDrag(fig, ax)
md  = getappdata(fig,'md');
pt  = get(0,'PointerLocation');
dp  = pt - md.pt;
md.pt = pt; setappdata(fig,'md',md);

switch md.mode
    case 'pan'   
        k = 1.0; 
        camdolly(ax, -k*dp(1), -k*dp(2), 0, 'movetarget','pixels');

    otherwise     
        camorbit(ax, -dp(1)*0.2, -dp(2)*0.2, 'camera');
end
drawnow limitrate
end