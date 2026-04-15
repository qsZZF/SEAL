function stopDrag(fig)
set(fig,'WindowButtonMotionFcn',[]);
if isappdata(fig,'md')
    rmappdata(fig,'md');
end
end
