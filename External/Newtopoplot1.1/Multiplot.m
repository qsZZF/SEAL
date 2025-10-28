function axesset = Multiplot(numplots)
axesset = {};
maxrownum = 6;
maxcolnum = 6;
maxnum = maxrownum*maxcolnum;
SIZEBOX = 200*max(1,1+log10(maxnum/numplots));
colnum = ceil(sqrt(numplots));
colnum = max(colnum,maxcolnum);
rownum = ceil(numplots/colnum);
rownum = max(rownum,maxrownum);

nbfig = ceil(numplots/maxnum);
for figind=1:nbfig
    curplotnum = (numplots>=figind*maxnum)*maxnum+((numplots<figind*maxnum))*mod(numplots,maxnum);
    rowcols = [round(sqrt(curplotnum)),ceil(sqrt(curplotnum))];
%     if figind>1
%         rowcols = [rownum colnum];
%     end
    curfig = figure;
    pos = get(curfig,'Position');
    posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
    posy = max(0,pos(2)+pos(4)-SIZEBOX*rowcols(1));
    set(curfig,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
    for i=1:curplotnum
        curaxes = subplot(rowcols(1),rowcols(2),i,'Parent',curfig);
        axesset{end+1} = curaxes;
    end
end

end