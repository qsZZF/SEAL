function [Valpha, cdata] = seal_bulidsource(SourceData, vc, rgb, colorMap, ax)
    if ~isempty(rgb)
        cortexcolor = rgb;
    else
        cortexcolor = [0.85 0.85 0.85];
    end
    cdata  = repmat(cortexcolor, length(SourceData), 1);
    roiIdx = 1:length(vc);
    CLim   = [0, 1];
    alphaValue = 1;

    if nargin < 4 || isempty(colorMap)
        colorMap = CBar('active');
        if ~isempty(SourceData) && min(SourceData) < 0
            colorMap = CBar('pnactive');
        end
    end
    if nargin < 5, ax = []; end

    if ~isempty(SourceData)
        if min(SourceData) < 0
            CLim = [-1 1];
        end
        thresh = otsu(abs(SourceData)/max(abs(SourceData)+eps));
        CLim   = CLim * max(abs(SourceData));
        cmin   = min(CLim); cmax = max(CLim);
        tlen   = size(colorMap,1);

        if cmin >= 0
            colorMap(1:floor(tlen*thresh),:) = repmat(cortexcolor, floor(tlen*thresh), 1);
        else
            lo = floor(tlen*(0.5 - thresh/2)) + 1;
            hi = floor(tlen*(0.5 + thresh/2));
            colorMap(lo:hi,:) = repmat(cortexcolor, hi-lo+1, 1);
        end

        s2plot = roiIdx;
        if cmax == cmin
            idx = ones(size(SourceData(s2plot)));
        else
            idx = floor((SourceData(s2plot)-cmin)/(cmax-cmin)*(tlen-1)) + 1;
        end
        idx = min(max(idx,1), tlen);
        idx(~isfinite(idx)) = 1;
        cdata(s2plot,:) = colorMap(idx,:);

        % 关键:只作用于指定的 axes,避免污染全局 gca
        if ~isempty(ax) && isvalid(ax)
            clim(ax, [cmin cmax]);
            colormap(ax, colorMap);
        else
            caxis([cmin cmax]);
            colormap(colorMap);
        end
    end

    Valpha = ones(size(vc,1), 1) * alphaValue;
end


function cbars = CBar(type)
if nargin<1
    type = 'active';
end
switch type
    case 'active'
        cbars = [1 1 0;1 0.3 0;0.7 0 0];
    case 'iactive'
        cbars = [0.7 0 0;1 0.3 0;1 1 0];
    case 'invactive'
        cbars = [0.5 0 0;0.8 0.4 0;1 1 1];
    case 'mjet'
        cbars = [0.8 0 0;0.8 0.4 0;0.8 0.8 0;0.4 0.8 0.4;0 0.8 0.8;0 0.4 0.8;0 0 0.8];
    case 'pnactive'
        cbars = [1 1 0.5;1 0.6 0;1 0 0;0 0 0;0 0 1;0 0.6 1; 0.5 1 1];
end


cbars = lineInterp(cbars,256,8);
cbars = cbars(end:-1:1,:);

end


function Cbar = lineInterp(cbar,cnum,pernum)
np = size(cbar,1);
nump = max(round(cnum/(np-1)),pernum);
for ic = 1:3
    tcolor = cbar(:,ic);
    ncolor = [];
    for itc = 1:np-1
        ncolor = [ncolor,linspace(tcolor(itc),tcolor(itc+1),nump)];
    end
    Cbar(:,ic) = ncolor';
end
while size(Cbar,1)>cnum
    Cbar(round(end/2),:)=[];
end
end
