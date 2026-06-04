function [Valpha, cdata] = seal_bulidsource(SourceData, vc, rgb, colorMap, ax, varargin)
    if ~isempty(rgb)
        cortexcolor = rgb;
    else
        cortexcolor = [0.85 0.85 0.85];
    end
    SourceData = double(SourceData(:));
    cdata  = repmat(cortexcolor, length(SourceData), 1);
    roiIdx = 1:length(vc);
    CLim   = [0, 1];
    alphaValue = 1;
    dataThreshold = [];
    isExplicitThreshold = false;

    for iArg = 1:2:length(varargin)
        Param = lower(varargin{iArg});
        Value = varargin{iArg + 1};
        switch Param
            case {'clim', 'datalimit', 'datalimits'}
                CLim = double(Value(:))';
            case {'threshold', 'thresh', 'datathreshold'}
                dataThreshold = double(Value);
                isExplicitThreshold = true;
            case {'valuemode', 'displaymode'}
                if ischar(Value) || isstring(Value)
                    mode = lower(char(Value));
                    if any(strcmp(mode, {'abs', 'absolute', 'amplitude', 'brainstorm'}))
                        SourceData = abs(SourceData);
                    end
                end
            case 'cortexcolor'
                cortexcolor = Value;
                cdata = repmat(cortexcolor, length(SourceData), 1);
        end
    end

    if nargin < 4 || isempty(colorMap)
        colorMap = CBar('active');
        if ~isempty(SourceData) && min(SourceData) < 0
            colorMap = CBar('pnactive');
        end
    end
    if nargin < 5, ax = []; end

    if ~isempty(SourceData)
        if isempty(CLim) || numel(CLim) ~= 2
            CLim = [0 1];
        end
        if isequal(CLim, [0 1])
            if min(SourceData) < 0
                CLim = [-1 1] * max(abs(SourceData));
            else
                CLim = CLim * max(abs(SourceData));
            end
        end
        cmin   = min(CLim); cmax = max(CLim);
        if cmax <= cmin
            cmax = cmin + eps;
            CLim = [cmin cmax];
        end
        tlen   = size(colorMap,1);

        plotData = SourceData;
        if isExplicitThreshold
            dataThreshold = max(0, min(1, dataThreshold));
            if cmin < 0 && abs(cmin + cmax) <= eps(max(abs(CLim)))
                threshValue = dataThreshold * max(abs(CLim));
                plotData(abs(plotData) < threshValue) = 0;
            else
                threshValue = cmin + (cmax - cmin) * dataThreshold;
                plotData(plotData < threshValue) = 0;
                plotData(plotData > cmax) = cmax;
            end
        else
            thresh = otsu(abs(SourceData)/max(abs(SourceData)+eps));
            if cmin >= 0
                colorMap(1:floor(tlen*thresh),:) = repmat(cortexcolor, floor(tlen*thresh), 1);
            else
                lo = floor(tlen*(0.5 - thresh/2)) + 1;
                hi = floor(tlen*(0.5 + thresh/2));
                colorMap(lo:hi,:) = repmat(cortexcolor, hi-lo+1, 1);
            end
        end

        s2plot = roiIdx;
        if isExplicitThreshold
            s2plot = s2plot(plotData(s2plot) ~= 0 & isfinite(plotData(s2plot)));
        end
        if cmax == cmin
            idx = ones(size(plotData(s2plot)));
        else
            idx = floor((plotData(s2plot)-cmin)/(cmax-cmin)*(tlen-1)) + 1;
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
    case {'brainstorm', 'source'}
        cbars = seal_brainstorm_colormap(256);
        return;
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
