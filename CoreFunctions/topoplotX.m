% eeglab topoplot
function [handle,pltchans,epos] = topoplotX(Values,chanlocs,varargin)

%% Set defaults
headrad = 0.5;          % actual head radius - Don't change this!
GRID_SCALE = 67;        % plot map on a 67X67 grid
CIRCGRID   = 201;       % number of angles to use in drawing circles
HEADCOLOR = [0 0 0];    % default head color (black)
HLINEWIDTH = 1.7;       % default linewidth for head, nose, ears
BLANKINGRINGWIDTH = .035;% width of the blanking ring
HEADRINGWIDTH    = .007;% width of the cartoon head ring

CONTOURNUM = 0;
ELECTRODES = 'on';
plotrad = 0.55;
SHADING = 'interp';
ax2plot = gca;

nargs = nargin;
if nargs > 2
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
            case 'plotrad'
                plotrad = Value;
            case 'shading'
                SHADING = lower(Value);
                if ~any(strcmp(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
            case 'axes'
                ax2plot = Value;
        end
    end
end

isEmptyValues = isempty(Values);  %% ==== CHANGED: 记录是否为空，进入仅标注模式
if ~isEmptyValues
    Values = Values(:); % make Values a column vector
end

%% Read channel location
labels={chanlocs.labels};
Th=[chanlocs.theta];
Rd=[chanlocs.radius];

Th = pi/180*Th;                              % convert degrees to radians
allchansind = 1:length(Th);
plotchans = 1:length(chanlocs);

%% remove infinite and NaN values
inds = [];
if ~isEmptyValues
    inds = union(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
end
for chani=1:length(chanlocs)
    if isempty(chanlocs(chani).X); inds = [inds chani]; end
end

plotchans   = setdiff(plotchans,inds);

[x,y]       = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans   = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
Th          = Th(plotchans);
Rd          = Rd(plotchans);
x           = x(plotchans);
y           = y(plotchans);
labels_char = char(labels(plotchans)); % 仅用于 numbers 分支
if ~isEmptyValues
    Values      = Values(plotchans);
end
intrad      = min(1.0,max(Rd)*1.02);   % default: just outside the outermost electrode location

%% Find plotting channels（头内通道）
pltchans = find(Rd <= plotrad); % plot channels inside plotting circle

%% Squeeze channel locations to <= headrad
squeezefac = headrad/plotrad;
Rd    = Rd*squeezefac;      % squeeze electrode arc_lengths towards the vertex
x     = x*squeezefac;
y     = y*squeezefac;

%% ==== CHANGED: Values 为空 → 仅绘制头型 + 电极点 + 标签，直接返回
if isEmptyValues
    h = ax2plot; cla(h); hold(h,'on');
    % 头部外圈
    AXHEADFAC = 1.05;
    set(h,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC);
    circ = linspace(0,2*pi,CIRCGRID);
    rx = sin(circ); ry = cos(circ);

    % 头圈
    hwidth = HEADRINGWIDTH;
    hin  = headrad*(1- hwidth/2);
    headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
    heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
    patch(h,headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'hittest','off');

    % 耳朵 & 鼻子（沿用原比例）
    base  = headrad-.0046;
    basex = 0.18*headrad;   % nose width
    tip   = 1.15*headrad;
    tiphw = .04*headrad;    % nose tip half width
    tipr  = .01*headrad;    % nose tip rounding
    q = .04;                % ear lengthening
    EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % headrad = 0.5
    EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
    sf    = headrad/plotrad;  % 与原实现一致

    plot3(h,[basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf, ...
          2*ones(5,1),'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');                 % nose
    plot3(h,EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')    % left ear
    plot3(h,-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')   % right ear

    % 仅绘制头内通道
    x_in = x(pltchans); y_in = y(pltchans);
    % 强制显示标签
    plot3(h, y_in, x_in, ones(size(x_in)), '.', 'Color', [0 0 0], 'markersize', 5, 'linewidth', .5, 'hittest','off');
    text(h, double(y_in), double(x_in), {chanlocs(plotchans(pltchans)).labels}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'hittest','off');

    epos   = [x_in; y_in];
    handle = [];              % 没有 surface 热图
    axis(h,'equal'); axis(h,'off');
    return;
end
%% ==== CHANGED 结束（以下为原有的插值/着色路径） ====

%% 插值用通道（建议：可将条件改为 abs(x)<=intrad & abs(y)<=intrad 或半径内）
intchans = find(x <= intrad & y <= intrad); % 原代码逻辑保持

allx  = x; ally  = y;
allchansind = allchansind(pltchans);
intTh = Th(intchans);
intRd = Rd(intchans);
intx  = x(intchans);
inty  = y(intchans);
Th    = Th(pltchans);
Rd    = Rd(pltchans);
x     = x(pltchans);
y     = y(pltchans);

intValues = Values(intchans);
Values = Values(pltchans);
labels_char = labels_char(pltchans,:);

%% create grid
xmin = min(-headrad,min(intx)); xmax = max(headrad,max(intx));
ymin = min(-headrad,min(inty)); ymax = max(headrad,max(inty));
xi   = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
yi   = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)

[Xi,Yi,Zi] = griddata(inty,intx,intValues,yi',xi,'v4'); % interpolate data

%% Mask out data outside the head
mask = (sqrt(Xi.^2 + Yi.^2) <= headrad); % mask outside the plotting circle
Zi(mask == 0) = NaN;                  % mask non-plotting voxels with NaNs
grid = plotrad; %#ok<NASGU>
delta = xi(2)-xi(1); % length of grid entry

%% Scale the axes and make the plot
h = ax2plot; % uses current axes
cla(h)
hold(h,'on');

AXHEADFAC = 1.05;     % do not leave room for external ears if head cartoon
set(h,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC);
unsh = (GRID_SCALE+1)/GRID_SCALE; % un-shrink the effects of 'interp' SHADING

if strcmp(SHADING,'interp')
    handle = surface(h,Xi*unsh,Yi*unsh,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
else
    handle = surface(h,Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
end
contour(h,Xi,Yi,Zi,CONTOURNUM,'k','hittest','off');
colormap(jet)

%% Plot filled ring to mask jagged grid boundary
hwidth = HEADRINGWIDTH;                   % width of head ring
hin  = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius

if strcmp(SHADING,'interp')
    rwidth = BLANKINGRINGWIDTH*1.3;       % width of blanking outer ring
else
    rwidth = BLANKINGRINGWIDTH;
end
rin    =  headrad*(1-rwidth/2);           % inner ring radius
if hin>rin
    rin = hin;                            % dont blank inside the head ring
end

circ = linspace(0,2*pi,CIRCGRID);
rx = sin(circ);
ry = cos(circ);
ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];
patch(h,ringx,ringy,0.01*ones(size(ringx)),get(gcf,'color'),'edgecolor','none','hittest','off'); hold on

%% Plot cartoon head, ears, nose
headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
patch(h,headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'hittest','off'); hold on

% Plot ears and nose
base  = headrad-.0046;
basex = 0.18*headrad;                   % nose width
tip   = 1.15*headrad;
tiphw = .04*headrad;                    % nose tip half width
tipr  = .01*headrad;                    % nose tip rounding
q = .04; % ear lengthening
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % headrad = 0.5
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
sf    = headrad/plotrad;

plot3(h,[basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf, ...
      2*ones(5,1),'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');
plot3(h,EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')
plot3(h,-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')

%% Mark electrode locations
if strcmp(ELECTRODES,'on')   % plot electrodes as spots
    plot3(h,y,x,ones(size(x)),'.','Color',[0 0 0],'markersize',5,'linewidth',.5,'hittest','off');
elseif strcmp(ELECTRODES,'labels')  % print electrode names (labels)
    text(h,double(y),double(x),{chanlocs(pltchans).labels},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'hittest','off')
elseif strcmp(ELECTRODES,'numbers')
    for i = 1:size(labels_char,1)
        ht = text(double(y(i)),double(x(i)),1,int2str(allchansind(i)),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0]);
        set(ht,'hittest','off');
    end
end

epos=[x; y];
axis off
axis equal
end
