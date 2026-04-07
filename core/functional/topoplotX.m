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
ax2plot = gca; % 默认获取当前坐标区
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

isEmptyValues = isempty(Values); 
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
labels_char = char(labels(plotchans));
if ~isEmptyValues
    Values      = Values(plotchans);
end
intrad      = min(1.0,max(Rd)*1.02);   % default: just outside the outermost electrode location

%% Find plotting channels
pltchans = find(Rd <= plotrad); % plot channels inside plotting circle

%% Squeeze channel locations to <= headrad
squeezefac = headrad/plotrad;
Rd    = Rd*squeezefac;      % squeeze electrode arc_lengths towards the vertex
x     = x*squeezefac;
y     = y*squeezefac;
h = ax2plot; % 绑定目标坐标区

if isEmptyValues
    cla(h); hold(h,'on');
    AXHEADFAC = 1.05;
    set(h,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC);
    circ = linspace(0,2*pi,CIRCGRID);
    rx = sin(circ); ry = cos(circ);
    hwidth = HEADRINGWIDTH;
    hin  = headrad*(1- hwidth/2);
    headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
    heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
    patch(h,headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'hittest','off');
    sf    = headrad/plotrad;
    offset = hwidth / sf;
    % Plot ears and nose
    base  = headrad-.0046+offset*7;
    basex = 0.18*headrad;                   % nose width
    tip   = 1.15*headrad+offset*7;
    tiphw = .04*headrad;                    % nose tip half width
    tipr  = .01*headrad;                    % nose tip rounding
    q = .04; % ear lengthening
    EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; 
    EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
    EarX = EarX + offset*8;
    plot3(h,[basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf, ...
        2*ones(5,1),'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');                 
    plot3(h,EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')    
    plot3(h,-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')   
    x_in = x(pltchans); y_in = y(pltchans);
    plot3(h, y_in, x_in, ones(size(x_in)), '.', 'Color', [0 0 0], 'markersize', 5, 'linewidth', .5, 'hittest','off');
    text(h, double(y_in), double(x_in), {chanlocs(plotchans(pltchans)).labels}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'hittest','off');
    epos   = [x_in; y_in];
    handle = [];            
    axis(h,'equal'); axis(h,'off');
    return;
end

intchans = find(x <= intrad & y <= intrad); 
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
xi   = linspace(xmin,xmax,GRID_SCALE);   
yi   = linspace(ymin,ymax,GRID_SCALE);   
[Xi,Yi,Zi] = griddata(inty,intx,intValues,yi',xi,'v4'); 

%% Mask out data outside the head
mask = (sqrt(Xi.^2 + Yi.^2) <= headrad); 
Zi(mask == 0) = NaN;                  
grid = plotrad; %#ok<NASGU>
delta = xi(2)-xi(1); 

%% Scale the axes and make the plot
cla(h)
hold(h,'on');
AXHEADFAC = 1.05;     
set(h,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC);
unsh = (GRID_SCALE+1)/GRID_SCALE; 
if strcmp(SHADING,'interp')
    handle = surface(h,Xi*unsh,Yi*unsh,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
else
    handle = surface(h,Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
end
contour(h,Xi,Yi,Zi,CONTOURNUM,'k','hittest','off');

% 【修复1】将 colormap 绑定到特定的坐标区 h，防止污染其他窗口
colormap(h, jet);

%% Plot filled ring to mask jagged grid boundary
hwidth = HEADRINGWIDTH;                   
hin  = squeezefac*headrad*(1- hwidth/2);  
if strcmp(SHADING,'interp')
    rwidth = BLANKINGRINGWIDTH*1.3;       
else
    rwidth = BLANKINGRINGWIDTH;
end
rin    =  headrad*(1-rwidth/2);           
if hin>rin
    rin = hin;                            
end
circ = linspace(0,2*pi,CIRCGRID);
rx = sin(circ);
ry = cos(circ);
ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];

% 【修复2】不再使用 get(gcf,'color')，直接使用白色 [1 1 1] 进行边缘遮罩，App中更安全
patch(h,ringx,ringy,0.01*ones(size(ringx)),[1 1 1],'edgecolor','none','hittest','off'); 

%% Plot cartoon head, ears, nose
headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
patch(h,headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'hittest','off'); 
sf    = headrad/plotrad;
offset = hwidth / sf; 

base  = headrad-.0046+offset;
basex = 0.18*headrad;                   
tip   = 1.15*headrad+offset;
tiphw = .04*headrad;                    
tipr  = .01*headrad;                    
q = .04; 
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; 
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
EarX = EarX + offset;
plot3(h,[basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf, ...
      2*ones(5,1),'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');
plot3(h,EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')
plot3(h,-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')

%% Mark electrode locations
if strcmp(ELECTRODES,'on')   
    plot3(h,y,x,ones(size(x)),'.','Color',[0 0 0],'markersize',5,'linewidth',.5,'hittest','off');
elseif strcmp(ELECTRODES,'labels')  
    text(h,double(y),double(x),{chanlocs(pltchans).labels},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'hittest','off')
elseif strcmp(ELECTRODES,'numbers')
    for i = 1:size(labels_char,1)
        ht = text(h, double(y(i)),double(x(i)),1,int2str(allchansind(i)),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0]);
        set(ht,'hittest','off');
    end
end
epos=[x; y];

% 【修复3】安全地关闭坐标轴显示
axis(h, 'off');
axis(h, 'equal');
end