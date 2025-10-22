<<<<<<< HEAD
<<<<<<< HEAD
function handles = Newtopoplot(Values,chanlocs,eegtopoSet,varargin)
%Newtopoplot plots  EEG channel locations topography/spatial weights topography/ spatial connectivity topography
%   handles =  Newtopoplot(Values,chanlocs,eegtopoSet,varargin);
%       Another two functions arrow & predint needed.
%   Input:
%      -Values(needed): []/ nchan X 1 vector/ nchanT X 2(3) matrix, if Values is
%               empty([]), then Newtopoplot plots EEG channel locations; if Values
%               is nchan by 1 vector(where nchan is the number of channels,length(chanlocs)),
%               then Newtopoplot plots spatial weights topography;
%               if Values is nchanT by 2 or 3 matrix, then Newtopoplot
%               plots connectivity topography, note that nchanT is the number of
%               pairs of electrodes to plot connectivity lines, and the first two
%               columns are the serial number of channels(e.g Values = [2,3]
%               means that channel 2 and 3 have connectivity); the third column(if
%               provided) is the strength of connectivity(e.g Values = [2,3,0.5])
%      -chanlocs(needed): Channel locations structure. Using EEGlab EEG.chanlocs or
%               location structure S containing:
%                   S.labels,S.theta,S.radiu,S.X,S.Y,S.Z
%      -eegtopoSet(needed but the programme load it by default):
%               Settings of topography, which has been provided as 'EEGtopoSet'.
%               Newtopoplot will load it if exists, or you need to provide it
%   Options:
%      -CLim(default [0,1]):color limits of topography, 2 by 1 vector
%      -electrodes(default 'on'): 'off': don't show electrodes; 'on': show
%               electordes; 'numbers': show electrodes and channel serial
%               number; 'text': show electrodes and channel labels
%      -plotchans(default is 1:length(chanlocs)): channel numbers to plot,vector
%      -plotrad(default is 0.55):radius of cartoon head, no need to change in general,float
%      -isInside(default is false): true:plot channels inside head(max(chanlocs.radius)<plotrad)
%                false: plot all channels
%      -headColor(default is [0 0 0]): color of cartoon head and ear; 3 by
%                1 vector.
%      -headLineWidth(default is 2.5): line width of cartoon head; float
%      -electrodeSize(default is 20): size of electrode. float
%      -electrodColor(default is [0 0 0]): color of electrode
%      -textColor(default is [0 0 0]): color of electrode numbers/labels.(see Options:electrodes)
%      -textSize(default is 10): size of electrode numbers/labels
%      -contournum(default is 0): number contour lines in spatial weight topography, int (not recommended)
%      -shading(default is 'interp'): 'interp'/'flat',
%      -colormap(default is jet):colormap, if you plot connectivity with
%               values(nchanT X 3), better specific it so that color of connectivity line is
%               consistent with colormap
%      -isDir(default is false): used when plot connectivity topography.
%               true: plot directed connectivity topography using arrow. The
%                   direction is from the first column to second column
%                   (e.g [2,3], the arrow starts from channel 2 point at channel 3, see Input:Values)
%               false: plot undirected connectivity topography
%      -LineWidth(defalut is 2): connectivity line width, float
%      -LineColor(defalut is [0 0 0]): connectivity line color, used when no connectivity strength
%               specified(nchanT X 2, see Input:Values); if strength has
%               been specified(nchanT X 3), then line color is calculated
%               using strength according to colormap(see Options:colormap)
%      -Axes: default is gca; plot topography on specific axis if provided.
%   Examples:
%       plot channel inside head locations with labels:
%           figure;Newtopoplot([],chanlocs,[],'electrodes','labels','isInside',true);
%       plot spatial weight topography
%           nchan = length(chanlocs);Values = rand(nchan,1);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','labels','isInside',false);
%       plot spatial weight topography in region of interest(ROI)
%           nchan = length(chanlocs); ROI = 1:10;% channel 1-10 as ROI
%           Values = nan(nchan,1); Values(ROI) = rand(ROI,1);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','labels','isInside',false);
%           Note: Values = nan(nchan,1), then Newtopoplot only plots values inside ROI;
%               if you want to plot values in the whole head, then specific
%               values outside ROI first, e.g Values = zeros(nchan,1);
%       plot undirected connectivity topography without strength values
%           nchan = length(chanlocs); Values = randi(nchan,12,2);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','on','LineWidth',3,'LineColor',[1 0 0]);
%       plot directed connectivity topography with strength values
%           nchan = length(chanlocs); Values = randi(nchan,12,3);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','on','LineWidth',3,'colormap','jet');

% load coordinates
if nargin<3||isempty(eegtopoSet)
    try
        path2file = which('Newtopoplot');
        path2set = fileparts(path2file);
        path2set = fullfile(path2set,'EEGtopoSet');
        load(path2set,'eegtopoSet');
    catch
        error('No topograph setting found');
    end
end


% default setting
CONTOURNUM = 0;
ELECTRODES = 'on';
isInside = false;
isDir = false;
plotrad = 0.55;
SHADING = 'interp';
headColor = [0 0 0];
electrodeColor = [0 0 0];
textColor = [0 0 0];
lineColor = [0 0 0];
LineWidth = 2;
CLim = [0 1];
colorMap = CBar;
electrodeSize = 20;
headLineWidth = 1.5;
plotchans = 1:length(chanlocs);
textSize = 10;
ax2plot = gca;

nargs = nargin;
if nargs > 3
    for iter = 1:2:length(varargin)
        Param = lower(varargin{iter});
        Value = varargin{iter+1};
        switch Param
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
            case 'plotrad'
                plotrad = Value;
            case 'shading'
                SHADING = Value;
                if ~any(strcmpi(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
            case 'headcolor'
                headColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'electrodecolor'
                electrodeColor = Value;
                if size(Value,2)~=3
                    error('Electrode color must be a 1 x 3 matrix')
                end
            case 'textcolor'
                textColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'linewidth'
                LineWidth = Value;
            case 'linecolor'
                lineColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'colormap'
                colorMap = Value;
                if size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
            case 'electrodesize'
                electrodeSize = Value;
            case 'headlinewidth'
                headLineWidth = Value;
            case 'clim'
                CLim = Value;
                if numel(CLim)~=2
                    error('Color Limit must be 1 x 2 vector');
                end
            case 'isinside'
                isInside = Value;
                if ~islogical(isInside)
                    error('isInside must be true or false');
                end
            case 'isdirection'
                isDir = Value;
                if ~islogical(isDir)
                    error('isInside must be true or false');
                end
            case 'plotchans'
                plotchans = Value;
            case 'textsize'
                textSize = Value;
            case 'axes'
                ax2plot = Value;
        end
    end
    
    
end
handles.axes = ax2plot;
colormap(colorMap);
%% Read channel location
Th=[chanlocs.theta];
Rd=[chanlocs.radius];
Rd = Rd*plotrad/0.55;
Th = pi/180*Th;                              % convert degrees to radians
allchansind = 1:length(Th);

%% remove infinite and NaN values
inds = [];
if ~isempty(Values)
    if size(Values,1)==1||size(Values,2)==1
        Values = Values(:);
    end
    inds = union(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
end
for chani=1:length(chanlocs)
    if isempty(chanlocs(chani).X); inds = [inds chani]; end
end
if isInside
    outs = find(Rd>=0.5495);
    plotchans = setdiff(plotchans,outs);
end
plotchans   = setdiff(plotchans,inds);

[x,y]       = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans   = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
x           = x(plotchans);
y           = y(plotchans);
Rd          = Rd(plotchans);
% project the coordinates of electrodes

yy = x*62.6; xx = y*62.6; yyt = xx;  yt = sqrt(344.^2-xx.^2).*sign(yy);
yid = yy>0; yid2 = yy<=0;
yyt(yid) = mean(predint(eegtopoSet.hhead,xx(yid)+409),2)-355.5;
yyt(yid2) = mean(predint(eegtopoSet.lhead,xx(yid2)+409),2)-355.5;
yy = yy.*(yyt./yt);
yy(isnan(yy)) = 0;
xx = xx*10;yy = yy*10;


k = 1;
%% Plot Head
% cla(ax2plot)
hold(ax2plot,'on');
plot(ax2plot,eegtopoSet.Chead(:,1)/k,eegtopoSet.Chead(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cleft(:,1)/k,eegtopoSet.Cleft(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cright(:,1)/k,eegtopoSet.Cright(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cnose(:,1)/k,eegtopoSet.Cnose(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);

if strcmp(ELECTRODES,'off')
    plot(ax2plot,xx/k,yy/k,'.','MarkerSize',0.001,'Color','none');
else
    plot(ax2plot,xx/k,yy/k,'.','MarkerSize',electrodeSize,'Color',electrodeColor);
    if strcmp(ELECTRODES,'labels')
        text(ax2plot,xx/k,yy/k+1.2*electrodeSize,{chanlocs(plotchans).labels},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',textColor,'hittest','off','FontName','Arial','FontSize',textSize);
    end
    if strcmp(ELECTRODES,'numbers')
        text(ax2plot,xx/k,yy/k+1.2*electrodeSize,num2str(allchansind(:)),'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',textColor,'hittest','off','FontName','Arial','FontSize',textSize);
        
    end
end
Rratio = max(1,max(Rd)/0.55);
set(ax2plot,'Xtick',[],'Ytick',[],'DataAspectRatio',[0.9,1,1],'xlim',[-450/k,450/k]*Rratio,...
    'Ylim',[-420/k,520/k]*Rratio,'clim',CLim,'Xcolor',[0.9 0.95 1],'Ycolor',[0.9 0.95 1],'Zcolor',[0.9 0.95 1],'visible','off');
% axis off
set(ax2plot.Parent,'Color',[1 1 1]);
if isempty(Values)
    title([num2str(length(plotchans)),' of ',num2str(length(Th)),' channels']);
end
if ~isempty(Values)
    
    if size(Values,2)==1
        % Plot topo weight/energy values
        Values = Values(plotchans,:);
        if ~isInside && max(Rd)>=0.5495
            % if some channel to plot is outside the head
            [mRd,mid] = max(Rd);
            ratio = mRd/0.5495;
            
            xi = round(-344*ratio*1.1):5:round(344*ratio*1.1);
            yi = round(-349*ratio*1.1)-5:5:round(446*ratio*1.1);
            
            [Xi,Yi,Zi] = griddata(xx,yy,Values,xi,yi','v4');
%             Zi = smooth2a(Zi,5,5);
            fXi = Xi(1,:);
            [Y_upLimit,Y_lowLimit] = outEllip(abs(xx(mid)),(yy(mid)),fXi);
            mask = Yi>repmat(Y_upLimit,size(Xi,1),1)|Yi<repmat(Y_lowLimit,size(Xi,1),1);
            Zi(mask==1) = NaN;
            [Xi,Yi,Zi] = sdownsample(Xi,Yi,Zi,k);
            [Y_upLimit,Y_lowLimit,fXi] = vdownsample(Y_upLimit,Y_lowLimit,fXi,k);
            ringx = [[fXi,fXi(end:-1:1),fXi(1)]*1.05,[fXi,fXi(end:-1:1),fXi(1)]*0.97];
            ringy = [[Y_upLimit*1.1,Y_lowLimit(end:-1:1)*1.1,Y_upLimit(1)*1.1],...
                [Y_upLimit*0.98,Y_lowLimit(end:-1:1)*0.96,Y_upLimit(1)*0.98]];
           
        else
            xi = [max(min(xx),-344),min(max(xx),344)];
            yi = [min(yy),max(yy)];
            x0 = mean(xi); y0 = mean(yi);
            xit = linspace(((xi(1)-x0)*sqrt(2)*1.3+x0),((xi(end)-x0)*sqrt(2)*1.3+x0),...
                max(150,round((xi(2)-xi(1))/1)));
            
            yit = linspace(((yi(1)-y0)*sqrt(2)*1.3+y0),((yi(end)-y0)*sqrt(2)*1.3+y0),...
                max(150,round((yi(2)-yi(1))/1)));
            [Xi,Yi,Zi] = griddata(xx,yy,Values,xit,yit','v4');
%             Zi = smooth2a(Zi,3,3);
            fXi = Xi(1,:);
            outId = abs(fXi)>344;
            [tup,tlow] = inEllip(x0,y0,1.1*(max(abs(xi-x0)))+x0,1.1*(max(abs(yi-y0)))+y0,fXi);
            Y_upLimit = mean(predint(eegtopoSet.hhead,fXi+409),2)-355.5;
            Y_lowLimit = mean(predint(eegtopoSet.lhead,fXi+409),2)-355.5;
            Y_lowLimit(outId) = 0;
            Y_upLimit(outId) = 0;
            Y_upLimit = min(Y_upLimit(:),tup(:));
            Y_lowLimit = max(Y_lowLimit(:),tlow(:));
            Y_upLimit = Y_upLimit';
            Y_lowLimit = Y_lowLimit';
            
            tup = Y_upLimit(~outId);
            tlow = Y_lowLimit(~outId);
            fXi = fXi(~outId);
            fXi = fXi-x0;
            fXi(fXi<0) = fXi(fXi<0)+0;
            [Xi,Yi,Zi] = sdownsample(Xi,Yi,Zi,k);
            [tup,tlow,fXi] = vdownsample(tup,tlow,fXi,k);
            ringx = [[Xi(1,:),Xi(1,end:-1:1),Xi(1,1)],[fXi,fXi(end:-1:1),fXi(1)]+x0];
            ringy = [[Yi(1,:),Yi(end,:),Yi(1,1)],...
                [y0+(tup-y0),y0+(tlow(end:-1:1)-y0),y0+(tup(1)-y0)]];
            
            
        end
        
        handles.hsuf = surface(ax2plot,Xi,Yi,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING,...
            'FaceLighting','gouraud');
        
        patch(ax2plot,ringx,ringy,zeros(size(ringx)),[1 1 1],'edgecolor','none');
        ax2plot.Children = [ax2plot.Children(3:end);ax2plot.Children(1);ax2plot.Children(2)];
        if CONTOURNUM>0
            contour(ax2plot,Xi,Yi,Zi,CONTOURNUM,'LineColor','k');
            ax2plot.Children = [ax2plot.Children(2:end-1);ax2plot.Children(1);ax2plot.Children(end)];
        end
    elseif size(Values,2)>1
        % plot connectivity map
        ChanPair = Values(:,1:2);
        selfIdx = (ChanPair(:,1)==ChanPair(:,2));
        ChanPair(selfIdx,:) = [];
        [vlist,tlist] = ismember(ChanPair,plotchans);
        vlist = and(vlist(:,1),vlist(:,2));
        ChanPair = tlist(vlist,:);
        Values = Values(vlist,:);
        if max(ChanPair(:))>length(chanlocs)|| min(ChanPair(:))<1
            error('Out of channel numbers')
        end
        if size(Values,2)==3
            lineC = Values(:,3);
            if length(lineC) == 1
                color2plot = repmat(lineColor,size(Values,1),1);
            else
                %             lineC = (lineC-min(lineC));
                cm = colorMap;
                lineC = lineC*abs(diff(CLim))/(max(lineC)-min(lineC));
                lineC = lineC+min(CLim)-min(lineC);
                lineC = round((lineC-min(lineC))/abs(diff(CLim))*(size(cm,1)-1))+1;
                color2plot = cm(lineC,:);
            end
        else
            color2plot = repmat(lineColor,size(Values,1),1);
        end
        if isDir
            Ddist = sqrt((xx(ChanPair(:,1))-xx(ChanPair(:,2))).^2+(yy(ChanPair(:,1))-yy(ChanPair(:,2))).^2);
            arrowLen = max(min(min(Ddist)/3,15),10);
            for ki = 1:size(ChanPair,1)
                set(gcf,'CurrentAxes',ax2plot);
                arrow([xx(ChanPair(ki,1)),yy(ChanPair(ki,1))],[xx(ChanPair(ki,2)),yy(ChanPair(ki,2))],...
                    [electrodeSize;electrodeSize],'width',LineWidth,'length',arrowLen,...
                    'TipAngle',15,'FaceColor',color2plot(ki,:),'EdgeColor','none');
            end
        else
            for ki = 1:size(ChanPair,1)
                line(ax2plot,xx(ChanPair(ki,:)),yy(ChanPair(ki,:)),'LineWidth',LineWidth,...
                    'Color',color2plot(ki,:));
            end
        end
    end
end
% pbaspect(ax2plot,[1 1 1]);
% daspect(ax2plot,[1 1 1]);

end



function matrixOut = smooth2a(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc)
%
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element "i" by the mean of the rectange centered on "i".  Any NaN
% elements are ignored in the averaging.  If element "i" is a NaN, then it
% will be preserved as NaN in the output.  At the edges of the matrix,
% where you cannot build a full rectangle, as much of the rectangle that
% fits on your matrix is used (similar to the default on Matlab's builtin
% function "smooth").
%
% "matrixIn": original matrix
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
%
% "matrixOut": smoothed version of original matrix
%
%
% 	Written by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
%
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
%
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se

%
% Initial error statements and definitions
%
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr;
if nargin < 3, N(2) = N(1); else N(2) = Nc; end

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the
% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by-
% (2*Nc+1) rectangle centered on element "i".
%
[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

%
% Setting all "NaN" elements of "matrixIn" to zero so that these will not
% affect the summation.  (If this isn't done, any sum that includes a NaN
% will also become NaN.)
%
A = isnan(matrixIn);
matrixIn(A) = 0;

%
% For each element, we have to count how many non-NaN elements went into
% the sums.  This is so we can divide by that number to get a mean.  We use
% the same matrices to do this (ie, "eL" and "eR").
%
nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;

%
% Actually taking the mean.
%
matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;
end

function [upY,lowY] = outEllip(x,y,Xi)
if y>0
    kup = sqrt((y)^2/446^2+(x)^2/344^2);
    kdown = kup;
end
if y<=0
    kdown = sqrt((y)^2/350^2+(x)^2/344^2);
    kup = kdown;
end

upY = (kup^2*446^2-446^2/344^2.*(Xi.^2));



lowY = (kdown^2*350^2-350^2/344^2.*(Xi.^2));
upY(upY<0) = 0;
lowY(lowY<0) = 0;
upY = sqrt(upY)*1.05;
lowY = -sqrt(lowY)*1.1;
end

function [upY,lowY] = inEllip(x0,y0,x,y,Xi)

r1 = 2*(y-y0)^2;
r2 = 2*(x-x0)^2;
yt = ((r1-r1/r2*(Xi-x0).^2));
yt(yt<0) = 0;
yt = sqrt(yt);
upY = (yt*1.1+y0);
lowY = (-yt*1.1+y0);

end

function cbars = CBar
cbars = [0.8 0 0;0.8 0.4 0;0.8 0.8 0;0.4 0.8 0.4;0 0.8 0.8;0 0.4 0.8;0 0 0.8];
% cbars = [68 104 188;116 139 209;174 188 225;220 218 229;251 235 220;254 221 150;232 160 108]/255;
% cbars = [251 245 247;244 226 240;217 204 222;237 184 204;241 184 201]/255;
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
if size(Cbar,1)>cnum
    Cbar(round(end/2),:)=[];
end
end

function [x,y,z] = sdownsample(x,y,z,k)
x = x/k; y = y/k;
t1 = size(x,1); t2 = size(x,2);
if k>=min(t1,t2)
   error('too large k value') 
end
x = x(:,round(1:(t2-1)/(t2/k-1):t2));
x = x(round(1:(t1-1)/(t1/k-1):t1),:);
y = y(:,round(1:(t2-1)/(t2/k-1):t2));
y = y(round(1:(t1-1)/(t1/k-1):t1),:);
z = z(:,round(1:(t2-1)/(t2/k-1):t2));
z = z(round(1:(t1-1)/(t1/k-1):t1),:);
end

function [ul,dl,fxi] = vdownsample(ul,dl,fxi,k)
ul = ul/k;
dl = dl/k;
fxi = fxi/k;
try
ul = ul(round(1:(length(ul)-1)/(length(ul)/k-1):length(ul)));
dl = dl(round(1:(length(dl)-1)/(length(dl)/k-1):length(dl)));
fxi = fxi(round(1:(length(fxi)-1)/(length(fxi)/k-1):length(fxi)));
catch
    error('too large k value');
end 
=======
function handles = Newtopoplot(Values,chanlocs,eegtopoSet,varargin)
%Newtopoplot plots  EEG channel locations topography/spatial weights topography/ spatial connectivity topography
%   handles =  Newtopoplot(Values,chanlocs,eegtopoSet,varargin);
%       Another two functions arrow & predint needed.
%   Input:
%      -Values(needed): []/ nchan X 1 vector/ nchanT X 2(3) matrix, if Values is
%               empty([]), then Newtopoplot plots EEG channel locations; if Values
%               is nchan by 1 vector(where nchan is the number of channels,length(chanlocs)),
%               then Newtopoplot plots spatial weights topography;
%               if Values is nchanT by 2 or 3 matrix, then Newtopoplot
%               plots connectivity topography, note that nchanT is the number of
%               pairs of electrodes to plot connectivity lines, and the first two
%               columns are the serial number of channels(e.g Values = [2,3]
%               means that channel 2 and 3 have connectivity); the third column(if
%               provided) is the strength of connectivity(e.g Values = [2,3,0.5])
%      -chanlocs(needed): Channel locations structure. Using EEGlab EEG.chanlocs or
%               location structure S containing:
%                   S.labels,S.theta,S.radiu,S.X,S.Y,S.Z
%      -eegtopoSet(needed but the programme load it by default):
%               Settings of topography, which has been provided as 'EEGtopoSet'.
%               Newtopoplot will load it if exists, or you need to provide it
%   Options:
%      -CLim(default [0,1]):color limits of topography, 2 by 1 vector
%      -electrodes(default 'on'): 'off': don't show electrodes; 'on': show
%               electordes; 'numbers': show electrodes and channel serial
%               number; 'text': show electrodes and channel labels
%      -plotchans(default is 1:length(chanlocs)): channel numbers to plot,vector
%      -plotrad(default is 0.55):radius of cartoon head, no need to change in general,float
%      -isInside(default is false): true:plot channels inside head(max(chanlocs.radius)<plotrad)
%                false: plot all channels
%      -headColor(default is [0 0 0]): color of cartoon head and ear; 3 by
%                1 vector.
%      -headLineWidth(default is 2.5): line width of cartoon head; float
%      -electrodeSize(default is 20): size of electrode. float
%      -electrodColor(default is [0 0 0]): color of electrode
%      -textColor(default is [0 0 0]): color of electrode numbers/labels.(see Options:electrodes)
%      -textSize(default is 10): size of electrode numbers/labels
%      -contournum(default is 0): number contour lines in spatial weight topography, int (not recommended)
%      -shading(default is 'interp'): 'interp'/'flat',
%      -colormap(default is jet):colormap, if you plot connectivity with
%               values(nchanT X 3), better specific it so that color of connectivity line is
%               consistent with colormap
%      -isDir(default is false): used when plot connectivity topography.
%               true: plot directed connectivity topography using arrow. The
%                   direction is from the first column to second column
%                   (e.g [2,3], the arrow starts from channel 2 point at channel 3, see Input:Values)
%               false: plot undirected connectivity topography
%      -LineWidth(defalut is 2): connectivity line width, float
%      -LineColor(defalut is [0 0 0]): connectivity line color, used when no connectivity strength
%               specified(nchanT X 2, see Input:Values); if strength has
%               been specified(nchanT X 3), then line color is calculated
%               using strength according to colormap(see Options:colormap)
%      -Axes: default is gca; plot topography on specific axis if provided.
%   Examples:
%       plot channel inside head locations with labels:
%           figure;Newtopoplot([],chanlocs,[],'electrodes','labels','isInside',true);
%       plot spatial weight topography
%           nchan = length(chanlocs);Values = rand(nchan,1);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','labels','isInside',false);
%       plot spatial weight topography in region of interest(ROI)
%           nchan = length(chanlocs); ROI = 1:10;% channel 1-10 as ROI
%           Values = nan(nchan,1); Values(ROI) = rand(ROI,1);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','labels','isInside',false);
%           Note: Values = nan(nchan,1), then Newtopoplot only plots values inside ROI;
%               if you want to plot values in the whole head, then specific
%               values outside ROI first, e.g Values = zeros(nchan,1);
%       plot undirected connectivity topography without strength values
%           nchan = length(chanlocs); Values = randi(nchan,12,2);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','on','LineWidth',3,'LineColor',[1 0 0]);
%       plot directed connectivity topography with strength values
%           nchan = length(chanlocs); Values = randi(nchan,12,3);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','on','LineWidth',3,'colormap','jet');

% load coordinates
if nargin<3||isempty(eegtopoSet)
    try
        path2file = which('Newtopoplot');
        path2set = fileparts(path2file);
        path2set = fullfile(path2set,'EEGtopoSet');
        load(path2set,'eegtopoSet');
    catch
        error('No topograph setting found');
    end
end


% default setting
CONTOURNUM = 0;
ELECTRODES = 'on';
isInside = false;
isDir = false;
plotrad = 0.55;
SHADING = 'interp';
headColor = [0 0 0];
electrodeColor = [0 0 0];
textColor = [0 0 0];
lineColor = [0 0 0];
LineWidth = 2;
CLim = [0 1];
colorMap = CBar;
electrodeSize = 20;
headLineWidth = 1.5;
plotchans = 1:length(chanlocs);
textSize = 10;
ax2plot = gca;

nargs = nargin;
if nargs > 3
    for iter = 1:2:length(varargin)
        Param = lower(varargin{iter});
        Value = varargin{iter+1};
        switch Param
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
            case 'plotrad'
                plotrad = Value;
            case 'shading'
                SHADING = Value;
                if ~any(strcmpi(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
            case 'headcolor'
                headColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'electrodecolor'
                electrodeColor = Value;
                if size(Value,2)~=3
                    error('Electrode color must be a 1 x 3 matrix')
                end
            case 'textcolor'
                textColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'linewidth'
                LineWidth = Value;
            case 'linecolor'
                lineColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'colormap'
                colorMap = Value;
                if size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
            case 'electrodesize'
                electrodeSize = Value;
            case 'headlinewidth'
                headLineWidth = Value;
            case 'clim'
                CLim = Value;
                if numel(CLim)~=2
                    error('Color Limit must be 1 x 2 vector');
                end
            case 'isinside'
                isInside = Value;
                if ~islogical(isInside)
                    error('isInside must be true or false');
                end
            case 'isdirection'
                isDir = Value;
                if ~islogical(isDir)
                    error('isInside must be true or false');
                end
            case 'plotchans'
                plotchans = Value;
            case 'textsize'
                textSize = Value;
            case 'axes'
                ax2plot = Value;
        end
    end
    
    
end
handles.axes = ax2plot;
colormap(colorMap);
%% Read channel location
Th=[chanlocs.theta];
Rd=[chanlocs.radius];
Rd = Rd*plotrad/0.55;
Th = pi/180*Th;                              % convert degrees to radians
allchansind = 1:length(Th);

%% remove infinite and NaN values
inds = [];
if ~isempty(Values)
    if size(Values,1)==1||size(Values,2)==1
        Values = Values(:);
    end
    inds = union(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
end
for chani=1:length(chanlocs)
    if isempty(chanlocs(chani).X); inds = [inds chani]; end
end
if isInside
    outs = find(Rd>=0.5495);
    plotchans = setdiff(plotchans,outs);
end
plotchans   = setdiff(plotchans,inds);

[x,y]       = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans   = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
x           = x(plotchans);
y           = y(plotchans);
Rd          = Rd(plotchans);
% project the coordinates of electrodes

yy = x*62.6; xx = y*62.6; yyt = xx;  yt = sqrt(344.^2-xx.^2).*sign(yy);
yid = yy>0; yid2 = yy<=0;
yyt(yid) = mean(predint(eegtopoSet.hhead,xx(yid)+409),2)-355.5;
yyt(yid2) = mean(predint(eegtopoSet.lhead,xx(yid2)+409),2)-355.5;
yy = yy.*(yyt./yt);
yy(isnan(yy)) = 0;
xx = xx*10;yy = yy*10;


k = 1;
%% Plot Head
% cla(ax2plot)
hold(ax2plot,'on');
plot(ax2plot,eegtopoSet.Chead(:,1)/k,eegtopoSet.Chead(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cleft(:,1)/k,eegtopoSet.Cleft(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cright(:,1)/k,eegtopoSet.Cright(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cnose(:,1)/k,eegtopoSet.Cnose(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);

if strcmp(ELECTRODES,'off')
    plot(ax2plot,xx/k,yy/k,'.','MarkerSize',0.001,'Color','none');
else
    plot(ax2plot,xx/k,yy/k,'.','MarkerSize',electrodeSize,'Color',electrodeColor);
    if strcmp(ELECTRODES,'labels')
        text(ax2plot,xx/k,yy/k+1.2*electrodeSize,{chanlocs(plotchans).labels},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',textColor,'hittest','off','FontName','Arial','FontSize',textSize);
    end
    if strcmp(ELECTRODES,'numbers')
        text(ax2plot,xx/k,yy/k+1.2*electrodeSize,num2str(allchansind(:)),'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',textColor,'hittest','off','FontName','Arial','FontSize',textSize);
        
    end
end
Rratio = max(1,max(Rd)/0.55);
set(ax2plot,'Xtick',[],'Ytick',[],'DataAspectRatio',[0.9,1,1],'xlim',[-450/k,450/k]*Rratio,...
    'Ylim',[-420/k,520/k]*Rratio,'clim',CLim,'Xcolor',[0.9 0.95 1],'Ycolor',[0.9 0.95 1],'Zcolor',[0.9 0.95 1],'visible','off');
% axis off
set(ax2plot.Parent,'Color',[1 1 1]);
if isempty(Values)
    title([num2str(length(plotchans)),' of ',num2str(length(Th)),' channels']);
end
if ~isempty(Values)
    
    if size(Values,2)==1
        % Plot topo weight/energy values
        Values = Values(plotchans,:);
        if ~isInside && max(Rd)>=0.5495
            % if some channel to plot is outside the head
            [mRd,mid] = max(Rd);
            ratio = mRd/0.5495;
            
            xi = round(-344*ratio*1.1):5:round(344*ratio*1.1);
            yi = round(-349*ratio*1.1)-5:5:round(446*ratio*1.1);
            
            [Xi,Yi,Zi] = griddata(xx,yy,Values,xi,yi','v4');
%             Zi = smooth2a(Zi,5,5);
            fXi = Xi(1,:);
            [Y_upLimit,Y_lowLimit] = outEllip(abs(xx(mid)),(yy(mid)),fXi);
            mask = Yi>repmat(Y_upLimit,size(Xi,1),1)|Yi<repmat(Y_lowLimit,size(Xi,1),1);
            Zi(mask==1) = NaN;
            [Xi,Yi,Zi] = sdownsample(Xi,Yi,Zi,k);
            [Y_upLimit,Y_lowLimit,fXi] = vdownsample(Y_upLimit,Y_lowLimit,fXi,k);
            ringx = [[fXi,fXi(end:-1:1),fXi(1)]*1.05,[fXi,fXi(end:-1:1),fXi(1)]*0.97];
            ringy = [[Y_upLimit*1.1,Y_lowLimit(end:-1:1)*1.1,Y_upLimit(1)*1.1],...
                [Y_upLimit*0.98,Y_lowLimit(end:-1:1)*0.96,Y_upLimit(1)*0.98]];
           
        else
            xi = [max(min(xx),-344),min(max(xx),344)];
            yi = [min(yy),max(yy)];
            x0 = mean(xi); y0 = mean(yi);
            xit = linspace(((xi(1)-x0)*sqrt(2)*1.3+x0),((xi(end)-x0)*sqrt(2)*1.3+x0),...
                max(150,round((xi(2)-xi(1))/1)));
            
            yit = linspace(((yi(1)-y0)*sqrt(2)*1.3+y0),((yi(end)-y0)*sqrt(2)*1.3+y0),...
                max(150,round((yi(2)-yi(1))/1)));
            [Xi,Yi,Zi] = griddata(xx,yy,Values,xit,yit','v4');
%             Zi = smooth2a(Zi,3,3);
            fXi = Xi(1,:);
            outId = abs(fXi)>344;
            [tup,tlow] = inEllip(x0,y0,1.1*(max(abs(xi-x0)))+x0,1.1*(max(abs(yi-y0)))+y0,fXi);
            Y_upLimit = mean(predint(eegtopoSet.hhead,fXi+409),2)-355.5;
            Y_lowLimit = mean(predint(eegtopoSet.lhead,fXi+409),2)-355.5;
            Y_lowLimit(outId) = 0;
            Y_upLimit(outId) = 0;
            Y_upLimit = min(Y_upLimit(:),tup(:));
            Y_lowLimit = max(Y_lowLimit(:),tlow(:));
            Y_upLimit = Y_upLimit';
            Y_lowLimit = Y_lowLimit';
            
            tup = Y_upLimit(~outId);
            tlow = Y_lowLimit(~outId);
            fXi = fXi(~outId);
            fXi = fXi-x0;
            fXi(fXi<0) = fXi(fXi<0)+0;
            [Xi,Yi,Zi] = sdownsample(Xi,Yi,Zi,k);
            [tup,tlow,fXi] = vdownsample(tup,tlow,fXi,k);
            ringx = [[Xi(1,:),Xi(1,end:-1:1),Xi(1,1)],[fXi,fXi(end:-1:1),fXi(1)]+x0];
            ringy = [[Yi(1,:),Yi(end,:),Yi(1,1)],...
                [y0+(tup-y0),y0+(tlow(end:-1:1)-y0),y0+(tup(1)-y0)]];
            
            
        end
        
        handles.hsuf = surface(ax2plot,Xi,Yi,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING,...
            'FaceLighting','gouraud');
        
        patch(ax2plot,ringx,ringy,zeros(size(ringx)),[1 1 1],'edgecolor','none');
        ax2plot.Children = [ax2plot.Children(3:end);ax2plot.Children(1);ax2plot.Children(2)];
        if CONTOURNUM>0
            contour(ax2plot,Xi,Yi,Zi,CONTOURNUM,'LineColor','k');
            ax2plot.Children = [ax2plot.Children(2:end-1);ax2plot.Children(1);ax2plot.Children(end)];
        end
    elseif size(Values,2)>1
        % plot connectivity map
        ChanPair = Values(:,1:2);
        selfIdx = (ChanPair(:,1)==ChanPair(:,2));
        ChanPair(selfIdx,:) = [];
        [vlist,tlist] = ismember(ChanPair,plotchans);
        vlist = and(vlist(:,1),vlist(:,2));
        ChanPair = tlist(vlist,:);
        Values = Values(vlist,:);
        if max(ChanPair(:))>length(chanlocs)|| min(ChanPair(:))<1
            error('Out of channel numbers')
        end
        if size(Values,2)==3
            lineC = Values(:,3);
            if length(lineC) == 1
                color2plot = repmat(lineColor,size(Values,1),1);
            else
                %             lineC = (lineC-min(lineC));
                cm = colorMap;
                lineC = lineC*abs(diff(CLim))/(max(lineC)-min(lineC));
                lineC = lineC+min(CLim)-min(lineC);
                lineC = round((lineC-min(lineC))/abs(diff(CLim))*(size(cm,1)-1))+1;
                color2plot = cm(lineC,:);
            end
        else
            color2plot = repmat(lineColor,size(Values,1),1);
        end
        if isDir
            Ddist = sqrt((xx(ChanPair(:,1))-xx(ChanPair(:,2))).^2+(yy(ChanPair(:,1))-yy(ChanPair(:,2))).^2);
            arrowLen = max(min(min(Ddist)/3,15),10);
            for ki = 1:size(ChanPair,1)
                set(gcf,'CurrentAxes',ax2plot);
                arrow([xx(ChanPair(ki,1)),yy(ChanPair(ki,1))],[xx(ChanPair(ki,2)),yy(ChanPair(ki,2))],...
                    [electrodeSize;electrodeSize],'width',LineWidth,'length',arrowLen,...
                    'TipAngle',15,'FaceColor',color2plot(ki,:),'EdgeColor','none');
            end
        else
            for ki = 1:size(ChanPair,1)
                line(ax2plot,xx(ChanPair(ki,:)),yy(ChanPair(ki,:)),'LineWidth',LineWidth,...
                    'Color',color2plot(ki,:));
            end
        end
    end
end
% pbaspect(ax2plot,[1 1 1]);
% daspect(ax2plot,[1 1 1]);

end



function matrixOut = smooth2a(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc)
%
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element "i" by the mean of the rectange centered on "i".  Any NaN
% elements are ignored in the averaging.  If element "i" is a NaN, then it
% will be preserved as NaN in the output.  At the edges of the matrix,
% where you cannot build a full rectangle, as much of the rectangle that
% fits on your matrix is used (similar to the default on Matlab's builtin
% function "smooth").
%
% "matrixIn": original matrix
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
%
% "matrixOut": smoothed version of original matrix
%
%
% 	Written by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
%
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
%
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se

%
% Initial error statements and definitions
%
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr;
if nargin < 3, N(2) = N(1); else N(2) = Nc; end

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the
% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by-
% (2*Nc+1) rectangle centered on element "i".
%
[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

%
% Setting all "NaN" elements of "matrixIn" to zero so that these will not
% affect the summation.  (If this isn't done, any sum that includes a NaN
% will also become NaN.)
%
A = isnan(matrixIn);
matrixIn(A) = 0;

%
% For each element, we have to count how many non-NaN elements went into
% the sums.  This is so we can divide by that number to get a mean.  We use
% the same matrices to do this (ie, "eL" and "eR").
%
nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;

%
% Actually taking the mean.
%
matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;
end

function [upY,lowY] = outEllip(x,y,Xi)
if y>0
    kup = sqrt((y)^2/446^2+(x)^2/344^2);
    kdown = kup;
end
if y<=0
    kdown = sqrt((y)^2/350^2+(x)^2/344^2);
    kup = kdown;
end

upY = (kup^2*446^2-446^2/344^2.*(Xi.^2));



lowY = (kdown^2*350^2-350^2/344^2.*(Xi.^2));
upY(upY<0) = 0;
lowY(lowY<0) = 0;
upY = sqrt(upY)*1.05;
lowY = -sqrt(lowY)*1.1;
end

function [upY,lowY] = inEllip(x0,y0,x,y,Xi)

r1 = 2*(y-y0)^2;
r2 = 2*(x-x0)^2;
yt = ((r1-r1/r2*(Xi-x0).^2));
yt(yt<0) = 0;
yt = sqrt(yt);
upY = (yt*1.1+y0);
lowY = (-yt*1.1+y0);

end

function cbars = CBar
cbars = [0.8 0 0;0.8 0.4 0;0.8 0.8 0;0.4 0.8 0.4;0 0.8 0.8;0 0.4 0.8;0 0 0.8];
% cbars = [68 104 188;116 139 209;174 188 225;220 218 229;251 235 220;254 221 150;232 160 108]/255;
% cbars = [251 245 247;244 226 240;217 204 222;237 184 204;241 184 201]/255;
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
if size(Cbar,1)>cnum
    Cbar(round(end/2),:)=[];
end
end

function [x,y,z] = sdownsample(x,y,z,k)
x = x/k; y = y/k;
t1 = size(x,1); t2 = size(x,2);
if k>=min(t1,t2)
   error('too large k value') 
end
x = x(:,round(1:(t2-1)/(t2/k-1):t2));
x = x(round(1:(t1-1)/(t1/k-1):t1),:);
y = y(:,round(1:(t2-1)/(t2/k-1):t2));
y = y(round(1:(t1-1)/(t1/k-1):t1),:);
z = z(:,round(1:(t2-1)/(t2/k-1):t2));
z = z(round(1:(t1-1)/(t1/k-1):t1),:);
end

function [ul,dl,fxi] = vdownsample(ul,dl,fxi,k)
ul = ul/k;
dl = dl/k;
fxi = fxi/k;
try
ul = ul(round(1:(length(ul)-1)/(length(ul)/k-1):length(ul)));
dl = dl(round(1:(length(dl)-1)/(length(dl)/k-1):length(dl)));
fxi = fxi(round(1:(length(fxi)-1)/(length(fxi)/k-1):length(fxi)));
catch
    error('too large k value');
end 
>>>>>>> origin
=======
function handles = Newtopoplot(Values,chanlocs,eegtopoSet,varargin)
%Newtopoplot plots  EEG channel locations topography/spatial weights topography/ spatial connectivity topography
%   handles =  Newtopoplot(Values,chanlocs,eegtopoSet,varargin);
%       Another two functions arrow & predint needed.
%   Input:
%      -Values(needed): []/ nchan X 1 vector/ nchanT X 2(3) matrix, if Values is
%               empty([]), then Newtopoplot plots EEG channel locations; if Values
%               is nchan by 1 vector(where nchan is the number of channels,length(chanlocs)),
%               then Newtopoplot plots spatial weights topography;
%               if Values is nchanT by 2 or 3 matrix, then Newtopoplot
%               plots connectivity topography, note that nchanT is the number of
%               pairs of electrodes to plot connectivity lines, and the first two
%               columns are the serial number of channels(e.g Values = [2,3]
%               means that channel 2 and 3 have connectivity); the third column(if
%               provided) is the strength of connectivity(e.g Values = [2,3,0.5])
%      -chanlocs(needed): Channel locations structure. Using EEGlab EEG.chanlocs or
%               location structure S containing:
%                   S.labels,S.theta,S.radiu,S.X,S.Y,S.Z
%      -eegtopoSet(needed but the programme load it by default):
%               Settings of topography, which has been provided as 'EEGtopoSet'.
%               Newtopoplot will load it if exists, or you need to provide it
%   Options:
%      -CLim(default [0,1]):color limits of topography, 2 by 1 vector
%      -electrodes(default 'on'): 'off': don't show electrodes; 'on': show
%               electordes; 'numbers': show electrodes and channel serial
%               number; 'text': show electrodes and channel labels
%      -plotchans(default is 1:length(chanlocs)): channel numbers to plot,vector
%      -plotrad(default is 0.55):radius of cartoon head, no need to change in general,float
%      -isInside(default is false): true:plot channels inside head(max(chanlocs.radius)<plotrad)
%                false: plot all channels
%      -headColor(default is [0 0 0]): color of cartoon head and ear; 3 by
%                1 vector.
%      -headLineWidth(default is 2.5): line width of cartoon head; float
%      -electrodeSize(default is 20): size of electrode. float
%      -electrodColor(default is [0 0 0]): color of electrode
%      -textColor(default is [0 0 0]): color of electrode numbers/labels.(see Options:electrodes)
%      -textSize(default is 10): size of electrode numbers/labels
%      -contournum(default is 0): number contour lines in spatial weight topography, int (not recommended)
%      -shading(default is 'interp'): 'interp'/'flat',
%      -colormap(default is jet):colormap, if you plot connectivity with
%               values(nchanT X 3), better specific it so that color of connectivity line is
%               consistent with colormap
%      -isDir(default is false): used when plot connectivity topography.
%               true: plot directed connectivity topography using arrow. The
%                   direction is from the first column to second column
%                   (e.g [2,3], the arrow starts from channel 2 point at channel 3, see Input:Values)
%               false: plot undirected connectivity topography
%      -LineWidth(defalut is 2): connectivity line width, float
%      -LineColor(defalut is [0 0 0]): connectivity line color, used when no connectivity strength
%               specified(nchanT X 2, see Input:Values); if strength has
%               been specified(nchanT X 3), then line color is calculated
%               using strength according to colormap(see Options:colormap)
%      -Axes: default is gca; plot topography on specific axis if provided.
%   Examples:
%       plot channel inside head locations with labels:
%           figure;Newtopoplot([],chanlocs,[],'electrodes','labels','isInside',true);
%       plot spatial weight topography
%           nchan = length(chanlocs);Values = rand(nchan,1);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','labels','isInside',false);
%       plot spatial weight topography in region of interest(ROI)
%           nchan = length(chanlocs); ROI = 1:10;% channel 1-10 as ROI
%           Values = nan(nchan,1); Values(ROI) = rand(ROI,1);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','labels','isInside',false);
%           Note: Values = nan(nchan,1), then Newtopoplot only plots values inside ROI;
%               if you want to plot values in the whole head, then specific
%               values outside ROI first, e.g Values = zeros(nchan,1);
%       plot undirected connectivity topography without strength values
%           nchan = length(chanlocs); Values = randi(nchan,12,2);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','on','LineWidth',3,'LineColor',[1 0 0]);
%       plot directed connectivity topography with strength values
%           nchan = length(chanlocs); Values = randi(nchan,12,3);
%           figure;Newtopoplot(Values,chanlocs,[],'electrodes','on','LineWidth',3,'colormap','jet');

% load coordinates
if nargin<3||isempty(eegtopoSet)
    try
        path2file = which('Newtopoplot');
        path2set = fileparts(path2file);
        path2set = fullfile(path2set,'EEGtopoSet');
        load(path2set,'eegtopoSet');
    catch
        error('No topograph setting found');
    end
end


% default setting
CONTOURNUM = 0;
ELECTRODES = 'on';
isInside = false;
isDir = false;
plotrad = 0.55;
SHADING = 'interp';
headColor = [0 0 0];
electrodeColor = [0 0 0];
textColor = [0 0 0];
lineColor = [0 0 0];
LineWidth = 2;
CLim = [0 1];
colorMap = CBar;
electrodeSize = 20;
headLineWidth = 1.5;
plotchans = 1:length(chanlocs);
textSize = 10;
ax2plot = gca;

nargs = nargin;
if nargs > 3
    for iter = 1:2:length(varargin)
        Param = lower(varargin{iter});
        Value = varargin{iter+1};
        switch Param
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
            case 'plotrad'
                plotrad = Value;
            case 'shading'
                SHADING = Value;
                if ~any(strcmpi(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
            case 'headcolor'
                headColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'electrodecolor'
                electrodeColor = Value;
                if size(Value,2)~=3
                    error('Electrode color must be a 1 x 3 matrix')
                end
            case 'textcolor'
                textColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'linewidth'
                LineWidth = Value;
            case 'linecolor'
                lineColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'colormap'
                colorMap = Value;
                if size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
            case 'electrodesize'
                electrodeSize = Value;
            case 'headlinewidth'
                headLineWidth = Value;
            case 'clim'
                CLim = Value;
                if numel(CLim)~=2
                    error('Color Limit must be 1 x 2 vector');
                end
            case 'isinside'
                isInside = Value;
                if ~islogical(isInside)
                    error('isInside must be true or false');
                end
            case 'isdirection'
                isDir = Value;
                if ~islogical(isDir)
                    error('isInside must be true or false');
                end
            case 'plotchans'
                plotchans = Value;
            case 'textsize'
                textSize = Value;
            case 'axes'
                ax2plot = Value;
        end
    end
    
    
end
handles.axes = ax2plot;
colormap(colorMap);
%% Read channel location
Th=[chanlocs.theta];
Rd=[chanlocs.radius];
Rd = Rd*plotrad/0.55;
Th = pi/180*Th;                              % convert degrees to radians
allchansind = 1:length(Th);

%% remove infinite and NaN values
inds = [];
if ~isempty(Values)
    if size(Values,1)==1||size(Values,2)==1
        Values = Values(:);
    end
    inds = union(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
end
for chani=1:length(chanlocs)
    if isempty(chanlocs(chani).X); inds = [inds chani]; end
end
if isInside
    outs = find(Rd>=0.5495);
    plotchans = setdiff(plotchans,outs);
end
plotchans   = setdiff(plotchans,inds);

[x,y]       = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans   = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
x           = x(plotchans);
y           = y(plotchans);
Rd          = Rd(plotchans);
% project the coordinates of electrodes

yy = x*62.6; xx = y*62.6; yyt = xx;  yt = sqrt(344.^2-xx.^2).*sign(yy);
yid = yy>0; yid2 = yy<=0;
yyt(yid) = mean(predint(eegtopoSet.hhead,xx(yid)+409),2)-355.5;
yyt(yid2) = mean(predint(eegtopoSet.lhead,xx(yid2)+409),2)-355.5;
yy = yy.*(yyt./yt);
yy(isnan(yy)) = 0;
xx = xx*10;yy = yy*10;


k = 1;
%% Plot Head
% cla(ax2plot)
hold(ax2plot,'on');
plot(ax2plot,eegtopoSet.Chead(:,1)/k,eegtopoSet.Chead(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cleft(:,1)/k,eegtopoSet.Cleft(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cright(:,1)/k,eegtopoSet.Cright(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);
plot(ax2plot,eegtopoSet.Cnose(:,1)/k,eegtopoSet.Cnose(:,2)/k,'Color',headColor,'LineWidth',headLineWidth);

if strcmp(ELECTRODES,'off')
    plot(ax2plot,xx/k,yy/k,'.','MarkerSize',0.001,'Color','none');
else
    plot(ax2plot,xx/k,yy/k,'.','MarkerSize',electrodeSize,'Color',electrodeColor);
    if strcmp(ELECTRODES,'labels')
        text(ax2plot,xx/k,yy/k+1.2*electrodeSize,{chanlocs(plotchans).labels},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',textColor,'hittest','off','FontName','Arial','FontSize',textSize);
    end
    if strcmp(ELECTRODES,'numbers')
        text(ax2plot,xx/k,yy/k+1.2*electrodeSize,num2str(allchansind(:)),'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',textColor,'hittest','off','FontName','Arial','FontSize',textSize);
        
    end
end
Rratio = max(1,max(Rd)/0.55);
set(ax2plot,'Xtick',[],'Ytick',[],'DataAspectRatio',[0.9,1,1],'xlim',[-450/k,450/k]*Rratio,...
    'Ylim',[-420/k,520/k]*Rratio,'clim',CLim,'Xcolor',[0.9 0.95 1],'Ycolor',[0.9 0.95 1],'Zcolor',[0.9 0.95 1],'visible','off');
% axis off
set(ax2plot.Parent,'Color',[1 1 1]);
if isempty(Values)
    title([num2str(length(plotchans)),' of ',num2str(length(Th)),' channels']);
end
if ~isempty(Values)
    
    if size(Values,2)==1
        % Plot topo weight/energy values
        Values = Values(plotchans,:);
        if ~isInside && max(Rd)>=0.5495
            % if some channel to plot is outside the head
            [mRd,mid] = max(Rd);
            ratio = mRd/0.5495;
            
            xi = round(-344*ratio*1.1):5:round(344*ratio*1.1);
            yi = round(-349*ratio*1.1)-5:5:round(446*ratio*1.1);
            
            [Xi,Yi,Zi] = griddata(xx,yy,Values,xi,yi','v4');
%             Zi = smooth2a(Zi,5,5);
            fXi = Xi(1,:);
            [Y_upLimit,Y_lowLimit] = outEllip(abs(xx(mid)),(yy(mid)),fXi);
            mask = Yi>repmat(Y_upLimit,size(Xi,1),1)|Yi<repmat(Y_lowLimit,size(Xi,1),1);
            Zi(mask==1) = NaN;
            [Xi,Yi,Zi] = sdownsample(Xi,Yi,Zi,k);
            [Y_upLimit,Y_lowLimit,fXi] = vdownsample(Y_upLimit,Y_lowLimit,fXi,k);
            ringx = [[fXi,fXi(end:-1:1),fXi(1)]*1.05,[fXi,fXi(end:-1:1),fXi(1)]*0.97];
            ringy = [[Y_upLimit*1.1,Y_lowLimit(end:-1:1)*1.1,Y_upLimit(1)*1.1],...
                [Y_upLimit*0.98,Y_lowLimit(end:-1:1)*0.96,Y_upLimit(1)*0.98]];
           
        else
            xi = [max(min(xx),-344),min(max(xx),344)];
            yi = [min(yy),max(yy)];
            x0 = mean(xi); y0 = mean(yi);
            xit = linspace(((xi(1)-x0)*sqrt(2)*1.3+x0),((xi(end)-x0)*sqrt(2)*1.3+x0),...
                max(150,round((xi(2)-xi(1))/1)));
            
            yit = linspace(((yi(1)-y0)*sqrt(2)*1.3+y0),((yi(end)-y0)*sqrt(2)*1.3+y0),...
                max(150,round((yi(2)-yi(1))/1)));
            [Xi,Yi,Zi] = griddata(xx,yy,Values,xit,yit','v4');
%             Zi = smooth2a(Zi,3,3);
            fXi = Xi(1,:);
            outId = abs(fXi)>344;
            [tup,tlow] = inEllip(x0,y0,1.1*(max(abs(xi-x0)))+x0,1.1*(max(abs(yi-y0)))+y0,fXi);
            Y_upLimit = mean(predint(eegtopoSet.hhead,fXi+409),2)-355.5;
            Y_lowLimit = mean(predint(eegtopoSet.lhead,fXi+409),2)-355.5;
            Y_lowLimit(outId) = 0;
            Y_upLimit(outId) = 0;
            Y_upLimit = min(Y_upLimit(:),tup(:));
            Y_lowLimit = max(Y_lowLimit(:),tlow(:));
            Y_upLimit = Y_upLimit';
            Y_lowLimit = Y_lowLimit';
            
            tup = Y_upLimit(~outId);
            tlow = Y_lowLimit(~outId);
            fXi = fXi(~outId);
            fXi = fXi-x0;
            fXi(fXi<0) = fXi(fXi<0)+0;
            [Xi,Yi,Zi] = sdownsample(Xi,Yi,Zi,k);
            [tup,tlow,fXi] = vdownsample(tup,tlow,fXi,k);
            ringx = [[Xi(1,:),Xi(1,end:-1:1),Xi(1,1)],[fXi,fXi(end:-1:1),fXi(1)]+x0];
            ringy = [[Yi(1,:),Yi(end,:),Yi(1,1)],...
                [y0+(tup-y0),y0+(tlow(end:-1:1)-y0),y0+(tup(1)-y0)]];
            
            
        end
        
        handles.hsuf = surface(ax2plot,Xi,Yi,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING,...
            'FaceLighting','gouraud');
        
        patch(ax2plot,ringx,ringy,zeros(size(ringx)),[1 1 1],'edgecolor','none');
        ax2plot.Children = [ax2plot.Children(3:end);ax2plot.Children(1);ax2plot.Children(2)];
        if CONTOURNUM>0
            contour(ax2plot,Xi,Yi,Zi,CONTOURNUM,'LineColor','k');
            ax2plot.Children = [ax2plot.Children(2:end-1);ax2plot.Children(1);ax2plot.Children(end)];
        end
    elseif size(Values,2)>1
        % plot connectivity map
        ChanPair = Values(:,1:2);
        selfIdx = (ChanPair(:,1)==ChanPair(:,2));
        ChanPair(selfIdx,:) = [];
        [vlist,tlist] = ismember(ChanPair,plotchans);
        vlist = and(vlist(:,1),vlist(:,2));
        ChanPair = tlist(vlist,:);
        Values = Values(vlist,:);
        if max(ChanPair(:))>length(chanlocs)|| min(ChanPair(:))<1
            error('Out of channel numbers')
        end
        if size(Values,2)==3
            lineC = Values(:,3);
            if length(lineC) == 1
                color2plot = repmat(lineColor,size(Values,1),1);
            else
                %             lineC = (lineC-min(lineC));
                cm = colorMap;
                lineC = lineC*abs(diff(CLim))/(max(lineC)-min(lineC));
                lineC = lineC+min(CLim)-min(lineC);
                lineC = round((lineC-min(lineC))/abs(diff(CLim))*(size(cm,1)-1))+1;
                color2plot = cm(lineC,:);
            end
        else
            color2plot = repmat(lineColor,size(Values,1),1);
        end
        if isDir
            Ddist = sqrt((xx(ChanPair(:,1))-xx(ChanPair(:,2))).^2+(yy(ChanPair(:,1))-yy(ChanPair(:,2))).^2);
            arrowLen = max(min(min(Ddist)/3,15),10);
            for ki = 1:size(ChanPair,1)
                set(gcf,'CurrentAxes',ax2plot);
                arrow([xx(ChanPair(ki,1)),yy(ChanPair(ki,1))],[xx(ChanPair(ki,2)),yy(ChanPair(ki,2))],...
                    [electrodeSize;electrodeSize],'width',LineWidth,'length',arrowLen,...
                    'TipAngle',15,'FaceColor',color2plot(ki,:),'EdgeColor','none');
            end
        else
            for ki = 1:size(ChanPair,1)
                line(ax2plot,xx(ChanPair(ki,:)),yy(ChanPair(ki,:)),'LineWidth',LineWidth,...
                    'Color',color2plot(ki,:));
            end
        end
    end
end
% pbaspect(ax2plot,[1 1 1]);
% daspect(ax2plot,[1 1 1]);

end



function matrixOut = smooth2a(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc)
%
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element "i" by the mean of the rectange centered on "i".  Any NaN
% elements are ignored in the averaging.  If element "i" is a NaN, then it
% will be preserved as NaN in the output.  At the edges of the matrix,
% where you cannot build a full rectangle, as much of the rectangle that
% fits on your matrix is used (similar to the default on Matlab's builtin
% function "smooth").
%
% "matrixIn": original matrix
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
%
% "matrixOut": smoothed version of original matrix
%
%
% 	Written by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
%
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
%
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se

%
% Initial error statements and definitions
%
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr;
if nargin < 3, N(2) = N(1); else N(2) = Nc; end

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the
% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by-
% (2*Nc+1) rectangle centered on element "i".
%
[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

%
% Setting all "NaN" elements of "matrixIn" to zero so that these will not
% affect the summation.  (If this isn't done, any sum that includes a NaN
% will also become NaN.)
%
A = isnan(matrixIn);
matrixIn(A) = 0;

%
% For each element, we have to count how many non-NaN elements went into
% the sums.  This is so we can divide by that number to get a mean.  We use
% the same matrices to do this (ie, "eL" and "eR").
%
nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;

%
% Actually taking the mean.
%
matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;
end

function [upY,lowY] = outEllip(x,y,Xi)
if y>0
    kup = sqrt((y)^2/446^2+(x)^2/344^2);
    kdown = kup;
end
if y<=0
    kdown = sqrt((y)^2/350^2+(x)^2/344^2);
    kup = kdown;
end

upY = (kup^2*446^2-446^2/344^2.*(Xi.^2));



lowY = (kdown^2*350^2-350^2/344^2.*(Xi.^2));
upY(upY<0) = 0;
lowY(lowY<0) = 0;
upY = sqrt(upY)*1.05;
lowY = -sqrt(lowY)*1.1;
end

function [upY,lowY] = inEllip(x0,y0,x,y,Xi)

r1 = 2*(y-y0)^2;
r2 = 2*(x-x0)^2;
yt = ((r1-r1/r2*(Xi-x0).^2));
yt(yt<0) = 0;
yt = sqrt(yt);
upY = (yt*1.1+y0);
lowY = (-yt*1.1+y0);

end

function cbars = CBar
cbars = [0.8 0 0;0.8 0.4 0;0.8 0.8 0;0.4 0.8 0.4;0 0.8 0.8;0 0.4 0.8;0 0 0.8];
% cbars = [68 104 188;116 139 209;174 188 225;220 218 229;251 235 220;254 221 150;232 160 108]/255;
% cbars = [251 245 247;244 226 240;217 204 222;237 184 204;241 184 201]/255;
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
if size(Cbar,1)>cnum
    Cbar(round(end/2),:)=[];
end
end

function [x,y,z] = sdownsample(x,y,z,k)
x = x/k; y = y/k;
t1 = size(x,1); t2 = size(x,2);
if k>=min(t1,t2)
   error('too large k value') 
end
x = x(:,round(1:(t2-1)/(t2/k-1):t2));
x = x(round(1:(t1-1)/(t1/k-1):t1),:);
y = y(:,round(1:(t2-1)/(t2/k-1):t2));
y = y(round(1:(t1-1)/(t1/k-1):t1),:);
z = z(:,round(1:(t2-1)/(t2/k-1):t2));
z = z(round(1:(t1-1)/(t1/k-1):t1),:);
end

function [ul,dl,fxi] = vdownsample(ul,dl,fxi,k)
ul = ul/k;
dl = dl/k;
fxi = fxi/k;
try
ul = ul(round(1:(length(ul)-1)/(length(ul)/k-1):length(ul)));
dl = dl(round(1:(length(dl)-1)/(length(dl)/k-1):length(dl)));
fxi = fxi(round(1:(length(fxi)-1)/(length(fxi)/k-1):length(fxi)));
catch
    error('too large k value');
end 
>>>>>>> 3afa99beaec32453db712d4621fd05217f53fcfd
end