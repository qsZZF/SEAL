% pop_Newtopoplot() - Plot topographs. This function pops up a GUI to call Newtopoplot(); 
% Turn to Newtopoplot() for more information
% Usage
%   >> pop_Newtopoplot(EEG)

function pop_Newtopoplot(EEG,nogui,varargin)

if nargin<1
    help Newtopoplot;
   return 
end

gheadline = [1,1,1,1.2];
gpopmen = [1.5 0.8 1 0.8 1 0.8];

geometry = {gheadline gpopmen [1.5 0.8 1 0.8 1.2 0.6] [1] gheadline gpopmen gpopmen [1] gheadline gpopmen gpopmen [1] [1 1]};
uilist = {...
    { 'Style', 'text', 'string' 'Topo Setting' }...
    { 'Style', 'checkbox', 'string' 'Numbers' 'value' 0 'tag' 'head_Numbers' }...
    { 'Style', 'checkbox', 'string' 'Labels' 'value' 0 'tag' 'head_Labels' }...
    { 'Style', 'checkbox', 'string' 'Show Inside Only'  'tag' 'head_isInside' } ...
    ...
    { 'Style', 'text', 'string' 'Channel To Plot'  }...
    { 'Style', 'edit', 'string' ''   'tag' 'head_chan2plot' }...
    { 'Style', 'text', 'string' 'Line Color' }...
    { 'Style', 'edit', 'string' '0 0 0 '  'tag' 'head_lineColor' }...
    { 'Style', 'text', 'string' 'Line Width' }...
    { 'Style', 'edit', 'string' '1.5'  'tag' 'head_lineWidth' }...
    ...
    { 'Style', 'text', 'string' 'Electrode Color'  }...
    { 'Style', 'edit', 'string' '0 0 0 '   'tag' 'head_eleColor' }...
    { 'Style', 'text', 'string' 'Electrode Size' }...
    { 'Style', 'edit', 'string' '20'  'tag' 'head_eleSize' }...
    { 'Style', 'checkbox', 'string' 'No Electrode '  'tag' 'head_noEle' } {} ...
    ...
    {}...
    ...
    { 'Style', 'text', 'string' 'Topo Weights' }...
    { 'Style', 'checkbox', 'string' 'Subplot' 'value' 1 'tag' 'topo_isSubplot' } ...
    { 'Style', 'checkbox', 'string' 'Show Inside Only' 'value' 1 'tag' 'topo_isInside' } {} ...
    ...
    { 'Style', 'text', 'string' 'EEG Time Range'  }...
    { 'Style', 'edit', 'string' ''   'tag' 'topo_EEGrange' }...
    { 'Style', 'text', 'string' 'ICA Weights' }...
    { 'Style', 'edit', 'string' ''  'tag' 'topo_ICAw' }...
    { 'Style', 'text', 'string' 'Other Data' }...
    { 'Style', 'edit', 'string' ''  'tag' 'topo_data' }...
    ...
    { 'Style', 'text', 'string' 'Color Limit'  }...
    { 'Style', 'edit', 'string' ''   'tag' 'topo_clim' }...
    { 'Style', 'text', 'string' 'Color Bar' }...
    { 'Style', 'edit', 'string' ''  'tag' 'topo_cbar' } {} {}...
    ...
    {}...
    ...
    { 'Style', 'text', 'string' 'Connectivity' }...
    { 'Style', 'checkbox', 'string' 'Subplot' 'value' 1 'tag' 'conn_isSubplot' } ...
    { 'Style', 'checkbox', 'string' 'Show Inside Only' 'value' 1 'tag' 'conn_isInside' } {} ...
    ...
    { 'Style', 'text', 'string' 'Other Data' }...
    { 'Style', 'edit', 'string' ''  'tag' 'conn_data' }...
    { 'Style', 'text', 'string' 'Line Width'  }...
    { 'Style', 'edit', 'string' ''   'tag' 'conn_lineWidth' }...
    { 'Style', 'text', 'string' 'Line Color' }...
    { 'Style', 'edit', 'string' '0 0 0'  'tag' 'conn_lineColor' }...
    ...
    { 'Style', 'text', 'string' 'Color Limit'  }...
    { 'Style', 'edit', 'string' ''   'tag' 'conn_clim' }...
    { 'Style', 'text', 'string' 'Color Bar' }...
    { 'Style', 'edit', 'string' ''  'tag' 'conn_cbar' }...
    { 'Style', 'checkbox', 'string' 'Direction' 'value' 0 'tag' 'conn_isdir' }...
    { 'Style', 'checkbox', 'string' 'Strength' 'value' 0 'tag' 'conn_isStrength' }...
    ...
    {}...
    ...    
    { 'Style', 'checkbox', 'string' 'Plot Topo Weights And Connect In One Fig' 'value' 0 'tag' 'com_plot' }...
    { 'Style', 'checkbox', 'string' 'Plot Gif' 'value' 0 'tag' 'plot_gif' }...    
    };

chanlocs = EEG.chanlocs;
% Default Settings
ELECTRODES = 'on';
isInside = false;
headColor = [0 0 0];
electrodeColor = [0 0 0];
lineColor = [0 0 0];
LineWidth = 2;
electrodeSize = 20;
headLineWidth = 1.5;
plotchans = 1:length(chanlocs);

istopo = 0;
isconn = 0;
com_plot = 0;

if nargin<2
    [ ~, ~, ~, structout ] = inputgui( geometry, uilist, ...
        'pophelp(''pop_Newtopoplot'');', 'Plot topograph');
    if isempty(structout)
        return;
    end
    
    % Check GUI Input
    % Topographic Setting
    if structout.head_noEle
        ELECTRODES = 'off';
    elseif structout.head_Numbers
        ELECTRODES = 'numbers';
    elseif structout.head_Labels
        ELECTRODES = 'labels';
    end
    if structout.head_isInside
        isInside = true;
    end
    if ~isempty(structout.head_lineColor)
        headColor = sscanf(structout.head_lineColor,'%f')';
    end
    if ~isempty(structout.head_chan2plot)
        tplotchans = sscanf(structout.head_chan2plot,'%c');
        eval(['plotchans = plotchans(',tplotchans,');']);
    end
    if ~isempty(structout.head_lineWidth)
        headLineWidth = sscanf(structout.head_lineWidth,'%f');
    end
    if ~isempty(structout.head_eleColor)
        electrodeColor = sscanf(structout.head_eleColor,'%f')';
    end
    if ~isempty(structout.head_eleSize)
        electrodeSize = sscanf(structout.head_eleSize,'%f');
    end
    % Topo Weights Setting
    topo_isInside = boolean(structout.topo_isInside);
    topo_isSubplot = boolean(structout.topo_isSubplot);
    topoValue = [];
    try
        if ~isempty(structout.topo_data)
            topoValue = evalin('base',structout.topo_data);
        elseif ~isempty(structout.topo_ICAw)
            topoRange = sscanf(structout.topo_ICAw,'%c');
            topoValue = eval(['EEG.icawinv(:,['  topoRange ']);']);
        elseif ~isempty(structout.topo_EEGrange)
            topoRange = sscanf(structout.topo_EEGrange,'%c');
            topoValue = eval(['EEG.data(:,[' topoRange ']);']);
        end
    catch e
        disp(e.message);
    end
    if ~isempty(topoValue)
        topoValue = double(topoValue);
        istopo = 1;
    end
    if ~isempty(structout.topo_clim)
        topo_clim = sscanf(structout.topo_clim,'%f')';
    else
        topo_clim = [0 1];
        topoValue = (topoValue-min(abs(topoValue)))./(max(abs(topoValue))-min(abs(topoValue)));
    end
    if ~isempty(structout.topo_cbar)
        topo_cbar = structout.topo_cbar;
    else
        topo_cbar = CBar;
    end
    
    % Connect Setting
    conn_isInside = boolean(structout.conn_isInside);
    conn_isSubplot = boolean(structout.conn_isSubplot);
    conn_isdir = boolean(structout.conn_isdir);
    conn_isStrength = boolean(structout.conn_isStrength);
    connValue = [];
    if ~isempty(structout.conn_data)
        try
            connValue = evalin('base',structout.conn_data);
        catch e
            disp(e.message)
        end
    end
    if ~isempty(connValue)
        connValue = double(connValue);
        isconn = 1;
    end
    if ~isempty(structout.conn_lineWidth)
        LineWidth = sscanf(structout.conn_lineWidth,'%f');
    end
    if ~isempty(structout.conn_lineColor)
        lineColor = sscanf(structout.conn_lineColor,'%f')';
    end
    if ~isempty(structout.com_plot)
        com_plot = boolean(structout.com_plot);
    end
    if ~isempty(structout.conn_clim)
        topo_clim = sscanf(structout.conn_clim,'%f')';
    else
        conn_clim = [0 1];
    end
    if ~isempty(structout.conn_cbar)
        conn_cbar = structout.conn_cbar;
    else
        conn_cbar = CBar;
    end
else
    switch nogui
        case 1
            ELECTRODES = 'labels';
        case 2
            ELECTRODES = 'numbers';
    end
end
% Plot
if com_plot
   if ~(istopo&&isconn)
       warning('lack of topograph/connectivity data');
       com_plot = false;
   end
end
if com_plot&&istopo&&isconn
    
    fnums = min(size(topoValue,2),size(connValue,3));
    if conn_isStrength
        if size(connValue,2)<3
            warning('Lack Of Connectivity Strength') ;
        end
    else
        connValue = connValue(:,1:2,:);
    end
    if topo_isSubplot||conn_isSubplot
        topoAxes = Multiplot(fnums);
        for inum = 1:fnums           
            Newtopoplot(connValue(:,:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',conn_isInside,...
                'clim',conn_clim,'colormap',conn_cbar,'axes',topoAxes{inum},'isdirection',conn_isdir,'linecolor',lineColor,...
                'linewidth',LineWidth);
            hold on
            Newtopoplot(topoValue(:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',topo_isInside,...
                'clim',topo_clim,'colormap',topo_cbar,'axes',topoAxes{inum});
            
        end
    else
        for inum = 1:fnums
            figure;         
            Newtopoplot(connValue(:,:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',conn_isInside,...
                'clim',conn_clim,'colormap',conn_cbar,'isdirection',conn_isdir,'linecolor',lineColor,...
                'linewidth',LineWidth);
            hold on
             Newtopoplot(topoValue(:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',topo_isInside,...
                'clim',topo_clim,'colormap',topo_cbar);
        end
    end
end
if istopo&&~com_plot
    fnums = size(topoValue,2);
    if topo_isSubplot
        topoAxes = Multiplot(fnums);
        for inum = 1:fnums
            Newtopoplot(topoValue(:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',topo_isInside,...
                'clim',topo_clim,'colormap',topo_cbar,'axes',topoAxes{inum});
        end
    else
        for inum = 1:fnums
            figure;
            Newtopoplot(topoValue(:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',topo_isInside,...
                'clim',topo_clim,'colormap',topo_cbar);
        end
    end
end
if isconn&&~com_plot
    if conn_isStrength
        if size(connValue,2)<3
            warning('Lack Of Connectivity Strength') ;
        end
    else
        connValue = connValue(:,1:2,:);
    end
    fnums = size(connValue,3);
    if conn_isSubplot
        connAxes = Multiplot(fnums);
        for inum = 1:fnums
            Newtopoplot(connValue(:,:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',conn_isInside,...
                'clim',conn_clim,'colormap',conn_cbar,'axes',connAxes{inum},'isdirection',conn_isdir,'linecolor',lineColor,...
                'linewidth',LineWidth);
        end
    else
        for inum = 1:fnums
            figure;
            Newtopoplot(connValue(:,:,inum),chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
                'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',conn_isInside,...
                'clim',conn_clim,'colormap',conn_cbar,'isdirection',conn_isdir,'linecolor',lineColor,...
                'linewidth',LineWidth);
        end
    end
end
if ~(istopo||isconn)
    figure;
    Newtopoplot([],chanlocs,[],'electrodes',ELECTRODES,'headcolor',headColor,'electrodecolor',electrodeColor,...
        'electrodesize',electrodeSize,'headlinewidth',headLineWidth,'plotchans',plotchans,'isinside',isInside);
end




function cbars = CBar
cbars = [0.8 0 0.1;0.8 0.5 0;1 1 0;0.5 1 0.5;0 1 1;0 0.5 1;0.1 0.4 0.7];
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


end