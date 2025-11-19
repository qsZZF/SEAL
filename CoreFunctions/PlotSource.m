function handles = PlotSource(SourceData,BrainModel,varargin)
%PlotSource plot cortex domain activity
% handles = PlotSource(SourceData,BrainModel,varargin)
%   Input:
%      -Values(needed): []/ nchan X 1 vector. If Values is
%               empty([]), then PlotSource plots EEG channel locations; if Values
%               is nchan by 1 vector(where nchan is the number of channels,length(chanlocs)), 
%               then topo3d plots spatial weights topography;
%      -BrainModel(needed but the programme load it by default): 
%               Settings of cortex
%   Options:
%      -CLim(default [0,1]):color limits of topography, 2 by 1 vector
%      -roiIdx (default is 1:length(vertices)): region of interst to plot,vector 
%      -cortexcolor(default is [0.75 0.75 0.75]/[0.9 0.7 0.7]): color of cortex; 3 by
%                1 vector.

%      -colormap:colormap, 
%      -Axes: default is gca; plot topography on specific axis if provided.

if nargin<2||isempty(BrainModel)
    try
        path2file = which('PlotSource');
        path2set = fileparts(path2file);
        path2set = fullfile(path2set,'BrainModel');
        load(path2set,'BrainModel');
    catch
        error('No brain model found');
    end
end
if nargin<1
    SourceData = [];
end

if isfield(BrainModel,'vertices')
    vc=BrainModel.vertices;
end
if isfield(BrainModel,'Vertices')
    vc=BrainModel.Vertices;
end
if isfield(BrainModel,'vert')
    vc=BrainModel.vert;
end
if isfield(BrainModel,'Vert')
    vc=BrainModel.Vert;
end
if isfield(BrainModel,'Faces')
    tri=BrainModel.Faces;
end
if isfield(BrainModel,'Face')
    tri=BrainModel.Face;
end
if isfield(BrainModel,'faces')
    tri=BrainModel.faces;
end
if isfield(BrainModel,'face')
    tri=BrainModel.face;
end
% vc = vc-mean(vc);
% default setting
smooth = [];

roiIdx = 1:length(vc);
hemisphere = [];
ax2plot = [];
CLim = [0,1];
figPosition = []; % 默认不指定位置
colorMap = CBar('active');
if ~isempty(SourceData)
    if min(SourceData)<0
        CLim = [-1 1];
        colorMap = CBar('pnactive');
    end
thresh = otsu(abs(SourceData)/max(abs(SourceData)));
CLim = CLim*max(abs(SourceData));
end
% colorMap = jet;
cortexcolor=[0.9, 0.7, 0.7];
if ~isempty(SourceData)
cortexcolor=[.85 .85 .85];
end
alphaValue = 1;
nargs = nargin;
if nargs > 3
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'smooth'
                smooth = Value;
                if Value<=0 || Value>100
                   warning('smooth must be in (0,100]');
                   smooth = 30;
                end
            case 'cortexcolor'
                cortexcolor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'colormap'
                colorMap = Value;
                if strcmpi(Value,'active')||strcmpi(Value,'pnactive')||strcmpi(Value,'mjet')||strcmpi(Value,'invactive')||strcmpi(Value,'iactive')
                    colorMap = CBar(Value);
                elseif size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
            case 'clim'
                CLim = Value;
                if numel(CLim)~=2
                    error('Color Limit must be 1 x 2 vector');
                end
            case 'roiidx'
                roiIdx = Value;
            case 'axes'
                ax2plot = Value;
            case 'hemisphere'
                hemisphere =  Value;
            case 'thresh'               
                if thresh>1||thresh<0
                    warning('thresh must be in (0 1),using 0.05 instead');
                else
                    thresh = Value;
                end
            case 'alpha'
                alphaValue = Value;
            case 'position'
                figPosition = Value;
                if numel(figPosition) ~= 4
                     error('Position must be [left, bottom, width, height] ');
                end
        end
    end
end

if isempty(ax2plot)
    handles.h = figure('color',[0 0 0],...
        'MenuBar','none',...
        'ToolBar','none',...
        'DockControls','off');

    % 设置窗口位置（如果提供了）
    if ~isempty(figPosition)
        set(handles.h, 'Position', figPosition);
    end
    handles.axes = axes('Parent',handles.h,'color','k');
else
    set(ax2plot,'color','k');
    handles.axes = ax2plot;
    handles.h = ax2plot.Parent;
end

if ~isempty(smooth)
    SmoothValue=smooth/100;
    VertConn=BrainModel.VertConn;
    if ~isempty(SmoothValue)
        iVertices=1:size(vc,1);
%         Smoothing factor
        SurfSmoothIterations = ceil(300 * SmoothValue * length(iVertices) / 100000);
%         Calculate smoothed vertices locations
        loc_sm=vc;
        loc_sm(iVertices,:) = tess_smooth(vc(iVertices,:), SmoothValue, SurfSmoothIterations, VertConn(iVertices,iVertices), 1);
%         Apply smoothed locations
        vc=loc_sm;
        
    end
end
    
handles.hp=patch(handles.axes,'vertices',vc,'faces',tri,...
    'FaceColor',cortexcolor,'edgecolor','none','edgealpha',0.2,'facelighting','gouraud'...
    ,'specularstrength',0.2,'ambientstrength',0.5,'diffusestrength',0.5,...
    'BackfaceLighting','lit','AmbientStrength',0.5,'SpecularExponent',1,...
    'SpecularColorReflectance',0.5,'EdgeLighting','gouraud','facealpha',1);

material dull;
%camlight('headlight','infinite');
set(handles.axes,'Xcolor',[0 0 0],'ycolor',[0 0 0],'zcolor',[0 0 0],'CameraPosition', [0 0 0],'CameraViewAngle',6);
view([ -90 90 ]);
axis equal;
axlim=[min(vc);max(vc)]*1.2;
axis(reshape(axlim,1,[]));

light('Position',[-100,0,-100]*mean(sum(vc.^2,2)));
light('Position',[100,0,100]*mean(sum(vc.^2,2)));
set(handles.axes, 'CameraUpVector', [0, 1, 0]);

if ~isempty(SourceData)
    [Valpha,cdata]=bulidsource(SourceData,vc);
    set(handles.hp,'FaceAlpha','flat','AlphaDataMapping','none',...
        'FaceVertexAlphaData',Valpha, ...
        'FaceVertexCData',cdata,...
        'facecolor', 'interp');
end

bind3DInteraction(handles.h, handles.axes);

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



