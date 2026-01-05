function CreateProcessTemplate(app, savePath , filename ,varargin)
% ensure the file is ended by .m
if ~endsWith(filename, '.m')
    filename = filename + ".m";
end

% get the function's name
funcName = erase(filename, '.m');

% check the savepath
if ~exist(savePath, 'dir')
    error("保存路径不存在: %s", savePath);
end
addpath(savePath);
fullFile = fullfile(savePath, filename);

% open the file
fid = fopen(fullFile, 'w');
if fid == -1
    error("无法创建文件: %s", fullFile);
end
nargs=nargin;
if nargs > 3
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'protocolpath'
                protocolPath=Value;
                processmode='single_protocol';
        end
    end
end




description1 = 'This is a automatically generated file which can be easily applied to your new protocols';
description2 =  ' with the same settings as you have used before';
fprintf(fid, 'function %s(app, protocolnode)\n', funcName);
fprintf(fid, '%% %s\n', upper(funcName));
fprintf(fid, '%% description: %s\n', description1);
fprintf(fid, '%%     %s\n', description2);

if isequal(processmode, 'single_protocol')
    ProcessinginfoFile = fullfile(protocolPath,'ProcessingInfo.mat');
    Processinginfo=load(ProcessinginfoFile);
    disp(Processinginfo);
    fieldNames = fieldnames(Processinginfo);

    %import the files the precess needed
    fprintf(fid, 'infoFile = fullfile(''%s'',''protocol_info.mat'');\n',protocolPath);
    fprintf(fid, 'global protocolInfo;\n');
    fprintf(fid, 's = load(infoFile); protocolInfo = s.protocolInfo;\n');
    fprintf(fid, 'CortexFileList = dir(fullfile(''%s'',''*cortex*.mat''));\n',protocolPath);
    fprintf(fid, 'LeadfieldFileList =dir( fullfile(''%s'',''*leadfield*.mat''));\n',protocolPath);
    fprintf(fid, 'DatafileList =dir( fullfile(''%s'',''*data*.mat''));\n',protocolPath);
    fprintf(fid, 'Data = load(fullfile(DatafileList(1).folder, DatafileList(1).name));\n');
    fprintf(fid, 'if isfield(Data,''data'')\n');
    fprintf(fid, '   Data=Data.data;\n');
    fprintf(fid, 'elseif isfield(Data,''Data'')\n');
    fprintf(fid, '   Data=Data.Data;\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'if ndims(Data)==3\n');
    fprintf(fid, '   Data= mean(Data, 3);\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'Cortex = load(fullfile(CortexFileList(1).folder, CortexFileList(1).name));\n');
    fprintf(fid, 'LeadField =load(fullfile(LeadfieldFileList(1).folder, LeadfieldFileList(1).name));\n');
    fprintf(fid, '\n');
    fprintf(fid, 'if isfield(LeadField,''GridAtlas'')\n');
    fprintf(fid, '    LeadField_fixed.Gain=...\n');
    fprintf(fid, '      bst_gain_orient(LeadField.Gain, ...\n');
    fprintf(fid, '      LeadField.GridOrient, ...\n');
    fprintf(fid, '      LeadField.Atlas);\n');
    fprintf(fid, 'else\n');
    fprintf(fid, '    LeadField_fixed.Gain=...\n');
    fprintf(fid, '      bst_gain_orient(LeadField.Gain, ...\n');
    fprintf(fid, '      LeadField.GridOrient, ...\n');
    fprintf(fid, '      []);\n');
    fprintf(fid, 'end\n');

    fprintf(fid, '\n');
    fprintf(fid, 'NoiseCovariance=eye(size(Data,1),size(Data,1));\n');%协方差矩阵导入功能？
    fprintf(fid, 'ResultSavePath = fullfile(protocolnode.NodeData.protocolPath, ''ResultData'');\n');
    fprintf(fid, 'if ~exist(ResultSavePath, ''dir'')\n');
    fprintf(fid, '    mkdir(ResultSavePath);\n');
    fprintf(fid, '    disp([''文件夹已创建: '' fullfile(pwd, ResultSavePath)]);\n');
    fprintf(fid, 'else\n');
    fprintf(fid, '    disp([''文件夹已存在: '' fullfile(pwd, ResultSavePath)]);\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'fileList = dir(fullfile(ResultSavePath, ''*.mat''));\n');
    fprintf(fid, 'global saveindex resultstruct; \n');
    fprintf(fid, 'resultstruct=struct();\n');
    fprintf(fid, 'counters = zeros(50,1);%%counter of each Algrithm\n');
    fprintf(fid, 'saveindex= length(fileList);%%counter of exsited resultfiles\n');
    fprintf(fid, 'counters = getIndex(fileList,counters);\n');
    preprocess=false;

        
      if  double(any(ismember(fieldNames,'filters')))
            for j = 1:size(Processinginfo.filters,2)
                filtertypes = Processinginfo.filters(j).filterTypes;
                switch  filtertypes
                    case 'highpass'
                        fprintf(fid, '\n');
                        fprintf(fid, 'highpass = designFIRHighpass(protocolInfo.metadata.SamplingRate,%f,%d,%d); \n', ...
                            Processinginfo.filters(j).FcArray, ...
                            Processinginfo.filters(j).RsArray, ...
                            Processinginfo.filters(j).NArray);
                        fprintf(fid, '[nChannels, ~] = size(Data);\n');
                        fprintf(fid, 'for i = 1:nChannels\n');
                        fprintf(fid, '    tempData = Data(i, :);\n');
                        fprintf(fid, '    Data(i, :)=filtfilt(highpass,1,tempData);\n');
                        fprintf(fid, 'end\n');
                        continue
                    case 'lowpass'
                        fprintf(fid, '\n');
                        fprintf(fid, 'lowpass = designFIRLowpass(protocolInfo.metadata.SamplingRate,%f,%d,%d); \n', ...
                            Processinginfo.filters(j).FcArray, ...
                            Processinginfo.filters(j).RsArray, ...
                            Processinginfo.filters(j).NArray);
                        fprintf(fid, '[nChannels, ~] = size(Data);\n');
                        fprintf(fid, 'for i = 1:nChannels\n');
                        fprintf(fid, '    tempData = Data(i, :);\n');
                        fprintf(fid, '    Data(i, :)=filtfilt(lowpass,1,tempData);\n');
                        fprintf(fid, 'end\n');
                        continue
                    case 'notch'
                        fprintf(fid, '\n');
                        fprintf(fid, 'notch = designFIRBandstop(protocolInfo.metadata.SamplingRate,%f,%f,%d,%d); \n', ...
                            Processinginfo.filters(j).FcArray(1), ...
                            Processinginfo.filters(j).FcArray(2), ...
                            Processinginfo.filters(j).RsArray, ...
                            Processinginfo.filters(j).NArray);
                        fprintf(fid, '[nChannels, ~] = size(Data);\n');
                        fprintf(fid, 'for i = 1:nChannels\n');
                        fprintf(fid, '    tempData = Data(i, :);\n');
                        fprintf(fid, '    Data(i, :)=filtfilt(notch,1,tempData);\n');
                        fprintf(fid, 'end\n');
                        continue
                end
            end
                preprocess=true;
      end
        %     case 'rereference'

        if  preprocess==true && double(any(ismember(fieldNames,'ESIAlgrithms')))
            fprintf(fid,'MicrostateMaps=[];\n');
            for j = 1:size(Processinginfo.ESIAlgrithms,2)
                Algrithms = Processinginfo.ESIAlgrithms{1,j};
                if isfield(Algrithms,'LearnNoise')
                    if Algrithms.LearnNoise
                        fprintf(fid,'LearnNoise=true;\n');
                    else
                        fprintf(fid, 'LearnNoise=false;\n');
                    end
                end

                switch Algrithms.Algrithm
                    case 'MNE'
                        fprintf(fid, '[Source_MNE,~]=seal_MNE(Data,LeadField.Gain,''RegularizationParameter'',%f,...\n', ...
                            Algrithms.RegularizationParameter);
                        fprintf(fid, '  ''DepthWeightExp'',%f,...\n',Algrithms.DepthWeightExp);
                        fprintf(fid, '  ''SourceCovariance'',''%s'',...\n',Algrithms.SourceCovariance);
                        fprintf(fid, '  ''NumOrientations'',%d,...\n',Algrithms.NumOrientations);
                        fprintf(fid, '  ''NoiseCovariance'',NoiseCovariance);\n');
                        fprintf(fid, 'counters(1)=counters(1)+1;\n');
                        fprintf(fid, 'saveindex = saveindex +1;\n');
                        fprintf(fid, 'type = ''MNE'';\n');
                        fprintf(fid, 'temp = reshape(Source_MNE, 3, [], size(Source_MNE, 2));\n');
                        fprintf(fid, 'Source_MNE = squeeze(sum(temp.^2, 1)); \n');
                        fprintf(fid, 'saveresult(app,Source_MNE,protocolnode,type,counters(1),ResultSavePath)\n');
                        fprintf(fid, '\n');
                    case 'BlockChampagne'
                        fprintf(fid, '[Source_BlockChampagne, ~] = seal_BlockChampagne(Data,LeadField_fixed.Gain,Cortex, ...\n');
                        fprintf(fid, '  ''NoiseCovariance'',NoiseCovariance,...\n');
                        fprintf(fid, '  ''LearnNoise'',LearnNoise ,...\n');
                        fprintf(fid, '  ''NeighborOrder'',%d,...\n',Algrithms.NeighborOrder);
                        fprintf(fid, '  ''NumOrientations'',%d,...\n',Algrithms.NumOrientations);
                        fprintf(fid, '  ''NoiseObservationMatrix'', eye(size(Data,1)));\n');
                        fprintf(fid, 'counters(8)=counters(8)+1;\n');
                        fprintf(fid, 'saveindex = saveindex +1;\n');
                        fprintf(fid, 'type = ''BlockChampagne'';\n');
                        fprintf(fid, 'saveresult(app,Source_BlockChampagne,protocolnode,type,counters(8),ResultSavePath)\n');
                        fprintf(fid, '\n');
                    case 'uSTAR'
                        fprintf(fid, '[Source_uSTAR, ~] = seal_uSTAR(Data,LeadField_fixed.Gain,Cortex, ...\n');
                        fprintf(fid, '  ''RunAlgorithm'',''%s'',...\n',Algrithms.RunAlgorithm);
                        fprintf(fid, '  ''MicrostateMaps'',MicrostateMaps);\n');
                        fprintf(fid, 'counters(9)=counters(9)+1;\n');
                        fprintf(fid, 'saveindex = saveindex +1;\n');
                        fprintf(fid, 'type = ''uSTAR'';\n');
                        fprintf(fid, 'saveresult(app,Source_uSTAR,protocolnode,type,counters(9),ResultSavePath)\n');
                        fprintf(fid, '\n');
                    case 'STARTS'
                        fprintf(fid, '[Source_STARTS, ~] = seal_STARTS(Data,LeadField_fixed.Gain,Cortex, ...\n');
                        fprintf(fid, '  ''NoiseCovariance'',NoiseCovariance, ...\n');
                        fprintf(fid, '  ''NumTBFs'',%d,...\n',Algrithms.NumTBFs);
                        fprintf(fid, '  ''LearnNoise'',LearnNoise ,...\n');
                        fprintf(fid, '  ''NoiseObservationMatrix'', eye(size(Data,1)),...\n');
                        fprintf(fid, '  ''NeighborOrder'',%d,...\n',Algrithms.NeighborOrder);
                        fprintf(fid, '  ''NumOrientations'',%d);\n',Algrithms.NumOrientations);
                        fprintf(fid, 'counters(7)=counters(7)+1;\n');
                        fprintf(fid, 'saveindex = saveindex +1;\n');
                        fprintf(fid, 'type = ''STARTS'';\n');
                        fprintf(fid, 'saveresult(app,Source_STARTS,protocolnode,type,counters(7),ResultSavePath)\n');
                        fprintf(fid, '\n');
                end
            end
        end
    fprintf(fid, 'end\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    targetText = 'Applyprocess...';  % 菜单文本（必须完全一致）
    targetSeparator = 'on';          % 分隔符状态（与创建代码一致）
    menuExists = false;
    if  isvalid(app.protocolmenu) && ~isempty(app.protocolmenu.Children)
        % 遍历 app.protocolmenu 的所有子菜单
        for idx = 1:length(app.protocolmenu.Children)
            childMenu = app.protocolmenu.Children(idx);
            % 过滤：确保是菜单对象，且匹配文本和分隔符
            if isa(childMenu, 'matlab.ui.container.Menu') && ...
                    strcmp(childMenu.Text, targetText) && ...
                    strcmp(childMenu.Separator, targetSeparator)
                % 找到匹配的已创建菜单，标记为存在
                menuExists = true;
                parentMenu = childMenu;
                break;
            end
        end
    end

    % 若未创建，则执行创建代码
    if ~menuExists
        parentMenu = uimenu(app.protocolmenu, 'Text', targetText, 'Separator', targetSeparator);
    end

    targetFuncName = funcName;

    % 初始化：假设子菜单未创建
    subMenuExists = false;

    if isvalid(parentMenu) && ~isempty(parentMenu.Children)
        % 遍历 parentMenu 的所有子菜单
        for idx = 1:length(parentMenu.Children)
            child = parentMenu.Children(idx);
            % 过滤：仅匹配「Menu 对象」且「Text 与 funcName 完全一致」
            if isa(child, 'matlab.ui.container.Menu') && strcmp(child.Text, targetFuncName)
                subMenuExists = true;
                break;
            end
        end
    end

    if ~subMenuExists
        uimenu(parentMenu, 'Text', targetFuncName, 'MenuSelectedFcn', @(s,e) applyprocess(app, filename));
    end

    fprintf(fid, '\n');
    fprintf(fid, 'function saveresult(app,results,protocolnode,type,Num,SavePath)\n');

    fprintf(fid, '    children = protocolnode.Children; \n');

    % 1. 如果没有子节点，直接返回空
    fprintf(fid, '    if isempty(children)\n');
    fprintf(fid, '        DataNode = [];\n');
    fprintf(fid, '        return;\n');
    fprintf(fid, '    end\n');

    % 2. 获取所有子节点的 Text 属性，并转换为 string 数组
    fprintf(fid, '    childTexts = string({children.Text});\n');

    % 3. 找出以 'Data' 开头的索引 (返回逻辑数组 [1 0 1 ...])
    fprintf(fid, '    isMatch = startsWith(childTexts, ''Data'');\n');

    % 4. 使用逻辑索引提取对应的节点对象
    fprintf(fid, '    DataNode = children(isMatch);\n');
    fprintf(fid, '    global saveindex resultstruct; \n');
    fprintf(fid, '    resultstruct(saveindex).filepath= sprintf(''%%s/%%s_Results_%%03d.mat'', SavePath, type, Num);\n');
    fprintf(fid, '    ID=sprintf(''%%s_Results_%%03d'', type, Num);\n');
    fprintf(fid, '    resultstruct(saveindex).ID=[''Result: '',ID];\n');
    fprintf(fid, '    resultstruct(saveindex).Source = DataNode.Text;\n');
    fprintf(fid, '    save(resultstruct(saveindex).filepath, ''results'', ''-v7.3'');\n');

    fprintf(fid, '    global protocolInfo;\n');
    fprintf(fid, '    protocolInfo.ResultData=resultstruct;\n');
    fprintf(fid, '    infoFilePath = fullfile(protocolnode.NodeData.protocolPath,''protocol_info.mat'');\n');
    fprintf(fid, '    save(infoFilePath, ''protocolInfo'');\n');
    fprintf(fid, '    node = uitreenode(DataNode,''Text'',ID,''NodeData'',resultstruct(saveindex).filepath);\n');
    fprintf(fid, '    node.ContextMenu = app.ResultDataFileMenu;\n');
    fprintf(fid, '    drawnow;\n');
    fprintf(fid, 'end\n');

    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, 'function counters = getIndex(filelist,counters)\n');
    fprintf(fid, '    for i=1:length(filelist)\n');
    fprintf(fid, '        if contains(filelist(i).name, ''MNE'')\n');
    fprintf(fid, '            counters(1)=counters(1)+1;\n');
    fprintf(fid, '        elseif contains(filelist(i).name, ''STARTS'')\n');
    fprintf(fid, '            counters(7)=counters(7)+1;\n');
    fprintf(fid, '        elseif contains(filelist(i).name, ''BlockChampagne'')\n');
    fprintf(fid, '            counters(8)=counters(8)+1;\n');
    fprintf(fid, '        elseif contains(filelist(i).name, ''uSTAR'')\n');
    fprintf(fid, '            counters(9)=counters(9)+1;\n');
    fprintf(fid, '        end\n');
    fprintf(fid, '    end\n');
    fprintf(fid, 'end\n');
end


fprintf(fid, '\n');
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');

fclose(fid);
end

function applyprocess(app,filename)
protocolnode=app.FileTree.SelectedNodes;
[~, name, ~] = fileparts(filename);  % 拆分路径、文件名、扩展名
funcHandle = str2func(name);
funcHandle(app,protocolnode);
end



