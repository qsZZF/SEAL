function CreateProcessTemplate(savePath , filename ,varargin)
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
fprintf(fid, 'function output = %s(input)\n', funcName);
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
    fprintf(fid, 's = load(infoFile); protocolInfo = s.protocolInfo;\n');
    fprintf(fid, 'CortexFileList = dir(fullfile(''%s'',''*cortex*.mat''));\n',protocolPath);
    fprintf(fid, 'LeadfieldFileList =dir( fullfile(''%s'',''*leadfield*.mat''));\n',protocolPath);
    fprintf(fid, 'DatafileList =dir( fullfile(''%s'',''*data*.mat''));\n',protocolPath);
    fprintf(fid, 'Data = load(fullfile(DatafileList(1).folder, DatafileList(1).name));\n');
    fprintf(fid, 'Cortex = load(fullfile(CortexFileList(1).folder, CortexFileList(1).name));\n');
    fprintf(fid, 'Leadfield =load(fullfile(LeadfieldFileList(1).folder, LeadfieldFileList(1).name));\n');
    fprintf(fid, '\n');
    fprintf(fid, 'NoiseCovariance=eye(size(Data,1),size(Data,1));\n');%协方差矩阵导入功能？

    for i = 1:length(fieldNames)
        field = fieldNames{i};
        switch field
            case 'filters'
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
            case 'rereference'

            case 'ESIAlgrithms'
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
                            fprintf(fid, '[Source_MNE,~]=seal_MNE(Data,Leadfield.Gain,''RegularizationParameter'',%f,...\n', ...
                                Algrithms.RegularizationParameter);
                            fprintf(fid, '''DepthWeightExp'',%f,...\n',Algrithms.DepthWeightExp);
                            fprintf(fid, '''SourceCovariance'',%s,...\n',Algrithms.SourceCovariance);
                            fprintf(fid, '''NumOrientations'',%d,...\n',Algrithms.NumOrientations);
                            fprintf(fid, '''NoiseCovariance'',NoiseCovariance);\n');
                        case 'BlockChampagne'
                            fprintf(fid, '[Source_BlockChampagne, ~] = seal_BlockChampagne(Data,Leadfield.Gain,Cortex, ...\n');
                            fprintf(fid, '''NoiseCovariance'',NoiseCovariance,...\n');
                            fprintf(fid, '''LearnNoise'',LearnNoise ,...\n');
                            fprintf(fid, '''NeighborOrder'',%d,...\n',Algrithms.NeighborOrder);
                            fprintf(fid, '''NumOrientations'',%d,...\n',Algrithms.NumOrientations);
                            fprintf(fid, '''NoiseObservationMatrix'', eye(size(app.Data,1)));\n');
                        case 'uSTAR'
                            fprintf(fid, '[Source_uSTAR, ~] = seal_uSTAR(Data,Leadfield.Gain,Cortex, ...\n');
                            fprintf(fid,'''RunAlgorithm'',%c,...\n',Algrithms.RunAlgorithm);

                            fprintf(fid,'''MicrostateMaps'',MicrostateMaps);\n');
                        case 'STARTS'
                            fprintf(fid, '[Source_STARTS, ~] = seal_STARTS(Data,Leadfield.Gain,Cortex, ...\n');
                            fprintf(fid, '''NoiseCovariance'',NoiseCovariance, ...\n');
                            fprintf(fid, '''NumTBFs'',%d,...\n',Algrithms.NumTBFs);
                            fprintf(fid, '''LearnNoise'',LearnNoise ,...\n');
                            fprintf(fid, '''NoiseObservationMatrix'', eye(size(app.Data,1)),...\n');
                            fprintf(fid, '''NeighborOrder'',%d,...\n',Algrithms.NeighborOrder);
                            fprintf(fid, '''NumOrientations'',%d);\n',Algrithms.NumOrientations);
                    end
                end
        end
    end
end










fprintf(fid, '\n');
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
fprintf(fid, '\n');    fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'end\n');

fclose(fid);
end
