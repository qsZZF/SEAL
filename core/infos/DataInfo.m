classdef DataInfo < handle
    %DATAINFO 信息管理类
    % 负责数据的元数据管理和存储

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        size int64
        srate int64
        chanlocs 
        isSimu = 0
        % 数据地址
        dataPath string
        resultPath string
        event 
        unit
        NoiseCovMat = []
        % 元数据
        metadata struct        % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime
    end
    
    properties (Dependent)
        dataCount int32     % 数量
    end
    
    methods
        function obj = DataInfo(name, dataPath, resultPath, desc, metadata)
            %DATAINFO 构造函数
            % 输入:
            %   name - 集合名称
            %   dataPath - 原始数据路径
            %   desc - 描述
            %   metadata - 元数据

            if nargin >= 1
                obj.name = name;
            end
            if nargin >= 2
                obj.dataPath = dataPath;
            else
                obj.dataPath = {};
            end
            if nargin >= 3
                obj.resultPath = resultPath;
            else
                obj.resultPath = {};
            end
            if nargin >= 4
                obj.desc = desc;
            end
            if nargin >= 5
                obj.metadata = metadata;
            else
                obj.metadata = struct();
            end
            
            % 初始化时间戳
            obj.createdDate = datetime('now');
            obj.modifiedDate = obj.createdDate;
            obj.dataPath = dataPath;
        end
        
        %% 文件操作方法
        function save(obj, path)
            %SAVE 保存信息

            checkDir(path);
            obj.updateModifiedDate();
            saveData(path, obj);
        end
    end
    
    methods (Access = private)
        function updateModifiedDate(obj)
            %UPDATEMODIFIEDDATE 更新修改时间
            obj.modifiedDate = datetime('now');
        end
    end
    
    %% 静态方法
    methods (Static)
        function data = openExisting(dataPath)
            %OPENEXISTING 打开现有信息
            % 输入:
            %   dataPath - 文件路径
            % 输出:
            %   data - 加载的DataInfo对象

            dataFile = dataPath;
            
            if ~isfile(dataFile)
                error('SEAL:DataInfo:FileNotFound', ...
                    'No such file exists: %s', dataFile);
            end
            
            % 加载数据
            data = loadData(dataFile);
            
            % 验证加载的对象类型
            if ~isa(data, 'DataInfo')
                error('SEAL:DataInfo:InvalidFileFormat', ...
                    'File does not contain a valid DataInfo object: %s', dataFile);
            end
        end

       function data = fromData(dataPath, varargin)
    %FROMDATA 智能导入数据 (兼容普通 .mat 与 Brainstorm 格式)
    %
    % 输入:
    %   dataPath - 数据文件路径
    %   varargin - 可选参数:
    %     'Name' - 名称（默认文件名）
    %     'Description' - 描述
    %     'Metadata' - 额外元数据
    
    % ==========================================
    % 0. 参数解析
    % ==========================================
    p = inputParser;
    addRequired(p, 'dataPath', @(x) isfile(x) && endsWith(x, '.mat', 'IgnoreCase', true));
    addParameter(p, 'Name', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Description', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Metadata', struct(), @isstruct);
    parse(p, dataPath, varargin{:});
    
    if isempty(p.Results.Name)
        [~, dataName, ~] = fileparts(dataPath);
    else
        dataName = p.Results.Name;
    end
    
    metadata = p.Results.Metadata;
    
    try
        % ==========================================
        % 1. 嗅探器：先看文件顶层有什么变量
        % ==========================================
        fileVars = whos('-file', dataPath);
        varNames = {fileVars.name};
        
        isBrainstorm = ismember('F', varNames) && ismember('Time', varNames);
        
        if isBrainstorm
            % ==========================================
            % 路线 A：Brainstorm 极速通道（延迟加载）
            % ==========================================
            metadata.originalFormat = 'Brainstorm';
            
            % 不用 load，直接从 whos 拿 size
            idxF = strcmp(varNames, 'F');
            dataSize = fileVars(idxF).size;
            
            % 只加载轻量级变量
            varsToLoad = {'Time', 'Comment', 'Events', 'ChannelFlag', 'ChannelLocs'};
            varsExist = varsToLoad(ismember(varsToLoad, varNames));
            lightData = load(dataPath, varsExist{:}, '-mat');
            
            % 采样率从 Time 向量计算
            if isfield(lightData, 'Time') && length(lightData.Time) > 1
                dataSrate = 1 / (lightData.Time(2) - lightData.Time(1));
            else
                dataSrate = 1000;
            end
            
            % 元数据
            if isfield(lightData, 'Events')
                metadata.events = lightData.Events;
            end
            if isfield(lightData, 'ChannelFlag')
                metadata.badChannels = lightData.ChannelFlag;
            end
            
            % Brainstorm 的通道位置（如果存在）
            if isfield(lightData, 'ChannelLocs')
                dataChanlocs = lightData.ChannelLocs;
            else
                dataChanlocs = [];
            end
            
            % 名字优先用 Comment
            if isempty(p.Results.Name) && isfield(lightData, 'Comment')
                dataName = lightData.Comment;
            end
            
            metadata.targetField = 'F';

            try
                matF = matfile(dataPath);
                fSize = size(matF, 'F');
                sampleCols = min(1000, fSize(2));
                sampleData = matF.F(:, 1:sampleCols);
                dataUnit = DataInfo.detectDataUnit(sampleData);
            catch
                % matfile 失败（比如文件不是 v7.3 格式），回退到 Brainstorm 默认
                dataUnit = "V";
            end
            
        else
            % ==========================================
            % 路线 B：通用 .mat 通道
            % ==========================================
            metadata.originalFormat = 'GenericMat';
            
            loadedData = load(dataPath, '-mat');
            
            % 拆俄罗斯套娃
            outerFields = fieldnames(loadedData);
            if length(outerFields) == 1 && isstruct(loadedData.(outerFields{1}))
                loadedData = loadedData.(outerFields{1});
            end
            
            availableFields = fieldnames(loadedData);
            
            % 字典
            dataAliases     = {'data', 'Data', 'DataFile', 'matrix', 'signal', 'EEG'};
            srateAliases    = {'srate', 'SamplingRate', 'fs', 'sfreq', 'sRate'};
            chanlocsAliases = {'chanlocs', 'Channel', 'channels', 'chaninfo', ...
                               'electrodes', 'locs', 'labels'};
            
            matchedDataField     = findFirstMatch(dataAliases,     availableFields);
            matchedSrateField    = findFirstMatch(srateAliases,    availableFields);
            matchedChanlocsField = findFirstMatch(chanlocsAliases, availableFields);
            
            % 数据字段：必须匹配，否则弹窗
            if isempty(matchedDataField)
                [indx, tf] = listdlg('ListString', availableFields, ...
                    'SelectionMode', 'single', ...
                    'PromptString', {'未识别到标准的【数据矩阵】字段！', ...
                                     '请在下方选择代表脑电信号的变量:'}, ...
                    'Name', 'SEAL 数据导入向导 (1/2)', ...
                    'ListSize', [250, 150]);
                if tf
                    matchedDataField = availableFields{indx};
                else
                    error('SEAL:Import:UserCancelled', '用户取消了数据字段匹配。');
                end
            end
            
            % 采样率字段：必须匹配，否则弹窗
            if isempty(matchedSrateField)
                [indx, tf] = listdlg('ListString', availableFields, ...
                    'SelectionMode', 'single', ...
                    'PromptString', {'未识别到标准的【采样率】字段！', ...
                                     '请在下方选择代表采样率的变量:'}, ...
                    'Name', 'SEAL 数据导入向导 (2/2)', ...
                    'ListSize', [250, 150]);
                if tf
                    matchedSrateField = availableFields{indx};
                else
                    error('SEAL:Import:UserCancelled', '用户取消了采样率字段匹配。');
                end
            end
            
            % 提取尺寸和采样率
            dataMatrix = loadedData.(matchedDataField);
            dataSize = size(dataMatrix);
            tempSrate = loadedData.(matchedSrateField);
            dataSrate = tempSrate(1);   % 防御：可能是数组
            
            % 单位检测
            dataUnit = DataInfo.detectDataUnit(dataMatrix);
            
            % 通道位置（可选）
            if ~isempty(matchedChanlocsField)
                dataChanlocs = loadedData.(matchedChanlocsField);
            else
                dataChanlocs = [];
            end
            
            metadata.targetField = matchedDataField;
        end
        
        % ==========================================
        % 2. 统一出口：创建 DataInfo
        % ==========================================
        data = DataInfo(...
            dataName, ...
            dataPath, ...
            "", ...
            p.Results.Description, ...
            metadata);
        
        data.size     = dataSize;
        data.srate    = round(dataSrate);
        data.unit     = dataUnit;
        data.chanlocs = dataChanlocs;
        
    catch ME
        error('SEAL:DataInfo:ImportFailed', ...
            '数据提取失败: %s\n详细原因: %s', dataPath, ME.message);
    end
end

% 辅助函数：在候选字段中找第一个匹配的别名
function matched = findFirstMatch(aliases, availableFields)
    matched = '';
    for i = 1:length(aliases)
        if ismember(aliases{i}, availableFields)
            matched = aliases{i};
            return;
        end
    end
end

        function data = fromEEGLAB(dataPath, varargin)
            %FROMEEGLAB 从 EEGLAB .set 文件直接创建 DataInfo
            % 输入:
            %   dataPath - .set 数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 名称（可选，默认为文件名）
            %     'Description' - 描述
            %     'Metadata' - 额外元数据
            % 输出:
            %   data - 组装好的 DataInfo 对象
            
            % ==========================================
            % 1. 参数解析
            % ==========================================
            p = inputParser;
            addRequired(p, 'dataPath', @(x) isfile(x) && endsWith(x, '.set', 'IgnoreCase', true));
            addParameter(p, 'Name', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Description', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Metadata', struct(), @isstruct);
            
            parse(p, dataPath, varargin{:});
            
            % 如果未提供名称，使用文件名
            if isempty(p.Results.Name)
                [~, name, ~] = fileparts(dataPath);
                dataName = name;
            else
                dataName = p.Results.Name;
            end
            
            try
                % ==========================================
                % 2. 加载 .set 文件 (本质为 .mat 文件)
                % ==========================================
                tempStruct = load(dataPath, '-mat');
                
                if ~isfield(tempStruct, 'EEG')
                    error('文件内未找到 "EEG" 结构体，这可能不是标准的 EEGLAB .set 文件。');
                end
                
                EEG = tempStruct.EEG;
                metadata = p.Results.Metadata;
                % ==========================================
                % 3. 核心：判断并读取外部 .fdt 数据文件
                % ==========================================
                if ischar(EEG.data) || isstring(EEG.data)
                    fdtFilename = char(EEG.data);
                    [setDir, ~, ~] = fileparts(dataPath);
                    
                    % 拼接出 .fdt 文件的完整路径
                    fdtPath = fullfile(setDir, fdtFilename);
                    
                    if ~isfile(fdtPath)
                        error('SEAL:DataInfo:MissingFDT', ...
                            '数据依赖外部文件但未找到：%s\n请确保 .set 和 .fdt 文件在同一文件夹下。', fdtFilename);
                    end
                    
                    
                % ==========================================
                % 4. 提取与封装元数据 (Metadata)
                % ==========================================
                
                metadata.fdtPath = fdtPath;
                end
                metadata.originalFormat = 'EEGLAB';
                metadata.nbchan = EEG.nbchan;
                metadata.pnts   = EEG.pnts;
                metadata.trials = EEG.trials;
                
                data = DataInfo(...
                    dataName, ...
                    dataPath, ...
                    "", ...
                    p.Results.Description, ...
                    metadata);
                
                % 核心属性赋值
                if ischar(EEG.data) || isstring(EEG.data)
                    % 数据在外部 .fdt 文件中，用结构体参数计算 size
                    if EEG.trials > 1
                        data.size = [EEG.nbchan, EEG.pnts, EEG.trials];
                    else
                        data.size = [EEG.nbchan, EEG.pnts];
                    end
                else
                    % 数据直接内嵌，直接取 size
                    data.size = size(EEG.data);
                end
                data.srate = EEG.srate;
                data.unit  = "μV"; 
                 % 提取通道标签（如果存在）
                if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'labels')
                    data.chanlocs = EEG.chanlocs;
                end
                
                % 提取事件/Marker信息（如果存在）
                if isfield(EEG, 'event') && ~isempty(EEG.event)
                    data.event = EEG.event;
                end
   
            catch ME
                % 统一抛出异常，方便上层（如 SessionNode）捕获并弹窗
                error('SEAL:DataInfo:FromEEGLABFailed', ...
                    'EEGLAB 数据解析失败: %s\n详细原因: %s', dataPath, ME.message);
            end
        end



        function unit = detectDataUnit(dataMatrix)
            %DETECTDATAUNIT 根据数据数量级猜测单位
            %
            % 判断依据：典型 EEG 信号的绝对值范围
            %   伏特 (V)  : 10^-6 ~ 10^-4 → 最大值 < 0.01
            %   毫伏 (mV) : 10^-3 ~ 10^-1 → 最大值 0.01 ~ 10
            %   微伏 (μV) : 1 ~ 100       → 最大值 10 ~ 1000

            % 用 99 分位数而不是 max，避开尖峰干扰
            magnitude = prctile(abs(dataMatrix(:)), 99);

            if magnitude < 0.01
                unit = "V";
            elseif magnitude < 10
                unit = "mV";
            else
                unit = "μV";
            end
        end


    end
end