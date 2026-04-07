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
        % 数据地址
        dataPath string
        resultPath string
        event 
        unit
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
            %FROMDATA 从数据文件直接创建DataInfo
            % 输入:
            %   dataPath - 数据文件路径
            %   varargin - 可选参数:
            %     'Name' - 名称（可选，默认为文件名）
            %     'Description' - 描述
            %     'Metadata' - 额外元数据
            % 输出:
            %   data - DataInfo对象
            
            % 参数解析
            p = inputParser;
            addRequired(p, 'dataPath', @(x) isfile(x) && endsWith(x, '.mat'));
            addParameter(p, 'Name', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Description', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Metadata', struct(), @isstruct);
            addParameter(p, 'AutoExtract', true, @islogical);
            
            parse(p, dataPath, varargin{:});
            
            % 如果未提供名称，使用文件名
            if isempty(p.Results.Name)
                [~, name, ~] = fileparts(dataPath);
                dataName = name;
            else
                dataName = p.Results.Name;
            end
            
            % 从数据文件中提取元数据
            metadata = p.Results.Metadata;
            
           % 创建DataInfo对象
            try
                data = DataInfo(...
                    dataName, ...
                    dataPath, ...
                    "", ...
                    p.Results.Description, ...
                    metadata);
                
                % 1. 获取文件中的所有变量名 (使用 load 获取 struct 的 fieldnames)
                loadedData = load(dataPath); % 如果你之前的 loadData 返回的是 struct，就继续用 loadData
                outerFields = fieldnames(loadedData);
                if length(outerFields) == 1 && isstruct(loadedData.(outerFields{1}))
                    % 把内层的结构体提取出来，作为真正的 loadedData
                    loadedData = loadedData.(outerFields{1});
                end
                
                availableFields = fieldnames(loadedData);
                
                % ==========================================
                % 2. 定义自动匹配字典 (Auto-match Dictionaries)
                % ==========================================
                dataAliases = {'data', 'Data', 'DataFile', 'matrix', 'signal', 'EEG'};
                srateAliases = {'srate', 'SamplingRate', 'fs', 'sfreq', 'sRate'};
                chanlocsAliases = {'chanlocs', 'Channel', 'channels', 'chaninfo', 'electrodes', 'locs', 'labels'};

                matchedDataField = '';
                matchedSrateField = '';
                matchedchanlocsField = '';
                % --- 尝试自动匹配 Data ---
                for i = 1:length(dataAliases)
                    if ismember(dataAliases{i}, availableFields)
                        matchedDataField = dataAliases{i};
                        break;
                    end
                end
                
                % --- 尝试自动匹配 Srate ---
                for i = 1:length(srateAliases)
                    if ismember(srateAliases{i}, availableFields)
                        matchedSrateField = srateAliases{i};
                        break;
                    end
                end
                
                for i = 1:length(chanlocsAliases)
                    if ismember(chanlocsAliases{i}, availableFields)
                        matchedchanlocsField = chanlocsAliases{i};
                        break;
                    end
                end

                % 找不到数据矩阵时触发弹窗
                if isempty(matchedDataField)
                    [indx, tf] = listdlg('ListString', availableFields, ...
                        'SelectionMode', 'single', ...
                        'PromptString', {'未识别到标准的【数据矩阵】字段！', '请在下方选择代表脑电信号的变量:'}, ...
                        'Name', 'SEAL 数据导入向导 (1/2)', ...
                        'ListSize', [250, 150]);
                    
                    if tf % tf 为 true 说明用户点了确定
                        matchedDataField = availableFields{indx};
                    else
                        error('SEAL:Import:UserCancelled', '用户取消了数据字段匹配，导入中止。');
                    end
                end
                
                % 找不到采样率时触发弹窗
                if isempty(matchedSrateField)
                    [indx, tf] = listdlg('ListString', availableFields, ...
                        'SelectionMode', 'single', ...
                        'PromptString', {'未识别到标准的【采样率】字段！', '请在下方选择代表采样率的变量:'}, ...
                        'Name', 'SEAL 数据导入向导 (2/2)', ...
                        'ListSize', [250, 150]);
                    
                    if tf
                        matchedSrateField = availableFields{indx};
                    else
                        % 对于采样率，你也可以选择不中断，而是给个默认值
                        % warning('用户跳过了采样率匹配，默认设置为 1000 Hz');
                        % matchedSrateField = 'default';
                        error('SEAL:Import:UserCancelled', '用户取消了采样率字段匹配，导入中止。');
                    end
                end
                
  
                data.size = size(loadedData.(matchedDataField));
                
                % 防御：如果上一步采样率用户没选（或采取了给默认值的策略）
                if strcmp(matchedSrateField, 'default')
                    data.srate = 1000;
                else
                    data.srate = loadedData.(matchedSrateField);
                end
                data.unit = DataInfo.detectDataUnit(loadedData.(matchedDataField));
                % [新增] 尝试为通道信息赋值
                if ~isempty(matchedchanlocsField)
                    % 如果成功匹配到了通道变量，将其存入你的数据对象中
                    data.chanlocs = loadedData.(matchedchanlocsField);
                else
                    data.chanlocs = [];                
                end

            catch ME
                % 精准抛出错误
                error('SEAL:DataInfo:ImportFailed', ...
                    '数据提取失败: %s\n详细原因: %s', dataPath, ME.message);
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