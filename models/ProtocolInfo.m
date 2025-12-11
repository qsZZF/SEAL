classdef ProtocolInfo < handle
    %PROTOCOLINFO 协议管理信息类
    % 负责协议的元数据管理、文件路径管理

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        type string = ""  % 协议类型
        desc string = ""  % 协议描述
        
        % 元数据
        metadata struct % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime
    end
    
    methods
        function obj = ProtocolInfo(name, type, desc)
            %PROTOCOLINFO 构造函数
            % 输入:
            %   name - 协议名称
            %   type - 协议类型
            %   desc - 协议描述

            if nargin >= 1
                obj.name = name;
            end
            if nargin >= 2
                obj.type = type;
            end
            if nargin >= 3
                obj.desc = desc;
            end
            
            % 初始化时间戳
            obj.createdDate = datetime('now');
            obj.modifiedDate = obj.createdDate;
            
            % 初始化空元信息
            obj.metadata = struct();
        end
        
        %% 文件操作方法
        function save(obj, path)
            %SAVE 保存协议信息
            % 输入:
            %   path - 保存路径
            
            if path == ""
                error('SEAL:ProtocolInfo:InvalidPath', ...
                    'No valid protocol path specified.');
            end
            
            if ~isfolder(path)
                mkdir(path);
            end
            
            % 保存协议元数据
            obj.updateModifiedDate();
            saveData(fullfile(path, strcat(obj.name, '.mat')), obj);
        end
        
        %% 验证方法
        function isValid = validateName(obj)
            %VALIDATENAME 验证协议名称
            isValid = ~isempty(obj.name) && length(obj.name) <= 50;
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
        function protocol = createNew(protocolName, protocolType, desc)
            %CREATENEW 创建新协议
            % 输入:
            %   protocolName - 协议名称
            %   protocolType - 协议类型
            %   desc - 协议描述
            % 输出:
            %   protocol - 新创建的ProtocolInfo对象
            
            protocol = ProtocolInfo(protocolName, protocolType, desc);
        end
        
        function protocol = openExisting(protocolPath)
            %OPENEXISTING 打开现有协议
            % 输入:
            %   protocolPath - 协议文件路径
            % 输出:
            %   protocol - 加载的ProtocolInfo对象

            protocolFile = protocolPath;
            
            if ~isfile(protocolFile)
                error('SEAL:ProtocolInfo:FileNotFound', ...
                    'No such file exists: %s', protocolFile);
            end
            
            % 加载协议数据
            protocol = loadData(protocolFile);
            
            % 验证加载的对象类型
            if ~isa(protocol, 'ProtocolInfo')
                error('SEAL:ProtocolInfo:InvalidFileFormat', ...
                    'File does not contain a valid ProtocolInfo object: %s', protocolFile);
            end
        end
        
        function isValid = validateProtocolName(name)
            %VALIDATEPROTOCOLNAME 验证协议名称格式
            % 输入:
            %   name - 协议名称
            % 输出:
            %   isValid - 是否有效
            
            isValid = ischar(name) || isstring(name);
            if isValid
                nameStr = char(name);
                isValid = ~isempty(nameStr) && length(nameStr) <= 50;
            end
        end
    end
end
