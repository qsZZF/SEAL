classdef SessionInfo < handle
    %SESSIONINFO 会话信息管理类
    % 负责会话的元数据管理，不直接管理数据节点

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        type string = ""  % 会话类型
        desc string = ""         % 会话描述
        
        % 元数据
        metadata struct                 % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime

    end
    
    methods
        function obj = SessionInfo(name, type, desc)
            %SESSIONINFO 构造函数
            % 输入:
            %   name - 会话名称
            %   type - 会话类型
            %   desc - 会话描述

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
            
            % 初始化元数据
            obj.metadata = struct();
        end
        
        %% 文件操作方法
        function save(obj, path)
            %SAVE 保存会话信息
            % 输入:
            %   path - 保存路径
            
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
        function session = createNew(name, varargin)
            %CREATENEW 创建新会话
            % 输入:
            %   name - 会话名称
            %   type - 会话类型
            %   desc - 会话描述
            %   varargin - 可选参数对:
            %     'SubjectID' - 被试ID
            %     'Experimenter' - 实验者
            %     'ExperimentDate' - 实验日期
            %     'ProtocolName' - 协议名称
            % 输出:
            %   session - 新创建的SessionInfo对象
            
            p = inputParser;
            addRequired(p, 'name', Validators.isText);
            
            parse(p, name);
            
            % 创建SessionInfo对象
            session = SessionInfo(...
                p.Results.name ...
                );
            
        end
        
        function session = openExisting(sessionPath)
            %OPENEXISTING 打开现有会话
            % 输入:
            %   sessionPath - 会话文件路径
            % 输出:
            %   session - 加载的SessionInfo对象

            sessionFile = sessionPath;
            
            if ~isfile(sessionFile)
                error('SEAL:SessionInfo:FileNotFound', ...
                    'No such file exists: %s', sessionFile);
            end
            
            % 加载会话数据
            session = loadData(sessionFile);
            
            % 验证加载的对象类型
            if ~isa(session, 'SessionInfo')
                error('SEAL:SessionInfo:InvalidFileFormat', ...
                    'File does not contain a valid SessionInfo object: %s', sessionFile);
            end
        end
       
    end
end