classdef ProjectInfo < handle
    %PROJECTINFO 项目管理信息类
    % 负责项目的元数据管理、文件路径管理

    properties (SetAccess = public)
        % 基本信息
        name string = "Untitled"
        desc string = ""
        
        % 元数据
        metadata struct % 用户自定义元数据
    end
    
    properties (SetAccess = private)
        % 只读信息
        createdDate datetime
        modifiedDate datetime
    end
    
    methods
        function obj = ProjectInfo(name, desc)
            %PROJECTINFO 构造函数
            % 输入:
            %   name - 项目名称
            %   path - 项目保存路径
            %   desc - 项目描述

            if nargin >= 1
                obj.name = name;
            end
            if nargin >= 2
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
            %SAVE 保存项目信息
            % 输入:
            %   savePath - 保存路径（可选）
            
            if path == ""
                error('SEAL:ProjectInfo:InvalidPath', ...
                    'No valid project path specified.');
            end
            
            if ~isfolder(path)
                mkdir(path);
            end
            
            % 保存项目元数据
            obj.updateModifiedDate();
            saveData(fullfile(path, strcat(obj.name, '.mat')), obj);
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
        function project = createNew(projectName, desc)
            %CREATENEW 创建新项目
            % 输入:
            %   projectName - 项目名称
            %   savePath - 保存路径
            % 输出:
            %   project - 新创建的ProjectInfo对象
            project = ProjectInfo(projectName, desc);
        end
        
        function project = openExisting(projectPath)
            %OPENEXISTING 打开现有项目
            % 输入:
            %   projectPath - 项目路径
            % 输出:
            %   project - 加载的ProjectInfo对象

            projectFile = projectPath;
            
            if ~isfile(projectFile)
                error('SEAL:ProjectInfo:FileNotFound', ...
                    'No such file exists: %s', projectFile);
            end
            
            % 加载项目数据
            project = loadData(projectFile);
        end
    end
end