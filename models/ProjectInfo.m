classdef ProjectInfo < handle
    %PROJECTINFO 项目管理信息类
    % 负责项目的元数据管理、文件路径管理
    
    properties (SetAccess = private)
        % 基本信息
        name string {mustBeNonzeroLengthText} = "Untitled"
        desc string = ""
        createdDate datetime
        modifiedDate datetime
        
        % 项目路径信息
        path string {mustBeFolder} = ""
        
        % 元数据
        metadata struct % 用户自定义元数据
    end
    
    methods
        function obj = ProjectInfo(name, path)
            %PROJECTINFO 构造函数
            % 输入:
            %   name - 项目名称
            %   path - 项目保存路径
            
            if nargin >= 1
                obj.name = name;
            end
            if nargin >= 2
                obj.path = path;
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
            
            if nargin > 1
                obj.path = path;
            end
            
            if obj.path == ""
                error('SEAL:ProjectInfo:InvalidPath', ...
                    'No valid project path specified.');
            end
            
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
            
            % 保存项目元数据
            obj.updateModifiedDate();
            save(fullfile(obj.path, [obj.name, '.mat']),'-struct', 'obj');
        end
        
        function load(obj, path)
            %LOAD 从文件加载项目信息
            % 输入:
            %   projectPath - 项目文件路径
            
            projectFile = path;
            
            if ~isfile(projectFile)
                error('SEAL:ProjectInfo:FileNotFound', ...
                    'No such file exists: %s', projectFile);
            end
            
            % 加载项目数据
            data = load(projectFile, 'obj');
            
            % 恢复属性
            obj.name = data.name;
            obj.description = data.description;
            obj.createdDate = data.createdDate;
            obj.modifiedDate = data.modifiedDate;
            obj.metadata = data.metadata;
            obj.path = path;
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
        function project = createNew(projectName, savePath)
            %CREATENEW 创建新项目
            % 输入:
            %   projectName - 项目名称
            %   savePath - 保存路径
            % 输出:
            %   project - 新创建的ProjectInfo对象
            
            project = ProjectInfo(projectName, savePath);
            
            % 创建项目目录结构
            if savePath ~= ""
                project.createDirectoryStructure();
            end
        end
        
        function project = openExisting(projectPath)
            %OPENEXISTING 打开现有项目
            % 输入:
            %   projectPath - 项目路径
            % 输出:
            %   project - 加载的ProjectInfo对象
            
            project = ProjectInfo();
            project.load(projectPath);
        end
    end
    
    methods (Access = private)
        function createDirectoryStructure(obj)
            %CREATEDIRECTORYSTRUCTURE 创建项目目录结构
            
            if obj.path == ""
                return;
            end
            
            % 创建主目录
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
            
            % 创建协议目录
            protocolsPath = fullfile(obj.path, 'Protocols');
            if ~isfolder(protocolsPath)
                mkdir(protocolsPath);
            end
        end
    end
end