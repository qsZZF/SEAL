classdef ProjectNode < BaseNode
    %PROJECTNODE 项目节点 - 整个数据树的根节点
    % 管理整个研究项目，包含多个协议
    
    properties (Constant)
        type = "ProjectNode"  % 节点类型标识符
    end
    
    properties
        % 项目特定属性
        projectInfo ProjectInfo
    end
    
    properties (Dependent)
        % 依赖属性
        name string
        protocolCount int32                 % 协议数量
        infoFile string
        protocolFolder string
    end
    
    methods
        function obj = ProjectNode()
            %PROJECTNODE 构造函数
            
            % 调用父类构造函数
            obj = obj@BaseNode([]);
        end

        function open(obj, projectPath)
            if(obj.isLoaded)
                return;
            end
            
            try
                % 检查输入路径是文件还是文件夹
                if isfile(projectPath)
                    % 如果是文件路径，提取文件夹路径
                    [projectPath, ~, ~] = fileparts(projectPath);
                end

                obj.path = projectPath;
                
                % 检查项目文件是否存在
                if ~isfile(obj.infoFile)
                    error('SEAL:ProjectNode:ProjectFileNotFound', ...
                        'No such project exists: %s', obj.infoFile);
                end
                
                % 加载项目元数据
                
                obj.projectInfo = ProjectInfo.openExisting(obj.infoFile);
                
                % 加载Protocol
                obj.openChildNodes();
                
                % 标记为已加载
                obj.isLoaded = false;
                
            catch ME
                error('SEAL:ProjectNode:LoadFailed', ...
                    'Failed to open: %s', ME.message);
            end
        end

        function save(obj)
            %SAVE 保存项目数据

            obj.createDirectoryStructure();
            obj.projectInfo.save(obj.infoFile);
            
            % 递归保存所有子节点
            for i = 1:obj.childCount
                obj.children(i).save();
            end
        end

        function load(obj)
            %LOAD 载入项目数据，对于Project而言无操作
            
            % 递归载入所有子节点
            for i = 1:obj.childCount
                obj.children(i).load();
            end
            obj.isLoaded = true;
        end

        function unload(obj)
            %UNLOAD 卸载项目数据，释放内存

            % 递归卸载所有子节点
            for i = 1:obj.childCount
                obj.children(i).unload();
            end
            obj.isLoaded = false;
        end

        function deleteFromDisk(obj)
            rmdir(obj.path, 's');
        end

        function rename(obj, newName)
%           等待进一步实现，这个功能是否要加入待定
            oldPath = obj.path;
            oldInfoFile = obj.infoFile;
            [parentPath, ~, ~] = obj.path;

            obj.projectInfo.name = newName;
            newPath = strcat(parentPath, newName);
            rename(oldInfoFile, obj.infoFile);
            rename(oldPath, newPath);
            obj.path = newPath;
        end
        
        %% 依赖属性get方法
        function count = get.protocolCount(obj)
            count = obj.childCount;
        end
        
        function name = get.name(obj)
            name = obj.projectInfo.name;
        end

        function filePath = get.infoFile(obj)
            [~, name, ~] = fileparts(obj.path);
            filePath = fullfile(obj.path, strcat(name, '.mat'));
        end

        function protocolFolder = get.protocolFolder(obj)
            protocolFolder = fullfile(obj.path, "Protocols");
        end
        
        %% 重写父类方法
        function addChild(obj, childNode)
            %ADDCHILD 重写添加子节点方法，确保只添加ProtocolNode
            
            if ~isa(childNode, 'ProtocolNode')
                error('SEAL:ProjectNode:InvalidChildType', ...
                    'ProjectNode只能添加ProtocolNode子节点');
            end
            
            % 调用父类方法
            addChild@BaseNode(obj, childNode);
        end
        
        function removeChild(obj, childNode)
            %REMOVECHILD 重写移除子节点方法
            
            % 调用父类方法
            removeChild@BaseNode(obj, childNode);
        end

        function protocol = createNewProtocol(obj, protocolName, varargin)
            protocol = ProtocolNode.createNew(protocolName, obj.protocolFolder, varargin{:});
            obj.addChild(protocol);
        end

        function protocol = openProtocol(obj, protocolFolder)
            protocol = ProtocolNode.openExisting(protocolFolder);
            % ProtocolInfo 相对轻量，且导入的Protocol不应该远程修改原始信息
            % 故此处强行设置当前project的path而非使用原始path
            % 在save时等价于自动复制并保存
            [~,protocolName,~] = fileparts(protocolFolder);
            protocol.path = fullfile(obj.protocolFolder, protocolName);
            obj.addChild(protocol);
        end

        function openChildNodes(obj)
            %OPENCHILDNODES 打开子节点
            % 这个方法会在打开协议时自动扫描并打开现有的子节点
            
            % 检查协议目录是否存在
            if ~isfolder(obj.protocolFolder)
                return;
            end
            
            % 获取所有协议目录
            protocols = dir(obj.protocolFolder);
            for i = 1:length(protocols)
                protocolDir = protocols(i);
                if protocolDir.isdir && ~startsWith(protocolDir.name, '.')
                    try
                        % 检查是否存在协议文件
                        folder = fullfile(obj.protocolFolder, protocolDir.name);
                        file = fullfile(folder, strcat(protocolDir.name, '.mat'));
                        if isfile(file)
                            % 创建并加载协议节点
                            protocolNode = ProtocolNode.openExisting(file);
                            protocolNode.path = folder;
                            obj.addChild(protocolNode);
                        end
                    catch ME
                        warning('SEAL:ProjectNode:ProtocolLoadFailed', ...
                            'Failed to open %s: %s', protocolDir.name, ME.message);
                    end
                end
            end
        end
    end

    methods (Static)
        function project = openExisting(projectPath)
            %OPENPEXISTING 打开现有项目
            % 输入:
            %   projectPath - 项目文件路径或项目文件夹路径
            
            if isfolder(projectPath)
                % 如果是文件夹，查找项目文件
                [~, projectName, ~] = fileparts(projectPath);
                projectFiles = dir(fullfile(projectPath, strcat(projectName, '.mat')));
                if isempty(projectFiles)
                    error('SEAL:ProjectNode:NoProjectFile', ...
                        'No such project in path: %s', projectPath);
                end
            else
                 [projectPath, ~, ~] = fileparts(projectPath);
            end
            
            % 创建项目节点
            project = ProjectNode();
            project.path = projectPath;
            
            % 加载项目数据
            project.open(project.path);
        end
        
        function project = createNew(projectName, projectPath, varargin)
            %CREATENEW 创建新项目
            % 输入:
            %   projectName - 项目名称
            %   projectPath - 项目存储路径
            %   varargin - 可选参数对
            
            p = inputParser;
            addRequired(p, 'projectName', Validators.isText);
            addRequired(p, 'projectPath', Validators.isText);
            addParameter(p, 'Description', '', Validators.isText);
            
            parse(p, projectName, projectPath, varargin{:});
            
            % 创建项目节点
            project = ProjectNode();
            project.path = fullfile(projectPath, projectName);
            project.projectInfo = ProjectInfo.createNew(p.Results.projectName, ...
                p.Results.Description);
            
            % 保存项目
            project.save();
            project.load();
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
            if ~isfolder(obj.protocolFolder)
                mkdir(obj.protocolFolder);
            end
        end
    end
end