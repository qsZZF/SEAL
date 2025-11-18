classdef (Abstract) BaseNode < handle & matlab.mixin.Heterogeneous
    %BASENODE 所有数据节点的抽象基类
    % 提供统一的节点接口、树形结构管理和事件通知机制
    % 注意，抽象基类永远不该被实例化！不要直接用这个类！
    % 如果你真的需要一个无类型的基类Node，请继承再实例化。
    
    properties (Abstract, Constant)
        % 节点类型标识符，由子类实现
        type
    end
    
    properties (SetAccess = public)
        path string
       
        % 树形结构关系
        parent BaseNode
        children BaseNode
        
        % 节点标识
        uuid string
        tags string
        isSelected logical
        isLoaded logical
        
        % 元数据
        metadata struct
        nodeInfo struct
    end
    
    properties (Dependent, SetAccess = public)
        % 依赖属性
        fullPath string
        isRoot logical
        isLeaf logical
        depth int16
        childCount double
        hasChildren logical
    end

    properties (Abstract, Dependent, SetAccess = public)
        % 依赖属性
        name string
    end
    
    events
        % 节点事件
        NodeAdded
        NodeRemoved
        NodeModified
        NodeSelected
        NodeLoaded
        NodeSaved
    end
    
    methods (Abstract, Access = protected)
        % 抽象方法 - 由子类实现
        % validateNode(obj)
    end

    methods (Abstract, Access = public)
        load(obj)
        unload(obj)
        open(obj, path)
        save(obj)
    end
    
    methods
        function obj = BaseNode(parentNode)
            %BASENODE 构造函数
            % 输入:
            %   nodeName - 节点名称
            %   parentNode - 父节点
            
            if nargin >= 2
                obj.parent = parentNode;
            end
            
            % 初始化基本属性
            obj.uuid = '123456';
            obj.children = BaseNode.empty;
            obj.tags = string.empty;
        end
        
        %% 依赖属性的get方法
        function fullPath = get.fullPath(obj)
            %GET.FULLPATH 获取完整路径
            if obj.isRoot
                fullPath = obj.name;
            else
                parentPath = obj.parent.fullPath;
                fullPath = strcat(parentPath, '/', obj.name);
            end
            fullPath = fullfile(obj.path, fullPath);
        end
        
        function root = get.isRoot(obj)
            %GET.ISROOT 判断是否为根节点
            root = isempty(obj.parent);
        end
        
        function leaf = get.isLeaf(obj)
            %GET.ISLEAF 判断是否为叶节点
            leaf = isempty(obj.children);
        end
        
        function depth = get.depth(obj)
            %GET.DEPTH 获取节点深度
            if obj.isRoot
                depth = 0;
            else
                depth = obj.parent.depth + 1;
            end
        end
        
        function count = get.childCount(obj)
            %GET.CHILDCOUNT 获取子节点数量
            count = length(obj.children);
        end
        
        function has = get.hasChildren(obj)
            %GET.HASCHILDREN 判断是否有子节点
            has = obj.childCount > 0;
        end
        
        %% 树形结构管理方法
        function addChild(obj, childNode)
            %ADDCHILD 添加子节点
            % 输入:
            %   childNode - 要添加的子节点
            
            validateattributes(childNode, {'BaseNode'}, {'scalar'});
            
            % 检查名称唯一性
            if any([obj.children.name] == childNode.name)
                error('SEAL:BaseNode:DuplicateChild', ...
                    'Child node named "%s" exists', childNode.name);
            end
            
            % 设置父子关系
            childNode.parent = obj;
            
            % 添加到子节点列表
            if isempty(obj.children)
                obj.children = childNode;
            else
                obj.children(end+1) = childNode;
            end
            
            % 触发事件
            obj.notify('NodeAdded', NodeEventData(childNode));
            obj.updateModifiedDate();
        end
        
        function removeChild(obj, childNode)
            %REMOVECHILD 移除子节点
            % 输入:
            %   childNode - 要移除的子节点或节点名称
            
            if ischar(childNode) || isstring(childNode)
                childNode = obj.getChild(childNode);
            end
            
            validateattributes(childNode, {'BaseNode'}, {'scalar'});
            
            % 检查是否为直接子节点
            if ~any(obj.children == childNode)
                error('SEAL:BaseNode:NotAChild', ...
                    '"%s" is not a direct descendant of "%s"', ...
                    childNode.name, obj.name);
            end
            
            % 移除父子关系
            childNode.parent = BaseNode.empty;
            obj.children(obj.children == childNode) = [];
            
            % 触发事件
            obj.notify('NodeRemoved', NodeEventData(childNode));
            obj.updateModifiedDate();
        end
        
        function child = getChild(obj, childName, recursive)
            %GETCHILD 获取子节点
            % 输入:
            %   childName - 子节点名称或索引
            %   recursive - 是否递归搜索，默认为false
            % 输出:
            %   child - 找到的子节点
            
            if nargin < 3
                recursive = false;
            end
            
            if isnumeric(childName)
                % 按索引获取
                if childName > 0 && childName <= obj.childCount
                    child = obj.children(childName);
                else
                    error('SEAL:BaseNode:IndexError', ...
                        'Index %d out of bound.', childName);
                end
            else
                % 按名称获取
                childNames = [obj.children.name];
                childIndex = find(childNames == childName, 1);
                
                if isempty(childIndex) && recursive
                    % 递归搜索
                    for i = 1:obj.childCount
                        child = obj.children(i).getChild(childName, true);
                        if ~isempty(child)
                            return;
                        end
                    end
                    child = BaseNode.empty;
                elseif ~isempty(childIndex)
                    child = obj.children(childIndex);
                else
                    child = BaseNode.empty;
                end
            end
        end
        
        function children = findChildren(obj, predicate, varargin)
            %FINDCHILDREN 查找满足条件的子节点
            % 输入:
            %   predicate - 判断函数或属性名
            %   varargin - 属性值（当predicate为属性名时）
            % 输出:
            %   children - 满足条件的子节点数组
            
            if ischar(predicate) || isstring(predicate)
                % 按属性值查找
                propName = predicate;
                if nargin < 3
                    error('SEAL:BaseNode:MissingValue', ...
                        'You need to specify a value to find!');
                end
                targetValue = varargin{1};
                
                predicate = @(node) isprop(node, propName) && ...
                    isequal(node.(propName), targetValue);
            end
            
            validateattributes(predicate, {'function_handle'}, {});
            
            children = BaseNode.empty;
            for i = 1:obj.childCount
                child = obj.children(i);
                if predicate(child)
                    children(end+1) = child; %#ok<AGROW> 强迫症看着黄线难受
                end
                
                % 递归搜索
                if child.hasChildren
                    grandchildren = child.findChildren(predicate);
                    children = [children, grandchildren]; %#ok<AGROW> 所以我把它给压制了
                end
            end
        end
        
        %% 节点状态管理
        function select(obj)
            %SELECT 选择节点
            obj.isSelected = true;
            obj.notify('NodeSelected', NodeEventData(obj));
        end
        
        function deselect(obj)
            %DESELECT 取消选择节点
            obj.isSelected = false;
        end
        
        function delete(obj)
            %DELETE 析构函数
            obj.unload();
            
            % 清理父子关系
            if ~obj.isRoot
                obj.parent.removeChild(obj);
            end
            
            % 删除所有子节点
            for i = obj.childCount:-1:1
                delete(obj.children(i));
            end
        end
        
        %% 工具方法
        function printTree(obj, level)
            %PRINTTREE 打印节点树结构
            % 输入:
            %   level - 当前层级（内部使用）
            
            if nargin < 2
                level = 0;
            end
            
            % 缩进
            indent = repmat('  ', 1, level);
            
            % 节点信息
            nodeStr = sprintf('%s├─ %s (%s)', indent, obj.name, obj.type);
            if obj.isSelected
                nodeStr = strcat(nodeStr, ' [SELECTED]');
            end
            if obj.isLoaded
                nodeStr = strcat(nodeStr, ' [LOADED]');
            end
            
            disp(nodeStr);
            
            % 打印子节点
            for i = 1:obj.childCount
                obj.children(i).printTree(level + 1);
            end
        end
    end
    
    methods (Access = protected)
        function updateModifiedDate(obj)
            %UPDATEMODIFIEDDATE 更新修改时间
            obj.modifiedDate = datetime('now');
        end
    end
end