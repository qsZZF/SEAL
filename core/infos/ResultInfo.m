classdef ResultInfo < handle
    %RESULTDATA 通用源成像结果数据管理类
    % 不绑定任何特定算法，支持动态存储任意参数组合
    
    properties (SetAccess = public)
        % 1. 基础身份信息
        name string = "Result_Untitled"
        algorithm string = "Unknown"  % 记录算法名称，如 'MNE', 'sLORETA', 'LCMV'
        parentDataName string         % 关联的原始数据名称
               
        datapath          
        
        % 3. 动态参数池 (核心改动)
        % 初始化为空结构体，无论什么算法的参数，都作为动态字段塞进来
        params struct = struct()      
        
        % 4. 辅助信息
        desc string = ""
        metadata struct = struct()    
    end
    
    properties (SetAccess = private)
        createdDate datetime
        modifiedDate datetime
    end
    
    methods
        function obj = ResultInfo(name, algorithm, datapath, paramsStruct)
            % 构造函数
            if nargin >= 1, obj.name = name; end
            if nargin >= 2, obj.algorithm = algorithm; end
            if nargin >= 3, obj.datapath = datapath; end
            % 直接接收一个封装好的结构体
            if nargin >= 4 && isstruct(paramsStruct)
                obj.params = paramsStruct; 
            end
            
            obj.createdDate = datetime('now');
            obj.modifiedDate = obj.createdDate;
        end
        
        %% 动态参数管理方法
        function setParameter(obj, paramName, paramValue)
            % 允许在对象创建后，动态添加或修改单个参数
            obj.params.(paramName) = paramValue;
            obj.updateModifiedDate();
        end
        
        function val = getParameter(obj, paramName)
            % 安全地获取参数，如果不存在则返回空
            if isfield(obj.params, paramName)
                val = obj.params.(paramName);
            else
                warning('Parameter "%s" does not exist in this result.', paramName);
                val = [];
            end
        end
        
        %% 数据保存
        function save(obj, savePath)
            obj.updateModifiedDate();
            resultObj = obj; 
            save(savePath, 'resultObj');
        end
    end
    
    methods (Access = private)
        function updateModifiedDate(obj)
            obj.modifiedDate = datetime('now');
        end
    end

    methods(Static)
        function data = fromData(datapath, infoStruct)
            %FROMDATA 从源计算结果直接创建 ResultInfo
            % 输入:
            %   datapath   - 计算出的源激活矩阵 (numeric) 或 保存路径 (string)
            %   infoStruct - 包含 'name', 'algorithm' 及算法参数的结构体
            
            % 1. 验证输入
            if ~isstruct(infoStruct)
                error('SEAL:ResultInfo:InvalidInput', '第二个参数必须是结构体。');
            end
            
            % 2. 提取必填基础信息，若未提供则赋予保底值
            if isfield(infoStruct, 'name')
                resName = infoStruct.name;
            else
                resName = "Untitled_Result";
            end
            
            if isfield(infoStruct, 'Algorithm')
                algo = infoStruct.Algorithm;
            else
                algo = "Unknown";
            end
            
            % 3. 调用基础构造函数创建对象
            try
                % 将整个 infoStruct 直接存入 params 中，保留所有算法的独特参数
                data = ResultInfo(resName, algo, datapath, infoStruct);
                
                % 提取并设置辅助属性（如果结构体中包含这些额外信息）
                if isfield(infoStruct, 'desc')
                    data.desc = infoStruct.desc;
                end
                if isfield(infoStruct, 'parentDataName')
                    data.parentDataName = infoStruct.parentDataName;
             
                end
                
            catch ME
                error('SEAL:ResultInfo:CreationFailed', ...
                    'Failed to create ResultInfo: %s', ME.message);
            end
        end

        
        function data = openExisting(filePath)
            %OPENEXISTING 从硬盘打开已保存的 ResultInfo .mat 文件
            
            if ~isfile(filePath)
                error('SEAL:ResultInfo:FileNotFound', ...
                    'No such file exists: %s', filePath);
            end
            
            % 加载 .mat 文件
            tmp = load(filePath);
            
            % 验证加载的对象类型 (兼容直接保存对象，或封装在 resultObj 变量中的情况)
            if isfield(tmp, 'resultObj') && isa(tmp.resultObj, 'ResultInfo')
                data = tmp.resultObj;
            elseif isa(tmp, 'ResultInfo')
                data = tmp;
            else
                error('SEAL:ResultInfo:InvalidFileFormat', ...
                    'File does not contain a valid ResultInfo object: %s', filePath);
            end
        end

    end

end