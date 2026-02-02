# SEAL 编码规范文档

## 目录

- [1. 概述](#1-概述)
- [2. 命名规范](#2-命名规范)
- [3. 文件组织](#3-文件组织)
- [4. 代码格式](#4-代码格式)
- [5. 注释规范](#5-注释规范)
- [6. MATLAB 特定规范](#6-matlab-特定规范)

---

## 1. 概述

本文档定义了 SEAL 项目的编码规范，旨在提高代码的可读性、可维护性和一致性。所有贡献者都应遵循这些规范。

### 原则

1. **一致性**: 代码风格应保持一致
2. **可读性**: 代码应易于理解
3. **可维护性**: 代码应易于修改和扩展
4. **文档化**: 代码应有充分的注释和文档
5. **测试**: 代码应有相应的测试

---

## 2. 命名规范

### 2.1 变量命名

#### 基本规则

- 使用驼峰式命名（camelCase）
- 变量名应具有描述性
- 避免使用单字母变量（循环变量除外）
- 避免使用缩写（除非是广泛认知的缩写）

#### 示例

```matlab
% 好的命名
userName = 'John';
dataFilePath = 'C:\data\file.mat';
samplingRate = 1000;
numberOfChannels = 64;

% 不好的命名
u = 'John';
path = 'C:\data\file.mat';
sr = 1000;
n = 64;
```

#### 布尔变量

布尔变量应以 `is`、`has`、`should`、`can` 等前缀开头：

```matlab
% 好的命名
isValid = true;
hasData = false;
shouldSave = true;
canLoad = false;

% 不好的命名
valid = true;
data = false;
save = true;
```

#### 常量

常量应使用全大写字母和下划线：

```matlab
% 好的命名

MAX_NUMBER_OF_CHANNELS = 64;
DEFAULT_SAMPLING_RATE = 1000;
MINIMUM_FREQUENCY = 0.1;
MAXIMUM_FREQUENCY = 100;

% 不好的命名
maxChannels = 64;
defaultSamplingRate = 1000;
```

### 2.2 函数命名

#### 基本规则

- 使用驼峰式命名（camelCase）
- 函数名应以动词开头
- 函数名应描述函数的功能
- 避免使用缩写（除非是广泛认知的缩写）

#### 示例

```matlab
% 好的命名
function data = loadData(filePath)
function result = calculateMean(inputArray)
function isValid = validateInput(inputValue)
function saveToFile(data, filePath)

% 不好的命名
function d = ld(fp)
function r = cm(ia)
function v = vi(iv)
function s = sf(d, fp)
```

#### 私有函数

私有函数应以下划线开头：

```matlab
% 私有函数
function result = _calculateInternalValue(input)
function isValid = _validateInternalState(state)
```

### 2.3 类命名

#### 基本规则

- 使用帕斯卡命名法（PascalCase）
- 类名应为名词
- 类名应描述类的功能

#### 示例

```matlab
% 好的命名
classdef ProjectNode
classdef DataProcessor
classdef SignalAnalyzer
classdef VisualizationManager

% 不好的命名
classdef projectnode
classdef data_processor
classdef Signal
classdef Manager
```

#### 抽象类

抽象类应以 `Abstract` 前缀开头：

```matlab
% 抽象类
classdef (Abstract) BaseNode
classdef (Abstract) AbstractProcessor
```

### 2.4 文件命名

#### 基本规则

- 文件名应与主函数或类名相同
- 使用驼峰式命名（camelCase）
- 文件扩展名：`.m` 用于 MATLAB 文件，`.mlapp` 用于 App Designer 应用

#### 示例

```matlab
% 函数文件
loadData.m
calculateMean.m
validateInput.m

% 类文件
ProjectNode.m
DataProcessor.m
SignalAnalyzer.m

% App Designer 应用
SEAL_GUI.mlapp
SEAL_Preprocess.mlapp
```

#### 测试文件

测试文件应以 `Test` 前缀开头：

```matlab
% 测试文件
TestProjectNode.m
TestDataProcessor.m
TestSignalAnalyzer.m
```

### 2.5 目录命名

#### 基本规则

- 使用驼峰式命名（camelCase）
- 目录名应具有描述性
- 避免使用空格和特殊字符

#### 示例

```
% 好的目录命名
Algorithms/
CoreFunctions/
DataProcessing/
Visualization/

% 不好的目录命名
algorithms/
core_functions/
data processing/
Visualization/
```

---

## 3. 文件组织

### 3.1 文件结构

每个 MATLAB 文件应按以下顺序组织：

1. 文件头注释
2. 函数/类定义
3. 属性定义
4. 方法定义
5. 私有方法定义
6. 静态方法定义

#### 示例

```matlab
% 文件头注释
%LOADDATA Load data from file
%   data = loadData(filePath) loads data from the specified file.
%
%   Inputs:
%       filePath - Path to the data file (string)
%
%   Outputs:
%       data - Loaded data (struct)
%
%   Example:
%       data = loadData('C:\data\file.mat');
%
%   See also: saveData, validateInput

function data = loadData(filePath)
    % 函数实现
end
```

### 3.2 类文件结构

类文件应按以下顺序组织：

1. 文件头注释
2. 类定义
3. 属性定义
   - 常量属性
   - 公共属性
   - 私有属性
   - 依赖属性
4. 事件定义
5. 方法定义
   - 构造函数
   - 公共方法
   - 私有方法
   - 静态方法

#### 示例

```matlab
% 文件头注释
%PROJECTNODE Project node class
%   ProjectNode represents a project in the SEAL system.
%
%   Properties:
%       name - Project name (string)
%       path - Project path (string)
%
%   Methods:
%       createNew - Create a new project
%       openExisting - Open an existing project
%       save - Save the project
%
%   See also: ProtocolNode, SessionNode

classdef ProjectNode < handle
    % 常量属性
    properties (Constant)
        type = "ProjectNode"
    end
    
    % 公共属性
    properties
        name string = "Untitled"
        path string = ""
    end
    
    % 私有属性
    properties (Access = private)
        createdDate datetime
        modifiedDate datetime
    end
    
    % 依赖属性
    properties (Dependent)
        isLoaded logical
    end
    
    % 事件
    events
        NodeModified
        NodeSaved
    end
    
    % 构造函数
    methods
        function obj = ProjectNode(name, path)
            % 构造函数实现
        end
    end
    
    % 公共方法
    methods
        function save(obj)
            % 保存方法实现
        end
    end
    
    % 私有方法
    methods (Access = private)
        function updateModifiedDate(obj)
            % 更新修改日期
        end
    end
    
    % 静态方法
    methods (Static)
        function project = createNew(name, path)
            % 创建新项目
        end
    end
end
```

### 3.3 目录结构

项目应按以下目录结构组织：

```
SEAL/
├── Algorithms/              # 算法库
│   ├── MinimumL2Norm/
│   └── SpatioTemporal/
├── CoreFunctions/           # 核心功能
│   ├── DataProcessing/
│   └── Utilities/
├── GUI/                    # 图形界面
├── core/                   # 核心系统
│   └── nodes/
├── models/                 # 数据模型
├── utilities/              # 工具函数
│   ├── IO/
│   └── Validators/
├── test/                   # 测试文件
└── docs/                   # 文档
```

---

## 4. 代码格式

### 4.1 缩进

- 使用 4 个空格缩进
- 不要使用 Tab 键

### 4.2 行长度

- 每行代码最多 120 个字符
- 如果超过 120 个字符，应换行

#### 示例

```matlab
% 好的格式
result = calculateMean(dataArray, weightingFactor, normalizationMethod, precision);

% 换行格式
result = calculateMean(dataArray, weightingFactor, ...
                    normalizationMethod, precision);
```

### 4.3 空行

- 函数之间使用 2 个空行
- 逻辑块之间使用 1 个空行
- 类定义内部的不同部分之间使用 1 个空行

#### 示例

```matlab
function result = calculateMean(inputArray)
    % 输入验证
    if isempty(inputArray)
        error('Input array cannot be empty');
    end
    
    % 计算均值
    result = sum(inputArray) / length(inputArray);
end


function result = calculateStandardDeviation(inputArray)
    % 输入验证
    if isempty(inputArray)
        error('Input array cannot be empty');
    end
    
    % 计算标准差
    meanValue = calculateMean(inputArray);
    result = sqrt(sum((inputArray - meanValue).^2) / length(inputArray));
end
```

### 4.4 括号

- 左括号后不换行
- 右括号前不换行
- 长表达式应换行

#### 示例

```matlab
% 好的格式
if isValid
    result = calculateValue(input);
end

% 长表达式
result = calculateValue(input1, input2, input3, ...
                      input4, input5, input6);
```

### 4.5 运算符

- 运算符前后应加空格
- 逗号后应加空格
- 分号后应加空格

#### 示例

```matlab
% 好的格式
result = a + b * c;
data = [1, 2, 3, 4, 5];

% 不好的格式
result = a+b*c;
data = [1,2,3,4,5];
```

---

## 5. 注释规范

### 5.1 文件头注释

每个文件都应有文件头注释，包含以下信息：

- 文件功能描述
- 输入参数说明
- 输出参数说明
- 示例
- 相关函数
- 作者
- 日期

#### 示例

```matlab
%LOADDATA Load data from file
%   data = loadData(filePath) loads data from the specified file.
%
%   Inputs:
%       filePath - Path to the data file (string)
%
%   Outputs:
%       data - Loaded data (struct)
%           .values - Data values (array)
%           .timestamps - Timestamps (datetime array)
%
%   Example:
%       data = loadData('C:\data\file.mat');
%       disp(data.values);
%
%   See also: saveData, validateInput
%
%   Author: Your Name
%   Date: 2024-01-01

function data = loadData(filePath)
    % 函数实现
end
```

### 5.2 函数注释

每个函数都应有注释，说明函数的功能、输入和输出。

#### 示例

```matlab
function result = calculateMean(inputArray)
    %CALCULATEMEAN Calculate the mean of an array
    %   result = calculateMean(inputArray) calculates the mean of the input array.
    %
    %   Inputs:
    %       inputArray - Input array (numeric array)
    %
    %   Outputs:
    %       result - Mean value (scalar)
    
    % 输入验证
    if isempty(inputArray)
        error('Input array cannot be empty');
    end
    
    % 计算均值
    result = sum(inputArray) / length(inputArray);
end
```

### 5.3 类注释

每个类都应有注释，说明类的功能、属性和方法。

#### 示例

```matlab
classdef ProjectNode < handle
    %PROJECTNODE Project node class
    %   ProjectNode represents a project in the SEAL system.
    %
    %   Properties:
    %       name - Project name (string)
    %       path - Project path (string)
    %       createdDate - Creation date (datetime)
    %       modifiedDate - Modification date (datetime)
    %
    %   Methods:
    %       createNew - Create a new project
    %       openExisting - Open an existing project
    %       save - Save the project
    %       load - Load the project
    %
    %   See also: ProtocolNode, SessionNode
end
```

### 5.4 行内注释

行内注释应解释代码的目的，而不是重复代码。

#### 示例

```matlab
% 好的注释
% 计算加权均值
result = sum(dataArray .* weights) / sum(weights);

% 不好的注释
% 乘以权重并求和，然后除以权重和
result = sum(dataArray .* weights) / sum(weights);
```

### 5.5 TODO 注释

使用 TODO 注释标记待完成的任务：

```matlab
% TODO: 实现数据验证
% TODO: 优化性能
% TODO: 添加错误处理
```

---

## 6. MATLAB 特定规范

### 6.1 输入验证

所有函数都应验证输入参数：

```matlab
function result = calculateMean(inputArray)
    % 输入验证
    if nargin < 1
        error('Not enough input arguments');
    end
    
    if isempty(inputArray)
        error('Input array cannot be empty');
    end
    
    if ~isnumeric(inputArray)
        error('Input must be numeric');
    end
    
    % 计算均值
    result = sum(inputArray) / length(inputArray);
end
```

### 6.2 输出验证

函数应验证输出参数：

```matlab
function [meanValue, stdValue] = calculateStatistics(inputArray)
    % 计算统计值
    meanValue = mean(inputArray);
    stdValue = std(inputArray);
    
    % 输出验证
    if nargout > 2
        error('Too many output arguments');
    end
end
```

### 6.3 使用 inputParser

对于复杂函数，使用 `inputParser`：

```matlab
function result = processData(varargin)
    % 参数解析
    p = inputParser;
    addRequired(p, 'data', @isnumeric);
    addParameter(p, 'Method', 'mean', @(x) ismember(x, {'mean', 'median', 'mode'}));
    addParameter(p, 'Normalize', false, @islogical);
    
    parse(p, varargin{:});
    
    % 获取参数
    data = p.Results.data;
    method = p.Results.Method;
    normalize = p.Results.Normalize;
    
    % 处理数据
    if normalize
        data = (data - mean(data)) / std(data);
    end
    
    switch method
        case 'mean'
            result = mean(data);
        case 'median'
            result = median(data);
        case 'mode'
            result = mode(data);
    end
end
```

### 6.4 避免全局变量

避免使用全局变量，使用类属性或函数参数代替。

### 6.5 使用持久变量

如果需要在函数调用之间保持状态，使用持久变量：

```matlab
function result = getCounter()
    persistent counter
    
    if isempty(counter)
        counter = 0;
    end
    
    counter = counter + 1;
    result = counter;
end
```
