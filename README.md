# SEAL

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

SEAL 是一个用于 MEG/EEG 源成像分析的 MATLAB 平台，提供完整的数据管理、预处理、源成像算法和可视化功能。

## 主要特性

- **节点系统**: 采用树形节点系统管理所有数据，所有节点继承自 `BaseNode` 抽象基类
- **数据与元数据分离**: 真实数据存储在文件系统中，内存中仅维护元数据和路径引用
- **按需加载**: 大数据集按需加载，避免内存溢出
- **可扩展架构**: 支持插件和算法扩展

## 功能特性

### 核心功能

- **数据管理**: 基于节点系统的数据管理，支持项目、协议、会话、数据的多级组织
- **数据导入**: 支持多种格式的 EEG/MEG 数据导入
- **预处理**: 提供完整的预处理流程，包括滤波、降采样、重参考、插值、拉普拉斯等
- **源成像**: 支持多种源成像算法，包括最小L2范数和时空算法
- **可视化**: 提供丰富的可视化功能，包括 EEG 显示、皮层显示、源成像结果可视化

### 算法库

#### 最小L2范数算法

- **MNE**: 最小范数估计
- **LORETA**: 低分辨率电磁层析成像
- **sLORETA**: 标准化LORETA
- **eLORETA**: 精确LORETA
- **LAURA**: 局部自回归平均
- **dSPM**: 动态统计参数映射

#### 时空算法

- **BlockChampagne**: Block Champagne算法
- **STARTS**: Spatio-Temporal Re-weighted Source Analysis
- **uSTAR**: Unconstrained Spatio-Temporal Re-weighted Source Analysis

### 预处理功能

- **滤波**: 高通、低通、带通、带阻滤波器
- **降采样**: 降低采样率
- **重参考**: 重新参考
- **插值**: 插值处理
- **拉普拉斯**: 拉普拉斯变换

## 快速开始

### 环境要求
- MATLAB R2022b 或更高版本
- MATLAB App Designer
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

### 安装

1. 克隆仓库到本地：

```bash
git https://github.com/qsZZF/SEAL.git
cd SEAL
```

2. 在 MATLAB 中打开项目：

```matlab
open('SEAL_MATLABToolbox.prj')
```

3. 添加到 MATLAB 路径
```matlab
addpath(genpath('SEAL'))
savepath
```

## 快速开始

### 使用 GUI

```matlab
% 启动主界面
SEAL_GUI
```

## 项目结构

```
SEAL/
├── algorithms/                     # 算法库
│   ├── MinimumL2Norm/              # 最小L2范数算法
│   └── Spatio Temporal/            # 时空算法
├── app/                           # 应用程序
│   ├── preprocess/                  # 预处理应用
│   ├── SEAL_GUI.mlapp               # 主界面
│   ├── SEAL_Preprocess.mlapp        # 预处理界面
│   ├── SEAL_showCortex.mlapp        # 皮层显示界面
│   └── SEAL_showEEG.mlapp           # EEG显示界面
├── core/                          # 核心系统
│   ├── functional/                  # 核心功能函数
│   │   ├── filter/                  # 滤波器设计
│   │   ├── sealImporter.m           # 数据导入
│   │   ├── PlotSource.m             # 源绘图
│   │   └── ...                      # 其他功能函数
│   ├── infos/                      # 数据模型
│   │   ├── ProjectInfo.m            # 项目信息
│   │   ├── ProtocolInfo.m           # 协议信息
│   │   └── ...                      # 其他Info类
│   └── nodes/                      # 数据节点系统
│       ├── BaseNode.m               # 基础节点类
│       ├── ProjectNode.m            # 项目节点
│       └── ...                      # 其他节点类
├── docs/                           # 文档
│   ├── coding_standards.md        # 编码规范
│   ├── data_dict.md               # 数据字典
│   ├── document.md                # 项目文档
│   ├── proj_structure.md          # 项目结构文档
│   ├── requirements_analysis.md    # 需求分析文档
│   └── storage_structure.md       # 存储结构文档
├── external/                       # 外部依赖
│   ├── Newtopoplot1.1/            # 拓扑图插件
│   ├── brainstorm/                 # Brainstorm 工具
│   └── eeglab/                    # EEGLAB 工具
├── sample/                         # 示例数据
│   ├── data/                       # 示例数据集
│   └── proj/                       # 示例项目
├── test/                           # 测试文件
├── utilities/                      # 工具函数
│   ├── IO/                         # 输入输出工具
│   ├── ProgressBar.m                # 进度条
│   └── ...                          # 其他工具函数
└── README.md                      # 项目说明
```

## 文档

- [项目文档](docs/document.md) - 完整的项目文档
- [项目结构](docs/proj_structure.md) - 项目结构说明
- [数据字典](docs/data_dict.md) - 数据字典文档
- [存储结构](docs/storage_structure.md) - 存储结构说明
- [需求分析](docs/requirements_analysis.md) - 需求分析文档
- [编码规范](docs/coding_standards.md) - 编码规范文档

## 开发指南

### 代码规范

- 变量和函数名：使用驼峰式命名，例如 myVariable，calculateSum
- 常量：使用大写字母和下划线，例如 MAX_SIZE
- 类名：使用帕斯卡命名法，即每个单词首字母大写，例如 MyClass，NodeData
- 文件名：通常与文件中的主要函数或类名相同，使用驼峰

详细的编码规范请参考 [编码规范](docs/coding_standards.md)。

### 设计模式

- **抽象工厂模式**: 节点创建使用静态工厂方法
- **组合模式**: 树形节点系统
- **模板方法模式**: BaseNode定义抽象方法，子类实现具体逻辑

### 扩展指南

#### 添加新节点

1. 继承 `BaseNode` 类
2. 实现抽象属性 `type` 和 `name`
3. 实现抽象方法 `load()`、`unload()`、`open()`、`save()`、`deleteFromDisk()`
4. 创建对应的 Info 类
5. 在父节点中添加子节点管理方法

#### 添加新算法

1. 在 `algorithms/` 目录下创建新算法文件
2. 实现算法逻辑
3. 在GUI中添加算法选择界面

#### 添加新预处理功能

1. 在 `core/functional/` 目录下创建新功能文件
2. 实现功能逻辑
3. 在GUI中添加预处理界面

## 测试

项目包含以下测试文件：

- **TestProjectNode.m**: 项目节点测试
- **TestProtocolNode.m**: 协议节点测试
- **TestChannelNode.m**: 通道节点测试
- **TestExternalProtocol.m**: 外部协议测试
- **TestSEALIntegration.m**: 集成测试

运行测试：

```matlab
% 运行所有测试
results = runtests('test');

% 显示测试结果
disp(results);
```

## 示例数据

项目包含示例数据，位于 `Database/` 目录：

- `Database/1/`: 示例数据集1
- `Database/2/`: 示例数据集2

## 外部依赖

### Newtopoplot1.1
拓扑图插件，提供EEG拓扑图绘制功能。

### Brainstorm
Brainstorm工具，提供部分辅助功能。

### EEGLAB
EEGLAB工具，提供ICA等功能。

## 贡献

欢迎贡献！然而写这个文档的人并不知道这是否合适，所以这里只是简单地写一下欢迎。

## 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

## 致谢

感谢以下开源项目：
- [EEGLAB](https://sccn.ucsd.edu/eeglab/)
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/)
- [Newtopoplot](https://www.mathworks.com/matlabcentral/fileexchange/17287-newtopoplot)

## 更新日志

### 2026/01/31
- 完成节点系统重构
- 建立完整的数据字典
- 更新项目结构文档

### 2025/11/01
- 增强鲁棒性
- 改进交互逻辑
- 修复大量路径和覆写异常

---

**注意**: 本项目仍在积极开发中，API 可能会发生变化。
