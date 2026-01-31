# SEAL - Source Estimation and Analysis Library

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

SEAL 是一个用于 MEG/EEG 源成像分析的 MATLAB 平台，提供完整的数据管理、预处理、源成像算法和可视化功能。

## 主要特性

- **节点系统**: 采用树形节点系统管理所有数据，支持父子关系和按需加载
- **数据与元数据分离**: 真实数据存储在文件系统中，内存中仅维护元数据和路径引用
- **按需加载**: 大数据集按需加载，避免内存溢出
- **可扩展架构**: 支持插件和算法扩展
- **图形化界面**: 基于 MATLAB App Designer 的用户友好界面
- **多种源成像算法**: 集成多种最小L2范数和时空算法

## 功能模块

### 数据管理
- 项目、协议、会话、数据的树形结构管理
- 通道、皮层、导联场信息管理
- 数据导入和导出

### 预处理
- 时域滤波（高通、低通、带通、带阻）
- 降采样和重参考
- 空间滤波（Laplacian、ICA）
- 插值处理

### 源成像算法

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

### 可视化
- EEG信号显示（时域、频域、拓扑图）
- 3D皮层模型显示
- 源估计结果可视化

## 安装

### 环境要求
- MATLAB R2022b 或更高版本
- MATLAB App Designer
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

### 安装步骤

1. 克隆仓库
```bash
git https://github.com/qsZZF/SEAL.git
cd SEAL
```

2. 在 MATLAB 中打开项目
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
├── Algorithms/              # 算法库
├── CoreFunctions/           # 核心功能函数
├── Database/               # 示例数据库
├── External/               # 外部依赖
├── GUI/                   # 图形用户界面
├── core/                  # 核心系统
│   └── nodes/             # 数据节点系统
├── models/                # 数据模型
├── utilities/              # 工具函数
├── test/                  # 测试文件
└── docs/                  # 文档
```

## 文档

- [数据字典](docs/data_dict.md) - 详细的数据结构说明
- [项目结构](docs/proj_structure.md) - 项目架构和目录结构
- [项目文档](docs/document.md) - 完整的项目文档
- [需求分析](docs/requirements_analysis.md) - 功能和需求说明
- [存储结构](docs/storage_structure.md) - 数据存储结构说明

## 开发指南

### 代码规范

- 变量和函数名：使用驼峰式命名，例如 `myVariable`，`calculateSum`
- 常量：使用大写字母和下划线，例如 `MAX_SIZE`
- 类名：使用帕斯卡命名法，即每个单词首字母大写，例如 `MyClass`，`NodeData`
- 文件名：通常与文件中的主要函数或类名相同，使用驼峰

### 添加新算法

1. 在 `Algorithms/` 目录下创建新算法文件
2. 实现算法逻辑
3. 在GUI中添加算法选择界面

### 添加新节点

1. 继承 `BaseNode` 类
2. 实现抽象属性 `type` 和 `name`
3. 实现抽象方法 `load()`、`unload()`、`open()`、`save()`、`deleteFromDisk()`
4. 创建对应的 Info 类
5. 在父节点中添加子节点管理方法

## 测试

运行项目测试：

```matlab
% 运行所有测试
runtests('test')
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
