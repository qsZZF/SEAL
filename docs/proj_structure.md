# SEAL 项目结构说明书

## 项目概述

SEAL 是一个用于 MEG/EEG 源成像分析的 MATLAB 平台，采用模块化架构设计，具有以下特点：

- **节点系统**: 采用树形节点系统管理所有数据，所有节点继承自 `BaseNode` 抽象基类
- **数据与元数据分离**: 真实数据存储在文件系统中，内存中仅维护元数据和路径引用
- **按需加载**: 大数据集按需加载，避免内存溢出
- **可扩展架构**: 支持插件和算法扩展

## 项目目录结构

```
SEAL/
├── algorithms/                     # 算法库
│   ├── MinimumL2Norm/              # 最小L2范数算法
│   │   ├── otsu.m
│   │   ├── seal_LAURA.m
│   │   ├── seal_LORETA.m
│   │   ├── seal_MNE.m
│   │   ├── seal_dSPM.m
│   │   ├── seal_eLORETA.m
│   │   └── seal_sLORETA.m
│   └── Spatio Temporal/            # 时空算法
│       ├── seal_BlockChampagne.m
│       ├── seal_STARTS.m
│       └── seal_uSTAR.m
├── app/                           # 应用程序
│   ├── preprocess/                  # 预处理应用
│   │   ├── Filters.mlapp
│   │   ├── SEAL_Preprocess.mlapp
│   │   ├── bandpass.mlapp
│   │   ├── bandsuppress.mlapp
│   │   ├── highpass.mlapp
│   │   └── lowpass.mlapp
│   ├── SEAL_DataCheck.mlapp          # 数据检查界面
│   ├── SEAL_GUI.mlapp               # 主界面
│   ├── SEAL_Newprotocol.mlapp        # 新建协议界面
│   ├── SEAL_selectAlgorithm.mlapp     # 算法选择界面
│   ├── SEAL_showCortex.mlapp        # 皮层显示界面
│   └── SEAL_showEEG.mlapp           # EEG显示界面
├── core/                          # 核心系统
│   ├── functional/                  # 核心功能函数
│   │   ├── filter/                  # 滤波器设计
│   │   │   ├── designFIRBandpass.m
│   │   │   ├── designFIRBandstop.m
│   │   │   ├── designFIRHighpass.m
│   │   │   └── designFIRLowpass.m
│   │   ├── CreateProcessTemplate.m
│   │   ├── PlotSource.m
│   │   ├── appendDataToMat.m
│   │   ├── applyDownsampleResample.m
│   │   ├── applyRereference.m
│   │   ├── bind3DInteraction.m
│   │   ├── bulidsource.m
│   │   ├── clamp.m
│   │   ├── doDrag.m
│   │   ├── getSmoothedVertices.m
│   │   ├── onKey.m
│   │   ├── onScroll.m
│   │   ├── parseChanList.m
│   │   ├── parseRefSpec.m
│   │   ├── runInterpolation.m
│   │   ├── runLaplacian.m
│   │   ├── sealImporter.m
│   │   ├── startDrag.m
│   │   ├── stopDrag.m
│   │   ├── tessSmooth.m
│   │   └── topoplotX.m
│   ├── infos/                      # 数据模型
│   │   ├── ChanlocsInfo.m           # 通道信息
│   │   ├── CortexInfo.m             # 皮层信息
│   │   ├── DataInfo.m               # 数据信息
│   │   ├── LeadfieldInfo.m          # 导联场信息
│   │   ├── ProjectInfo.m            # 项目信息
│   │   ├── ProtocolInfo.m           # 协议信息
│   │   └── SessionInfo.m            # 会话信息
│   └── nodes/                      # 数据节点系统
│       ├── BaseNode.m               # 基础节点类
│       ├── ChanlocsNode.m           # 通道节点
│       ├── CortexNode.m             # 皮层节点
│       ├── DataNode.m               # 数据节点
│       ├── LeadfieldNode.m          # 导联场节点
│       ├── ProjectNode.m            # 项目节点
│       ├── ProtocolNode.m           # 协议节点
│       └── SessionNode.m            # 会话节点
├── docImg/                        # 文档图片
│   ├── S.png
│   ├── alphaframe.png
│   ├── circle1.png
│   ├── circle2.png
│   ├── circle3.png
│   ├── image.png
│   ├── point.png
│   └── solidframe.png
├── docs/                           # 文档
│   ├── coding_standards.md        # 编码规范
│   ├── data_dict.md               # 数据字典
│   ├── document.md                # 项目文档
│   ├── proj_structure.md          # 项目结构文档
│   ├── requirements_analysis.md    # 需求分析文档
│   └── storage_structure.md       # 存储结构文档
├── external/                       # 外部依赖
│   ├── Newtopoplot1.1/            # 拓扑图插件
│   │   ├── EEGtopoSet.mat
│   │   ├── Multiplot.m
│   │   ├── Newtopoplot.m
│   │   ├── ReadMe.md
│   │   ├── ReadMe.txt
│   │   ├── arrow.m
│   │   ├── eegplugin_Newtopoplot.m
│   │   └── pop_Newtopoplot.m
│   ├── brainstorm/                 # Brainstorm 工具
│   │   ├── blk_diag.m
│   │   └── bst_gain_orient.m
│   └── eeglab/                    # EEGLAB 工具
│       ├── binica.m
│       ├── binica.sc
│       └── runica.m
├── sample/                         # 示例数据
│   ├── data/                       # 示例数据集
│   │   ├── 1/
│   │   │   ├── Chanlocs_60.mat
│   │   │   ├── Cortex_stroke.mat
│   │   │   ├── ERPset.mat
│   │   │   ├── LeadField_stroke.mat
│   │   │   └── PatientVisMS.mat
│   │   └── 2/
│   │       ├── sub-01_run-01_Cortex.mat
│   │       ├── sub-01_run-01_Data.mat
│   │       └── sub-01_run-01_Forward.mat
│   └── proj/                       # 示例项目
│       └── SampleProj/
│           ├── SampleProj.mat
│           └── Protocols/
│               └── Protocol 1/
│                   ├── Protocol 1.mat
│                   ├── cortex_info.mat
│                   ├── leadfield_info.mat
│                   ├── Results/
│                   │   └── mne_result.mat
│                   └── Sessions/
│                       └── Session 1/
│                           ├── ERPset.mat
│                           └── session_info.mat
├── test/                           # 测试文件
│   ├── TestChannelNode.m
│   ├── TestExternalProtocol.m
│   ├── TestProjectNode.m
│   ├── TestProtocolNode.m
│   └── TestSEALIntegration.m
├── utilities/                      # 工具函数
│   ├── IO/                         # 输入输出工具
│   │   ├── checkDir.m
│   │   ├── loadData.m
│   │   ├── rename.m
│   │   └── saveData.m
│   ├── ProgressBar.m                # 进度条
│   ├── Validators.m                # 验证器
│   ├── findAvailableName.m          # 查找可用名称
│   └── getmat.m                    # 获取MAT文件
├── .gitattributes                   # Git 属性配置
├── .gitignore                      # Git 忽略配置
├── ProgessBarDemo.m               # 进度条演示
├── README.md                      # 项目说明
└── SEAL_MATLABToolbox.prj         # MATLAB项目文件
```

## 架构层次说明

### 核心系统层

```
SEAL/core/
```

**职责**: 核心系统，包含节点系统、数据模型和功能函数

#### 节点系统

```
core/nodes/
```

- **BaseNode.m** - 所有节点的抽象基类，提供统一的节点接口和树形结构管理
- **ProjectNode.m** - 项目节点，整个数据树的根节点
- **ProtocolNode.m** - 协议节点，管理研究协议
- **SessionNode.m** - 会话节点，管理单个会话
- **DataNode.m** - 数据节点，管理数据文件
- **ChanlocsNode.m** - 通道节点，管理通道信息
- **CortexNode.m** - 皮层节点，管理皮层信息
- **LeadfieldNode.m** - 导联场节点，管理导联场信息

#### 数据模型

```
core/infos/
```

- **ProjectInfo.m** - 项目信息管理类
- **ProtocolInfo.m** - 协议信息管理类
- **SessionInfo.m** - 会话信息管理类
- **DataInfo.m** - 数据信息管理类
- **ChanlocsInfo.m** - 通道信息管理类
- **CortexInfo.m** - 皮层信息管理类
- **LeadfieldInfo.m** - 导联场信息管理类

#### 核心功能

```
core/functional/
```

- **filter/** - 滤波器设计
- **sealImporter.m** - 数据导入
- **PlotSource.m** - 源绘图
- **tessSmooth.m** - 网格平滑
- **applyDownsampleResample.m** - 降采样和重采样
- **applyRereference.m** - 重参考
- **runInterpolation.m** - 插值
- **runLaplacian.m** - 拉普拉斯变换

### 应用层

```
SEAL/app/
```

**职责**: 用户界面和交互

- **SEAL_GUI.mlapp** - 主界面
- **SEAL_PrePreprocess.mlapp** - 预处理界面
- **SEAL_showCortex.mlapp** - 皮层显示界面
- **SEAL_showEEG.mlapp** - EEG显示界面
- **SEAL_Newprotocol.mlapp** - 新建协议界面
- **SEAL_selectAlgorithm.mlapp** - 算法选择界面
- **SEAL_DataCheck.mlapp** - 数据检查界面
- **preprocess/** - 预处理子应用

### 算法层

```
SEAL/algorithms/
```

**职责**: 纯算法实现，无状态计算

- **MinimumL2Norm/** - 最小L2范数算法
  - MNE, LORETA, sLORETA, eLORETA, LAURA, dSPM
- **Spatio Temporal/** - 时空算法
  - BlockChampagne, STARTS, uSTAR

### 工具层

```
SEAL/utilities/
```

**职责**: 工具函数

- **IO/** - 输入输出工具
- **Validators.m** - 验证器
- **ProgressBar.m** - 进度条
- **findAvailableName.m** - 查找可用名称
- **getmat.m** - 获取MAT文件

### 外部依赖层

```
SEAL/external/
```

**职责**: 第三方依赖代码

- **Newtopoplot1.1/** - 拓扑图插件
- **brainstorm/** - Brainstorm工具
- **eeglab/** - EEGLAB工具

### 测试层

```
SEAL/test/
```

**职责**: 测试文件

- **TestProjectNode.m** - 项目节点测试
- **TestProtocolNode.m** - 协议节点测试
- **TestChannelNode.m** - 通道节点测试
- **TestExternalProtocol.m** - 外部协议测试
- **TestSEALIntegration.m** - 集成测试

### 示例层

```
SEAL/sample/
```

**职责**: 示例数据和项目

- **data/** - 示例数据集
- **proj/** - 示例项目

### 文档层

```
SEAL/docs/
```

**职责**: 项目文档

- **coding_standards.md** - 编码规范
- **data_dict.md** - 数据字典
- **document.md** - 项目文档
- **proj_structure.md** - 项目结构文档
- **requirements_analysis.md** - 需求分析文档
- **storage_structure.md** - 存储结构文档

## 数据流架构

```
用户操作 → app → core/nodes → algorithms → 数据存储
↑          ↓         ↓         ↓         ↓
可视化 ←── core/functional ←── 结果数据 ←── 处理流程 ←── 配置
```

## 设计模式

- **抽象工厂模式**: 节点创建使用静态工厂方法
- **组合模式**: 树形节点系统
- **模板方法模式**: BaseNode定义抽象方法，子类实现具体逻辑

## 数据存储设计

### 文件系统结构

```
ProjectName/
├── ProjectName.mat              # 项目信息文件
└── Protocols/                   # 部议目录
    └── ProtocolName/            # 部议目录
        ├── ProtocolName.mat      # 部议信息文件
        ├── chanlocs_info.mat    # 通道信息文件
        ├── cortex_info.mat      # 皮层信息文件
        ├── leadfield_info.mat   # 导联场信息文件
        └── Sessions/           # 会话目录
            └── SessionName/     # 会话目录
                ├── session_info.mat  # 会话信息文件
                └── DataName.mat     # 数据文件
```

## 节点系统设计

### 节点类型

1. **ProjectNode**: 项目节点，根节点
2. **ProtocolNode**: 部议节点，包含通道、皮层、导联场和会话
3. **SessionNode**: 会话节点，包含数据
4. **DataNode**: 数据节点，管理数据文件
5. **ChanlocsNode**: 通道节点，管理通道信息
6. **CortexNode**: 皮层节点，管理皮层信息
7. **LeadfieldNode**: 导联场节点，管理导联场信息

### 节点生命周期

1. **创建创建新节点
2. **打开**: 使用 `openExisting()` 打开现有节点
3. **加载**: 调用 `load()` 方法加载数据到内存
4. **保存**: 调用 `save()` 方法保存数据到磁盘
5. **卸载**: 调用 `unload()` 方法释放内存
6. **删除**: 调用 `deleteFromDisk()` 方法从磁盘删除

### 节点特性

- **树形结构**: 支持父子关系管理
- **按需加载**: 大数据集按需加载，避免内存溢出
- **事件通知**: 支持节点事件通知
- **双向绑定**: 父子节点双向绑定
- **元数据管理**: 支持元数据管理
