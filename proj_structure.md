# SEAL 项目结构说明书

## 项目概述


- **数据与元数据分离**: 真实数据存储在文件系统中，内存中仅维护元数据和路径引用
- **按需加载**: 大数据集按需加载，避免内存溢出
- **处理管道**: 借鉴Brainstorm的处理流水线思想
- **可扩展架构**: 支持插件和算法扩展

## 计划项目目录结构

```
SEAL/
├── app/                           # 应用入口和主界面
│   ├── SEAL_App.mlapp             # 主应用程序
│   ├── components/...             # 可复用界面组件
│   └── layouts/...                # 界面布局管理
├── core/                          # 核心系统
│   ├── mediator/                  # 中介者核心
│   │   ├── SealMediator.m         # 中央协调器
│   │   ├── EventManager.m         # 事件管理器
│   │   └── CommandDispatcher.m    # 命令分发器
│   ├── nodes/                     # 数据节点系统
│   │   ├── NodeManager.m          # 节点管理器
│   │   ├── BaseNode.m             # 基础节点类
│   │   ├── ProjectNode.m          # 项目节点
│   │   ├── ProtocolNode.m         # 协议节点
│   │   ├── ChannelNode.m          # 通道节点
│   │   ├── CortexNode.m           # 皮层节点
│   │   ├── LeadfieldNode.m        # 导联场节点
│   │   ├── SessionNode.m          # 会话节点
│   │   ├── DataNode.m             # 数据节点
│   │   ├── ResultNode.m           # 结果节点
│   │   └── NodeFactory.m          # 节点工厂
│   └── workspace/                 # 工作区管理（该部分我们得讨论一下）
│       ├── Workspace.m            # 工作区上下文
│       ├── SelectionManager.m     # 选择管理器
│       └── HistoryManager.m       # 操作历史
├── commands/                      # 命令系统
│   ├── BaseCommand.m              # 基础命令类
│   ├── CommandManager.m           # 命令管理器
│   ├── processing/                # 处理命令（下方可以继续分包，这里只举例）
│   │   ├── ImportCommand.m        # 导入命令
│   │   ├── PreprocessCommand.m    # 预处理命令
│   │   └── AnalyzeCommand.m       # 分析命令
│   ├── visualization/             # 可视化命令
│   │   ├── PlotCommand.m          # 绘图命令
│   │   ├── ViewCommand.m          # 视图命令
│   │   └── ExportPlotCommand.m    # 导出命令
│   └── data_management/           # 数据管理命令
│       ├── CreateNodeCommand.m    # 创建节点
│       ├── DeleteNodeCommand.m    # 删除节点
│       └── UpdateNodeCommand.m    # 更新节点
├── algorithms/                    # 算法库（无状态）
│   ├── registry/                  # 算法注册表
│   │   └── AlgorithmRegistry.m    # 算法注册器
│   ├── import/                    # 数据导入
│   │   ├── EEGImporter.m          # EEG导入
│   │   ├── MEGImporter.m          # MEG导入
│   │   └── FormatAdapter.m        # 格式适配器
│   ├── preprocessing/             # 预处理
│   │   ├── FilterAlgorithm.m      # 滤波算法
│   │   ├── ICA.m                  # ICA去除
│   │   └── SegmentAlgorithm.m     # 数据分段
│   ├── analysis/                  # 分析算法
│   │   ├── SourceLocalization.m   # 源定位
│   │   ├── ConnectivityAnalysis.m # 连接性分析
│   │   └── TimeFrequencyAnalysis.m # 时频分析
│   └── visualization/             # 可视化算法
│       ├── TopoPlotter.m          # 拓扑图
│       ├── TimeSeriesPlotter.m    # 时序图
│       └── BrainViewer.m          # 3D脑图
├── models/                        # 数据模型
│   ├── ProjectInfo.m              # 项目信息
│   ├── ProtocolInfo.m             # 协议信息
│   ├── SessionInfo.m              # 会话信息
│   ├── DataInfo.m                 # 数据数据
│   ├── ChannelData.m              # 通道数据
│   ├── CortexData.m               # 皮层数据
│   ├── LeadfieldData.m            # 导联场数据
│   ├── EEGData.m                  # EEG数据
│   ├── MEGData.m                  # MEG数据
│   └── ResultData.m               # 结果数据
└── utils/                         # 工具函数
    ├── validation/                # 验证工具
    ├── logging/                   # 日志工具
    ├── helpers/                   # 辅助函数
    └── constants.m                # 常量定义
```

## 架构层次说明

### 应用层

```
SEAL/app/
```

**职责**：用户界面和交互
- **SEAL_App.mlapp** - 主应用程序入口，MATLAB App Designer应用
- **components/** - 可复用UI组件（按钮、面板、图表等）
- **layouts/** - 界面布局管理，负责窗口排列和视图组织

### 核心业务层

```
SEAL/core/
```

**职责**：系统核心逻辑和数据流管理

#### **中介者 (mediator/)**

中央协调器，模块间通信枢纽

#### **数据节点 (nodes/)**

- **NodeManager.m** - 节点生命周期管理
- **BaseNode.m** - 所有节点的抽象基类
- **具体节点类型**：ProtocolNode（协议）、ChannelNode（通道）、CortexNode（皮层）等
- **NodeFactory.m** - 节点创建工厂，实现对象创建模式

#### **工作区管理 (workspace/)**

待讨论

### 命令层

```
SEAL/commands/
```

**职责**：用户操作封装和执行
- **命令模式实现**，将操作封装为对象
- **CommandManager.m** - 命令执行、撤销、重做管理
- **分类命令**：processing/(处理)、visualization/(可视化)、data_management/(数据管理)

### 算法层

```
SEAL/algorithms/
```

**职责**：纯算法实现，无状态计算
- **registry/** - 算法注册
- **import/** - 数据导入和格式适配
- **preprocessing/** - 信号预处理算法
- **analysis/** - 分析算法
- **visualization/** - 可视化渲染算法

### 数据层
```
SEAL/models/
```
**职责**：数据结构定义和持久化
- **EEGData.m** - EEG数据模型
- **MEGData.m** - MEG数据模型  
- **ResultData.m** - 分析结果数据模型

## 数据流架构

```
用户操作 → 命令层 → 中介者 → 算法层 → 数据节点 → 可视化反馈
↑          ↓         ↓         ↓         ↓         ↓
应用层 ←── 工作区 ←── 事件系统 ←── 结果数据 ←── 处理流程 ←── 配置
```

## 设计模式
主要是中介者模式+命令模式

## 数据存储设计

### 件系统结构
```
Projects/
├── Project1/
│   ├── project_info.mat
│   ├── Protocols/
│   │   ├── Protocol1/
│   │   │   ├── protocol_info.mat
│   │   │   ├── Singleton/
│   │   │   │   ├── channel_info.mat
│   │   │   │   ├── cortex_info.mat
│   │   │   │   └── leadfield_info.mat
│   │   │   └── Sessions/
│   │   │       └── Session1/
│   │   │           ├── session_info.mat
│   │   │           ├── OriginalData/
│   │   │           │   └── data_info.mat
│   │   │           └── FilteredData1/
│   │   │               └── data_info.mat
```
