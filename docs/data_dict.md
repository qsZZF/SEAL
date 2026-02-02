# SEAL 数据字典

## 节点系统概述

SEAL 采用树形节点系统来管理所有数据，所有节点继承自 `BaseNode` 抽象基类。节点系统支持树形结构管理、按需加载、事件通知等功能。

## BaseNode (抽象基类)

所有节点的基类，提供统一的节点接口和树形结构管理。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `path` | string | 节点存储路径 | 必填 |
| `parent` | BaseNode | 父节点 | 树结构双向绑定 |
| `children` | BaseNode数组 | 子节点列表 |  |
| `uuid` | string | 节点唯一标识符 | 自动生成 |
| `tags` | string数组 | 标签数组 |  |
| `isSelected` | logical | 是否被选中 |  |
| `isLoaded` | logical | 是否已加载 |  |
| `metadata` | struct | 元数据 |  |
| `nodeInfo` | struct | 节点信息 |  |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `isRoot` | logical | 是否为根节点 |
| `isLeaf` | logical | 是否为叶节点 |
| `childCount` | double | 子节点数量 |
| `hasChildren` | logical | 是否有子节点 |

### 抽象属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `type` | string | 节点类型标识符（由子类实现） |
| `name` | string | 节点名称（由子类实现） |

### 方法

| 方法名 | 描述 |
|--------|------|
| `addChild(childNode)` | 添加子节点 |
| `removeChild(childNode)` | 移除子节点 |
| `getChild(childName, recursive)` | 获取子节点 |
| `findChildren(predicate, varargin)` | 查找满足条件的子节点 |
| `select()` | 选择节点 |
| `deselect()` | 取消选择节点 |
| `load()` | 加载节点数据（抽象方法） |
| `unload()` | 卸载节点数据（抽象方法） |
| `open(path)` | 打开节点（抽象方法） |
| `save()` | 保存节点（抽象方法） |
| `deleteFromDisk()` | 从磁盘删除节点（抽象方法） |
| `printTree(level)` | 打印节点树结构 |

---

## ProjectNode (项目节点)

项目的根节点，节点类型为 "ProjectNode"。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `projectInfo` | ProjectInfo | 项目信息对象 | 必填 |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 项目名称 |
| `protocolCount` | int32 | 协议数量 |
| `infoFile` | string | 项目信息文件路径 |
| `protocolFolder` | string | 协议文件夹路径 |

### 方法

| 方法名 | 描述 |
|--------|------|
| `createNewProtocol(protocolName, varargin)` | 创建新协议 |
| `openProtocol(protocolFolder)` | 打开现有协议 |
| `openChildNodes()` | 打开所有子节点 |

### 静态方法

| 方法名 | 描述 |
|--------|------|
| `createNew(projectName, projectPath, varargin)` | 创建新项目 |
| `openExisting(projectPath)` | 打开现有项目 |

---

## ProtocolNode (协议节点)

协议节点，节点类型为 "ProtocolNode"。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `protocolInfo` | ProtocolInfo | 协议信息对象 | 必填 |
| `chanlocsNode` | ChanlocsNode | 通道节点 |  |
| `cortexNode` | CortexNode | 皮层节点 |  |
| `leadfieldNode` | LeadfieldNode | 导联场节点 |  |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 协议名称 |
| `protocolType` | string | 协议类型 |
| `desc` | string | 协议描述 |
| `sessionCount` | int32 | 会话数量 |
| `chanlocsInfoFile` | string | 通道信息文件路径 |
| `cortexInfoFile` | string | 皮层信息文件路径 |
| `leadfieldInfoFile` | string | 导联场信息文件路径 |
| `infoFile` | string | 协议信息文件路径 |
| `sessionFolder` | string | 会话文件夹路径 |
| `headModelType` | string | 头模型类型 |
| `orientation` | array | 方向 |

### 方法

| 方法名 | 描述 |
|--------|------|
| `addSession(sessionNode)` | 添加会话节点 |
| `removeSession(sessionNode)` | 移除会话节点 |
| `setChanlocsNode(chanlocsNode)` | 设置通道节点 |
| `setCortexNode(cortexNode)` | 设置皮层节点 |
| `setLeadfieldNode(leadfieldNode)` | 设置导联场节点 |
| `openChanlocsFromData(dataPath)` | 从数据打开通道节点 |
| `openCortexFromData(dataPath)` | 从数据打开皮层节点 |
| `openLeadfieldFromData(dataPath)` | 从数据打开导联场节点 |
| `openDataFromData(dataPath)` | 从数据打开数据节点 |
| `getAllSessions()` | 获取所有会话节点 |

### 静态方法

| 方法名 | 描述 |
|--------|------|
| `createNew(protocolName, path, varargin)` | 创建新协议 |
| `openExisting(srcPath)` | 打开现有协议 |

---

## SessionNode (会话节点)

会话节点，节点类型为 "SessionNode"。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `sessionInfo` | SessionInfo | 会话信息对象 | 必填 |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 会话名称 |
| `sessionType` | string | 会话类型 |
| `desc` | string | 会话描述 |
| `dataCount` | int32 | 数据节点数量 |
| `infoFile` | string | 会话信息文件路径 |

### 方法

| 方法名 | 描述 |
|--------|------|
| `addDataNode(dataNode)` | 添加数据节点 |
| `removeDataNode(dataNode)` | 移除数据节点 |
| `getDataNodeByName(name)` | 根据名称获取数据节点 |
| `openDataFromData(dataPath)` | 从数据打开数据节点 |

### 静态方法

| 方法名 | 描述 |
|--------|------|
| `createNew(name, path, varargin)` | 创建新会话 |
| `openExisting(srcPath)` | 打开现有会话 |

---

## DataNode (数据节点)

数据节点，节点类型为 "DataNode"。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `dataInfo` | DataInfo | 数据信息对象 | 必填 |
| `data_cache` | any | 数据缓存 |  |
| `res_cache` | any | 结果缓存 |  |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 数据名称 |
| `data` | any | 数据（自动加载） |
| `result` | any | 结果（自动加载） |
| `infoFile` | string | 信息文件路径 |
| `size` | int64 | 数据大小 |
| `srate` | int64 | 采样率 |

### 方法

| 方法名 | 描述 |
|--------|------|
| `setResultPath(path)` | 设置结果路径 |

### 静态方法

| 方法名 | 描述 |
|--------|------|
| `fromInfo(dataInfo)` | 从信息创建数据节点 |
| `fromData(path)` | 从数据文件创建数据节点 |
| `openExisting(srcPath)` | 打开现有数据节点 |

---

## ChanlocsNode (通道节点)

通道节点，节点类型为 "ChanlocsNode"。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `chanlocsInfo` | ChanlocsInfo | 通道信息对象 | 必填 |
| `cache` | any | 通道数据缓存 |  |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 通道名称 |
| `data` | any | 通道数据（自动加载） |
| `infoFile` | string | 信息文件路径 |

### 静态方法

| 方法名 | 描述 |
|--------|------|
| `fromInfo(chanlocsInfo)` | 从信息创建通道节点 |
| `fromData(path)` | 从数据文件创建通道节点 |
| `openExisting(srcPath)` | 打开现有通道节点 |

---

## CortexNode (皮层节点)

皮层节点，节点类型为 "CortexNode"。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `cortexInfo` | CortexInfo | 皮层信息对象 | 必填 |
| `cache` | any | 皮层数据缓存 |  |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 皮层名称 |
| `data` | any | 皮层数据（自动加载） |
| `infoFile` | string | 信息文件路径 |

### 静态方法

| 方法名 | 描述 |
|--------|------|
| `fromInfo(cortexInfo)` | 从信息创建皮层节点 |
| `fromData(path)` | 从数据文件创建皮层节点 |
| `openExisting(srcPath)` | 打开现有皮层节点 |

---

## LeadfieldNode (导联场节点)

导联场节点，节点类型为 "LeadfieldNode"。

### 属性

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `leadfieldInfo` | LeadfieldInfo | 导联场信息对象 | 必填填 |
| `cache` | any | 导联场数据缓存 |  |

### 依赖属性

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 导联场名称 |
| `data` | any | 导联场数据（自动加载） |
| `infoFile` | string | 信息文件路径 |
| `headModelType` | string | 头模型类型 |
| `orientation` | array | 方向 |

### 静态方法

| 方法名 | 描述 |
|--------|------|
| `fromInfo(leadfieldInfo)` | 从信息创建导联场节点 |
| `fromData(path)` | 从数据文件创建导联场节点 |
| `openExisting(srcPath)` | 打开现有导联场节点 |

---

## Info 类

### ProjectInfo

项目信息管理类。

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 项目名称 |
| `desc` | string | 项目描述 |
| `metadata` | struct | 用户自定义元数据 |
| `createdDate` | datetime | 创建日期 |
| `modifiedDate` | datetime | 修改日期 |

### ProtocolInfo

协议信息管理类。

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 协议名称 |
| `type` | string | 协议类型 |
| `desc` | string | 协议描述 |
| `metadata` | struct | 用户自定义元数据 |
| `createdDate` | datetime | 创建日期 |
| `modifiedDate` | datetime | 修改日期 |

### SessionInfo

会话信息管理类。

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 会话名称 |
| `type` | string | 会话类型 |
| `desc` | string | 会话描述 |
| `metadata` | struct | 用户自定义元数据 |
| `createdDate` | datetime | 创建日期 |
| `modifiedDate` | datetime | 修改日期 |

### DataInfo

数据信息管理类。

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 数据名称 |
| `desc` | string | 描述 |
| `size` | int64 | 数据大小 |
| `srate` | int64 | 采样率 |
| `dataPath` | string | 数据文件路径 |
| `resultPath` | string | 结果文件路径 |
| `metadata` | struct | 用户自定义元数据 |
| `createdDate` | datetime | 创建日期 |
| `modifiedDate` | datetime | 修改日期 |

### ChanlocsInfo

通道信息管理类。

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 通道名称 |
| `desc` | string | 描述 |
| `dataPath` | string | 通道数据文件路径 |
| `metadata` | struct | 用户自定义元数据 |
| `createdDate` | datetime | 创建日期 |
| `modifiedDate` | datetime | 修改日期 |

### CortexInfo

皮层信息管理类。

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 皮层名称 |
| `desc` | string | 描述 |
| `dataPath` | string | 皮层数据文件路径 |
| `metadata` | struct | 用户自定义元数据 |
| `createdDate` | datetime | 创建日期 |
| `modifiedDate` | datetime | 修改日期 |

### LeadfieldInfo

导联场信息管理类。

| 字段名 | 数据类型 | 描述 |
|--------|----------|------|
| `name` | string | 导联场名称 |
| `desc` | string | 描述 |
| `headModelType` | string | 头模型类型 |
| `orientation` | array | 方向 |
| `dataPath` | string | 导联场数据文件路径 |
| `metadata` | struct | 用户自定义元数据 |
| `createdDate` | datetime | 创建日期 |
| `modifiedDate` | datetime | 修改日期 |

---

## 数据层级关系

```
ProjectNode (项目)
  └── ProtocolNode (协议)
      ├── ProtocolInfo
      ├── ChanlocsNode (通道)
      ├── CortexNode (皮层)
      ├── LeadfieldNode (导联场)
      └── SessionNode (会话)
          ├── SessionInfo
          └── DataNode (数据)
              ├── DataInfo
              ├── data (数据缓存)
              └── result (结果缓存)
```

## 节点生命周期

1. **创建**: 使用静态工厂方法 `createNew()` 创建新节点
2. **打开**: 使用静态工厂方法 `openExisting()` 打开现有节点
3. **加载**: 调用 `load()` 方法加载数据到内存
4. **保存**: 调用 `save()` 方法保存数据到磁盘
5. **卸载**: 调用 `unload()` 方法释放内存
6. **删除**: 调用 `deleteFromDisk()` 方法从磁盘删除
