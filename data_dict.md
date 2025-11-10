# SEAL 数据字典

## 项目与协议数据结构

### Project
| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `projectName` | string | 项目名称 | 必填，最大长度50字符 |
| `description` | string | 项目描述 | 最大长度200字符 |
| `createdDate` | datetime | 创建日期 | 自动生成 |
| `modifiedDate` | datetime | 最后修改日期 | 自动更新 |
| `filePath` | string | 项目文件存储路径 | 必填 |
| `metadata` | struct | 用户自定义元数据 |  |
| `protocolList` | Protocol数组 | 子协议列表 | 必填，内容可空|

### Protocol
| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `project` | string | 所属项目 | 树结构双向绑定 |
| `protocolType` | string | 协议类型 | 'real_data', 'simulation' |
| `protocolName` | string | 协议名称 | 最大长度50字符 |
| `description` | string | 协议描述 |  |
| `created_date` | datetime | 创建日期 | 自动生成 |
| `dataList`| Data数组 | 所有Data，包括经过预处理的数组|必填，内容可空|
| `cortex` | Cortex| 皮层模型 ||
| `leadField` | LeadField| 引导场 ||
| `channel`|

当新建Protocol时，在Project文件夹下创建子文件夹，而非新增路径指针，该行为在大型软件工程上会因为不必要的拷贝被批斗，但是对于该项目而言，复制开销不大，为方便访问这么做问题不大。

## 核心数据对象

### Data
| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `protocol` | string | 所属协议 | 双向绑定 |
| `dataType` | string | 数据类型 | 考虑到未来有sim，保留 |
| `dataName` | string | 数据文件名 | 必填 |
| `preprocessedList`| preprocessed数组 | 所有Result|必填，内容可空|
| `metadata` | struct | 数据元数据 | 可选 |

Data很大，要存指针

### Channel

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `channelNames` | cell | 通道名称列表 | 必填 |
| `channelLocations` | struct | 通道位置信息 | 可选 |

### PreData

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `eventInfo` | struct 数组 | 事件信息，填入预处理啥的 | 可选 |
| `result`| Result | 所有Result|必填，内容可空|

### Cortex
| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `protocol` | string | 所属协议 |  |
| `name` | string | 模型名称 |  |
| `vertices` | matrix | 顶点坐标(N×3) |  |
| `faces` | matrix | 面片索引(M×3) |  |
| `numVertices` | int32 | 顶点数量 |  |
| `scouts` | struct | 脑区标记信息 | |
| `coordinateSystem` | string | 坐标系 |  |

### LeadField
| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `protocol` | string | 所属协议 |  |
| `name` | string | 导联场名 |  |
| `leadfieldMat` | matrix | 导联场矩阵 |  |
| `numSources` | int32 | 源数量 | 必填 |
| `orientationConstraint` | string | 方向约束 | 应当是个enum |
| `gainMatrix` | matrix | 增益矩阵 |  |

## 处理结果数据

### Result
| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `resultName` | string | 结果名 | 必填，唯一 |
| `data` | Data | 数据 | 必填，外键 |
| `algorithmName` | string | 算法名称 |  |
| `algorithmParams` | struct | 算法参数 |  |
| `sourceActivity` | matrix | 源活动数据 |  |
| `timeVector` | vector | 时间向量 | 必填 |

### PreprocessingParameter
| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `process_type` | enum | 处理类型 | 'filter', 'ica', 'resample', 'rereference' |
| `parameters` | struct | 处理参数 | 必填 |
