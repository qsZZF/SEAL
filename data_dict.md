# SEAL 数据字典

## 项目与协议数据结构

### ProjectNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `name` | string | 项目名称 | 必填，最大长度50字符 |
| `path` | string | 项目文件存储路径 | 必填 |
| `desc` | string | 项目描述 | 最大长度200字符 |
| `createdDate` | datetime | 创建日期 | 自动生成 |
| `modifiedDate` | datetime | 最后修改日期 | 自动更新 |
| `metadata` | struct | 用户自定义元数据 |  |
| `protocolList` | Protocol数组 | 子协议列表 | 必填，内容可空|

ProjectNode应当是整棵树的root

### ProtocolNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `project` | string | 所属项目 | 树结构双向绑定 |
| `type` | string | 协议类型 |  |
| `name` | string | 协议名称 | 最大长度50字符 |
| `desc` | string | 协议描述 |  |
| `createdDate` | datetime | 创建日期 | 自动生成 |
| `sessionList`| Session数组 | 所有Session||
| `cortex` | Cortex| 皮层模型 ||
| `leadfield` | LeadField| 引导场对象| |
| `channel`| Channel | 通道对象 ||

当新建Protocol时，在Project文件夹下创建子文件夹。

### SessionNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `protocol` | string | 所属协议 | 双向绑定 |
| `type` | string | 数据类型 | 考虑到未来有sim，保留 |
| `name` | string | 数据文件名 | 必填 |
| `dataList`| DateNode数组 | 所有DateNode||
| `metadata` | struct | 元数据 | 可选 |

### ChannelNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `name` | cell | 通道名称列表 | 必填 |
| `path` | string | 通道数据储存的位置 | 必填 |
| `metadata` | struct | 元数据 | 可选 |

### DataNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `path` | string | 数据储存的位置 | 必填 |
| `eventInfo` | struct 数组 | 事件信息，填入预处理啥的 | 可选 |
| `resultList`| Result | 所有Result|必填，内容可空|
| `metadata` | struct | 元数据 | 可选 |

### CortexNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `protocol` | string | 所属协议 |  |
| `name` | string | 模型名称 |  |
| `path` | string | 数据储存的位置 | 必填 |
| `metadata` | struct | 元数据 | 可选 |

### LeadfieldNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `protocol` | string | 所属协议 |  |
| `name` | string | 导联场名 |  |
| `path` | matrix | 导联场矩阵 |  |
| `metadata` | struct | 元数据 | 可选 |

### ResultNode

| 字段名 | 数据类型 | 描述 | 约束 |
|--------|----------|------|------|
| `name` | string | 结果名 | 必填，唯一 |
| `data` | Data | 数据 | 必填，树结构双向绑定 |
| `metadata` | struct | 元数据 | 可选 |
