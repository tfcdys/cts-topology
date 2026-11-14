# 自适应四叉 H-tree 算法：伪代码与说明

本文档描述 `h_tree.cpp` 中 `HTreeBuilder` 实际实现的 H-tree 构造算法。
算法的目标不是复刻规则阵列上的等长经典 H-tree，而是在任意平面 sink 分布上，
构造一棵**四叉递归、共享 H 臂、可曼哈顿布线**的时钟树。

实现位置：

- `HTreeBuilder`：[../h_tree.cpp#L29-L266](../h_tree.cpp#L29-L266)
- `choose_split`：[../h_tree.cpp#L73-L116](../h_tree.cpp#L73-L116)
- `quadrisect`：[../h_tree.cpp#L161-L199](../h_tree.cpp#L161-L199)
- `build_group`：[../h_tree.cpp#L200-L235](../h_tree.cpp#L200-L235)
- `connect_side`：[../h_tree.cpp#L237-L265](../h_tree.cpp#L237-L265)

---

## 1. 输入、输出与约定

### 1.1 输入

平面 sink 集合

\[
S=\{s_i=(x_i,y_i,c_i)\mid i=1,\dots,n\}
\]

其中：

- `x_i, y_i`：sink 的二维坐标，单位 um；
- `c_i`：sink 电容负载，单位 fF；
- `n`：sink 数量。

### 1.2 输出

一棵有根树

\[
T=(V,E,r)
\]

并满足：

1. 每个 sink 对应一个叶子节点，且恰好出现一次；
2. 存在唯一根节点 `r`，模拟时钟源；
3. 除根外每个节点恰有一个父节点；
4. 所有节点从根可达；
5. `|E| = |V| - 1`；
6. 每条边带曼哈顿实现标签 `V`、`HV` 或 `VH`，用于可视化 L 形走线。

节点类型：

| kind | 含义 |
|---|---|
| `root` | 时钟树根节点，位于全部 sink 的包围盒中心 |
| `internal` | 递归子树的内部节点，位于该组 sink 的包围盒中心 |
| `bend` | H 臂上的拐点/分支点，位于父节点正上/下方 |
| `sink` | 真实 sink 叶子，保留原始坐标、`sink_id`、`cap_ff` |

### 1.3 构造/合并顺序

H-tree 不采用 DME 的 cluster 两两合并。它是自顶向下的区域细分：

```text
BUILD_GROUP(ALL)
  └─ QUADRISECT -> LB, LT, RB, RT
       ├─ BUILD_GROUP(LB)     // 深度优先
       ├─ BUILD_GROUP(LT)
       ├─ BUILD_GROUP(RB)
       └─ BUILD_GROUP(RT)
     CONNECT_SIDE(TOP:    LT, RT)
     CONNECT_SIDE(BOTTOM: LB, RB)
```

- 父节点先创建；
- 四个子组按 `LB → LT → RB → RT` 的固定顺序递归；
- 子树全部构建完成后，先连接上侧 H 臂，再连接下侧 H 臂；
- 排序和分组使用 `(x,y,sink_id)` / `(y,x,sink_id)` 全序，结果确定。

### 1.4 关键几何符号

- `BBOX(G)`：点集 `G` 的轴对齐包围盒；
- `P(G)`：包围盒中心 `((xmin+xmax)/2, (ymin+ymax)/2)`；
- `ε`：浮点零阈值，实现中取 `1e-12`；
- `HV`：先水平后垂直的曼哈顿走线；
- `VH`：先垂直后水平的曼哈顿走线；
- `V`：纯垂直走线。

---

## 2. 总体伪代码

```text
ALGORITHM BUILD_HTREE(S)

Input : sink set S = {s_1, ..., s_n}, n >= 1
Output: rooted H-tree T = (V, E, root)

1  ids ← [1, 2, ..., n]
2  T.root ← BUILD_GROUP(ids, is_root = true)
3  assert number of sink nodes in T = n
4  assert |E| = |V| - 1
5  return T
```

`BUILD_GROUP` 是核心递归过程：

```text
FUNCTION BUILD_GROUP(ids, is_root)

Input : ids  —— 当前子树包含的 sink 下标列表
        is_root —— 是否为整棵树的根
Output: 当前子树的根节点 id

1  if |ids| = 1 then
2      i ← ids[0]
3      return NEW_SINK_NODE(s_i)          // 叶子保留原始坐标/电容
4  end if
5
6  b ← BBOX(ids)
7  cx ← (b.xmin + b.xmax) / 2
8  cy ← (b.ymin + b.ymax) / 2
9  p ← NEW_NODE((cx, cy), kind = is_root ? ROOT : INTERNAL)
10
11 (G_LB, G_LT, G_RB, G_RT) ← QUADRISECT(ids)
12 Top ← empty list
13 Bottom ← empty list
14
15 for each non-empty G_q in [G_LB, G_LT, G_RB, G_RT] do
16     child ← BUILD_GROUP(G_q, is_root = false)
17     if q ∈ {LT, RT} then
18         Top.append(child)
19     else
20         Bottom.append(child)
21     end if
22 end for
23
24 CONNECT_SIDE(p, (cx, cy), Top, side = TOP)
25 CONNECT_SIDE(p, (cx, cy), Bottom, side = BOTTOM)
26 return p
```

注意：

- 第 9 行先创建父节点，再递归创建四个子节点；
- 第 15–22 行只处理非空子集，因此实际度数可能是 1、2、3 或 4；
- `Top` 包含左上、右上两个子区；`Bottom` 包含左下、右下两个子区；
- 每个非空子区都会得到一次且仅一次递归调用。

---

## 3. 四叉划分：QUADRISECT

算法使用“先按 X 二分、再在每个半区按 Y 二分”的方式得到四个几何子集。
默认情况下，切分点取数量中位数，因此四个子集数量均衡；开启
`--cap-aware` 后，切分点在 40%–60% 数量范围内选择总电容最均衡的位置。

与直接按全局中位线切四象限相比，按“先 X 后半区 Y”的方式对非均匀分布
更稳健，不会因为某个象限过空而退化成极深的链。

```text
FUNCTION QUADRISECT(ids)

Input : ids —— 当前组的 sink 下标
Output: (G_LB, G_LT, G_RB, G_RT) 四个子集

1  A ← SORT(ids, key = (x, y, sink_id))
2  m ← CHOOSE_SPLIT(A)
3  L ← A[0 : m]                 // 左半区
4  R ← A[m : |A|]               // 右半区
5
6  for each side in [L, R] do
7      A_side ← SORT(side, key = (y, x, sink_id))
8      k ← CHOOSE_SPLIT(A_side)
9      bottom ← A_side[0 : k]
10     top    ← A_side[k : |A_side|]
11 end for
12
13 G_LB ← L.bottom
14 G_LT ← L.top
15 G_RB ← R.bottom
16 G_RT ← R.top
17 return (G_LB, G_LT, G_RB, G_RT)
```

切分点选择：

```text
FUNCTION CHOOSE_SPLIT(sorted)

Input : sorted —— 已按某个坐标轴排序的 sink 下标列表
Output: 低侧元素个数 m

1  n ← |sorted|
2  if cap_aware = false or n < 5 then
3      return floor(n / 2)          // 默认：数量中位数
4  end if
5
6  lo ← ceil(0.40 * n)
7  hi ← floor(0.60 * n)
8  total ← Σ_{i ∈ sorted} EFFECTIVE_CAP(i)
9  if total <= 0 then return floor(n / 2)
10
11 best_score ← +∞
12 best_m ← floor(n / 2)
13 running ← 0
14 for i ← 0 to n-1 do
15     running ← running + EFFECTIVE_CAP(sorted[i])
16     m ← i + 1
17     if m < lo or m > hi then continue
18     cap_imbalance ← |2 * running - total| / total
19     count_dev ← |m - n/2| / (n/2)
20     score ← cap_imbalance + 0.10 * count_dev
21     if score < best_score then
22         best_score ← score
23         best_m ← m
24     end if
25 end for
26 return best_m
```

其中：

```text
EFFECTIVE_CAP(i) = cap_i     if cap_aware = true
                 = 1.0       if cap_aware = false
```

默认关闭电容感知时，所有 sink 的有效电容完全相同，`CHOOSE_SPLIT` 退化为
数量中位数，因此不会产生任何电容驱动的偏置。

性质：

- 四个子集两两不相交，并集等于 `ids`；
- 默认模式下每个子集大小不超过 `ceil(n / 4)`；开启电容感知后，由于每级
  切分限制在 40%–60%，每个子集大小不超过约 `0.36n`，仍严格小于 `n`，
  递归必然终止；
- 采用 `(x,y,sink_id)` 和 `(y,x,sink_id)` 的全序比较，结果确定、可复现；
- 左半区所有点的 X 坐标不大于右半区；同一半区内，下半区 Y 不大于上半区。

`Top` 与 `Bottom` 的划分为：

```text
Top    = {G_LT, G_RT}
Bottom = {G_LB, G_RB}
```

---

## 4. H 臂连接：CONNECT_SIDE

`CONNECT_SIDE` 把同一侧（上侧或下侧）的 1 或 2 个子树根连接到父节点，
并在需要时生成共享的水平 H 臂。

```text
FUNCTION CONNECT_SIDE(p, P, C, side)

Input : p    —— 父节点
        P    —— 父节点坐标 (P.x, P.y)
        C    —— 同侧子节点列表，最多 2 个
        side —— TOP 或 BOTTOM
Output: 无；直接向 T 中添加节点和边

1  if C is empty then
2      return
3  end if
4
5  ybar ← (1 / |C|) * Σ_{c ∈ C} c.y
6  if side = TOP then
7      ybar ← max(P.y, ybar)
8  else
9      ybar ← min(P.y, ybar)
10 end if
11
12 if |ybar − P.y| < ε then
13     // H 臂退化到父节点高度，不需要额外 bend 节点
14     for each c in C do
15         ADD_EDGE(p → c, route = HV)
16     end for
17     return
18 end if
19
20 b ← NEW_BEND_NODE((P.x, ybar))
21 ADD_EDGE(p → b, route = V)
22
23 for each c in C do
24     ADD_EDGE(b → c, route = HV)
25 end for
```

### 4.1 几何解释

对两个上侧子节点，`CONNECT_SIDE` 形成如下结构：

```text
          C_TL                              C_TR
            │                                 │
            └──────────── B_T ────────────────┘      y = ybar_T
                          │
                          P                          root/internal
                          │
            ┌──────────── B_B ────────────────┐      y = ybar_B
            │                                 │
          C_BL                              C_BR
```

- `P → B_T`、`P → B_B` 是纯垂直的 H-tree 主干；
- `B_T → C_TL` 与 `B_T → C_TR` 使用 `HV` 走线：
  - 先沿 `y = ybar_T` 水平走到子节点的 X 坐标；
  - 再垂直走到子节点实际 Y 坐标；
- 左右两条水平段共享 `B_T` 并形成一条连续的顶部 H 臂；
- 子节点自身又是下一层 H-tree 的根，从而递归形成嵌套 H 结构。

当某一侧只有一个子节点时，只要 `ybar` 与父节点高度不同，仍会生成
“竖干 + 单侧 H 臂”，使整棵树保持正交 H 结构，而不是退化为斜线；
若 `ybar == P.y`，则直接使用 `HV` 边，避免零长度竖干。

### 4.2 为什么使用 `ybar`

- 若直接使用父节点的 `P.y` 作为 H 臂高度，臂可能离子节点太远，视觉上
  不再像 H-tree；
- 若直接使用某个子节点的 Y，另一子节点会产生额外竖直绕行；
- 使用同侧子节点 Y 的平均值作为 `ybar`，再夹紧到父节点同侧，可以同时：
  - 保持 H 臂水平；
  - 让两个子节点共享水平臂；
  - 避免 H 臂越过父节点跑到错误的一侧。

### 4.3 边的曼哈顿实现

拓扑 JSON 中的边不一定在欧式平面上是直线。`route` 标签定义实际走线：

| route | 实际路径 |
|---|---|
| `V` | `(x1,y1) → (x1,y2)` |
| `H` | `(x1,y1) → (x2,y1)` |
| `HV` | `(x1,y1) → (x2,y1) → (x2,y2)` |
| `VH` | `(x1,y1) → (x1,y2) → (x2,y2)` |

因此，一条 `b → c` 的 `HV` 边在物理上等价于一个水平段加一个垂直段；
节点 `c` 仍然是该段末端，不需要为每个拐角单独建图节点。

---

## 5. 递归构造示例

以四个 sink 为例，假设四叉划分得到：

```text
G_LB = {s1}      G_LT = {s2}
G_RB = {s3}      G_RT = {s4}
```

构造过程：

```text
BUILD_GROUP({s1,s2,s3,s4})
├─ 创建根 p，位于四个点的包围盒中心
├─ 划分得到四个单点组
├─ 递归创建 4 个 sink 叶子
├─ CONNECT_SIDE(p, Top={s2_node, s4_node}, TOP)
│  ├─ 创建 B_T
│  ├─ 添加 p → B_T（V）
│  ├─ 添加 B_T → s2_node（HV）
│  └─ 添加 B_T → s4_node（HV）
└─ CONNECT_SIDE(p, Bottom={s1_node, s3_node}, BOTTOM)
   ├─ 创建 B_B
   ├─ 添加 p → B_B（V）
   ├─ 添加 B_B → s1_node（HV）
   └─ 添加 B_B → s3_node（HV）
```

假设上、下两侧的 `ybar` 都不等于父节点高度，则得到 `7` 个节点、`6` 条边，
满足 `|E| = |V| − 1`。若某一侧的 `ybar` 恰好等于父节点高度，该侧不创建
bend 节点，节点数和边数相应各少 1，但仍然是一棵树。

对于更多 sink，每个子节点继续作为下一层的 `internal` 节点，重复同样的
四叉划分和 H 臂连接，直到每个分支只剩一个 sink。

---

## 6. 正确性说明

### 6.1 覆盖性

**命题：** 每个输入 sink 恰好对应输出树中的一个 sink 叶子。

证明要点：

1. `QUADRISECT` 产生的四个子集互不相交，并集等于当前 `ids`；
2. 每次递归只把子集传给下一层，不会重复处理；
3. 当 `|ids| = 1` 时创建且只创建一个 sink 节点；
4. 因此每个 sink 下标最终恰好到达一个叶子。

实现中 `build()` 还会检查 `sink_nodes == n`，否则报内部错误。

### 6.2 连通性

**命题：** 根可以到达所有节点。

证明要点：

- 对每个非空子组，`BUILD_GROUP` 返回其子树根；
- `CONNECT_SIDE` 通过 `p → b → child` 或 `p → child` 把子树根接到父节点；
- 递归结构天然形成从根出发的路径；
- 因此每个叶子都有一条从根到它的有向路径。

### 6.3 无环性与树结构

**命题：** 输出是一棵有根树，且 `|E| = |V| − 1`。

证明要点：

- 每个子树的根只被其父层连接一次；
- 边方向始终从父层指向子层，不会产生回边；
- 新节点通过一条边加入，因此每加入一个非根节点恰好增加一条边；
- 所以最终 `|E| = |V| − 1`，加上连通性即为一棵树。

实现中 `build()` 会断言上述边数关系。

### 6.4 正交可布线性

**命题：** 每条边都存在对应的曼哈顿实现。

证明要点：

- `V` 边是纯垂直段；
- `HV` 边可分解为水平段加垂直段；
- H 臂使用 `HV` 时，水平段落在 `y = ybar` 上，垂直段落在子节点 X 上；
- 所有段均为轴对齐线段，因此拓扑可按曼哈顿规则布线。

---

## 7. 复杂度

设当前组大小为 `m`。

- 四叉划分需要对 `m` 个点排序，主排序代价 `O(m log m)`；
- 四个子组大小至多 `ceil(m/4)`；
- 递归深度至多 `O(log n)`；
- 总复杂度：

默认模式：

\[
T(n)=4T(n/4)+O(n\log n)=O(n\log^2 n)
\]

开启电容感知后，每个子问题规模不超过约 `0.36n`，递推仍然满足
`O(n log^2 n)` 的上界。

- 空间：输出本身 `O(n)`，递归栈 `O(log n)`，中间排序/分组向量合计 `O(n)`；
- 节点数：每个内部组最多增加 1 个 internal 节点和 2 个 bend 节点，因此
  总节点数为 `O(n)`，边数同样为 `O(n)`。

本机实测：100000 sink 的 H-tree 构造约 0.3 秒，拓扑 JSON 约 30 MB。

---

## 8. 与经典 H-tree 的区别

经典 H-tree 通常假设规则阵列或规则网格，每一层使用等长 H 臂，从而让所有
叶子到根的线长尽量相等。本实现面向任意平面 sink 分布，因此做了如下调整：

1. **几何自适应**：父节点取当前组包围盒中心，而不是固定网格中心；
2. **数量均衡四叉划分**：按 X/Y 中位数划分，避免空象限导致退化；
3. **H 臂高度自适应**：`ybar` 取同侧子节点 Y 均值，并夹紧到父节点同侧；
4. **允许非满四叉**：空子组直接跳过，度数可以是 1–4；
5. **不插入延时补偿**：不通过蛇形绕线强制等长，因此 max/min 路径差是
   拓扑本身的 skew，而不是强制为零；
6. **电容感知默认关闭**：`cap_ff` 不参与默认划分，所有 sink 使用统一的
   有效权重 `1.0`；原始 `cap_ff` 仍保留在拓扑中用于追溯。开启
   `--cap-aware` 后，切分点才使用电容加权，并限制在 40%–60% 数量范围内。

这些调整使算法能够处理不规则平面 sink 集合，同时保留 H-tree 的递归分层、
共享水平臂和曼哈顿可布线特征。

---

## 9. 实现映射

| 伪代码 | C++ 实现 |
|---|---|
| `BUILD_HTREE` | `HTreeBuilder::build()` |
| `BUILD_GROUP` | `HTreeBuilder::build_group()` |
| `QUADRISECT` | `HTreeBuilder::quadrisect()` |
| `CHOOSE_SPLIT` | `HTreeBuilder::choose_split()` |
| `EFFECTIVE_CAP` | `HTreeBuilder::cap_weight()` |
| `CONNECT_SIDE` | `HTreeBuilder::connect_side()` |
| `NEW_SINK_NODE` / `NEW_NODE` | `HTreeBuilder::add_node()` |
| `ADD_EDGE` | `HTreeBuilder::add_edge()` |
| `BBOX` | `bounds_of()` |

相关输出字段：

```json
{
  "capacitance_aware": false,
  "capacitance_mode": "uniform",
  "root": 0,
  "nodes": [
    {"id": 0, "x": 500.0, "y": 500.0, "kind": "root"},
    {"id": 1, "x": 500.0, "y": 700.0, "kind": "bend"}
  ],
  "edges": [
    {"source": 0, "target": 1, "route": "V"}
  ]
}
```

评估器按边两端坐标的曼哈顿距离计算：

```text
len(u,v) = |u.x - v.x| + |u.y - v.y|
```

因此 `HV`、`VH`、`V`、`H` 在评估器中都统一折算为曼哈顿长度。
