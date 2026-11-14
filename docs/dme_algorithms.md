# ZST-DME 与 BST-DME 实现说明

本项目实现了两种基于 **Deferred-Merge Embedding (DME)** 的时钟树拓扑，
每种算法是独立源文件和独立可执行文件：

| 算法 | 源文件 | 可执行文件 | 功能 |
|---|---|---|---|
| `zst_dme` | [zst_dme.cpp](../zst_dme.cpp) | `build/zst_dme` | Zero-Skew Tree DME，零 skew 嵌入 |
| `bst_dme` | [bst_dme.cpp](../bst_dme.cpp) | `build/bst_dme` | Bounded-Skew Tree DME，有界 skew 嵌入 |
| 公共部分 | [dme_common.cpp](../dme_common.cpp) | — | Manhattan TRR 几何、NNG 拓扑、输出 |

公共的 case 读取、JSON 输出和 CLI 位于
[cts_common.cpp](../cts_common.cpp)。

算法概念参考：

- ZSTDME: <https://github.com/geziangfinn/ZSTDME>
- BST-DME: <https://github.com/abk-openroad/BST-DME>
- J. Cong, A. B. Kahng, C.-K. Koh, C.-W. A. Tsao,
  “Bounded-Skew Clock and Steiner Routing”.

本实现是独立编写的简化版本，不复用参考仓库源码。与参考实现相比：

- 使用 **线性延迟模型**：`delay == Manhattan 走线长度`；
- 两种算法共用 **最近邻聚类（NNG）** 生成二叉拓扑；
- 不做 buffer insertion、障碍规避、Elmore 延迟优化；
- 电容感知默认关闭，所有 sink 的延迟权重相同；
- 输出统一为本项目的 rooted-tree JSON schema。

因此本实现适合比较拓扑结构、覆盖性、曼哈顿线长和根到 sink 路径 skew，
不等价于参考仓库的完整工业流程。

---

## 1. 总体流程

```text
BUILD_DME(S, variant):
  1. T_topo ← BUILD_NEAREST_NEIGHBOR_TOPOLOGY(S)
  2. depth ← TREE_DEPTH(T_topo)
  3. for each node in bottom-up order:
         compute merging segment / TRR and delay
  4. choose root point on root merging segment
  5. for each child top-down:
         choose child point on its merging segment
  6. emit rooted tree with Manhattan-routed edges
```

### 1.1 最近邻拓扑与合并顺序

两种 DME 使用**完全相同**的合并顺序：

```text
active ← {sink_0, ..., sink_{n-1}}
next_id ← n

while |active| > 1:
    (a,b) ← argmin_{i<j} L1(rep_i, rep_j)
            距离相同时按 (node_id_i, node_id_j) 字典序打破平局
    parent ← next_id++
    left[parent] ← a
    right[parent] ← b
    rep[parent] ← midpoint(rep_a, rep_b)
    active ← active \ {a,b} ∪ {parent}

root ← 最后一个 parent
```

- 每轮只合并一对，共 `n-1` 轮；
- 父 cluster 代表点取两个子代表的几何中点；
- 叶节点 id 为 `0..n-1`，内部节点 id 为 `n..2n-2`，根最后创建；
- 父 id 总大于子 id，因此自底向上只需按内部节点 id 递增遍历；
- `zst_dme` 与 `bst_dme` 的区别不在合并顺序，而在合并半径和 skew 预算。

该步骤是确定性的。测试规模 16 / 100 时，朴素 `O(n^3)` 实现已足够快。

### 1.2 二叉拓扑的输出

对于 `n` 个 sink：

- 叶子节点：`0 .. n-1`，对应真实 sink；
- 内部节点：`n .. 2n-2`，最后一个为根；
- 节点数：`2n-1`；
- 边数：`2n-2`；
- 所有边方向为父 → 子，`route = "HV"`。

---

## 2. DME 几何

### 2.1 合并段

每个子树用一个 merge segment（MS）表示可行的根位置集合。叶子 MS 退化为
sink 点：

```text
MS_leaf = {sink_position}
```

内部节点的 MS 由两个子 MS 和两个延迟半径决定。

### 2.2 Manhattan TRR

点 `p` 到合并段 `MS` 的曼哈顿距离记为：

```text
d(p, MS) = min_{q ∈ MS} (|p.x-q.x| + |p.y-q.y|)
```

以 `MS` 为核心、半径 `r` 的 TRR 为：

```text
TRR(MS, r) = { p | d(p, MS) <= r }
```

其边界由斜率为 `0 / ∞ / +1 / -1` 的线段组成。实现中通过“核心段两个端点
分别与 L1 菱形四个顶点做 Minkowski 和，再求凸包”生成边界多边形，避免为
不同斜率单独写分支。

### 2.3 合并段计算

给定子段 `MS_L, MS_R`、半径 `r_L, r_R`：

```text
L = min_{p∈MS_L, q∈MS_R} d(p,q)
```

在零 skew 条件下：

```text
r_L = (delay_R - delay_L + L) / 2
r_R = (delay_L - delay_R + L) / 2
r_L + r_R = L
```

父合并段取：

```text
MS_parent = boundary(TRR(MS_L, r_L)) ∩ filled(TRR(MS_R, r_R))
```

实现上把 `MS_L` 的边界多边形边逐条裁剪到 `MS_R` 的填充多边形中，再取
最长结果；反向再做一次以增强数值鲁棒性。若几何退化，使用最近点 fallback。

### 2.4 自顶向下嵌入

根点取根合并段的中点。对每个子节点：

- 若子节点是 sink，位置固定为 sink 坐标；
- 否则在子合并段上寻找与父点距离等于 `r_child` 的点；
- 通过 L1 菱形边界与子合并段求交，选择靠近子段中点的交点。

理论上，任取满足距离约束的点都保持零 skew；实现选择靠近段中点的点以减少
数值极端。若由于半径夹紧导致某些叶子的路径长度仍短于目标值，ZST-DME 会在
叶子父边上插入确定性 U 形 snaking，使所有根到 sink 路径长度相等；这会增加
少量 bend 节点和线长。

---

## 3. ZST-DME

ZST-DME 使用零 skew 半径：

```text
e_L = (delay_R - delay_L + L) / 2
e_R = (delay_L - delay_R + L) / 2
```

若某一半径为负，则夹紧到 0，并把另一半径设为延迟差；随后计算父合并段。

自顶向下精确选择子节点位置，使每条父子边的实际曼哈顿长度等于计划半径。
对测试中的规则/随机个案，实际根到所有 sink 路径长度相等。

命令：

```bash
./build/zst_dme --case sinks_100
# 或
./run_cts.sh --case sinks_100 --topology zst_dme
```

---

## 4. BST-DME

BST-DME 允许每个 merge 的左右延迟存在有限差异。设全局 skew 上界为 `B`，
拓扑深度为 `D`，每个 merge 使用局部预算：

```text
b = B / D
delta = b / 2
e_L ← e_L - delta
e_R ← e_R + delta
```

约束：

- `e_L + e_R = L` 保持不变；
- `e_L, e_R >= 0`；
- 允许的左右延迟差由 `delta` 控制，因此整体 skew 目标不超过 `B`；
- 当前实现固定把预算向右子节点倾斜（负号方向），以保证结果完全确定。

这里的倾斜方向是一个确定性启发式，不代表最优分配。

该做法移动合并段，但不增加局部两边长度之和；和参考仓库完整的 TRR 面积
优化相比，这是一个简化但有界的实现。

默认 `B = 100um`，可通过 `--skew-bound` 修改：

```bash
./build/bst_dme --case sinks_100 --skew-bound 100
# 或
./run_cts.sh --case sinks_100 --topology bst_dme --skew-bound 100
```

---

## 5. 复杂度

- NNG 拓扑：朴素 `O(n^3)`（`n=100` 足够快）；
- DME 自底向上：`O(n)` 次 merge，每次只处理常数条多边形边；
- DME 自顶向下：`O(n)`；
- 内存：`O(n)`。

参考测试只覆盖 16 和 100 sink。当前 NNG 实现限制最多 2000 个 sink；
若要扩展到 10000+，建议将 NNG 替换为空间网格/优先队列版本。

---

## 6. 测试结果

所有测试均使用默认关闭电容感知，`sinks_grid16` 的 16 个 sink 电容均为
20 fF。BST-DME 使用 `--skew-bound 100`。

### 6.1 sinks_grid16

| 拓扑 | sinks | nodes | edges | Manhattan WL (um) | max path (um) | min path (um) | max-min (um) |
|---|---:|---:|---:|---:|---:|---:|---:|
| htree | 16 | 31 | 30 | 4416.666667 | 833.333333 | 750.000000 | 83.333333 |
| zst_dme | 16 | 31 | 30 | 4291.666667 | 750.000000 | 750.000000 | 0.000000 |
| bst_dme | 16 | 31 | 30 | 4301.666667 | 790.000000 | 730.000000 | 60.000000 |

### 6.2 sinks_100

| 拓扑 | sinks | nodes | edges | Manhattan WL (um) | max path (um) | min path (um) | max-min (um) |
|---|---:|---:|---:|---:|---:|---:|---:|
| htree | 100 | 195 | 194 | 14711.230118 | 1267.637265 | 714.660247 | 552.977017 |
| zst_dme | 100 | 199 | 198 | 7880.606106 | 775.274461 | 775.274461 | 0.000000 |
| bst_dme | 100 | 199 | 198 | 7913.688927 | 791.941127 | 743.743640 | 48.197487 |

BST-DME 的实测 skew 为 `48.197487um`，低于设定上界 `100um`。评估器会读取
拓扑中的 `dme.skew_bound_um`，检查 `max_path - min_path` 是否满足该上界；
本页所有 DME case 均通过该检查。

---

## 7. 评估与可视化

两种 DME 输出与本项目其他拓扑共用同一格式，因此可直接复用：

```bash
.venv/bin/python scripts/evaluate_topology.py \
  --topology results/zst_dme_sinks_100.json \
  --case-file data/cases/sinks_100.csv

.venv/bin/python scripts/visualize_topology.py \
  --input results/zst_dme_sinks_100.json --mode both --outdir results
```

DME 的合并段/TRR 只用于嵌入阶段；最终 JSON 中每条边是父到子的抽象边，
可视化时按 `HV` 规则绘制曼哈顿 L 形走线，评估器按曼哈顿长度累计。
