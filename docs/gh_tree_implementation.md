# GH-tree 实现与实验报告

GH-tree 参考 TritonCTS：

- <https://github.com/The-OpenROAD-Project/TritonCTS>
- K. Han, A. B. Kahng, J. Li, "Optimal Generalized H-Tree Topology and
  Buffering for High-Performance and Low-Power Clock Distribution",
  IEEE TCAD, 2018.

参考仓库中的 `src/ghtree` 包含 generalized H-tree、K-means clustering、
buffer/latency 优化和布线。本项目的 GH-tree 是第一版**纯拓扑实现**，复用其
“四区域聚类 + 递归广义 H-tree”思想，不包含：

- buffer insertion；
- RC/Elmore 延迟优化；
- LUT 驱动的 buffer 选择；
- blockage 规避；
- TritonCTS 的 Verilog/DEF/LEF 解析。

## 1. 新增文件

| 文件 | 说明 |
|---|---|
| [gh_tree.cpp](../gh_tree.cpp) | GH-tree 完整实现 + `main()` |
| [gh_tree_tests.cpp](../gh_tree_tests.cpp) | 0–5、1000 random/collinear/duplicate、确定性测试 |
| [Makefile](../Makefile) | `build/gh_tree`、`build/gh_tree_tests` |
| [run_cts.sh](../run_cts.sh) | `--topology gh_tree`；GH 专用输出策略 |

## 2. 算法

### 2.1 四区域 K-means 聚类

每个递归组：

1. 先按 X/Y 中位数得到四个初始区域，保证非空；
2. 计算四个区域的加权 centroid；
3. 用 Manhattan 距离做 8 次 K-means 迭代；
4. 若出现空 cluster，或最大 cluster 超过 `0.75N`，回退到初始四叉划分；
5. 每个 cluster 递归构造子树。

默认关闭电容感知；开启 `--cap-aware` 后，centroid 使用 sink 电容作为权重。

### 2.2 广义 H 连接

当前组的根节点位于包围盒中心。四个 cluster 的根节点按 Y 坐标分为：

- 上侧 group；
- 下侧 group。

每侧共享一条 H 臂：

- 根 → 侧 bend 为竖直线；
- bend → 子根为 `HV` 曼哈顿连接；
- 一侧可有 1–4 个子根，因此支持非标准四叉分布；
- 子 cluster 继续递归。

这与经典 H-tree 的固定四象限 H 骨架相似，但区域由 K-means 聚类得到，因此称为
generalized H-tree。

### 2.3 复杂度

- 每层 K-means：`O(N * 4 * 8) = O(N)`；
- 每层排序：`O(N log N)`；
- 树高期望 `O(log N)`，总体 `O(N log² N)`；
- 1000 sink 实测约 1 ms。

## 3. 输出策略

按需求，GH-tree 的评估结果**只输出**：

- `results/gh_tree_<case>_manhattan.png`
- `results/gh_tree_<case>_report.txt`

不输出：

- 欧式连线图；
- `.eval.json`；
- 其他评估图片。

`run_cts.sh` 会自动将 GH-tree 的 mode 固定为 `manhattan`，并且不传
`--json-out`。报告.txt 由评估器生成，包含：

```text
Algorithm runtime    : x.xxx ms
Sink coverage
Complete clock tree
Total Manhattan wirelength
Max root-to-sink path
Min root-to-sink path
Max - min path difference
Average root-to-sink path
```

## 4. 运行

```bash
make
make test

./build/gh_tree --case sinks_grid16
./build/gh_tree --case sinks_100
./build/gh_tree --case sinks_1000

./run_cts.sh --case sinks_1000 --topology gh_tree
```

## 5. 实验结果

所有结果 PASS（覆盖率、完整树、单根、无环、全节点可达）。

| Case | Sinks | Nodes/Edges | Manhattan WL (um) | Max path (um) | Min path (um) | Max-min (um) | Runtime (ms) |
|---|---:|---:|---:|---:|---:|---:|---:|
| sinks_grid16 | 16 | 31 / 30 | 4416.666667 | 833.333333 | 750.000000 | 83.333333 | 0.011 |
| sinks_100 | 100 | 247 / 246 | 9728.441502 | 1121.131350 | 480.547590 | 640.583760 | 0.076 |
| sinks_1000 | 1000 | 2485 / 2484 | 40033.837041 | 1024.868992 | 788.417678 | 236.451314 | 0.981 |

运行时间来自五组实验的中位数；单次报告.txt 记录的是该次拓扑构建器的实际耗时。

### 5.1 对比

1000 sink：

| 算法 | WL (um) | Max path | Skew |
|---|---:|---:|---:|
| H-tree | 52255.853 | 1175.999 | 311.564 |
| GH-tree | 40033.837 | 1024.869 | 236.451 |
| PT-DME | 43088.005 | 1934.805 | 1125.250 |
| HP-DME | 56978.416 | 1254.971 | 439.999 |
| ZST-DME | 34536.295 | 924.139 | 0.000 |
| BST-DME | 34075.794 | 949.139 | 83.638 |

GH-tree 在 H-tree/PT-DME/HP-DME 中具有最低的 1000-sink 线长、最大路径和 skew；
ZST-DME/BST-DME 在 skew 上仍更优，但运行时间约 90ms。

## 6. 测试

```bash
make test
```

GH-tree 测试：

- 0/1/2/3/4/5 节点；
- 1000 random；
- 1000 collinear；
- 1000 duplicate；
- 重复构建确定性。

日志：

- `results/gh_tree_tests.log`
- `results/gh_tree_tests_time.txt`

结果目录不提交，运行 `make test` 可在本地生成。

## 7. 参考与限制

- 参考 TritonCTS 的 generalized H-tree/K-means 思路；
- 未复制参考仓库代码；
- 不实现 buffer、RC、slew、blockage；
- K-means 为确定性 capacitated 初始化 + 有限迭代；
- 全共线/重复坐标走确定性回退；
- 输出只保留 GH 要求的 Manhattan 图和报告.txt。
