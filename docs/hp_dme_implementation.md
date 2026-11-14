# HP-DME 实现与 PT-DME 对比报告

本文档对应 [hp_dme_opti.md](hp_dme_opti.md)。HP-DME 在 PT-DME 的
Delaunay + Greedy MIS 框架上增加：

- 累计路径区间 `path_min/path_max`；
- 路径感知统一代价 `J`；
- 几何离群点检测；
- 离群点与选中三角形的一对一匹配；
- 四终端共享 H 臂候选；
- 路径感知三元/二元合并；
- 确定性的并行候选计算和主线程拓扑提交。

优化前 PT-DME 的结果已完整保留在
`results/baseline_pt_dme/`（本地生成，不提交）。

## 1. 新增/修改文件

| 文件 | 说明 |
|---|---|
| [hp_dme.cpp](../hp_dme.cpp) | HP-DME 完整实现 + `main()` |
| [hp_dme.hpp](../hp_dme.hpp) | `HpDmeOptions`、`HpDmeResult`、`build_hp_dme()` |
| [hp_dme_tests.cpp](../hp_dme_tests.cpp) | 0–5、1000 random/skewed/collinear/duplicate/outlier、并行确定性测试 |
| [Makefile](../Makefile) | `build/hp_dme`、`build/hp_dme_tests`、`make test` |
| [run_cts.sh](../run_cts.sh) | `--topology hp_dme` |
| [scripts/run_experiments.py](../scripts/run_experiments.py) | 加入 HP-DME 并生成 PT/HP 对比报告 |
| [docs/hp_dme_implementation.md](hp_dme_implementation.md) | 本报告 |

复用的现有模块：

- `dme_common`：Manhattan / TRR / 二元合并基元；
- `third_party/delaunator.hpp`：Delaunay；
- `cts_common`：树结构、JSON、CLI、运行时字段。

## 2. 数据结构

```cpp
struct HpActiveNode {
    int id;
    Point pos;
    double capacitance;
    double path_min;   // 当前子树根到 sink 的最小累计路径
    double path_max;   // 当前子树根到 sink 的最大累计路径
    int tree_root;
    int wait_rounds;
};
```

`MergeCandidate` 保存实际局部拓扑：

```cpp
struct LocalNode { Point pos; int parent; };

struct MergeCandidate {
    Point root;
    std::vector<LocalNode> local_nodes;  // [0] = 新根
    std::vector<int> children;           // active_nodes 索引
    std::vector<int> attachments;        // 每个 child 挂接的 local node
    std::vector<double> branch_lengths;  // root -> child 实际路径长度
    double path_min, path_max;
    double wire_length, skew, max_path, cost;
    int serial;
    bool valid;
};
```

主线程提交候选时才分配全局拓扑 ID，因此并行候选不会修改共享拓扑。

## 3. 路径感知代价

对候选 merge：

\[
H=\max_i(path\_max_i + l_i),\quad
L=\min_i(path\_min_i + l_i),\quad
S=H-L
\]

其中 `l_i` 是候选实际局部拓扑从新根到第 i 个子树根的长度，而不是简单假设
`l_i=d1(R,P_i)`。

参考尺度：

```text
H_ref = max(1, max_i path_max_i)
S_ref = max(1, max_i (path_max_i - path_min_i))
W_ref = max(1, bbox_width + bbox_height)
```

统一代价：

\[
J = \alpha H/H_{ref} + \beta S/S_{ref} + \gamma W/W_{ref}
\]

默认：

```text
alpha = 1.0
beta  = 1.0
gamma = 0.2
lambda = 0.2
outlier_ratio_threshold = 2.5
max_wait_rounds = 4
```

这些参数可通过 `HpDmeOptions` 修改。

## 4. 路径感知三元合并

保留 PT-DME 的三元候选路径集合：

1. 最近点对 A-B；
2. A-B 的 HV 折线中点 M；
3. M-C 的 HV/VH 候选线段；
4. 每条线段采样端点、32 个等距点，以及子节点 x/y 对齐的临界点；
5. 每个候选构造“根 + 三条直接 Manhattan 分支”的实际局部拓扑；
6. 计算 H/L/S/W/J；
7. 按 J、max_path、skew、wire_length、根 x/y、候选序号排序。

三元候选同时用于三角形排序代价和最终提交。

## 5. 离群点检测与匹配

最近邻距离：

\[
d_{nn}(i)=\min_{j\neq i}d_1(P_i,P_j)
\]

中位数 `D_med`，离群度 `O_i=d_nn/(D_med+eps)`，阈值 `2.5`。
只对未被 MIS 覆盖的节点检测。离群点排序：

1. 离群度降序；
2. `wait_rounds` 降序；
3. 节点 ID 升序。

匹配：

- 离群点到三角形区域的曼哈顿距离；
- 对每个候选匹配构造四终端 H-tree，计算 `J4`；
- 匹配代价 `M = J4 + lambda * D / D_ref`；
- 按离群点优先级串行贪心一对一匹配；
- 平局按三角形顶点 ID 字典序。

匹配前，所有 `(outlier, selected triangle)` 候选的 `D` 和 `J4` 由线程池并行
预计算，主线程只做确定性贪心分配。

## 6. 四终端 H-tree

四终端输入为三角形 `(A,B,C)` 加离群点 `P`。

枚举三种二二分组：

```text
(A,B) | (C,P)
(A,C) | (B,P)
(A,P) | (B,C)
```

每种分组枚举两种上下方向，构造共享 H 臂：

- 两对终点分别形成顶部/底部 H 臂；
- 选择 trunk x 的有限候选：bbox 中心、均值、中位、各终端 x；
- 选择 root y 的有限候选：bbox 中心、均值、y_top、y_bottom、臂中点；
- H 臂连接点 `(clamp(x0, pair_x_range), y_pair)`；
- 根到 H 臂连接点是实际局部边；
- 连接点到四个终端的边也是实际局部边；
- 相同坐标的局部节点会去重；
- `W` 是全部局部边长度之和，共享边只计一次。

三种分组 × 两种方向 × 有限根候选全部计算 J，按最终平局规则选择。

## 7. 并行与确定性

并行阶段：

- 所有 Delaunay 三角形的路径感知三元候选；
- 所有离群点-选中三角形匹配单元的四终端候选。

线程只读 `active_nodes`，写自己的 `Triangle` 或 `MatchCell` 槽位。
主线程负责：

- MIS；
- 离群点一对一匹配；
- 新活动节点 ID；
- 局部拓扑节点/边 ID；
- carry-over；
- 路径范围传播。

已验证 1 线程与 4 线程生成逐节点、逐边一致的拓扑。

## 8. 退化处理

- 重复坐标：按 ID 升序做路径感知零长度二元合并；
- 全部共线：主坐标排序，每三个连续节点三元合并，剩余 carry-over；
- MIS 为空 / Delaunay 无三角形：最近点对路径感知二元兜底；
- 长等待：未匹配且 `wait_rounds >= 4` 的离群点，在存在至少两个未使用节点时，
  选择一个最近未使用节点做二元兜底；
- 最终两个节点：路径感知二元收尾。

## 9. 正确性验证

实现内部断言：

- 每轮 `N_next = N - 2k3 - 3k4 - fallback`；
- 重复坐标合并数量关系；
- 共线三元合并数量关系；
- selected 三角形两两不共享顶点；
- MIS 极大性；
- 最终 `|E|=|V|-1`；
- 根可达所有节点；
- 最终 DFS 路径 `min/max` 与传播的 `path_min/path_max` 一致（误差 1e-6）。

## 10. 测试

```bash
make test
```

HP-DME 测试结果：

- `small_0..5`：PASS
- `A_random_1000`：PASS
- `B_skewed_1000`：PASS
- `C_collinear_1000`：PASS
- `D_duplicates_1000`：PASS
- `E_outliers`：PASS
- `parallel_determinism`：PASS

测试日志：

- `results/hp_dme_tests.log`
- `results/hp_dme_tests_time.txt`

本机整组测试：

```text
Wall time: 0.05 s
Maximum resident set size: 8556 KB
```

## 11. PT-DME vs HP-DME

对比使用优化前基线
`results/baseline_pt_dme/experiment_summary_before_hp.csv`。

| Case | PT WL | HP WL | WL 变化 | PT Max | HP Max | Max 变化 | PT Skew | HP Skew | Skew 变化 | PT Runtime | HP Runtime |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sinks_grid16 | 5833.333 | 6666.667 | +14.29% | 1302.083 | 756.348 | -41.91% | 854.167 | 96.029 | -88.76% | 0.855 | 0.750 |
| sinks_100 | 11471.135 | 15814.219 | +37.86% | 1646.126 | 1084.992 | -34.09% | 1005.140 | 614.621 | -38.85% | 2.238 | 2.223 |
| sinks_1000 | 43088.005 | 56978.416 | +32.24% | 1934.805 | 1254.971 | -35.14% | 1125.250 | 439.999 | -60.90% | 5.401 | 22.824 |

结论：

- 最大路径降低 34%–42%；
- 路径差降低 39%–89%；
- 总线长增加 14%–38%；
- 1000 sink 的 HP 运行时间约 22.8ms，仍低于 ZST-DME/BST-DME 的约 90ms；
- 离群点匹配率在 100/1000 标准 case 上为 100%，平均等待轮数为 0；
- 四终端合并数：100 sink 为 10，1000 sink 为 42。

完整对比文件：

- `results/hp_vs_pt_comparison.md`
- `results/hp_vs_pt_comparison.csv`
- `results/experiment_summary.md`

## 12. 近似与未实现部分

- 三元根使用有限采样 + 临界点，不保证连续空间全局最优；
- 四终端 H-tree 使用有限 root/trunk 候选，不保证全局最优；
- 四终端 H-tree 的 pair 上下方向、H 臂高度采用有限候选；
- 二元收尾沿 L1 路径采样 65 点，不保证全局最优；
- 离群点匹配是确定性贪心，不是全局最优二分图匹配；
- 路径代价只使用线长 proxy，不包含真实 RC/Elmore 延迟；
- 不做 buffer 插入、slew、障碍规避；
- 并行只覆盖候选几何计算；拓扑提交、ID 分配和 carry-over 由主线程串行；
- 长等待兜底只处理每轮一个最高优先级离群点。

## 13. 关键接口

```cpp
HpDmeResult build_hp_dme(
    const std::vector<Sink>& sinks,
    const HpDmeOptions& options = {});
```

主要内部函数：

```cpp
MergeCandidate optimize_ternary_merge(...);
MergeCandidate optimize_four_terminal_htree(...);
MergeCandidate optimize_binary_merge(...);
std::vector<int> detect_outliers(...);
double point_to_triangle_distance(...);
```

所有最终拓扑仍使用项目统一的 rooted-tree JSON，可直接复用评估器和可视化。
