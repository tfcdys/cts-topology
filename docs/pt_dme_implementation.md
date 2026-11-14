# PT-DME 实现说明与测试报告

本文档对应设计说明 [pt-dme.md](pt-dme.md)，记录当前实现、依赖、正确性机制、
并行方式、测试结果和已知限制。

## 1. 新增/修改文件

| 文件 | 作用 |
|---|---|
| `pt_dme.cpp` | PT-DME 完整实现 + `main()` |
| `pt_dme.hpp` | `PtDmeOptions` / `PtDmeResult` / `build_pt_dme()` |
| `pt_dme_tests.cpp` | A-E 回归测试（0–5、1000 random、skewed、collinear、duplicates） |
| `third_party/delaunator.hpp` | Delaunay 三角剖分，vendored 第三方头文件 |
| `third_party/DELAUNATOR_LICENSE` | Delaunator-C++ MIT 许可证 |
| `cts_common.hpp/.cpp` | 公共树结构、JSON 输出、CLI；新增 `Topology::metadata_json` |
| `Makefile` | 新增 `build/pt_dme`、`build/pt_dme_tests`、`make test` |
| `run_cts.sh` | 新增 `--topology pt_dme` |

可执行文件：

```bash
./build/pt_dme --case <case|file.csv> [--out results/pt_dme_<case>.json]
make test
```

## 2. 核心数据结构和接口

```cpp
struct Point { double x, y; };

struct ActiveNode {
    int id;                 // 活动节点 id，用于确定性 tie-break
    Point pos;              // 子树代表位置
    double capacitance;
    double delay_min, delay_max;
    int tree_root;          // 最终拓扑中的子树根节点 id
    int wait_rounds;
};

struct Triangle {
    std::array<int, 3> v;   // active_nodes 索引
    Point root;             // 三元根
    double cost;            // F(T)
};
```

公开接口：

```cpp
PtDmeResult build_pt_dme(const std::vector<Sink>& sinks,
                         const PtDmeOptions& options = {});
```

`PtDmeResult` 包含最终 `Topology`、每轮统计、运行时间和线程数。统计写入
拓扑 JSON 的 `algorithm_metadata` 字段。

## 3. Delaunay 实现与依赖

使用 **Delaunator-C++**：

- 仓库：<https://github.com/delfrrr/delaunator-cpp>
- 许可证：MIT
- 头文件：[third_party/delaunator.hpp](../third_party/delaunator.hpp)
- 许可证：[third_party/DELAUNATOR_LICENSE](../third_party/DELAUNATOR_LICENSE)

流程：

1. 先删除重复坐标（合并为子树），保证输入点唯一；
2. 将 `ActiveNode.pos` 转成 `[x0,y0,x1,y1,...]`；
3. 调用 `delaunator::Delaunator`；
4. 将三角形顶点映射回活动节点索引；
5. 过滤零面积退化三角形；
6. 按 `ActiveNode.id` 规范化三个顶点顺序；
7. 三角形整体按顶点 ID 字典序排序，得到确定性结果。

如果 Delaunay 抛异常或没有有效三角形，进入二元备用合并，不会丢失节点。

## 4. 三元根节点几何

严格按设计文档第 5 节：

1. 在 A/B/C 中按曼哈顿距离选择最近点对；
   - 距离相同按两个 ActiveNode ID 排序后的字典序；
2. 固定 A-B 使用 **HV** 路径；
3. 按 L 型路径累计长度求中点 M，不使用端点坐标平均；
4. 构造 M-C 的 HV 和 VH 两条候选 L 路径，展开为水平/垂直候选线段集合；
5. 对每条候选线段参数化 `P(t)`，最小化
   `F(t)=max(d1(P,A), d1(P,B), d1(P,C))`；
6. 每条线段固定执行 80 次三分搜索，再评估端点；
7. 候选点比较顺序：
   - F 最小；
   - 距 M 最近；
   - x 较小；
   - y 较小。

分支连接固定使用 HV，从 R 连到每个子树的 `tree_root`。

## 5. 极大独立集

流程：

```text
triangles = Delaunay(...)
计算每个三角形 cost 和 root
按 (cost, -max_wait, 顶点ID字典序) 排序
used = [false] * active_count
selected = []
for t in triangles:
    if !used[t.a] && !used[t.b] && !used[t.c]:
        select(t)
        used[t.a]=used[t.b]=used[t.c]=true
```

实现中同时断言：

- selected 中任意两个三角形不共享顶点；
- 任一未入选三角形至少与一个 selected 三角形共享顶点，
  即该集合是极大独立集，而不是要求最大基数。

## 6. Carry-over

每轮：

1. 先按 selected 顺序创建三元根节点；
2. 未被任何 selected 三角形使用的旧 ActiveNode 原样保留；
3. 保留其 `tree_root`、`pos`、电容、延迟区间、身份；
4. `wait_rounds++`；
5. 新三元根 `wait_rounds = 0`；
6. 断言 `N_next = N - 2k`（k 为三元合并数）。

二元备用合并时断言 `N_next = N - 1`。重复坐标预处理的合并数也单独断言。

## 7. 并行执行与线程安全

- 并行阶段 1：对 Delaunay 产生的所有三角形并行计算三元根 `root` 和 `cost`；
- 并行阶段 2：MIS 确定后，对 selected 中每个三角形并行生成局部
  `LocalMergeResult`（根、电容汇总、延迟区间）；
- 每个任务只读公共 `active_nodes`，只写自己的结果槽位，没有共享写；
- 主线程负责 Delaunay、排序、MIS、拓扑节点 ID 分配和边写入；
- 结果按三角形索引写回，线程完成顺序不影响结果；
- selected 顺序、节点 ID 和边顺序完全由确定性排序决定；
- 线程数默认为 `std::thread::hardware_concurrency()`，并限制为不超过
  三角形数量。

因此串行/并行执行产生相同的拓扑。

## 8. 重复坐标与共线退化

### 8.1 重复坐标

- 按坐标排序后分组；
- 同组按 ActiveNode ID 升序；
- 依次执行零长度二元合并；
- 每次合并保留两个子树拓扑；
- 不添加随机扰动；
- 合并后的 ActiveNode 进入下一阶段。

### 8.2 共线

- 用叉积判断所有点是否共线；
- 选择范围较大的轴作为主轴；
- 排序键：`(主轴, 次轴, node_id)`；
- 每三个连续节点组成一个不重叠三元组；
- 剩余 1–2 个节点 carry-over；
- 断言 `N_next = N - 2k`。

## 9. 二元 DME 收尾和备用合并

以下情况调用二元合并：

- 最终只剩两个活动节点；
- MIS 为空；
- Delaunay 无有效三角形；
- 重复坐标零长度合并。

二元合并复用 `dme_common` 的点核合并几何：

- 两子树代表点作为点核；
- 半径各取两者曼哈顿距离的一半；
- 通过 `dme_merge_segment` 求合并段并取其中点作为根；
- 根连接两个子树的 `tree_root`。

第一版只使用位置几何，不传播完整 TRR 或真实 RC 延迟。

## 10. 回归测试

```bash
make test
```

结果保存在：

- `results/pt_dme_tests.log`
- `results/pt_dme_tests_time.txt`

测试覆盖：

| 测试 | 内容 | 结果 |
|---|---|---|
| small_0..5 | 0/1/2/3/4/5 节点终止与二元收尾 | PASS |
| A_random_1000 | 1000 随机 sink | PASS |
| B_skewed_1000 | 900 个密集 + 100 个远离 sink | PASS |
| C_collinear_1000 | 1000 个共线 sink | PASS |
| D_duplicates_1000 | 1000 sink、8 个重复坐标组 | PASS |
| parallel_determinism | 300 sink，1 线程与 4 线程拓扑逐节点/逐边一致 | PASS |

每个测试检查：

- sink 不丢失、不重复；
- 单根、无环、全节点可达；
- 边数 = 节点数 - 1；
- sink 为叶子；
- 每轮活动节点数量关系；
- MIS 顶点独立性和极大性（构建过程断言）。

## 11. 1000 random sink 结果

测试 A 固定 seed `20261114`：

```text
sinks=1000
nodes=1500
edges=1499
Manhattan WL=50091.8 um
max path=1554.59 um
min path=709.937 um
max-min=844.65 um
threads=12
runtime=7.05124 ms
max_wait_rounds=3
```

轮次记录：

| round | active before | triangles | selected | ternary | carry | active after |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1000 | 1978 | 272 | 272 | 184 | 456 |
| 1 | 456 | 891 | 126 | 126 | 78 | 204 |
| 2 | 204 | 393 | 55 | 55 | 39 | 94 |
| 3 | 94 | 175 | 26 | 26 | 16 | 42 |
| 4 | 42 | 71 | 11 | 11 | 9 | 20 |
| 5 | 20 | 29 | 5 | 5 | 5 | 10 |
| 6 | 10 | 11 | 3 | 3 | 1 | 4 |
| 7 | 4 | 2 | 1 | 1 | 1 | 2 |
| final | 2 | — | — | — | — | 1（二元收尾） |

## 12. 其他 1000 sink 测试

| 测试 | WL (um) | max path | min path | max-min | 轮数 |
|---|---:|---:|---:|---:|---:|
| B_skewed_1000 | 8867.94 | 938.767 | 294.986 | 643.78 | 8 |
| C_collinear_1000 | 4550.00 | 644.25 | 280.25 | 364.00 | 6 |
| D_duplicates_1000 | 990.00 | 315.00 | 135.00 | 180.00 | 2 |

## 13. 16/100/1000 标准 case 结果

| Case | sinks | nodes/edges | WL (um) | max path | min path | max-min |
|---|---:|---:|---:|---:|---:|---:|
| sinks_grid16 | 16 | 24 / 23 | 5833.333333 | 1302.083333 | 447.916667 | 854.166667 |
| sinks_100 | 100 | 150 / 149 | 11471.134675 | 1646.126291 | 640.986648 | 1005.139642 |
| sinks_1000 | 1000 | 1500 / 1499 | 43088.004856 | 1934.805124 | 809.555211 | 1125.249913 |

全部通过覆盖率和完整树检查。PT-DME 第一版只做几何拓扑，不优化 skew，
因此路径差较大是预期行为。

## 14. 运行时间与内存

`make test` 整体（含 0–5、A–D 四组 1000-sink、并行确定性测试）在本机：

```text
User time: 0.01 s
System time: 0.01 s
Wall time: 0.02 s
Maximum resident set size: 4400 KB
```

标准 1000-sink CLI 构建（含 Delaunay、拓扑、JSON 写出、评估和可视化）：

```text
elapsed=0.71 s
maxrss=103964 KB
```

其中拓扑构建本身约 5 ms；主要时间在 Python 评估和 Matplotlib 可视化。

## 15. 已知限制

- 第一版不传播完整 DME merging region / TRR 区间；
- 不使用真实 RC、buffer、slew 优化；
- 三元根只在文档定义的候选路径集合上最优，不保证全局最优；
- 每轮重建 Delaunay，尚未做增量优化；
- 二元收尾使用点核几何；
- 拓扑 skew 不保证最优。

这些限制与设计文档第 14 节一致。
