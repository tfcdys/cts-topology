# SALT 与 SLLT 实现报告

参考：

- <https://github.com/OSCC-Project/iEDA/tree/master/src/operation/iCTS>
- iEDA 中的 `src/third_party/salt`
- <https://github.com/chengengjie/salt>
- <https://github.com/OpenROAD/OpenROAD-flow-scripts> 中使用的 FLUTE/FLUTE3

两种拓扑：

| 拓扑 | 含义 | 可执行文件 |
|---|---|---|
| `salt` | SALT shallow-light Steiner tree | `build/salt` |
| `sllt` | BST-DME 与 SALT 融合 | `build/sllt` |

## 1. 新增文件

| 文件 | 说明 |
|---|---|
| [salt.cpp](../salt.cpp) | SALT 独立入口 |
| [sllt.cpp](../sllt.cpp) | SLLT 独立入口 |
| [salt_common.hpp](../salt_common.hpp) / [salt_common.cpp](../salt_common.cpp) | SALT/SLLT builder、树转换、SLLT 融合 |
| [salt_tests.cpp](../salt_tests.cpp) | SALT/SLLT 回归测试 |
| [third_party/salt/](../third_party/salt/) | vendored SALT（MIT） |
| [third_party/flute3/](../third_party/flute3/) | vendored FLUTE3（BSD） |
| [third_party/SALT_LICENSE](../third_party/SALT_LICENSE) | SALT 许可证 |
| [third_party/FLUTE_LICENSE](../third_party/FLUTE_LICENSE) | FLUTE 许可证 |

## 2. SALT

SALT（Steiner shAllow-Light Tree）在最小线长 Steiner 树和最短路径树之间折中。
实现直接复用 SALT 的 `SaltBuilder`：

1. FLUTE 构造 rectilinear SMT；
2. SALT 对 SMT 做 shallow-light DFS；
3. RSA 替换 source 子树；
4. `Refine::cancelIntersect`、`flip`、`uShift`、`substitute` 后处理；
5. 输出带 Steiner 点的树。

本项目额外处理：

- 坐标从 `double um` 量化到 `int nm`；
- source 取 sink 包围盒中心；
- 重复 sink 坐标加 1nm 级偏移，避免 FLUTE 合并；
- SALT 树可能穿过 sink pin，转换时把 sink 保持为叶节点，用同坐标内部 clone 继续子树；
- 极端情况下 SALT/FLUTE 仍可能折叠重复 pin，转换后把缺失 sink 重新挂到最近的非 sink 路由节点。

运行：

```bash
./build/salt --case sinks_1000
./run_cts.sh --case sinks_1000 --topology salt
```

## 3. SLLT

SLLT 参考 iCTS 的 `bstSaltTree`：

1. 用项目已有 BST-DME 构造 bounded-skew 树；
2. 将 BST 树转换成 SALT `Tree`/`Net`；
3. 用 iCTS `TreeSaltBuilder` 的 SALT/CBS 思路做 shallow-light 精化；
4. 输出 SALT 树拓扑。

与纯 SALT 的区别是初始树来自 BST-DME，而不是 FLUTE SMT。SLLT 的
shallow-light 精化会改变 BST 的原始路径长度，因此不保证保持 BST 的 skew bound；
这是 shallow-light 融合的取舍。

运行：

```bash
./build/sllt --case sinks_1000
./run_cts.sh --case sinks_1000 --topology sllt
```

## 4. 输出策略

按实验要求，`salt` 和 `sllt` 的评估结果**只输出**：

- `results/<topology>_<case>_manhattan.png`
- `results/<topology>_<case>_report.txt`

`run_cts.sh` 会自动：

- 强制 Manhattan 模式；
- 不生成 `.eval.json`；
- 不生成欧式图；
- 在 report.txt 中记录 `Algorithm runtime`。

## 5. 实验结果

| 拓扑 | Case | Sinks | Nodes/Edges | WL (um) | Max path | Min path | Max-min | Runtime (ms) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SALT | sinks_grid16 | 16 | 30/29 | 4166.669000 | 750.000667 | 416.666667 | 333.334000 | 25.923 |
| SALT | sinks_100 | 100 | 197/196 | 6501.058730 | 871.419044 | 17.974083 | 853.444961 | 27.768 |
| SALT | sinks_1000 | 1000 | 1997/1996 | 27634.600534 | 988.707000 | 16.341702 | 972.365298 | 63.093 |
| SLLT | sinks_grid16 | 16 | 30/29 | 4250.002000 | 770.000667 | 396.666667 | 373.334000 | 0.166 |
| SLLT | sinks_100 | 100 | 197/196 | 6293.011028 | 775.275509 | 8.339579 | 766.935929 | 1.358 |
| SLLT | sinks_1000 | 1000 | 1997/1996 | 27828.323578 | 929.139000 | 34.875828 | 894.263172 | 109.219 |

所有结果 PASS：覆盖率、完整树、单根、无环、全节点可达、sink 为叶节点。
表中 runtime 为 5 次实验的中位数；report.txt 记录单次实际耗时。

说明：

- SALT 的线长明显低于 H-tree/PT-DME/HP-DME，但 shallow-light 目标不保证低 skew；
- SALT 的 runtime 包含一次 FLUTE LUT 读取（约 25ms）；
- SLLT 不调用 FLUTE，因此小规模更快，但 1000 sink 的 SALT substitute 较慢；
- SLLT 由 BST 种子和 SALT/CBS 精化产生，路径指标与纯 SALT 有差异，但
  shallow-light 精化仍不保证 BST 的 skew bound。

## 6. 测试

```bash
make test
```

覆盖：

- 0/1/2/3/4/5 节点；
- 1000 random；
- 1000 collinear；
- 1000 duplicate；
- SALT/SLLT 重复构建确定性。

日志：

- `results/salt_tests.log`
- `results/salt_tests_time.txt`

本机：

```text
Wall time: 0.43 s
Maximum resident set size: 49644 KB
```

## 7. 依赖与限制

- SALT：MIT；
- FLUTE3：BSD 3-Clause；
- SALT refine 使用 Boost.Geometry 头文件；
- FLUTE LUT 路径固定为 `third_party/flute3/etc/POWV9.dat` 和 `POST9.dat`，
  因此二进制需要在项目根目录运行；`run_cts.sh` 会先 `cd` 到项目根目录；
- 本实现只做拓扑，不做 buffer/RC/slew/blockage；
- SALT 的重复 pin 回退是工程性容错，不改变原算法主体。
