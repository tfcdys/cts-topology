# HP-DME：路径感知异构并行 DME 算法实现任务书

## 1. 实现目标

在现有 PT-DME（Parallel Triangular DME）基础上，实现 HP-DME（Heterogeneous Parallel DME）。

本次改进针对 PT-DME 的主要问题：

- 三元合并对几何距离优化较好，但没有充分考虑已有子树的累计路径长度。
- 几何离群点可能连续多轮无法参与合并，最终在算法后期才被合并。
- 离群点晚合并时，已有子树的路径长度差异可能很大，导致最终时钟树的最大路径和路径差值明显增大。
- 只优化当前合并的几何线长，不能保证最终树的路径平衡。

本次实现采用：

**Delaunay + Greedy MIS + Outlier Matching + Four-terminal H-tree + Path-aware Merge Cost**

核心原则：

1. 保留原 PT-DME 的 Delaunay 和 MIS 框架。
2. 利用离群点与三角形组合，增加四终端 H-tree 合并方式。
3. 三元合并和四终端合并统一考虑子树累计路径长度。
4. 保留所有未参与合并的节点，不丢失任何原始 sink。
5. 所有并行合并任务必须满足输入节点互不重叠。
6. 所有决策必须确定性，保证相同输入得到相同拓扑。

---

## 2. 数据结构修改

### 2.1 ActiveNode

在现有 ActiveNode 基础上增加：

```cpp
struct ActiveNode {
    int id;
    Point pos;

    double capacitance;

    // 当前子树根节点到所有原始 sink 的最短、最长累计路径长度
    double path_min;
    double path_max;

    // 当前子树根节点
    int tree_root;

    // 连续等待未参与合并的轮数
    int wait_rounds;
};
```

其中：

- `pos`：当前活动节点代表位置。
- `path_min`：当前子树根到所有原始 sink 的最小累计路径长度。
- `path_max`：当前子树根到所有原始 sink 的最大累计路径长度。
- `wait_rounds`：该活动节点连续多少轮未参与合并。

对于初始 sink：

```cpp
path_min = 0.0;
path_max = 0.0;
wait_rounds = 0;
```

注意：

`path_min` 和 `path_max` 表示从当前活动节点代表位置出发，沿已经构造好的子树到达原始 sink 的累计路径长度。

它们不能在新一轮合并时清零。

这两个量是后续路径平衡优化的核心数据。

### 2.2 MergeCandidate

统一三元合并和四终端合并的结果：

```cpp
struct MergeCandidate {
    Point root;

    double path_min;
    double path_max;

    double wire_length;
    double skew;
    double max_path;
    double cost;

    std::vector<int> child_ids;

    // 保存实际生成的拓扑结构
    TreeTopology topology;

    // 从新根到每个子树代表节点的实际连接路径长度
    std::vector<double> branch_lengths;
};
```

`branch_lengths[i]` 必须对应 `child_ids[i]`。

必须保存真实拓扑结构，不能只保存根坐标和子节点 ID。

---

## 3. 统一的路径感知代价函数

这是 HP-DME 的核心改进。

对于一个待合并的候选拓扑，设：

- 输入子树数量为 k。
- 第 i 个输入子树的路径范围为 `[path_min_i, path_max_i]`。
- 新根到第 i 个子树代表节点的实际拓扑路径长度为 `l_i`。
- 新根位置为 R。

则合并后：

\[
H(R)=\max_{1\leq i\leq k}
\left(path\_max_i+l_i\right)
\]

\[
L(R)=\min_{1\leq i\leq k}
\left(path\_min_i+l_i\right)
\]

\[
S(R)=H(R)-L(R)
\]

其中：

- H：合并后子树的最大 sink 路径长度。
- L：合并后子树的最小 sink 路径长度。
- S：合并后子树的路径差值。

新增拓扑的布线长度为：

\[
W(R)=\text{新增拓扑所有不同物理边的曼哈顿长度之和}
\]

共享边只计算一次。

统一代价函数：

\[
J(R)=
\alpha\frac{H(R)}{H_{ref}}
+
\beta\frac{S(R)}{S_{ref}}
+
\gamma\frac{W(R)}{W_{ref}}
\]

其中：

- `H_ref`：当前输入子树路径长度的参考尺度。
- `S_ref`：当前输入子树路径差值的参考尺度。
- `W_ref`：当前候选新增拓扑布线长度的参考尺度。

参考尺度必须大于零，可以使用对应输入量的最大值与一个很小的正数中的较大值。

第一版权重固定为：

```cpp
alpha = 1.0;
beta  = 1.0;
gamma = 0.2;
```

这些权重是第一版实验参数，不代表理论最优值。应支持通过配置修改。

### 3.1 代价函数实现要求

每个候选拓扑都必须按照以下步骤计算：

1. 生成实际拓扑。
2. 从实际拓扑计算每个输入子树对应的连接路径长度 `l_i`。
3. 计算 H。
4. 计算 L。
5. 计算 S。
6. 计算新增拓扑的实际边长 W。
7. 计算 J。

不能直接假设：

\[
l_i=d_1(R,P_i)
\]

只有当实际生成的连接路径确实是从 R 到 P_i 的最短 Manhattan 路径时，才可以这样计算。

如果 H-tree 存在共享干线、转折点或中间分支，必须沿实际拓扑累计路径长度。

### 3.2 新子树路径范围传播

合并完成后：

\[
path\_max_{new}=H(R)
\]

\[
path\_min_{new}=L(R)
\]

即：

```cpp
new_node.path_max = candidate.max_path;
new_node.path_min =
    candidate.max_path - candidate.skew;
```

也可以直接保存候选计算得到的 L 和 H，避免浮点数重复计算。

特别注意：

不能仅使用新根到子树代表点的距离作为新节点路径长度，也不能丢弃原子树的 `path_min/path_max`。

---

## 4. 每轮 Delaunay 和 Greedy MIS

每轮对所有活动节点执行 Delaunay 三角剖分。

得到三角形集合：

\[
T=\{T_1,T_2,\ldots,T_m\}
\]

每个三角形包含三个不同的活动节点。

两个三角形只要共享任意一个顶点，就视为冲突。

采用原 PT-DME 的确定性贪心 MIS 选择方式。

这里的 MIS 指 maximal independent set，不要求是 maximum independent set。

### 4.1 三角形排序

每个三角形定义排序键：

1. `cost(T)` 升序。
2. 三角形三个顶点中的最大 `wait_rounds` 降序。
3. 将三个顶点 ID 升序排序后，按字典序升序。

第一版的 `cost(T)` 使用三元路径感知合并的最优代价：

\[
cost(T)=\min_{R\in\Gamma_T}J_3(R)
\]

其中：

- \(\Gamma_T\) 是现有三元根点构造器产生的合法候选位置集合。
- \(J_3\) 是本任务定义的统一路径感知代价函数。

如果现有实现暂时不能在 MIS 阶段计算完整代价，可以先使用原 PT-DME 的几何代价排序，但必须保留接口，后续切换为路径感知代价。

### 4.2 贪心 MIS

按照排序后的三角形依次遍历：

- 如果该三角形的三个顶点均未被占用，则选中该三角形。
- 将三个顶点标记为已占用。
- 如果有任意顶点已经被占用，则跳过。

最终得到 Selected 三角形集合。

Selected 中的三角形必须两两顶点不相交。

未被 Selected 覆盖的活动节点集合记为：

\[
U=V-\bigcup_{T\in Selected}V(T)
\]

注意：

U 中的节点不全是几何离群点。

只有满足后续离群点检测条件的节点，才可以作为四终端 H-tree 的额外输入。

其他节点正常 carry-over。

---

## 5. 几何离群点检测

对每个活动节点计算最近邻曼哈顿距离：

\[
d_{nn}(i)=
\min_{j\neq i}d_1(P_i,P_j)
\]

其中：

\[
d_1(P_i,P_j)=
|x_i-x_j|+|y_i-y_j|
\]

计算所有活动节点最近邻距离的中位数：

\[
D_{med}=\operatorname{median}
\{d_{nn}(1),\ldots,d_{nn}(N)\}
\]

定义离群度：

\[
O_i=
\frac{d_{nn}(i)}
{D_{med}+\epsilon}
\]

第一版阈值：

```cpp
outlier_ratio_threshold = 2.5;
```

仅对 U 中的节点执行离群点标记：

\[
Outliers=
\{i\in U\mid O_i>2.5\}
\]

其中 epsilon 是用于避免除零的极小正数。

如果所有活动节点位置完全相同，或者最近邻距离中位数为零，应采用确定性退化处理，不得出现 NaN 或无穷大。

离群点按照以下优先级排序：

1. `outlier_ratio` 降序。
2. `wait_rounds` 降序。
3. 节点 ID 升序。

说明：

这里的离群点定义是几何稀疏程度，不等同于电气离群点。第一版先采用几何定义，不引入额外的统计分类器。

---

## 6. 离群点与 Selected 三角形匹配

这是 HP-DME 的第二个核心改进。

对每个离群点 P，从尚未分配的 Selected 三角形中寻找一个匹配三角形。

若三角形：

\[
T=(A,B,C)
\]

则将：

\[
Q=(A,B,C,P)
\]

作为四终端合并输入。

### 6.1 匹配限制

必须满足：

1. 一个离群点最多匹配一个三角形。
2. 一个 Selected 三角形最多匹配一个离群点。
3. 被匹配的三角形不能再执行普通三元合并。
4. 离群点只能来自 U。
5. 已被其他任务使用的活动节点不能重复分配。
6. 没有合适三角形的离群点继续 carry-over。

### 6.2 点到三角形的几何距离

计算离群点 P 到三角形几何区域的曼哈顿距离：

\[
D(P,T)=
\min_{X\in Triangle(T)}d_1(P,X)
\]

如果 P 在三角形区域内部，则：

\[
D(P,T)=0
\]

否则计算 P 到三条三角形边界线段的最小曼哈顿距离。

该距离仅作为匹配的空间代价，不应单独决定匹配结果。

### 6.3 四终端匹配代价

对每个候选匹配 `(P,T)`，调用四终端 H-tree 优化器，计算：

\[
J_4(P,T)=
\min_{R,\mathcal H\in\mathcal C_4(T,P)}
J(R,\mathcal H)
\]

其中：

- \(\mathcal C_4(T,P)\) 是该四个输入节点能够构造出的合法 H-tree 候选集合。
- J 使用第 3 节定义的统一路径感知代价。
- 每个候选必须计算实际拓扑路径长度。

定义匹配代价：

\[
M(P,T)=J_4(P,T)+\lambda
\frac{D(P,T)}{D_{ref}}
\]

第一版：

```cpp
lambda = 0.2;
```

`D_ref` 为当前轮匹配距离的正参考尺度，可使用候选匹配距离的中位数；若为零，则使用一个正的数值下限。

### 6.4 确定性贪心匹配

按照离群点优先级遍历：

对于当前离群点 P：

1. 遍历所有尚未分配的 Selected 三角形。
2. 计算每个候选 `(P,T)` 的匹配代价 M。
3. 选择 M 最小的三角形。
4. 如果代价相同，按三角形排序后的顶点 ID 字典序选择。
5. 将 P 和 T 标记为已分配。
6. 继续处理下一个离群点。

注意：

这里是确定性贪心一对一匹配，不要求求解全局最优二分图匹配。

如果候选四终端拓扑全部无效，则该离群点不匹配，保留为 carry-over 节点。

---

## 7. 四终端 H-tree 构造与优化

对匹配成功的四个输入节点：

\[
A,B,C,P
\]

构造合法的四终端 H-tree。

### 7.1 候选分组

枚举三种不重复的二二分组：

\[
(A,B)|(C,P)
\]

\[
(A,C)|(B,P)
\]

\[
(A,P)|(B,C)
\]

每种分组都必须调用统一的四终端 H-tree 几何构造器，生成合法候选拓扑。

候选构造器必须输出：

- 根节点 R。
- 中间转折点或分支节点。
- 完整物理边集合。
- 四个输入节点对应的实际连接路径长度。
- 新增拓扑的曼哈顿边长总和。

如果工程已有 HTreeBuilder，优先复用其合法拓扑生成逻辑。

如果已有构造器只支持固定的根位置或固定分组方式，需要扩展它，使其能够生成上述三种分组的候选。

不能把三种分组本身当成完整的 H-tree 几何定义。每个分组仍然必须生成实际可连接、无断边的树形拓扑。

### 7.2 候选位置优化

对于每一种合法 H-tree 拓扑，按照工程 HTreeBuilder 支持的合法几何参数搜索根位置及转折点。

要求：

- 所有连接均为水平或垂直线段。
- 不能产生断开的子树。
- 不能丢失任意输入终端。
- 不能将同一物理边重复计入布线长度。
- 必须保存真实的物理边和分支关系。

如果候选位置的搜索空间为连续区间，第一版可以采用确定性的有限候选搜索：

1. 几何中位位置。
2. BBOX 中心。
3. 输入节点坐标的 X/Y 坐标组合候选。
4. 现有 HTreeBuilder 产生的其他合法候选位置。

对所有合法候选计算 J，选择代价最小者。

第一版不要求解决连续空间的全局最优化问题，但必须保证候选生成确定性，且实际代价计算正确。

### 7.3 最终选择

在三种分组产生的所有合法候选中：

\[
H^*=
\arg\min_{H\in\mathcal C_4}
J(H)
\]

平局规则：

1. `cost` 最小。
2. `max_path` 最小。
3. `skew` 最小。
4. `wire_length` 最小。
5. 根节点 X 坐标最小。
6. 根节点 Y 坐标最小。
7. 拓扑结构 ID 字典序最小。

最终生成一个四终端 H-tree 合并节点。

---

## 8. 三元合并改造

普通 Selected 三角形中，没有匹配离群点的三角形继续执行三元合并。

但必须把原有纯几何根点代价替换为第 3 节定义的路径感知代价。

三元合并继续使用现有 PT-DME 的合法根点候选集合：

\[
\Gamma_T
\]

对每个候选根位置 R：

1. 构造三元实际拓扑。
2. 计算三个输入子树到新根的实际路径长度。
3. 计算新路径最大值 H。
4. 计算新路径最小值 L。
5. 计算路径差值 S。
6. 计算新增布线长度 W。
7. 计算 J。

选择：

\[
R^*=
\arg\min_{R\in\Gamma_T}J_3(R)
\]

平局规则与四终端合并保持一致。

必须保留原有三元拓扑合法性约束，不允许为了降低 J 生成不符合工程拓扑规则的连接。

---

## 9. 每轮合并任务生成

完成离群点匹配后，将本轮任务划分为：

- `FourTerminalTasks`：匹配成功的三角形与离群点。
- `TernaryTasks`：没有匹配离群点的 Selected 三角形。
- `CarryOverNodes`：所有未参与上述任务的活动节点。

设：

- 本轮三元合并数量为 k3。
- 本轮四终端合并数量为 k4。
- 本轮输入活动节点数为 N。

则下一轮活动节点数量必须满足：

\[
\boxed{
N_{next}=N-2k_3-3k_4
}
\]

解释：

- 三元合并消耗 3 个节点，生成 1 个节点，净减少 2 个。
- 四终端合并消耗 4 个节点，生成 1 个节点，净减少 3 个。

所有 carry-over 节点都必须完整保留。

对 carry-over 节点：

```cpp
node.wait_rounds += 1;
```

对参与合并的新节点：

```cpp
new_node.wait_rounds = 0;
```

不能因为离群点匹配失败而丢弃离群点。

---

## 10. 退化情况与兜底

### 10.1 MIS 为空

如果本轮 MIS 为空：

- 不得进入死循环。
- 使用原 PT-DME 的确定性最近邻二元 DME 作为备用。
- 二元合并消耗两个节点，生成一个节点。
- 本轮活动节点数量至少减少 1。

### 10.2 共线节点

如果所有活动节点共线：

- 按主坐标排序。
- 使用连续、不重叠的三元组执行三元合并。
- 剩余 1～2 个节点 carry-over。
- 仍然可以对满足匹配条件的节点尝试四终端合并，但第一版允许先使用上述确定性共线策略，不要求强行构造不合法的二维 Delaunay 三角形。

### 10.3 重复坐标

对于重复坐标：

- 不使用随机扰动。
- 先按节点 ID 确定性处理重复坐标。
- 可以使用零长度二元合并消除重复代表位置。
- 必须保留所有原始 sink 的拓扑归属。

### 10.4 长时间等待

为避免离群点无限 carry-over，增加：

```cpp
max_wait_rounds = 4;
```

当节点等待轮数超过阈值时：

- 提高其离群点处理优先级。
- 如果没有可行四终端合并，则触发确定性二元 DME 兜底。
- 兜底后必须记录该节点的等待轮数和对应合并代价。

不能因为节点等待时间过长而直接丢弃节点。

---

## 11. 并行实现要求

MIS 选中的三角形之间没有共享顶点。

离群点来自 U，且每个离群点和三角形均是一对一匹配。

因此，匹配完成后的四终端任务和普通三元任务应当互不共享输入节点。

可以并行执行：

- `optimize_ternary_merge()`
- `optimize_four_terminal_htree()`

推荐接口：

```cpp
MergeCandidate optimize_ternary_merge(
    const ActiveNode& a,
    const ActiveNode& b,
    const ActiveNode& c
);

MergeCandidate optimize_four_terminal_htree(
    const ActiveNode& a,
    const ActiveNode& b,
    const ActiveNode& c,
    const ActiveNode& p
);
```

每个线程只能生成自己的本地 MergeCandidate。

禁止多个线程同时修改：

- 全局节点 ID。
- 全局拓扑结构。
- 全局活动节点数组。
- 全局边集合。

所有任务完成后，由主线程按照确定性顺序：

1. 分配新节点 ID。
2. 提交拓扑。
3. 写入新节点。
4. 生成下一轮活动节点数组。

保证并行执行不会改变最终拓扑结果。

---

## 12. HP-DME 主流程

实现以下主流程：

```cpp
while (active_nodes.size() > 2) {

    // 1. 处理重复坐标和退化情况
    handle_duplicate_positions(active_nodes);

    // 2. Delaunay triangulation
    auto triangles = delaunay_triangulation(active_nodes);

    // 3. 生成确定性排序
    sort_triangles_deterministically(triangles);

    // 4. Greedy maximal independent set
    auto selected = greedy_mis(triangles);

    // 5. 获取未被 MIS 覆盖的节点
    auto uncovered = get_uncovered_nodes(
        active_nodes, selected
    );

    // 6. 识别几何离群点
    auto outliers = detect_outliers(
        uncovered, active_nodes
    );

    // 7. 离群点和 Selected 三角形一对一匹配
    auto matches = match_outliers_to_triangles(
        outliers, selected
    );

    // 8. 生成四终端任务
    auto four_tasks = build_four_terminal_tasks(matches);

    // 9. 未匹配的 Selected 三角形生成普通三元任务
    auto ternary_tasks = build_ternary_tasks(
        selected, matches
    );

    // 10. 剩余节点 carry-over
    auto carry_over = build_carry_over_nodes(
        active_nodes,
        four_tasks,
        ternary_tasks
    );

    // 11. 并行执行所有合并任务
    auto results = execute_parallel(
        ternary_tasks,
        four_tasks
    );

    // 12. 主线程统一提交结果
    active_nodes = commit_results_deterministically(
        results, carry_over
    );

    // 13. 检查节点数量
    assert(
        active_nodes.size() ==
        previous_size - 2 * ternary_tasks.size()
                       - 3 * four_tasks.size()
    );

    // 14. 检查 sink 覆盖、拓扑连通性和无环性
    validate_tree_invariants();
}

// 最后两个节点执行普通二元 DME
if (active_nodes.size() == 2) {
    binary_dme(active_nodes[0], active_nodes[1]);
}

// 一个节点则直接结束
```

注意：

如果本轮 MIS 为空，必须执行二元备用合并，不得让主循环原地不动。

对于共线和重复坐标等特殊情况，允许进入专门的确定性退化处理流程，但必须保证活动节点数持续减少，或者最终进入正常终止状态。

---

## 13. 代码模块划分

建议将实现拆分为以下模块：

```cpp
// 1. 路径代价计算
MergeCost evaluate_merge_cost(
    const TreeTopology& topology,
    const std::vector<ActiveNode>& children
);

// 2. 三元合并优化
MergeCandidate optimize_ternary_merge(
    const ActiveNode& a,
    const ActiveNode& b,
    const ActiveNode& c
);

// 3. 四终端 H-tree 优化
MergeCandidate optimize_four_terminal_htree(
    const ActiveNode& a,
    const ActiveNode& b,
    const ActiveNode& c,
    const ActiveNode& p
);

// 4. 离群点检测
std::vector<int> detect_outliers(
    const std::vector<ActiveNode>& uncovered,
    const std::vector<ActiveNode>& all_nodes
);

// 5. 点到三角形距离
double manhattan_distance_to_triangle(
    const Point& p,
    const Triangle& triangle
);

// 6. 离群点匹配
std::vector<OutlierMatch> match_outliers_to_triangles(
    const std::vector<int>& outliers,
    const std::vector<Triangle>& selected
);

// 7. 主流程
void run_hp_dme(...);
```

允许根据现有工程代码调整接口，但必须保持上述算法语义。

不要为了实现新算法重复开发已经存在的 Delaunay、DME、HTreeBuilder 或拓扑管理模块。

优先复用原有代码，并通过新增接口扩展功能。

---

## 14. 正确性验证

必须验证以下不变量。

### 14.1 Sink 覆盖

最终树中每个原始 sink 必须且只能出现一次。

\[
\forall s\in Sinks,\quad count(s)=1
\]

不能遗漏 sink，也不能重复连接 sink。

### 14.2 拓扑连通性

最终拓扑必须是一棵连通树。

对于 N 个原始 sink 和 N 个叶节点：

\[
|E|=|V|-1
\]

这里的 V、E 是最终物理拓扑的完整节点和边，不是每轮 ActiveNode 数量。

### 14.3 无环

最终拓扑必须无环。

### 14.4 路径统计正确性

对最终树执行一次独立的 DFS 路径遍历，直接计算：

- 所有 sink 的根到 sink 路径长度。
- 最大路径长度。
- 最小路径长度。
- 路径差值。

将 DFS 计算结果与 HP-DME 在每轮传播的 `path_min/path_max` 进行对比。

两者应在浮点误差范围内一致。

### 14.5 节点数量

每轮必须满足：

\[
N_{next}=N-2k_3-3k_4
\]

如果使用二元兜底，还需要额外减去对应的净节点数量。

---

## 15. 实验对比

必须和以下算法在相同输入上对比：

- H-tree
- ZST-DME
- BST-DME
- 原 PT-DME
- 新 HP-DME

测试规模至少包括：

- 16 sinks。
- 100 sinks。
- 1000 sinks。

还必须包括：

- 均匀网格。
- 随机均匀分布。
- 密集区域与稀疏区域混合。
- 一个或多个明显离群点。
- 多个离群点。
- 全共线。
- 重复坐标。
- 少量节点。

### 15.1 必须统计的指标

| 指标 | 说明 |
|---|---|
| Total Wire Length | 总曼哈顿布线长度 |
| Max Path | 最大 root-to-sink 路径 |
| Min Path | 最小 root-to-sink 路径 |
| Path Skew | 最大路径减最小路径 |
| Runtime | 算法运行时间 |
| Node Count | 最终拓扑节点数 |
| Edge Count | 最终拓扑边数 |
| Ternary Count | 三元合并次数 |
| Four-terminal Count | 四终端 H-tree 合并次数 |
| Outlier Count | 离群点数量 |
| Matched Outlier Count | 成功匹配的离群点数量 |
| Match Rate | 离群点匹配率 |
| Average Wait Rounds | 离群点平均等待轮数 |
| Fallback Count | 二元兜底次数 |

### 15.2 核心验证目标

重点验证：

1. HP-DME 是否降低 PT-DME 的最大路径长度。
2. HP-DME 是否降低 PT-DME 的路径差值。
3. 离群点是否不再集中到算法最后阶段合并。
4. 总布线长度增加是否处于可接受范围。
5. 四终端 H-tree 是否确实改善离群点与三元子树的路径平衡。
6. 相同输入是否得到完全一致的拓扑。
7. 并行执行是否与串行执行产生一致的拓扑和代价结果。

不能只以最大路径改善作为算法成功的标准。如果最大路径显著下降，但总线长大幅增加，也需要明确报告这一代价。

---

## 16. 实现交付要求

完成后必须交付：

1. HP-DME 完整实现。
2. 所有新增数据结构和接口。
3. 路径感知代价函数实现。
4. 离群点检测和匹配实现。
5. 四终端 H-tree 构造与优化实现。
6. 三元合并路径感知改造。
7. 并行任务提交和确定性 ID 分配。
8. 退化输入处理。
9. 单元测试和完整测试报告。
10. 与 PT-DME 的实验对比结果。

必须在报告中明确说明：

- 实际修改了哪些文件。
- 新增了哪些函数。
- 复用了哪些现有模块。
- 哪些部分仍然使用近似搜索。
- 哪些部分没有实现全局最优保证。
- 是否通过所有拓扑正确性测试。

如果现有工程的拓扑结构与本任务书中的接口不同，应先检查现有实现，再做最小范围适配，不得为了套用示例接口而破坏原有拓扑数据结构。

**最终目标：在保持拓扑正确、算法确定性和可并行性的前提下，利用离群点四终端 H-tree 与累计路径感知代价，改善 PT-DME 的最大路径和路径差值，同时控制总布线长度。**