# 工程 profile

## 实测范围

CPU本地运行，Python/R环境见 `results/environment.json`、`results/R_environment.txt`。设置主要数值线程为1；宿主机有其他分析进程，数据未作严格独占机器测速。不报告亚秒差异的统计显著性。原生方法在suite中逐个运行，每个方法处理全部16样本并落盘；部分时段与独立空间/分析进程重叠。

| method      |   samples |   native_API_seconds |   peak_RSS_MiB |   permutations_per_sample |
|:------------|----------:|---------------------:|---------------:|--------------------------:|
| CellPhoneDB |        16 |               493.1  |          488.2 |                      1000 |
| FastCCC     |        16 |                63.04 |          894.6 |                         0 |
| LIANA       |        16 |                52.07 |          511.9 |                      1000 |
| CellChat    |        16 |              1027    |          811.8 |                       100 |

CPDB/LIANA各1000置换；FastCCC Mean/Minimum/Arithmetic；CellChat100 bootstrap、trim=0，全部已测LR特征，包含cofactor。CPDB官方v5与LIANA打包的cellphonedb资源不是同一个版本快照。所有原生API实跑，不以自实现近似代替。

## 原型阶段

| stage                           |   seconds |   peak_rss_mb_process_final |   n_features |
|:--------------------------------|----------:|----------------------------:|-------------:|
| load                            |   0.6327  |                       993.6 |  nan         |
| aggregation                     |   1.797   |                       993.6 |  nan         |
| communication_CellChat          |   0.441   |                       993.6 |    3.155e+04 |
| fit_CellChat                    |   0.1757  |                       993.6 | 9216         |
| communication_LIANA_CellPhoneDB |   0.5359  |                       993.6 |    2.995e+04 |
| fit_LIANA_CellPhoneDB           |   0.00934 |                       993.6 | 9367         |

其中CellChat资源一条读入→聚合→C→拟合链为 **3.046s**；上表峰值列为整个核查进程最终high-water mark，不是逐阶段额外RAM。原始主运行的IC、null和power阶段另见 `results/kang_run_03/profile.tsv`。构图和解算单独计时，并逐图确认最小值与冻结结果一致：

| database          | algorithm         | stage         |   seconds |   n_programs |
|:------------------|:------------------|:--------------|----------:|-------------:|
| CellChat          | exact_MVC         | graph_seconds |   0.06814 |          814 |
| CellChat          | exact_MVC         | solve_seconds |   0.2689  |          814 |
| CellChat          | exact_hitting_set | graph_seconds |   0.02951 |          100 |
| CellChat          | exact_hitting_set | solve_seconds |   0.09978 |          100 |
| LIANA_CellPhoneDB | exact_MVC         | graph_seconds |   0.05772 |          814 |
| LIANA_CellPhoneDB | exact_MVC         | solve_seconds |   0.2437  |          814 |

**速度比较边界：** 原生API包含specificity评估与输出，本原型只直接计算C后做样本层模型；它们不是相同统计输出。该现实multi-sample工作流有明显时间节省空间，但不能写成“相同准确性推断全面快若干倍”。原生工具也可跳过不需要的specificity流程。**本轮未满足“时间更少且整体RAM更低”的组合判据。** 原型核查进程峰值约994MiB，原生API进程约488–895MiB；两侧输入处理范围不同，不能把这组数当公平端到端比较，更不能声称已获得RAM收益。native输入已完成公用的library归一化和LR裁剪，prep内存不能被忽略。

2k/5k/10k/full的aggregation、固定450特征的100次样本拟合和100×199 cell permutation分别计时：`results/scaling_run_01/profile.tsv`。post-aggregation fit不访问细胞矩阵；cell-number scaling不是把新的供体数伪造成更多细胞。50k/100k duplication为可选，未运行。

## 空间同口径计算

| stage                          |   seconds |   peak_rss_mb |   n_spots | mode        | n_lr   | edges   | storage_bytes   | copies   |
|:-------------------------------|----------:|--------------:|----------:|:------------|:-------|:--------|:----------------|:---------|
| sparse_graph                   |  0.004438 |          1065 |      2695 | contact     | —      | 15656.0 | 198656.0        | —        |
| dense_distance_comparator      |  0.0518   |          1065 |      2695 | contact     | —      | —       | 58104200.0      | —        |
| dense_communication_comparator |  2.338    |          1065 |      2695 | contact     | 1864.0 | —       | —               | —        |
| sparse_communication           |  0.2389   |          1065 |      2695 | contact     | 1864.0 | —       | —               | —        |
| sparse_graph                   |  0.006665 |          1065 |      2695 | exponential | —      | 46202.0 | 565208.0        | —        |
| dense_distance_comparator      |  0.04469  |          1065 |      2695 | exponential | —      | —       | 58104200.0      | —        |
| dense_communication_comparator |  1.883    |          1065 |      2695 | exponential | 1864.0 | —       | —               | —        |
| sparse_communication           |  0.2492   |          1065 |      2695 | exponential | 1864.0 | —       | —               | —        |

sparse/dense aggregation使用相同LR向量、相同kernel和相同分组操作，仅矩阵表示不同；覆盖全部LR和全部group pairs。dense距离/权重只存在于此比较脚本，正式`cci_core.opportunity_graph`使用KD-tree+CSR。

Spatial CellChat原生20-bootstrap profiling：85.26s，2695 spots；它的群体距离/概率模型与本原型不同，因此只作现实参照，不作数值等价基线。

## 复杂度与资源限制

- scRNA：读取/归一化/聚合需处理输入非零元素；聚合后C成本主要随samples×celltype pairs×LR变化，sample层拟合与cell数无关。
- spatial：KD-tree建图通常O(N logN+|E|)，每个LR的机会聚合O(|E|)，全资源为O(K|E|)，再加输入/输出成本；并不是无视K的总O(|E|)。半径过大时|E|仍可逼近N²，不能作无条件线性保证。
- simple IC：二部图maximum matching/MVC；复杂体/同基因跨角色的全局靶点为ILP，最坏情形NP-hard。100个实测最优不代表任意规模可快解。
- CSR列反复切片的初版很慢；改为CSC按列读取、复用complex组件和按receiver复用空间传播。最终全输出与冻结数值/完整dense对照一致，见 `results/algebra_checks.json` 和 `results/spatial_run_04/numerical_equivalence.json`。旧run保留且不混入最终profile。

KT3：共享聚合与稀疏图的工程收益成立；未给出所有统计输出完全相同且整体RAM更低的全面优越性证据。对PLAN的影响：H4支持保留小型实现，不能挽救H2/H3创新失败。
