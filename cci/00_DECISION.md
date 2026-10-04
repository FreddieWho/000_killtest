# 决定：KILL 统一新方法主线

2026-09-28。完成本轮kill-test；不进入完整package/论文开发，保留小型工程原型。

- **统计问题存在，但一般设计不是新贡献。** 六类细胞、4951特征的pooled条件置换FPR为15.91%，配对样本检验1.69%；LIANA原生样本分数接同一模型为4.56%。前者是自实现条件差异reference，不是CPDB原生specificity p值失效。
- **IC跨库独立性失败。** 留一供体的count→IC单调预测R²：CellChatDB 0.736，LIANA/CellPhoneDB 0.918，后者ρ≈0.963；去重物理gene后仍约0.918，触发原KILL条件。跨库排名ρ=0.838，不等于完全随机。
- **工程有用，不能补足科学条件。** 原型读入→聚合→C→拟合3.046s；空间全LR dense对照通过，typed kernel移除9,188个支持通道。原型核查进程峰值约994MiB，本轮未满足“时间更少且整体RAM更低”的组合判据。

| Module | Evidence for gap | Collision risk | Statistical validity | Speed ceiling | Biology interpretability | Decision |
|---|---|---|---|---|---|---|
| replicate-aware inference | 已有一般设计 | 高 | 条件null保守/近nominal | 聚合可复用 | 非因果 | KILL：独立统计贡献 |
| scRNA opportunity operator | 仅工程价值 | 高 | paired/factorial已跑 | 聚合后与cell数无关 | mixing假设 | GO：保留工程 |
| spatial opportunity operator | 稀疏计算收益 | 高 | 全LR数值一致 | O(K·E) | spot邻接proxy | GO：保留工程 |
| IC100/IC90 | 跨库独立性失败 | 未完全排除 | MVC/ILP精确；IC90可选未做 | complex最坏NP-hard | 干预未验证 | KILL：通用新统计量 |
| unified package | 不满足统计+IC层GO | 尚有全文缺口 | 不足以升级 | 不能以速度补足 | 未验证 | KILL |

**边界：真实blockade验证NOT_VALIDATED**，仅完成Kang刺激sanity和合成图删除。Spatial CellChat 2.2.0仅20次bootstrap作profiling；scRNA CPDB/LIANA各1000、CellChat100。MOSANIC全文403及矩阵U项保留，不宣称首次或穷尽碰撞；本次KILL不依赖证明文献不存在。

证据：[碰撞](01_collision_matrix.md) · [统计](02_statistical_killtest.md) · [工程](03_engineering_profile.md) · [IC](04_intervention_complexity.md) · [空间](05_spatial_homology.md) · [benchmark](results/benchmark.tsv) · [复现](README.md)。

对PLAN：H2广义设计创新、H3通用独立性未获支持；H4工程同构获支持，H1不支持首次声明。
