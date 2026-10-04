# 决策记录

## D001 · 2026-09-28 · 先审计后编码
遵循原指令，先完成近期方法碰撞审计，再实现最小 Python 原型。环境已发现 numpy/scipy/anndata/scanpy/liana/statsmodels/networkx 和 R；安装成功与版本、实际运行可用性尚待核对。方法证据不能由包是否安装代替。

## D002 · 2026-09-28 · 比较同一统计问题
CellPhoneDB 原生 specificity 与条件差异 null 不同。报告分开 specificity 输出和样本级 score 差异；不将其原生 p 值的条件 null 拒绝率称为失准。LIANA interaction_pvalue 也不自动视为联合检验。依据见 01_collision_matrix.md。

## D003 · 2026-09-28 · 广义设计创新已碰撞，继续规定 kill-test
MultiNicheNet 支持 paired/multifactorial 设计，广义一般设计首创主张 KILL。按控制文件授权继续评估直接 C 估计对象、工程和 IC；不修改 PLAN。MOSANIC 全文与同义最小集合仍保留未知。

## D004 · 2026-09-28 · 复合体与分子身份
复合体先在每个细胞取必需亚基表达的最小值，再聚合到样本；两模态一致。先取亚基均值再取最小值不满足 all-to-all Ω 的同构，因此 kang_run_02 被 kang_run_03 替代。若同一基因同时是 ligand/receptor，另用去重基因 ILP 检查物理分子靶点数；角色分离 MVC 不自动等于全局分子干预数。

## D005 · 2026-09-28 · Visium 距离单位和注释匹配
官方 v1 名义 spot 直径55µm、中心距100µm，但此文件 JSON 直径配55µm得到近邻84.216µm，教程配65µm约99.5µm。采用观测近邻像素间距与官方100µm中心距标定，保留两种直径推算作为冲突记录；spatial_run_04 为规范版本。2695 spots 中匹配作者注释1072，另1个作者条码不在当前官方 filtered counts；1623 spots 保留 UNANNOTATED。

## D006 · 2026-09-28 · 明示统计范围
2k 起始规模及每组至少5细胞的共同集合只剩两种细胞类型、450个非恒定 simple LR 特征。该表仅用于同特征规模曲线；另完整运行六主要细胞类型、4951特征的全细胞 null，不能用450特征曲线替代六类型主实验。

## D007 · 2026-09-28 · 原生基线参数
CPDB与LIANA各样本1000次置换，FastCCC Mean/Minimum/Arithmetic。scRNA CellChat使用默认100次bootstrap、trim=0和全部已测LR特征；并补入全部已测cofactor，不能因先做LR子集而漏掉调节分子。独立 Spatial CellChat 仅作20次bootstrap profiling及语义参照，其 specificity 不用于条件差异FPR。

## D008 · 2026-09-28 · 本轮结束，统一新方法KILL
按原kill规则，替代资源count→IC的留一供体单调预测R²约0.918且近单调；去重物理gene后失败仍成立。现有一般设计路线也未形成可支持独立统计方法的优势。因此不进入统一package/论文开发，保留已运行的小型工程原型。该决定不是证明IC在所有图上无意义，也不声称所有文献碰撞已排除。真实blockade验证NOT_VALIDATED；整体RAM优势未获证明。证据与模块表见00_DECISION.md，逐项工作验收见results/completion_audit.tsv。

## D009 · 2026-09-30 · 优先消除可测量的工程内存阻塞
用户授权继续优化。使用全基因library总量归一化、按两库所需基因联合裁剪、分块读取；保持cell-level complex minimum与样本模型。比较独立进程中的相同输出，历史报告与冻结数值不覆盖。只有工程同口径证据可更新，H2/H3失败不因性能改写。

## D010 · 2026-09-30 · 推荐4096行分块，保留原默认接口
4096行块实测429MiB/8.81s，原版996MiB/11.60s；1024行块391MiB但13.43s。推荐4096行作为stream模式默认，原run_kang默认全基因路径仍保留供全基因统计需求使用。全部分数/模型和完整Kang结果通过1e-12等价检查；H4工程证据增强，D008科学KILL不变。见06_engineering_optimization.md。

## D011 · 2026-10-02 · 独立审阅保留停止决定与解释边界
独立重建1628个简单图并重算两库留一供体R²，未发现推翻D008的计算错误。总体R²不等于精确靶数可替代，条件随机化不等于独立生物校准，复合体完整独立性尚未检验。维持停止决定；本次只交付审阅和建议，不修改PLAN或历史报告，不启动后续科研分支。证据：reviews/20261002/independent_review.md、metrics.json。
