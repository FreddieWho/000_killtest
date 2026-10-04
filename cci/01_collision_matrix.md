# Collision audit · 2026-09-28

检索窗口：2024-01-01 至 2026-09-28；按任务补查窗口前的原始方法。来源为原始论文、作者仓库、官方文档及已安装 LIANA 1.6.1 源码。检索未发现不等于证明不存在；未核实的能力不写成“不支持”。尚无证据触发五条件同时满足的整体 collision KILL。

## 能力矩阵

Y=有直接来源；P=部分/外部工作流；D=可由输出汇总，非独立方法贡献；NR=本次所查来源未报告；U=未充分核实。design 指样本层设计；时间排序不自动等于重复测量推断。最后两列区分单靶 knockout 与覆盖全部 LR 通道的最小分子集合。

| Method | scRNA | spatial | biological replicates | arbitrary design matrix | paired/longitudinal | analytical/no cell permutation | count | strength | heterogeneity | perturbation/vulnerability | minimum intervention set |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CellPhoneDB v5 [1] | Y | P：microenvironment | P：外部 DEG | P：外部 DEG | P：外部 DEG | P：非统计/DEG；specificity 用置换 | D | Y | P | P：CellSign 下游 TF | NR |
| CellChat v2 [2] | Y | Y | P：样本比较/空间聚合 | NR | NR | NR：置换 | Y | Y | Y：网络模式 | P：比较条件 | NR |
| CellChat v3 / Spatial CellChat [3] | Y | Y：单细胞空间 | U：不能由版本名推断 | U | U | U | Y | Y | Y | U | NR |
| FastCCC [4] | Y | NR | P：reference/condition comparison | NR | NR | Y：FFT 分布 | D | Y | P：多个评分 | P：reference deviation | NR |
| LIANA+ [5] | Y | Y | Y：pseudobulk/by_sample | Y：基因级外部 DEA | Y：DEA 设计 | P：DEA/部分 LR 方法 | D | Y | Y：多视图/因子 | Y：因果子网；非最小断网 | NR |
| MultiNicheNet ≥2 [6] | Y | P：可补充空间信息 | Y：muscat pseudobulk | Y：协变量/contrasts | Y：paired/多因素教程 | Y：pseudobulk DE | D | Y：prioritization | Y：样本差异 | Y：下游 target/干预条件 | NR |
| Tensor-cell2cell v1/v2 [7] | Y | P：可接受空间 CCC 分数 | Y：context 轴 | NR：因子模型不等于设计检验 | P：time 轴；非重复测量误差模型 | Y：分解本身 | D | Y：输入分数 | Y：耦合程序 | P：跨 context | NR |
| dominoSignal/DCST [8] | Y | NR | Y：subject presence | NR：Fisher exact；GLM 为未来方向 | NR：不可视为配对检验 | Y：subject contingency | Y：linkage occurrence | P：上游；DCST 二值 | Y | Y：治疗队列差异 | NR |
| scDiffCom [9] | Y | NR | NR：原生置换单位为细胞 | NR | NR | NR：cell label/condition shuffle | Y | Y | P | P：差异条件 | NR |
| SpatialDM [10] | NR：需坐标 | Y | Y：跨样本 global z | P：covariate LRT | P：time/covariate；配对误差 U | Y：解析 z；置换可选 | D | P：空间共表达 | Y | P：差异空间信号 | NR |
| COMMOT [11] | P：需真实/外部空间坐标 | Y | P：下游比较 | NR | NR | P：OT 分数；显著性可置换 | D | Y | Y | P：下游基因关联 | NR |
| scRICH [12] | Y | Y | NR：metacells 不是供体重复 | NR | NR | P：动态模型 | D | Y | Y：核心目标 | Y：EGF 共表达实验 | NR |
| Renoir [13] | P：参考表达 | Y | NR | NR | NR | P：neighborhood activity | D | Y：ligand-target | Y：niches | P：靶点活性排序 | NR |
| MOSANIC [14] | NR：核心需空间图 | Y | NR | NR | NR | Y：训练后 readout，不等于解析检验 | D | Y：attention | Y | Y：hub / in-silico knockout | NR：公开接口未见 exact cover |
| scSeqCommDiff [15] | Y | P：空间验证 | Y：multi-sample 模式 | U：两条件模式已核实 | U | P：sample score statistics | D | Y | Y | P：下游响应 | NR |
| receptor-target hitting set [16]（2022 先例） | Y | NR | Y：患者优化，非推断单位 | NR | NR | Y：ILP | NR | NR | Y：肿瘤细胞覆盖 | Y：靶点组合 | Y：覆盖肿瘤细胞，**并非 LR 通道** |

## 对剩余 gap 的约束

1. **一般设计不是空白。** MultiNicheNet 官方文档包含样本级配对、batch、时间×条件 contrasts。LIANA+ 的 pseudobulk DEA 也可接一般设计。不能把“细胞不是 biological replicate”或调用现有线性模型写成首次贡献。[5,6]
2. **估计对象仍需区分。** LIANA+ `df_to_lr` 将基因统计拼接到 LR 表，官方教程明确 interaction_* 为两端统计量均值。已查看本机 `liana/multi/df_to_lr.py:182–184`：逐列 mean。平均 p 值不自动成为 `C=L^TΩR` 的联合 p 值。直接 C 的 beta/SE 与基因 DEA 不同，但不能只凭这一差异认定有重要新方法。[5]
3. **不混淆 null。** CellPhoneDB specificity p 值回答 LR 是否细胞类型特异，不能拿 condition-label null 下的拒绝率称为它的 type-I error。基线须标明 specificity、sample-score 二阶段差异，或自实现的 pooled-condition 检验。[1]
4. **统一和速度均有先例。** LIANA+/CellChat 已覆盖两模态，FastCCC 已消除 permutation，SpatialDM 有解析空间检验。`L^TΩR` 是统一记号和实现选择，单凭公式不构成创新。[2–5,10]
5. **IC 算法不是新算法。** 二部图 MVC、hitting set、干预靶点组合已有成熟先例。剩余待检验的是“LR 通道覆盖数”是否在真实 programs 上独立、有数据库稳定性和可干预解释。[16]
6. **MOSANIC 碰撞需保留边界。** 作者接口明确支持 hub、relay、单基因 in-silico knockout；尚未核实其全文是否给出同义的最小集合及真实 blockade 验证。其 DOI 页面本次浏览失败。不可宣称已经彻底排除干预方向碰撞。[14]

初步决策：A 的广义“样本层/一般设计首创”主张 KILL；直接 C 估计对象的实用性待测。IC 与完整组合暂 HOLD，继续任务规定的实证 kill-tests。此处未用未运行结果支持 GO。

## 原始来源

[1] CellPhoneDB v5 [官方方法说明](https://cellphonedb.readthedocs.io/en/stable/RESULTS-DOCUMENTATION.html)；[API](https://cellphonedb.readthedocs.io/en/stable/cellphonedb.src.core.methods.html)。DOI: [10.1038/s41596-024-01137-1](https://doi.org/10.1038/s41596-024-01137-1)。

[2] Jin et al. CellChat protocol，2024 online。[10.1038/s41596-024-01045-4](https://doi.org/10.1038/s41596-024-01045-4)；[官方 spatial tutorial](https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.Rmd)，已保存 `sources/CellChat_spatial_tutorial.Rmd`。

[3] [CellChat 官方 README 的 v3 声明](https://github.com/jinworks/CellChat/blob/main/README.md)；不将 v2 实测自动当成 v3 实测。

[4] Hou et al. FastCCC，2025。[10.1038/s41467-025-66272-z](https://doi.org/10.1038/s41467-025-66272-z)；[官方仓库](https://github.com/Svvord/FastCCC)。

[5] Dimitrov et al. LIANA+，2024。[10.1038/s41556-024-01469-w](https://doi.org/10.1038/s41556-024-01469-w)；[官方 targeted workflow](https://liana.readthedocs.io/en/stable/tutorials/notebooks/targeted.html)；[df_to_lr API](https://liana.readthedocs.io/en/stable/generated/liana.multi.df_to_lr.html)。

[6] [MultiNicheNet 官方仓库及 2024-05 v2 更新](https://github.com/saeyslab/multinichenetr)；[paired 官方教程](https://github.com/saeyslab/multinichenetr/blob/main/vignettes/paired_analysis_SCC.knit.md)；[原始预印本 10.1101/2023.06.13.544751](https://doi.org/10.1101/2023.06.13.544751)。

[7] Tensor-cell2cell v2，2026：[10.1093/bioinformatics/btaf667](https://doi.org/10.1093/bioinformatics/btaf667)，PMID 41719185；[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12937581/)。v1 为其引述的既有 tensor framework。

[8] Mitchell et al. dominoSignal，2026：[10.1093/bioinformatics/btag089](https://doi.org/10.1093/bioinformatics/btag089)；[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12998610/)。

[9] Lagger et al. scDiffCom，2023（窗口前但任务指定）：[10.1038/s43587-023-00514-x](https://doi.org/10.1038/s43587-023-00514-x)。

[10] SpatialDM，2023（窗口前）：[10.1038/s41467-023-39608-w](https://doi.org/10.1038/s41467-023-39608-w)。

[11] COMMOT，2023（窗口前）：[10.1038/s41592-022-01728-4](https://doi.org/10.1038/s41592-022-01728-4)。

[12] scRICH，2026：[10.1016/j.cels.2026.101703](https://doi.org/10.1016/j.cels.2026.101703)；[作者论文](https://www.sciencedirect.com/science/article/pii/S2405471226001857)。

[13] Renoir，2026：[10.1038/s41467-026-72388-7](https://doi.org/10.1038/s41467-026-72388-7)。

[14] MOSANIC，2026 预印本：[作者仓库](https://github.com/debraj-55555/MOSANIC)；[10.64898/2026.07.05.736582](https://doi.org/10.64898/2026.07.05.736582)。全文本次尚未取得。

[15] scSeqCommDiff，2025：[10.1093/nargab/lqaf084](https://doi.org/10.1093/nargab/lqaf084)。

[16] The landscape of receptor-mediated precision cancer combination therapy via a single-cell perspective，2022：[10.1038/s41467-022-29154-2](https://doi.org/10.1038/s41467-022-29154-2)。

检索词包含每个指定方法名，以及 differential/sample-aware CCC、minimum vertex cover、minimum hitting set、minimum intervention、network vulnerability、perturbation。精确 IC 词检索未找到同定义不是不存在证明；特别保留 [14] 全文缺口。广泛检索命中的综述仅用于发现来源，不作为方法能力依据。

对PLAN的影响：H1未获可宣称首次的完整排除证据；H2广义一般设计首创已被现有工作反驳，后续只评价直接C估计对象与工程价值。
