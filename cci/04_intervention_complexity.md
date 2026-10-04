# Intervention Complexity

## 定义与实际支持图

simple LR：在sample×sender×receiver内，双方detection≥0.1且score>0为active；每组至少5细胞。去重相同分子LR边，不把重复数据库条目当新通道。主分析为CellChatDB的Secreted Signaling+Cell-Cell Contact，与LIANA打包的CellPhoneDB资源分别映射到全部实测基因。

IC100为role-tagged二部图的minimum vertex cover；由maximum matching精确求解。代码逐图检查每条边被覆盖且cover大小等于matching大小。Count=active simple边数，strength=这些边的C总和。n_programs是相关图数量，**生物重复仍为8 donors**。

若同一gene既是L又是R，角色节点的两个target不必对应两个物理分子。另对这类图求去重gene最小hitting set，不能把普通MVC的多项式保证照搬到合并基因图。

complex采用每细胞必需亚基min后聚合，支持率为所有必需亚基在同一细胞出现的比例。每个active interaction的可击中集合为配体/受体全部必需亚基；任一个必需亚基被完全阻断即假设该通道失活。对含complex且通道数最多的100个sample programs运行exact ILP：**100/100 OPTIMAL，gap=0**，保存靶点证书。不把可选agonist/cofactor自动当作必需靶点。

IC90是原指令的可选项，本轮未实现；未用strength-weighted替代主定义。

## 经验独立性

仅nonempty programs；一元线性/单调模型都按**留一供体**预测，避免同供体program混入训练和评估。未设定事后排除失败数据库的规则。

| database          | predictor   |   n_programs |   pearson |   spearman |   heldout_donor_isotonic_R2 |   heldout_donor_loglinear_R2 |   similar_value_different_IC_program_fraction |
|:------------------|:------------|-------------:|----------:|-----------:|----------------------------:|-----------------------------:|----------------------------------------------:|
| CellChat          | count       |          797 |    0.8522 |     0.8454 |                      0.7356 |                       0.67   |                                        0.7716 |
| CellChat          | strength    |          797 |    0.3414 |     0.5982 |                      0.4173 |                       0.3009 |                                        0.8482 |
| LIANA_CellPhoneDB | count       |          750 |    0.9562 |     0.9633 |                      0.9185 |                       0.86   |                                        0.5653 |
| LIANA_CellPhoneDB | strength    |          750 |    0.7811 |     0.8212 |                      0.682  |                       0.6711 |                                        0.82   |

**按原阈值判定：LIANA/CellPhoneDB的count预测R²>0.9且近单调，失败。** CellChat资源中不等价，不能由此把跨库失败改写为全面PASS；也不能由一个资源的失败断言数学上所有图都等价于count。

去重分子gene的敏感性：

| database          | predictor   |   n_programs |   molecular_IC_heldout_donor_R2 |   spearman |   programs_with_role_overlap |   programs_with_changed_IC |
|:------------------|:------------|-------------:|--------------------------------:|-----------:|-----------------------------:|---------------------------:|
| CellChat          | count       |          797 |                          0.7356 |     0.8454 |                           23 |                          0 |
| CellChat          | strength    |          797 |                          0.4173 |     0.5982 |                           23 |                          0 |
| LIANA_CellPhoneDB | count       |          750 |                          0.9183 |     0.9639 |                           29 |                         22 |
| LIANA_CellPhoneDB | strength    |          750 |                          0.68   |     0.822  |                           29 |                         22 |

CellChat中23个nonempty programs有跨角色同gene，但未改变最小值；替代资源29个有重合、22个改变最小值。count预测约0.918的失败仍在。

相似count/strength但IC不同的具体图见 `results/kang_run_03/similar_value_cases.tsv`（每类最多展示30对，参与比例按全部匹配计算）；四象限各示例见 `results/kang_run_03/quadrants.tsv`。相似定义为差≤5%、IC差至少2且至少1.5倍，仅用于描述案例，不冒充显著性检验。

![IC quadrants](results/figures/IC_quadrants.png)

## 数据库敏感性

在747个两库均非空的同sample×sender×receiver上，IC Spearman=0.8383，count Spearman=0.5749，按各库median分高低的一致率=83.13%。排名未完全随机，但约16.9%高低类别翻转。资源的类别覆盖与节点/边数量不同，这属于此原生资源比较的实际敏感性，不是只改变一个因素的因果实验。

共享LR交集会产生相同支持图和IC，不能把这种构造性相同包装成独立稳定性证据。两库原始/映射表、程序配对表均已保存；未以少量挑选LR替代全映射扫描。

## 最低成本干预有效性

`results/ic_sanity_01/synthetic_intervention.json`给出两个count=12、total strength=12的图：star的IC=1，12条独立LR轴的IC=12；最优单靶分别删除100%和1/12通道。这是图删除恒等式，证明定义有结构含义；**不证明downstream biology**。低IC预测的是理想、正确选中的最小靶集，并不意味着任意单靶都有效。

检索了公开刺激/扰动资料：

| 数据 | 来源 | 本轮用途/边界 |
|---|---|---|
| Kang IFN-β | GSE96583；PMID29227470 | 已运行刺激响应方向检查，非内源LR blockade真值 |
| Immune Dictionary | [10.1038/s41586-023-06816-9](https://doi.org/10.1038/s41586-023-06816-9)，GSE202186 | 外源cytokine刺激；不能直接回答固定内源网络最小拆解靶数，未下载 |
| Perturb-CITE-seq | [10.1038/s41588-021-00779-1](https://doi.org/10.1038/s41588-021-00779-1)，SCP1064 | co-culture/CRISPR更相关，但尚未建立control网络—特定target—response的IC配准；未下载218k全集，未作验证结论 |

**真实独立干预验证：NOT_VALIDATED。** 按原指令允许的fallback，本轮仅synthetic proof与刺激sanity；不把“未找到即取的数据”说成数据不存在。固定LR图不含剂量、蛋白定位、补偿和反馈，IC不是可直接开药的靶点数量。

KT4：跨合理资源的通用独立性主张KILL。KT5：有相当排名一致性，但非完全稳定；生物学保持HOLD。对PLAN的影响：H3未满足整体GO所需条件，不能通过挑数据库或换成strength score保留原主张。
