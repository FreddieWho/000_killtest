# 统计 kill-test

## 估计对象与输入

Kang官方mirror：24,673 cells ×15,706 genes，8 donors ×2 conditions。矩阵为整数counts，逐细胞行和与nCount_RNA差为0。映射：patient=replicate，condition=label，sample=patient__condition，celltype=原cell_type。原始GEO标题/设计见 `sources/GSE96583.txt`，PMID 29227470；本次不用细胞重聚类。

按 sample×celltype 保存 raw sum、library size、log1p(CP10k) mean/variance、detection、cell number：`results/kang_run_03/sufficient_statistics.npz`、`groups.tsv`。library normalization在全部测量基因上完成，随后索引LR；LR裁剪不改变library size。simple C=mean(L)×mean(R)，complex先取每细胞必需亚基min再聚合。

真实设计 `~ condition + patient`：16行、9列、残差df=7。另以48个独立合成样本、200个响应跑通 `~ treatment*time + batch`，输出beta/SE/p/FDR；见 `factorial_fit.tsv`、`factorial_metadata.tsv`、`designs.json`。未运行可选LMM，不声称已验证随机效应或真实纵向队列。

## Null 构造与指标

枚举8供体的256种within-pair sign flips（双侧统计只有128种符号对称图样）。每个真实样本的所有细胞整块保留，因此供体异质性、组成、library size、不平衡和基因协方差保留。另固定随机种子20260928抽100种外层分配，每种做199次真正的细胞condition-label置换（按celltype保留分组大小）。后者是自实现的 pooled 条件差异reference，**不是CPDB specificity方法的复刻**。

这是对现有配对表达的**条件随机化**，不是独立生成的100/256个生物队列，也不是已经去除IFN机制的sharp-null生物数据。完整随机化p的校准具有构造性，不能单独拿来验证方法。FPR=未校正p<0.05比例；FDR=每次BH后FDP的均值，全null下等于至少一次发现的比例。所有报告保留同一比较中的特征分母；不将相关LR或256个排列当独立患者估计置信区间。

### 六种主要细胞类型，全24,673细胞

| method                            |   n_features |     FPR |   FDR |
|:----------------------------------|-------------:|--------:|------:|
| paired_sample_t                   |         4951 | 0.01688 |  0.01 |
| pooled_cell_condition_permutation |         4951 | 0.1591  |  0.28 |
| pooled_cell_delta                 |         4951 | 0.1102  |  1    |

使用 B/CD14/CD4/CD8/FCGR3A/NK 六类；分析全部基因可映射、条件标签不变筛选后非恒定的4951个simple跨类型LR特征。观测到pooled置换失准，样本层保守。细胞置换p的最小值是1/200，BH受到该离散分辨率限制；不把其较低FDR读作高精度保证。

### 原型全可用细胞类型、256次null

| method                                          |   n_features |     FPR |      FDR |
|:------------------------------------------------|-------------:|--------:|---------:|
| CellChat_paired_product                         |         9216 | 0.01593 | 0.007812 |
| CellChat_paired_product_exact_signflip          |         9216 | 0.01488 | 0        |
| LIANA_CellPhoneDB_paired_product                |         9367 | 0.01575 | 0.007812 |
| LIANA_CellPhoneDB_paired_product_exact_signflip |         9367 | 0.01455 | 0        |

大量低表达LR只有少数供体出现配对差值，离散性使总体保守。事后信息分层（描述，不替代主分母）中，8位供体都有非零差值的结果为：

| resource          |   nonzero_paired_donors |   n_features |   sample_t_FPR |   exact_signflip_FPR |
|:------------------|------------------------:|-------------:|---------------:|---------------------:|
| CellChat          |                       8 |         1735 |        0.04534 |              0.04688 |
| LIANA_CellPhoneDB |                       8 |         1823 |        0.0452  |              0.04688 |

### 原生方法每样本分数 + 同一配对OLS

| method      |   n_features |     FPR |      FDR |
|:------------|-------------:|--------:|---------:|
| CellChat    |        10364 | 0.01607 | 0.007812 |
| CellPhoneDB |         4929 | 0.03148 | 0.03906  |
| FastCCC     |         4929 | 0.03148 | 0.03906  |
| LIANA       |         5390 | 0.04558 | 0.05469  |

真正运行了四个原生包的16个样本；其后对缓存的样本分数做256种条件交换。重新交换标签不改变每样本原生输入，所以无需重复计算同一原生specificity null。原生分数、方法参数、每样本日志与receipt见 `results/native_*_01/`；配对beta/SE/p/FDR和每次null见native analysis目录。

CPDB/FastCCC为各自原生arithmetic score，LIANA为其原生lr_means，CellChat为native probability；非恒定且16样本均有输出的特征各自成family。资源、过滤、分数定义不同，不能仅按这些FPR对方法作统一优劣排名。它们表明已有样本工作流可直接接受一般统计模型。LIANA `interaction_pvalue` 均值未冒充联合p；CPDB/CellChat原生specificity p未用于条件差异FPR。

## Cell-number scaling

|   n_cells | method                            |   n_features |      FPR |   FDR |
|----------:|:----------------------------------|-------------:|---------:|------:|
|      2000 | paired_sample_t                   |          450 | 0.006289 |  0.01 |
|      2000 | pooled_cell_condition_permutation |          450 | 0.05351  |  0    |
|      2000 | pooled_cell_delta                 |          450 | 0.02849  |  0.35 |
|      5000 | paired_sample_t                   |          450 | 0.01062  |  0.01 |
|      5000 | pooled_cell_condition_permutation |          450 | 0.09089  |  0.08 |
|      5000 | pooled_cell_delta                 |          450 | 0.05396  |  0.68 |
|     10000 | paired_sample_t                   |          450 | 0.01562  |  0.01 |
|     10000 | pooled_cell_condition_permutation |          450 | 0.1302   |  0.22 |
|     10000 | pooled_cell_delta                 |          450 | 0.09296  |  0.92 |
|     24673 | paired_sample_t                   |          450 | 0.02196  |  0.01 |
|     24673 | pooled_cell_condition_permutation |          450 | 0.1841   |  0.39 |
|     24673 | pooled_cell_delta                 |          450 | 0.148    |  0.95 |

**缩减范围明示：** 从2k到full都至少每sample×celltype有5细胞的共同集合，仅剩CD14+ monocytes和CD4 T、450个特征。此表仅用于同特征规模曲线；六类型主结果单独完整执行如上。所有规模使用真实细胞的分层子样本，保存具体cell indices；未用复制细胞生成biology结果。

![QQ](results/figures/null_QQ.png)
![FPR scaling](results/figures/FPR_scaling.png)

QQ图针对六类型主实验；显示y轴上限8，delta-method更极端的尾部超出画幅，完整分位数保存在对应QQ TSV。离散p、稀疏零差值下不期待严格连续uniform；不能把QQ偏离全部解释为失准。

## Power：仅score层已知注入

在label-invariant非零供体差值方差的特征中选择约10%，分别注入0.5/1/2 donor-SD的正向score差异，再做BH。256种符号排列重复。原型的921个注入特征已核实全部非零。该power是统计量层的条件合成结果，不能等同真实LR扰动效应或方法重建biology的能力。

|   effect_in_donor_SD |     power |      FDR |   n_injected |
|---------------------:|----------:|---------:|-------------:|
|                  0.5 | 0.0004284 | 0.006777 |          921 |
|                  1   | 0.001039  | 0.005941 |          921 |
|                  2   | 0.8417    | 0.008823 |          921 |

原生分数的同类合成比较（各自feature universe，不能作严格共同真值排行榜）：

| method      |   effect_in_donor_SD |     power |      FDR |   n_injected |
|:------------|---------------------:|----------:|---------:|-------------:|
| CellChat    |                  0.5 | 0.0003582 | 0.006785 |         1036 |
| CellChat    |                  1   | 0.0007428 | 0.006099 |         1036 |
| CellChat    |                  2   | 0.8707    | 0.008555 |         1036 |
| CellPhoneDB |                  0.5 | 0.001588  | 0.0345   |          492 |
| CellPhoneDB |                  1   | 0.003867  | 0.03112  |          492 |
| CellPhoneDB |                  2   | 0.7169    | 0.01944  |          492 |
| FastCCC     |                  0.5 | 0.001588  | 0.0345   |          492 |
| FastCCC     |                  1   | 0.003867  | 0.03112  |          492 |
| FastCCC     |                  2   | 0.7169    | 0.01944  |          492 |
| LIANA       |                  0.5 | 0.02089   | 0.04745  |          539 |
| LIANA       |                  1   | 0.1402    | 0.04784  |          539 |
| LIANA       |                  2   | 0.9395    | 0.04702  |          539 |

## 刺激方向与决定

ISG15/IFIT1/IFIT3/MX1/OAS1的供体层方向检查见 `results/ic_sanity_01/IFN_response_sanity.tsv`：35/35 个celltype×marker均值正向。它只检查IFN刺激响应，不能验证endogenous sender→receiver blockade。

KT1：在该真实结构条件随机化下发现pooled失准及随细胞数增加恶化；样本层总体保守、有信息子集接近nominal。KT2：paired/factorial跑通，但LIANA/MultiNicheNet已有一般设计路线，**统计创新贡献KILL**。对PLAN的影响：支持关注实验单位，不支持把现成样本层建模包装为新方法。
