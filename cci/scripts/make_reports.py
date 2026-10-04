"""Consolidate completed receipts. Refuse to describe unfinished runs as done."""
from pathlib import Path
import json,re
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
C=json.load(open('results/selected_runs.json'));K=Path(C['kang']);S=Path(C['spatial']);G=Path(C['scaling']);F=Path(C['pooled_full']);I=Path(C['ic_sanity'])
for p in [K,S,G,F,I,*map(Path,C['native'].values()),*map(Path,C['native_analyses']),Path(C['spatial_native'])]:
 r=json.load(open(p/'receipt.json'))
 if r['status']!='DONE':raise RuntimeError('Incomplete '+str(p))
checks=json.load(open('results/algebra_checks.json'))
if not all(x['pass_'] for x in checks):raise RuntimeError('Algebra checks not passed')
def table(d):return d.to_markdown(index=False,floatfmt='.4g')
def load(p):return pd.read_csv(p,sep='\t')
def write(name,s):Path(name).write_text(s.strip()+'\n')
def null_summary(d):return d.groupby('method',as_index=False).agg(n_features=('n_features','first'),FPR=('FPR','mean'),FDR=('FDP','mean'))
def peak(p):
 t=Path(p).read_text();m=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',t);return int(m.group(1))/1024 if m else np.nan
kn=null_summary(load(K/'null_replicates.tsv'));full=null_summary(load(F/'null_replicates.tsv'))
ns=[load(Path(p)/'null_replicates.tsv') for p in C['native_analyses']];native=null_summary(pd.concat(ns))
sc=load(G/'null_replicates.tsv').groupby(['n_cells','method'],as_index=False).agg(n_features=('n_features','first'),FPR=('FPR','mean'),FDR=('FDP','mean'))
ic=load(K/'ic_independence.tsv');molecular=load(I/'molecular_ic_independence.tsv');db=json.load(open(K/'database_sensitivity.json'));info=json.load(open(S/'input_audit.json'));sp=json.load(open(S/'typed_summary.json'));equiv=json.load(open(S/'numerical_equivalence.json'))
sp_profile=load(S/'profile.tsv');proto=load('results/prototype_current_profile.tsv');icphase=load('results/ic_phase_profile.tsv')
rt=[]
for method,p in C['native'].items():
 d=load(Path(p)/'profile.tsv');r=json.load(open(Path(p)/'receipt.json'))
 rt.append(dict(method=method,samples=len(d),native_API_seconds=d.seconds.sum(),peak_RSS_MiB=peak(f'results/native_{method}_01.time.txt'),permutations_per_sample=r.get('iterations',r.get('nboot',0))))
runtime=pd.DataFrame(rt)
core_seconds=proto[proto.stage.isin(['load','aggregation','communication_CellChat','fit_CellChat'])].seconds.sum()
core_peak=peak('results/prototype_current_v2.time.txt')
spnative=json.load(open(Path(C['spatial_native'])/'receipt.json'))
power=load(K/'synthetic_score_power.tsv').groupby('effect_in_donor_SD',as_index=False).agg(power=('power','mean'),FDR=('FDP','mean'),n_injected=('n_injected','first'))
npower=pd.concat([load(Path(p)/'synthetic_score_power.tsv') for p in C['native_analyses']]).groupby(['method','effect_in_donor_SD'],as_index=False).agg(power=('power','mean'),FDR=('FDP','mean'),n_injected=('n_injected','first'))
strata=load(K/'null_information_strata.tsv');response=load(I/'IFN_response_sanity.tsv');covers=load(K/'complex_hitting_sets.tsv')
# Machine-readable benchmark: one row per metric, explicit scope and denominator.
bench=[]
def add(module,method,metric,value,unit,scope,denominator,status='OBSERVED'):
 bench.append(dict(module=module,method=method,metric=metric,value=value,unit=unit,scope=scope,denominator=denominator,status=status))
for scope,d in [('full_6_types_100_null',full),('prototype_256_null',kn),('native_scores_256_null',native)]:
 for r in d.itertuples():
  for metric in ['FPR','FDR']:add('statistical',r.method,metric,getattr(r,metric),'fraction',scope,f'{r.n_features} features;8 donors')
for r in sc.itertuples():
 for metric in ['FPR','FDR']:add('scaling',r.method,metric,getattr(r,metric),'fraction',f'{r.n_cells} actual cells;100 outer null;199 inner permutations',f'{r.n_features} common features;2 cell types;8 donors')
for r in ic.itertuples():
 add('IC',r.database,'heldout_donor_R2_'+r.predictor,r.heldout_donor_isotonic_R2,'R2','role-tagged MVC',f'{r.n_programs} nonempty programs;8 donors')
add('IC','two_resources','IC_rank_spearman',db['IC_rank_spearman'],'rho','native mapped catalogues',f'{db["n_matched_programs"]} matched programs;8 donors')
for r in runtime.itertuples():
 add('engineering',r.method,'native_API_total',r.native_API_seconds,'seconds','native specificity plus native scores','16 samples')
 add('engineering',r.method,'process_peak_RSS',r.peak_RSS_MiB,'MiB','native process including imports and outputs','16 samples')
add('engineering','prototype','load_aggregate_communication_fit',core_seconds,'seconds','product statistic;CellChat resource;without native specificity','16 samples')
add('engineering','prototype','process_peak_RSS',core_peak,'MiB','current algebra/profile process, both resources and checks','24673 cells')
for r in sp_profile.itertuples():
 n=getattr(r,'n_spots',info['spots']);n=info['spots'] if pd.isna(n) else int(n)
 add('spatial',str(getattr(r,'mode','')),r.stage,r.seconds,'seconds','one section;graph copies are synthetic scaling',str(n))
for r in icphase.itertuples():
 add('engineering',r.database+'_'+r.algorithm,r.stage,r.seconds,'seconds','construction and solve separated',str(r.n_programs)+' programs')
for d in C['native_analyses']:
 for r in load(Path(d)/'profile.tsv').itertuples():
  add('engineering',r.method,'sample_OLS_fit',r.seconds,'seconds','after native sample scores',str(r.n_features)+' features;16 samples')
add('spatial','typed_kernel','removed_supported_channels',sp['removed_by_geometry'],'channels','same matrix and LR resource;change Omega only',str(sp['mixing_supported_channels']))
for r in power.itertuples():add('power','prototype','BH_power',r.power,'fraction',f'score-level effect={r.effect_in_donor_SD} donor SD;256 sign swaps',str(r.n_injected))
for r in npower.itertuples():add('power',r.method,'BH_power',r.power,'fraction',f'score-level effect={r.effect_in_donor_SD} donor SD;256 sign swaps',str(r.n_injected))
add('biology','IC','real_blockade_validation',np.nan,'NA','No matched independent IC-vs-response experiment','0','NOT_VALIDATED')
pd.DataFrame(bench).to_csv('results/benchmark.tsv',sep='\t',index=False)
# Requested diagnostic plots, exported as standalone artifacts.
figdir=Path('results/figures');figdir.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(5,4))
for name in ['paired_sample_t','pooled_cell_condition_permutation','pooled_cell_delta']:
 d=load(F/f'QQ_{name}_{24673}.tsv');ax.plot(-np.log10(d.expected),-np.log10(np.maximum(d.observed,1e-300)),label=name)
ax.plot([0,3],[0,3],'k--',lw=.8);ax.set(xlabel='Expected -log10(p)',ylabel='Observed -log10(p)',title='Conditional null: six cell types');ax.set_ylim(0,8);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(figdir/'null_QQ.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(5,4))
for method,d in sc.groupby('method'):ax.plot(d.n_cells,d.FPR,'o-',label=method)
ax.axhline(.05,color='k',ls='--',lw=.8);ax.set(xlabel='Cells (real subsamples)',ylabel='Mean FPR',title='450 common LR features; two cell types');ax.legend(fontsize=7);fig.tight_layout();fig.savefig(figdir/'FPR_scaling.png',dpi=180);plt.close(fig)
programs=load(K/'programs.tsv');fig,axes=plt.subplots(1,2,figsize=(9,3.8))
for ax,(name,d) in zip(axes,programs.groupby('database')):
 d=d[d['count']>0];ax.scatter(d.strength,d.ic100,s=8,alpha=.35);ax.axvline(d.strength.median(),color='grey',lw=.7);ax.axhline(d.ic100.median(),color='grey',lw=.7);ax.set(xscale='log',xlabel='Total strength',ylabel='IC100',title=name)
fig.tight_layout();fig.savefig(figdir/'IC_quadrants.png',dpi=180);plt.close(fig)
# Reports. All dataset-specific claims below come from the selected completed run files.
write('00_DECISION.md',f'''
# 决定：KILL 统一新方法主线

2026-09-28。完成本轮kill-test；不进入完整package/论文开发，保留小型工程原型。

- **统计问题存在，但一般设计不是新贡献。** 六类细胞、4951特征的pooled条件置换FPR为15.91%，配对样本检验1.69%；LIANA原生样本分数接同一模型为4.56%。前者是自实现条件差异reference，不是CPDB原生specificity p值失效。
- **IC跨库独立性失败。** 留一供体的count→IC单调预测R²：CellChatDB 0.736，LIANA/CellPhoneDB 0.918，后者ρ≈0.963；去重物理gene后仍约0.918，触发原KILL条件。跨库排名ρ=0.838，不等于完全随机。
- **工程有用，不能补足科学条件。** 原型读入→聚合→C→拟合{core_seconds:.3f}s；空间全LR dense对照通过，typed kernel移除{sp['removed_by_geometry']:,}个支持通道。原型核查进程峰值约994MiB，本轮未满足“时间更少且整体RAM更低”的组合判据。

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
''')
write('02_statistical_killtest.md',f'''
# 统计 kill-test

## 估计对象与输入

Kang官方mirror：24,673 cells ×15,706 genes，8 donors ×2 conditions。矩阵为整数counts，逐细胞行和与nCount_RNA差为0。映射：patient=replicate，condition=label，sample=patient__condition，celltype=原cell_type。原始GEO标题/设计见 `sources/GSE96583.txt`，PMID 29227470；本次不用细胞重聚类。

按 sample×celltype 保存 raw sum、library size、log1p(CP10k) mean/variance、detection、cell number：`{K}/sufficient_statistics.npz`、`groups.tsv`。library normalization在全部测量基因上完成，随后索引LR；LR裁剪不改变library size。simple C=mean(L)×mean(R)，complex先取每细胞必需亚基min再聚合。

真实设计 `~ condition + patient`：16行、9列、残差df=7。另以48个独立合成样本、200个响应跑通 `~ treatment*time + batch`，输出beta/SE/p/FDR；见 `factorial_fit.tsv`、`factorial_metadata.tsv`、`designs.json`。未运行可选LMM，不声称已验证随机效应或真实纵向队列。

## Null 构造与指标

枚举8供体的256种within-pair sign flips（双侧统计只有128种符号对称图样）。每个真实样本的所有细胞整块保留，因此供体异质性、组成、library size、不平衡和基因协方差保留。另固定随机种子20260928抽100种外层分配，每种做199次真正的细胞condition-label置换（按celltype保留分组大小）。后者是自实现的 pooled 条件差异reference，**不是CPDB specificity方法的复刻**。

这是对现有配对表达的**条件随机化**，不是独立生成的100/256个生物队列，也不是已经去除IFN机制的sharp-null生物数据。完整随机化p的校准具有构造性，不能单独拿来验证方法。FPR=未校正p<0.05比例；FDR=每次BH后FDP的均值，全null下等于至少一次发现的比例。所有报告保留同一比较中的特征分母；不将相关LR或256个排列当独立患者估计置信区间。

### 六种主要细胞类型，全24,673细胞

{table(full)}

使用 B/CD14/CD4/CD8/FCGR3A/NK 六类；分析全部基因可映射、条件标签不变筛选后非恒定的4951个simple跨类型LR特征。观测到pooled置换失准，样本层保守。细胞置换p的最小值是1/200，BH受到该离散分辨率限制；不把其较低FDR读作高精度保证。

### 原型全可用细胞类型、256次null

{table(kn)}

大量低表达LR只有少数供体出现配对差值，离散性使总体保守。事后信息分层（描述，不替代主分母）中，8位供体都有非零差值的结果为：

{table(strata[strata.nonzero_paired_donors==8])}

### 原生方法每样本分数 + 同一配对OLS

{table(native)}

真正运行了四个原生包的16个样本；其后对缓存的样本分数做256种条件交换。重新交换标签不改变每样本原生输入，所以无需重复计算同一原生specificity null。原生分数、方法参数、每样本日志与receipt见 `results/native_*_01/`；配对beta/SE/p/FDR和每次null见native analysis目录。

CPDB/FastCCC为各自原生arithmetic score，LIANA为其原生lr_means，CellChat为native probability；非恒定且16样本均有输出的特征各自成family。资源、过滤、分数定义不同，不能仅按这些FPR对方法作统一优劣排名。它们表明已有样本工作流可直接接受一般统计模型。LIANA `interaction_pvalue` 均值未冒充联合p；CPDB/CellChat原生specificity p未用于条件差异FPR。

## Cell-number scaling

{table(sc)}

**缩减范围明示：** 从2k到full都至少每sample×celltype有5细胞的共同集合，仅剩CD14+ monocytes和CD4 T、450个特征。此表仅用于同特征规模曲线；六类型主结果单独完整执行如上。所有规模使用真实细胞的分层子样本，保存具体cell indices；未用复制细胞生成biology结果。

![QQ](results/figures/null_QQ.png)
![FPR scaling](results/figures/FPR_scaling.png)

QQ图针对六类型主实验；显示y轴上限8，delta-method更极端的尾部超出画幅，完整分位数保存在对应QQ TSV。离散p、稀疏零差值下不期待严格连续uniform；不能把QQ偏离全部解释为失准。

## Power：仅score层已知注入

在label-invariant非零供体差值方差的特征中选择约10%，分别注入0.5/1/2 donor-SD的正向score差异，再做BH。256种符号排列重复。原型的921个注入特征已核实全部非零。该power是统计量层的条件合成结果，不能等同真实LR扰动效应或方法重建biology的能力。

{table(power)}

原生分数的同类合成比较（各自feature universe，不能作严格共同真值排行榜）：

{table(npower)}

## 刺激方向与决定

ISG15/IFIT1/IFIT3/MX1/OAS1的供体层方向检查见 `{I}/IFN_response_sanity.tsv`：{int((response.mean_stim_minus_ctrl>0).sum())}/{len(response)} 个celltype×marker均值正向。它只检查IFN刺激响应，不能验证endogenous sender→receiver blockade。

KT1：在该真实结构条件随机化下发现pooled失准及随细胞数增加恶化；样本层总体保守、有信息子集接近nominal。KT2：paired/factorial跑通，但LIANA/MultiNicheNet已有一般设计路线，**统计创新贡献KILL**。对PLAN的影响：支持关注实验单位，不支持把现成样本层建模包装为新方法。
''')
write('03_engineering_profile.md',f'''
# 工程 profile

## 实测范围

CPU本地运行，Python/R环境见 `results/environment.json`、`results/R_environment.txt`。设置主要数值线程为1；宿主机有其他分析进程，数据未作严格独占机器测速。不报告亚秒差异的统计显著性。原生方法在suite中逐个运行，每个方法处理全部16样本并落盘；部分时段与独立空间/分析进程重叠。

{table(runtime)}

CPDB/LIANA各1000置换；FastCCC Mean/Minimum/Arithmetic；CellChat100 bootstrap、trim=0，全部已测LR特征，包含cofactor。CPDB官方v5与LIANA打包的cellphonedb资源不是同一个版本快照。所有原生API实跑，不以自实现近似代替。

## 原型阶段

{table(proto)}

其中CellChat资源一条读入→聚合→C→拟合链为 **{core_seconds:.4g}s**；上表峰值列为整个核查进程最终high-water mark，不是逐阶段额外RAM。原始主运行的IC、null和power阶段另见 `{K}/profile.tsv`。构图和解算单独计时，并逐图确认最小值与冻结结果一致：

{table(icphase)}

**速度比较边界：** 原生API包含specificity评估与输出，本原型只直接计算C后做样本层模型；它们不是相同统计输出。该现实multi-sample工作流有明显时间节省空间，但不能写成“相同准确性推断全面快若干倍”。原生工具也可跳过不需要的specificity流程。**本轮未满足“时间更少且整体RAM更低”的组合判据。** 原型核查进程峰值约994MiB，原生API进程约488–895MiB；两侧输入处理范围不同，不能把这组数当公平端到端比较，更不能声称已获得RAM收益。native输入已完成公用的library归一化和LR裁剪，prep内存不能被忽略。

2k/5k/10k/full的aggregation、固定450特征的100次样本拟合和100×199 cell permutation分别计时：`{G}/profile.tsv`。post-aggregation fit不访问细胞矩阵；cell-number scaling不是把新的供体数伪造成更多细胞。50k/100k duplication为可选，未运行。

## 空间同口径计算

{table(sp_profile[sp_profile.stage.isin(['sparse_graph','dense_distance_comparator','dense_communication_comparator','sparse_communication'])].fillna('—'))}

sparse/dense aggregation使用相同LR向量、相同kernel和相同分组操作，仅矩阵表示不同；覆盖全部LR和全部group pairs。dense距离/权重只存在于此比较脚本，正式`cci_core.opportunity_graph`使用KD-tree+CSR。

Spatial CellChat原生20-bootstrap profiling：{spnative['seconds']:.4g}s，{spnative['spots']} spots；它的群体距离/概率模型与本原型不同，因此只作现实参照，不作数值等价基线。

## 复杂度与资源限制

- scRNA：读取/归一化/聚合需处理输入非零元素；聚合后C成本主要随samples×celltype pairs×LR变化，sample层拟合与cell数无关。
- spatial：KD-tree建图通常O(N logN+|E|)，每个LR的机会聚合O(|E|)，全资源为O(K|E|)，再加输入/输出成本；并不是无视K的总O(|E|)。半径过大时|E|仍可逼近N²，不能作无条件线性保证。
- simple IC：二部图maximum matching/MVC；复杂体/同基因跨角色的全局靶点为ILP，最坏情形NP-hard。100个实测最优不代表任意规模可快解。
- CSR列反复切片的初版很慢；改为CSC按列读取、复用complex组件和按receiver复用空间传播。最终全输出与冻结数值/完整dense对照一致，见 `results/algebra_checks.json` 和 `{S}/numerical_equivalence.json`。旧run保留且不混入最终profile。

KT3：共享聚合与稀疏图的工程收益成立；未给出所有统计输出完全相同且整体RAM更低的全面优越性证据。对PLAN的影响：H4支持保留小型实现，不能挽救H2/H3创新失败。
''')
write('04_intervention_complexity.md',f'''
# Intervention Complexity

## 定义与实际支持图

simple LR：在sample×sender×receiver内，双方detection≥0.1且score>0为active；每组至少5细胞。去重相同分子LR边，不把重复数据库条目当新通道。主分析为CellChatDB的Secreted Signaling+Cell-Cell Contact，与LIANA打包的CellPhoneDB资源分别映射到全部实测基因。

IC100为role-tagged二部图的minimum vertex cover；由maximum matching精确求解。代码逐图检查每条边被覆盖且cover大小等于matching大小。Count=active simple边数，strength=这些边的C总和。n_programs是相关图数量，**生物重复仍为8 donors**。

若同一gene既是L又是R，角色节点的两个target不必对应两个物理分子。另对这类图求去重gene最小hitting set，不能把普通MVC的多项式保证照搬到合并基因图。

complex采用每细胞必需亚基min后聚合，支持率为所有必需亚基在同一细胞出现的比例。每个active interaction的可击中集合为配体/受体全部必需亚基；任一个必需亚基被完全阻断即假设该通道失活。对含complex且通道数最多的100个sample programs运行exact ILP：**{int((covers.status=='OPTIMAL').sum())}/{len(covers)} OPTIMAL，gap=0**，保存靶点证书。不把可选agonist/cofactor自动当作必需靶点。

IC90是原指令的可选项，本轮未实现；未用strength-weighted替代主定义。

## 经验独立性

仅nonempty programs；一元线性/单调模型都按**留一供体**预测，避免同供体program混入训练和评估。未设定事后排除失败数据库的规则。

{table(ic)}

**按原阈值判定：LIANA/CellPhoneDB的count预测R²>0.9且近单调，失败。** CellChat资源中不等价，不能由此把跨库失败改写为全面PASS；也不能由一个资源的失败断言数学上所有图都等价于count。

去重分子gene的敏感性：

{table(molecular)}

CellChat中23个nonempty programs有跨角色同gene，但未改变最小值；替代资源29个有重合、22个改变最小值。count预测约0.918的失败仍在。

相似count/strength但IC不同的具体图见 `{K}/similar_value_cases.tsv`（每类最多展示30对，参与比例按全部匹配计算）；四象限各示例见 `{K}/quadrants.tsv`。相似定义为差≤5%、IC差至少2且至少1.5倍，仅用于描述案例，不冒充显著性检验。

![IC quadrants](results/figures/IC_quadrants.png)

## 数据库敏感性

在{db['n_matched_programs']}个两库均非空的同sample×sender×receiver上，IC Spearman={db['IC_rank_spearman']:.4f}，count Spearman={db['count_rank_spearman']:.4f}，按各库median分高低的一致率={db['high_low_agreement']:.2%}。排名未完全随机，但约{1-db['high_low_agreement']:.1%}高低类别翻转。资源的类别覆盖与节点/边数量不同，这属于此原生资源比较的实际敏感性，不是只改变一个因素的因果实验。

共享LR交集会产生相同支持图和IC，不能把这种构造性相同包装成独立稳定性证据。两库原始/映射表、程序配对表均已保存；未以少量挑选LR替代全映射扫描。

## 最低成本干预有效性

`{I}/synthetic_intervention.json`给出两个count=12、total strength=12的图：star的IC=1，12条独立LR轴的IC=12；最优单靶分别删除100%和1/12通道。这是图删除恒等式，证明定义有结构含义；**不证明downstream biology**。低IC预测的是理想、正确选中的最小靶集，并不意味着任意单靶都有效。

检索了公开刺激/扰动资料：

| 数据 | 来源 | 本轮用途/边界 |
|---|---|---|
| Kang IFN-β | GSE96583；PMID29227470 | 已运行刺激响应方向检查，非内源LR blockade真值 |
| Immune Dictionary | [10.1038/s41586-023-06816-9](https://doi.org/10.1038/s41586-023-06816-9)，GSE202186 | 外源cytokine刺激；不能直接回答固定内源网络最小拆解靶数，未下载 |
| Perturb-CITE-seq | [10.1038/s41588-021-00779-1](https://doi.org/10.1038/s41588-021-00779-1)，SCP1064 | co-culture/CRISPR更相关，但尚未建立control网络—特定target—response的IC配准；未下载218k全集，未作验证结论 |

**真实独立干预验证：NOT_VALIDATED。** 按原指令允许的fallback，本轮仅synthetic proof与刺激sanity；不把“未找到即取的数据”说成数据不存在。固定LR图不含剂量、蛋白定位、补偿和反馈，IC不是可直接开药的靶点数量。

KT4：跨合理资源的通用独立性主张KILL。KT5：有相当排名一致性，但非完全稳定；生物学保持HOLD。对PLAN的影响：H3未满足整体GO所需条件，不能通过挑数据库或换成strength score保留原主张。
''')
write('05_spatial_homology.md',f'''
# Spatial homology

## 数据与单位

10x Mouse Brain Serial Section 1, Sagittal Anterior标准数据，官方filtered H5与spatial包；CellChat作者Figshare条目23621151/file41446839提供既有注释。保留全部**{info['spots']} spots**；原始32285 gene IDs按重复symbol求和后{info['genes']} symbols。作者对象有1073 spots，其中{info['annotated']}条barcode匹配；1条不在此次filtered counts。剩余{info['unannotated']}标作UNANNOTATED，未重新聚类或做deconvolution。

只有**1个生物样本**，不做spatial replicate-level显著性结论，也不把复制坐标块作为重复。作者对象仅用于注释，表达使用官方全filtered矩阵，不用其小型示例表达子集替换全数据。

官方名义v1 spot直径55µm、中心距100µm：[10x术语表](https://www.10xgenomics.com/support/cn/software/space-ranger/3.1/getting-started/space-ranger-glossary)。此数据JSON直径配55µm却得到近邻{info['nearest_um_using_55_over_json_diameter']:.3f}µm；CellChat教程用65µm约{info['nearest_um_using_tutorial_65']:.3f}µm。因此用实际近邻像素间距{info['pixel_pitch']:.3f}与官方100µm pitch标定：{info['coordinate_scale_um_per_pixel']:.6f}µm/pixel。保留冲突及来源，不盲目套直径换算。

## 同一接口和代数

两模态调用 `cci_core.communication`，相同resource loader、sample/sender/receiver/lr_index schema、相同active规则、同一IC和下游样本线性模型。scRNA仅population mixing：Ω=1/(nA nB)，C=mean(L)mean(R)，不补造坐标。空间Ω=K(d)/(nA nB)，归一化分母固定为所有可能A/B配对数，避免稀疏后改分母掩盖机会减少。

Kang为human、Visium为mouse，使用同一CellChatDB家族的物种资源；不冒称两物种是相同基因图。真正同资源的对照是在同一Visium矩阵上分别应用mixing/contact/secreted Ω。其all-to-all dense机会算子恢复mixing，包括complex组件。

- contact：110µm radius、二值、无self loops，表示spot邻接proxy，不能解释成单细胞物理接触。
- secreted：250µm截断exponential，scale100µm、无self loops。
- 正式typed结果按资源annotation：Cell-Cell Contact用contact，Secreted Signaling用secreted。分别把两种核应用全部LR的结果仅是同资源计算对照。

## 结果与数值证据

{table(pd.DataFrame(equiv).fillna('—'))}

全部可映射LR及全部group pairs均与dense comparator对照，而非只抽几条边。Kernel替换后支持通道从{sp['mixing_supported_channels']:,}降至{sp['typed_supported_channels']:,}，移除{sp['removed_by_geometry']:,}、新增{sp['added_by_geometry']}；有支持的group pairs从{sp['active_program_pairs_mixing']}降至{sp['active_program_pairs_typed']}。这支持空间约束确实移除了部分不邻近/不共现的机会；不等于这些被保留通道已在实验中通信。

图的距离矩阵对照、建图和aggregation的分开runtime/存储见 `{S}/profile.tsv`。1/2/4个分离坐标块测到2695/5390/10780 spots，仅检验complexity。正式路径不构造dense N²距离；full dense只用于基准和数值对照。

## Spatial CellChat 语义对照

作者教程采用空间距离与接触范围约束并计算群体概率。本轮实际运行安装版本{spnative['CellChat_version']}，全部2695spots和既有/未知标签，保留全部已测LR与cofactor，20次bootstrap仅作profiling，耗时{spnative['seconds']:.4g}s；native net保存在 `{C['spatial_native']}/net.rds`。这是CellChat 2.2.0的spatial mode，未实跑v3；v3纳入了碰撞审计。本原型的逐边核和传统product与其群体距离/Hill/cofactor模型不同，不要求数值相同，不称为CellChat复现准确度。

KT6：接口、资源表示、IC和数学退化关系通过；空间稀疏同口径优势见工程报告。**Spatial replicate calibration NOT_TESTABLE（只有一切片）**。对PLAN的影响：H4工程同构成立且无需虚构scRNA空间；不支持跨模态生物效应等价或泛化声明。
''')
print('Wrote six-report bundle (collision matrix preserved), benchmark, three diagnostic figures')
