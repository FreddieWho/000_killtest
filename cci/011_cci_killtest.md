# 011 后继项目：CCI/CCC Statistical Redesign Kill-Test Prompt

> 目标：在正式开发前，以最低成本验证一个候选方向是否同时满足：
> 1) 统计问题真实存在；2) 不是 FastCCC/LIANA+/MultiNicheNet/CellChat v3 等近期工作的重复；3) scRNA 与 spatial 可以由同一个数学对象自然覆盖；4) 计算复杂度从原理上更低；5) 新统计量 **Intervention Complexity** 提供 count/strength 之外真正独立、可实验验证的信息。
>
> 本文是 kill-test，不是完整论文施工。任何一个核心假设失败，都允许直接杀掉对应模块。

---

# 1. 不可漂移的科学问题

现有 CCC/CCI 工具长期主要输出：

- interaction count / number；
- interaction strength / probability / magnitude；
- significance / specificity。

本项目不以“再造一个更准的 strength”作为创新。

需要同时检验两个候选贡献：

## Contribution A — Replicate-aware, modality-homologous CCC inference

核心原则：

> **cells are observations; biological samples/patients are replicates.**

单细胞数不应该虚增实验自由度。目标是在 sample/patient 层做 inference，同时让 scRNA 与 spatial 共用同一套 communication algebra。

## Contribution B — Intervention Complexity (IC)

不再只描述“有多少 interaction / interaction 多强”，而回答：

> **这套 sender→receiver communication 要阻断多少个分子靶点才能被真正拆掉？**

这是结构性、干预导向统计量；不是另一个 continuous strength score。

---

# 2. 必做 collision audit（先于编码）

检索 2024-01-01 至今，至少覆盖：

- CellPhoneDB v5；
- CellChat / Spatial CellChat / CellChat v3；
- FastCCC；
- LIANA+；
- MultiNicheNet；
- Tensor-cell2cell v1/v2；
- dominoSignal；
- scDiffCom；
- SpatialDM；
- COMMOT；
- scRICH；
- Renoir；
- MOSANIC；
- 任意 differential CCC / sample-aware CCC / perturbation-aware CCC / network vulnerability / minimum-cut / hitting-set CCI 方法。

输出 `collision_matrix.md`：

| Method | scRNA | spatial | biological replicates | arbitrary design matrix | paired/longitudinal | analytical/no cell permutation | count | strength | heterogeneity | perturbation/vulnerability | minimum intervention set |
|---|---|---|---|---|---|---|---|---|---|---|---|

### Collision kill rule

若已有方法同时做到：

1. sample/patient 是 inference unit；
2. 支持一般 design / paired / longitudinal；
3. 同一数学框架覆盖 scRNA + spatial；
4. 明确定义并计算“破坏 sender→receiver LR 网络所需最小分子集合”或等价指标；
5. 已有真实 perturbation 验证；

则本方向整体 **KILL**。

如果只撞其中一部分，必须在报告中清楚重定义剩余 gap，不得夸大“首次”。

---

# 3. 核心数学原型：同一个 Opportunity Operator

只做最简单、最透明的原型。不要先上 GNN/Transformer。

设每个 sample 中：

- ligand expression/support vector：`L`
- receptor expression/support vector：`R`
- communication opportunity operator：`Ω`

基础 communication quantity：

\[
C = L^T \Omega R
\]

这只是统一代数对象，不要求把它包装成新的 strength。

## 3.1 scRNA mode

scRNA 没有真实 geometry。不得伪造空间关系。

`Ω_scRNA` 只表示 **population mixing / exchangeability**。在 cell-type 聚合后，本质退化为 sender/receiver sample-level sufficient statistics：

\[
C_{s,A\to B,k}=g(\bar L_{s,A},\bar R_{s,B})
\]

Kill-test 第一版 `g` 请直接采用最简单/传统形式（product 或 geometric mean），不要发明新 score。

## 3.2 spatial mode

空间数据：

\[
\Omega_{ij}=K(d_{ij})
\]

至少测试两类 kernel：

- contact：稀疏邻接 / radius graph；
- secreted：截断 exponential/Gaussian kernel。

必须构建 sparse graph，禁止 dense `N×N` distance matrix 作为正式实现。

复杂度目标：

\[
O(|E|)\quad \text{而不是}\quad O(N^2)
\]

### 模态统一的判据

“统一”不是指两头结果一样，而是：

- 同一接口；
- 同一 LR database；
- 同一 per-sample output schema；
- 同一 downstream replicate-level inference；
- 仅 `Ω` 的观测信息不同。

若 scRNA 必须引入大量猜测 spatial topology 才能工作，则禁止，改成 scRNA-only 主线 + spatial extension。

---

# 4. Replicate-aware inference：真正的第一 Kill-Test

## 4.1 问题

检验经典 pooled-cell / cell-label-permutation CCC 在 multi-donor data 上是否发生 pseudoreplication / type-I error inflation。

真正 experimental unit：

\[
N_{rep}=N_{patient/sample}
\]

不是：

\[
N_{cell}
\]

## 4.2 原型流程

1. 每个 `sample × celltype × gene` 建 sufficient statistics：
   - raw sum / library size；
   - mean；
   - variance；
   - detection rate；
   - number of cells。
2. 对每个 sample 得到 `sender × receiver × LR` communication statistic。
3. downstream inference 在 sample 层拟合：

简单版：

\[
C \sim condition
\]

paired 版：

\[
C \sim condition + patient
\]

可选 mixed model：

\[
C \sim treatment*time + (1|patient)
\]

4. 第一版允许用 limma/GLS/LMM，不要一开始发明复杂 likelihood。

## 4.3 Type-I error simulation

从真实数据生成 null，保持：

- donor heterogeneity；
- cell-type composition；
- library size；
- cell-number imbalance；
- gene–gene covariance。

但打乱/消除 condition effect。

至少做 100–500 null replicates（kill-test 够用，不需正式论文规模）。

比较：

- CellPhoneDB（可降低 permutation 数做 profiling；正式 type-I 用可承受设置）；
- CellChat；
- FastCCC；
- LIANA+/sample-aware baseline；
- 本原型。

核心指标：

- empirical FPR at 0.05；
- FDR；
- QQ plot / p-value uniformity；
- FPR 随 cells per donor 增加是否恶化；
- runtime。

### A-module PASS

若 pooled-cell 方法在 realistic multi-donor imbalance 下明显 inflation，而 sample-level 原型接近 nominal，并且现有 sample-aware 工具不能自然完成一般 design，则 PASS。

### A-module KILL

若 LIANA+/现有方法已能以同等简单、正确方式完成上述一般 inference，且本方法只剩速度差异，则该统计贡献 KILL，可保留工程部分。

---

# 5. 新统计量：Intervention Complexity (IC)

## 5.1 目标

对某个 sample、sender A、receiver B，得到当前支持的 ligand–receptor 图。

传统指标：

- Count：多少条 LR 通道；
- Strength：总体表达/通信强度。

我们新增：

> **需要阻断最少多少个 molecular nodes，才能破坏全部（或指定比例）当前支持的 LR 通道？**

## 5.2 最干净的主定义：IC100

先限制到 simple one-ligand ↔ one-receptor interactions，构建二部图：

\[
G=(L,R,E)
\]

每条 active LR 是一条 edge。

定义：

\[
IC_{100}=\min |S|,\quad S\subseteq L\cup R
\]

使每一条 active edge 至少有一个端点在 `S` 中。

这就是 **minimum vertex cover**。

二部图上可通过 maximum matching 精确、快速求解（Kőnig theorem）。

### 为什么它不是 count/strength 的重复

两个 network 可以有相同 edge count、相同 total strength，但：

- 网络 A 所有边共享一个 receptor → `IC100=1`；
- 网络 B 是多个独立 LR 轴 → `IC100` 很高。

它测的是 **molecular bottleneck / intervention redundancy**。

## 5.3 复杂体升级

CellPhoneDB/CellChat 中存在：

- heteromeric ligand；
- heteromeric receptor；
- cofactor。

对此不要强行塞进普通二部图。定义 hyperedge / target set：

- 若复合体任一必需 subunit 被阻断即可使 interaction 失活，则 interaction 对应一个“可击中节点集合”；
- 求覆盖全部 active interactions 的 minimum hitting set。

Kill-test 阶段：

- 先对 simple LR 用 exact MVC；
- 再对 top 100–500 complex-containing programs 用 ILP / OR-Tools exact hitting set；
- 若复杂度高再测试 greedy approximation。

## 5.4 可选 IC90

为了避免一个极弱 edge 迫使 IC100 增加，可以探索：

\[
IC_{90}
\]

即破坏 90% **supported channels** 所需最少靶点。

注意：主版本优先按 edge/support 计数，不按原 strength 加权，以免退化成 strength 的函数。

strength-weighted IC 只能作为 secondary analysis。

## 5.5 新统计量独立性 Kill-Test

在所有 sender→receiver programs 上计算：

- interaction count；
- total strength；
- IC100；
- IC90（若实现）。

要求：

1. IC 与 count 不应近似单调等价；
2. IC 与 strength 不应高相关到可被简单线性/单调函数替代；
3. 找出明确的四象限案例：
   - high strength / low IC；
   - high strength / high IC；
   - low strength / low IC；
   - low strength / high IC。

若 IC 基本只是 count 的 rescale，则 KILL。

---

# 6. Intervention validity：最低成本验证

Kill-test 阶段不要求自己做湿实验。

优先寻找公开 perturbation 数据或明确 ligand/receptor stimulation data，测试：

- low-IC program 是否更容易被单一 perturbation 显著打掉；
- high-IC program 是否单靶 perturbation 后仍保留更多 downstream response；
- 若没有干净数据，仅做 in-silico/synthetic proof，并明确“尚未验证 biology”，不要制造结论。

正式项目若继续，才设计 3–5 条 co-culture/blockade 实验。

---

# 7. 小型、经典、可快速下载的数据集

不要为了 kill-test 下载百万细胞 atlas。优先下面 3 个。第 1 和第 3 必做，第 2 可选。

## Dataset 1 — Kang et al. 2018 PBMC IFN-β（必做，scRNA replicate benchmark）

**用途**：

- paired donor design；
- pseudoreplication/type-I error；
- differential CCC；
- runtime；
- IFN stimulation 提供一个非常干净的 biology sanity check。

**规模**：约 24.7k cells、15.7k genes、8 donors × control/stim，共 16 sample units。LIANA/scverse 处理版首次下载约 40 MB。

**推荐获取方式 A（最方便）**：

```python
import liana as li
adata = li.ds.kang_2018()
```

返回对象包含 `sample`, `condition`, `patient`, `cell_abbr` 等字段。

**推荐获取方式 B（固定 h5ad mirror）**：

```python
import scanpy as sc
adata = sc.read(
    "kang_2018.h5ad",
    backup_url="https://ndownloader.figshare.com/files/34464122",
)
```

**原始来源**：GEO `GSE96583`。

GEO RAW bundle ~72.7 MB；不要下载 SRA FASTQ，kill-test 没必要。

**使用要求**：

- 不重新聚类；直接用已有 cell type / donor / condition；
- 第一轮只保留 4–6 个主要 cell types 可进一步提速；
- 构造 paired design：`~ condition + patient`；
- 用 donor-level label shuffle / within-pair sign flip 做 null。

**已知性质**：IFN-β 6h stimulation 是经典强 perturbation；适合检查 IFN-response 相关 signaling 是否至少方向合理，但不要把它当作真正 endogenous sender→receiver blockade ground truth。

---

## Dataset 2 — CellChat human skin NL/LS（可选，经典 compatibility benchmark）

**用途**：

- 与 CellChat 经典教程直接对照；
- 检查 count/strength 网络是否大体复现；
- 测 IC 与传统 count/strength 的独立性；
- 不把它当作主要 replicate-level statistical proof。

**数据来源**：CellChat 官方教程使用的 `data_humanSkin_CellChat.rda`。

```r
load(url("https://ndownloader.figshare.com/files/25950872"))
```

其中包含 normalized expression 和 metadata（condition、patient id、cell labels）。官方教程中 LS 子集约 5,011 cells × 17,328 genes。

**使用要求**：

- 第一轮只跑官方 CellChat DB 中 Secreted Signaling + Cell-Cell Contact；
- 不做复杂重新预处理；
- 先复现 CellChat tutorial 的主要 sender/receiver pattern，再计算 IC；
- 若数据的有效 biological replicate 结构不足，不得用它证明 type-I calibration。

---

## Dataset 3 — 10x Visium mouse brain sagittal anterior（必做，small spatial benchmark）

**用途**：

- 验证同一 `Ω` operator 的 spatial specialization；
- sparse radius/KNN graph；
- 比较 dense pairwise vs sparse computation；
- 与 Spatial CellChat / CellChat 官方 spatial tutorial 做语义对照。

**官方数据名**：10x Genomics `Mouse Brain Serial Section 1, Sagittal Anterior, Visium`（标准 1.0.0 数据集）。CellChat 官方 spatial tutorial 直接使用该数据。

只下载：

- filtered expression matrix / H5；
- `spatial/` coordinates + scale factors；
- 必要的 low-res image（若代码需要）。

**不要下载 BAM。**

CellChat tutorial reference：

`https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.Rmd`

10x dataset 页面可从该 tutorial 的链接进入。

**使用要求**：

- 直接采用教程已有 spot/cell-type annotation，或使用教程 processed object（若能直接取得）；
- kill-test 不做 deconvolution 方法比较；
- 构建半径 graph，目标 `|E| << N^2`；
- contact 与 secreted signaling 使用不同 kernel/radius；
- 记录 dense distance、sparse graph construction、communication aggregation 三部分 runtime/RAM。

---

# 8. Baseline 最小集合

不要第一轮装十几个工具。

必做：

1. **CellPhoneDB**：传统 permutation reference；
2. **FastCCC**：解析/FFT speed reference；
3. **CellChat / Spatial CellChat**：strength + spatial reference；
4. **LIANA+ 或其 sample-aware workflow**：multi-sample/differential reference。

可选，仅 collision audit 或必要时执行：

- MultiNicheNet；
- dominoSignal；
- scDiffCom；
- SpatialDM；
- COMMOT；
- scRICH；
- MOSANIC。

若安装成本高，先读方法/代码并仅运行最直接竞争者，不要让 kill-test 被环境配置吞掉。

---

# 9. 工程实现原型

目标不是做完整 package，而是证明复杂度与数据流成立。

建议 Python 原型即可；Rust/C++ 暂时不要作为前提。

## 9.1 数据结构

```text
Input AnnData / sparse matrix
       ↓
per-sample × celltype sufficient statistics
       ↓
LR resource indexed by integer gene IDs
       ↓
per-sample sender→receiver active-LR graph
       ↓
replicate-aware inference + IC
```

Spatial：

```text
coordinates
   ↓
radius/KNN sparse graph Ω
   ↓
sparse aggregation
```

避免 pandas 大量 row-wise loop；优先 sparse matrix / NumPy / Polars/Arrow（哪个最简单就用哪个）。

## 9.2 profiling

记录：

- read/load；
- aggregation；
- LR lookup；
- opportunity operator；
- statistical fit；
- IC graph construction；
- MVC / hitting-set solve；
- total peak RAM。

对 Kang 做 cell-number scaling：

- 2k；
- 5k；
- 10k；
- full ~25k；

可通过 resampling/duplication构造 50k/100k synthetic scaling，但不用于 biological result。

Spatial 做 spot/cell scaling时可复制坐标块或抽样，仅测 complexity。

---

# 10. 关键 Kill Tests

## KT1 — Statistical correctness

问：pooled-cell / cell-label methods 是否在 multi-donor null 下产生明显 FPR inflation？

PASS signal：

- baseline FPR 明显偏离 0.05，且随每 donor cells 增多恶化；
- sample-level method 接近 nominal。

若没有这个现象，不要强行写“wrong experimental unit”。

## KT2 — General design advantage

至少跑通：

```text
~ condition + patient
```

并在 synthetic metadata 上跑：

```text
~ treatment * time + batch
```

如果现有 LIANA+/其他工具已经同样自然完成 interaction-level beta/SE/p/FDR，则本贡献弱化。

## KT3 — Speed / complexity

目标不是击败 FastCCC 单样本 2×。

目标场景：multi-sample / repeated design。

要求证明：

- scRNA inference after aggregation 基本不再随 cell number 线性爆炸；
- spatial 使用 sparse `Ω`，不形成 dense `N²`；
- 相比“每 sample 单独 CellChat/CPDB + second-stage stats”，现实 multi-sample workload 至少有数量级潜力，或至少 >5× 且 RAM 更低。

## KT4 — IC is not count/strength

计算相关性和四象限案例。

最低要求：

- 存在大量 count 相似但 IC 显著不同的 network；
- 存在 strength 相似但 IC 显著不同的 network；
- IC 不应被简单 count/strength 单变量模型高精度预测。

若 `R² > 0.9` 且几乎单调，IC KILL。

## KT5 — Database sensitivity

至少比较两个 LR resource（例如 CellChatDB vs OmniPath/CellPhoneDB resource 的交集或可映射子集）。

问：IC 排名/高低是否完全由数据库边数决定？

若更换合理数据库后 IC 大幅随机翻转，必须降低 claim 或 KILL。

## KT6 — Spatial homology

在 Visium dataset 上验证：

- 同一接口、同一 LR graph、同一 IC；
- 只有 `Ω` 变成 spatial adjacency/kernel；
- sparse 实现具有明确速度优势；
- spatial constraint 应使部分远距离/不相邻 interaction 消失，这是 sanity check。

若为了让 spatial 有意义必须重写一套完全不同算法，则停止“统一双模态”主 claim，专注 scRNA。

---

# 11. 决策矩阵

最终输出：

| Module | Evidence for gap | Collision risk | Statistical validity | Speed ceiling | Biology interpretability | Decision |
|---|---|---|---|---|---|---|
| replicate-aware inference | | | | | | GO/HOLD/KILL |
| scRNA opportunity operator | | | | | | |
| spatial opportunity operator | | | | | | |
| IC100/IC90 | | | | | | |
| unified package | | | | | | |

### 整体 GO 条件

至少同时满足：

1. **统计层**：发现现有主流/常用 workflow 在 multi-replicate design 中存在真实空白，或我们的 general-design inference 明显更直接/正确；
2. **计算层**：multi-sample 或 spatial workload 有 >5× 的现实空间，最好来自消除 permutation/dense pair enumeration，而不是语言重写；
3. **IC 层**：IC 与 count/strength 提供明显正交信息，并对 LR database 不至于极度不稳定；
4. **碰撞层**：没有近期论文已经把 replicate-aware unified operator + minimum molecular intervention set 作为完整方法做掉。

### 整体 KILL 条件

满足任一即可考虑停止：

- 现有 sample-aware CCC 已完整支持一般统计设计，本方法只剩换 API；
- FastCCC/其他方法已经覆盖我们全部 speed claim；
- scRNA/spatial 统一只能靠大量额外假设，导致一头明显劣化；
- IC 本质等价于 interaction count 或 database degree；
- 新统计量没有明确可实验解释。

---

# 12. 交付物

只交付以下内容：

1. `00_DECISION.md` — 一页结论，明确 GO/HOLD/KILL；
2. `01_collision_matrix.md`；
3. `02_statistical_killtest.md` — type-I/FDR/power；
4. `03_engineering_profile.md` — runtime/RAM/complexity；
5. `04_intervention_complexity.md` — 定义、算法、独立性、database sensitivity；
6. `05_spatial_homology.md`；
7. `results/benchmark.tsv`；
8. 最小可重复代码。

不要写正式 paper，不要过度包装。如果某一假设失败，要直接写“失败，因此杀掉/降级”，不要为了保项目而换 claim。

---

# 13. 数据来源核验备忘（agent 必须重新核对）

- Kang et al. 2018 GEO: `GSE96583`; processed LIANA loader约 25k cells/8 donors×2 conditions。
- LIANA loader docs: https://liana.readthedocs.io/en/stable/generated/liana.ds.kang_2018.html
- Kang processed figshare mirror: https://ndownloader.figshare.com/files/34464122
- CellChat human skin tutorial data: https://ndownloader.figshare.com/files/25950872
- CellChat spatial tutorial: https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.Rmd
- 10x Visium mouse brain sagittal anterior：从上述官方 tutorial 指向的 10x 页面下载 filtered matrix + spatial only。

若链接失效，用 accession/dataset name 找官方镜像；不要换成完全不同的超大数据集。

