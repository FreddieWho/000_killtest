# 011 后继项目：四路线并行 Kill-Test Prompt

> 适用对象：可访问源码、可编译运行、可进行 profiling/benchmark 的 coding/research agent。
> 
> 本文不是施工方案，不允许直接开始“大重写”。唯一目标是在最低成本下回答：**这条路线是否存在足以支撑一篇有意义方法/系统论文的结构性性能空间？**

---

## 0. 总任务与硬约束

我们已经否决了“为了快而快”的优化路线。STARsolo/NanoFilter 的核心教训是：如果可修改阶段只占端到端运行的一小部分，即使局部 10–100×，最终也可能只剩 <10–20% 总收益，无法形成强 claim。

本轮只允许评估以下四条路线：

1. **GTDB-Tk / pplacer**：现代化 phylogenetic placement execution，保持 GTDB/pplacer 语义尽可能不变。
2. **RepeatModeler2**：寻找 de novo repeat modeling 中可删除的中间文件、重复 pass、单线程 sampling/调度和数据搬运。
3. **HLA typing**：重点审计 HLA*LA / HLA-HD 等仍被使用、但 graph/DP/IO 可能存在高成本关键路径的实现；不要碰已经被近期工作明显优化过的 OptiType 对齐部分。
4. **GATK GenomicsDBImport**：仅作为低优先级保底路线；只有发现 >3× 端到端或数量级 I/O/RAM 改善的架构空间才值得继续。

### 禁止事项

- 不准以“换 Rust/C++ 就会更快”为前提。
- 不准先写完整替代品再看有没有收益。
- 不准把 compiler flags、机器特定参数、NUMA pinning 当成主要创新。
- 不准把科学输出改掉后再用“差不多”掩盖性能收益。
- 不准用微型 toy 数据证明速度；至少要有一个能触发真实瓶颈的现实 workload。
- 不准在 profiling 未证明 Amdahl ceiling 前投入复杂重写。

### 统一继续阈值

除非某路线满足下面至少一条，否则默认 **KILL**：

- 现实 workload 的理论/实测端到端加速上限 **≥3×**；或
- RAM / temporary storage / file-handle 等资源门槛降低 **≥5×**，且该门槛当前确实限制普通用户运行；或
- 能删除一个完整的数据 pass / 大型中间文件阶段，并预计带来 **≥2× 端到端 + 显著资源下降**；或
- 存在一个很清楚的 production pain point（例如 100+ GB RAM 才能运行），即使速度收益略低，也能把任务从 HPC 变成 commodity workstation。

若只能稳定做到 10–50% 小幅加速，除非能证明极大规模部署具有显著成本意义，否则 KILL。

---

# Track A — GTDB-Tk / pplacer

## A1. 核心问题

不要问“pplacer 能否更快”，而要先问：

> **在 2026 年现代 GTDB-Tk 工作流里，pplacer 仍然占多少端到端时间和 RAM？新版 ANI/skani prescreen 是否已经把大多数 genome 从 placement critical path 移走？**

如果现代真实 workload 中只有少数 genome 进入 pplacer，则即使 pplacer 100×，也可能重演 STARsolo。

## A2. 必做 collision audit

检索 2025-01-01 至今：

- GTDB-Tk release/changelog；
- pplacer fork / rewrite；
- EPA-ng / SCAMPP / APPLES / atlas-place 等替代 placement；
- GTDB-specific placement optimization；
- 是否已有“pplacer-compatible / likelihood-equivalent modern engine”。

输出 `A_collision.md`：

| 方法 | 年份 | 是否 GTDB 实际可用 | 是否改变 placement algorithm | RAM | runtime | 与我们潜在路线冲突 |
|---|---:|---|---|---:|---:|---|

只要发现近两年已有团队公开实现“GTDB/pplacer 语义基本不变 + 数量级降 RAM/加速”，本路线默认 KILL，除非还有明显未覆盖瓶颈。

## A3. 现实 workload

至少建立两个 workload：

### A-W1：普通 reference-near genome
用于测 ANI/skani prescreen 后是否根本不进入 pplacer。

### A-W2：novel/environmental MAG / 需要 placement 的 genome
至少 50–200 个，确保进入 bacterial placement；允许从 GTDB 教程、公开 MAG 集合或官方测试数据中取小批量。

不要一开始跑 full GTDB 全树；先用官方 split-tree / 标准 classify_wf 路线重现真实行为。

## A4. Profiling

必须同时记录：

- wall time；
- CPU time；
- peak RSS；
- read/write bytes；
- major/minor page faults；
- pplacer 占 classify_wf 的 wall-time fraction；
- pplacer 内部热点函数；
- partial likelihood vector 内存；
- tree/reference state vs query-specific state；
- allocation/free；
- cache miss（若方便）；
- single-thread vs multi-thread scaling。

Linux 推荐：`/usr/bin/time -v`, `perf record/report`, `pidstat`, `iostat`; 必要时 heap profiler。

## A5. 架构假设，只允许做最小原型

仅当 profile 支持时才测试：

1. **reference-state sharing / memory layout**：reference tree likelihood/state 是否被低效复制或按对象布局造成 cache/RAM 膨胀；
2. **gap-aware / sparse likelihood**：alignment 中大量 gap 是否仍做不必要运算；
3. **SIMD vectorization**：高占比 likelihood kernel 是否适合 AVX2/AVX-512；
4. **query batching**：相同 reference state 对多个 query 是否可以批处理复用；
5. **mmap / compact serialization**：是否可让 reference state 可共享、按需读取；
6. **split-tree scheduling**：是否存在明显 load imbalance / 重复加载。

只做能验证 ceiling 的 prototype，不做完整 replacement。

## A6. 语义一致性

至少比较：

- placement edge；
- likelihood weight ratio / relevant placement statistics；
- GTDB-Tk 最终 taxonomy；
- borderline genomes（多个候选 edge 接近）必须单独列出。

目标不是“分类大体一致”，而是尽可能明确差异来源。若新策略必须改变算法近似才能快很多，要把路线重新归类为“新算法”，不要冒充 infra rewrite。

## A7. PASS / KILL

**PASS** 若满足任一：

- pplacer 在需要 placement 的现实 workload 占 wall time ≥40%，且存在预计 ≥3× pplacer speedup，最终 classify_wf ≥2×；
- peak RAM 可从当前数量级降低 ≥5×，足以从大内存节点降到常规 32–64 GB workstation；
- 能证明一个明确的 reference-state/data-layout 根因导致当前内存异常，并有小原型实证。

**KILL** 若：

- skani prescreen 后 pplacer 已很少进入关键路径；
- RAM 主要由无法压缩的必要 reference state 决定；
- EPA-ng/其他现有方法已经能满足同样语义和资源目标；
- 预计端到端 <2× 且无硬件门槛突破。

---

# Track B — RepeatModeler2

## B1. 核心问题

验证 RepeatModeler2 的长耗时到底来自：

- 真正不可避免的 similarity/all-vs-all computation；
- 还是 sampling、文本中间文件、批次切分、RECON/RepeatScout 之间反复 materialize、调度空转和低并行环节。

只有后者占比足够大，才值得做现代 execution engine。

## B2. collision audit

检索 2025–2026：

- RepeatModeler2 2.0.9 及近期 changelog/issues；
- EDTA、RepeatExplorer、RED、Earl Grey 等是否已经替代核心场景；
- 是否有新的 RepeatModeler-compatible fast rewrite；
- Earth BioGenome / VGP / genome annotation pipelines 当前是否仍推荐或实际调用 RepeatModeler2。

重点不是“有新工具”，而是**RepeatModeler2 是否仍是 de novo repeat library 的主流默认之一**。

## B3. 测试 genome

Kill-test 不要跑 2–3 Gb genome。

至少两档：

- **small**：~100–200 Mb eukaryotic genome，用于完整跑通 round 1–N；
- **medium**：~300–600 Mb genome，用于观察 I/O/temporary scaling；若完整运行预计 >3 h，可只跑到能稳定暴露瓶颈的 round，但必须说明外推方法。

优先使用官方 benchmark/教程中有已知行为的物种（例如 Drosophila / rice 的较小版本或公开 assembly），不下载巨大原始 reads，只使用 assembled genome FASTA。

## B4. profiling 分解

强制把总时间拆成：

\[
T=T_{sampling}+T_{similarity}+T_{parse}+T_{cluster}+T_{IO}+T_{scheduling}+T_{other}
\]

每一 round 都记录：

- wall time；
- CPU utilization；
- peak RSS；
- temporary file count / bytes；
- read/write bytes；
- 每个 external process 的时间；
- 单线程区段；
- 等待 I/O / child-process 的时间。

用 `strace -c` 或 eBPF/`perf` 识别 open/read/write/stat/rename 等系统调用是否异常突出。

## B5. 允许测试的最小架构原型

仅当 profile 支持：

1. parallel deterministic sampling；
2. binary/columnar intermediate 替代巨大文本 HSP/cluster 文件；
3. pipe/stream 直接连接 producer→consumer，避免写盘再读盘；
4. persistent parsed genome / memory-map；
5. round-to-round incremental state，避免重复扫描；
6. task-level scheduler 解决 CPU idle / imbalance；
7. 若 RMBlast 占绝大多数：只测试输入削减/结果 streaming，禁止重新发明 aligner。

## B6. 科学一致性

比较：

- final consensus families 数量；
- family sequence identity；
- class/family assignment（如有）；
- 使用同一 RepeatMasker 对 genome remask 后的 masked fraction / family coverage；
- 任何随机 sampling 需固定 seed 或分析 stochastic variance。

若架构变更使 repeat family discovery 本质改变，则必须停止把项目称作纯 infra。

## B7. PASS / KILL

**PASS**：

- `sampling + parse + I/O + scheduling` ≥40% wall time；或
- temporary data ≥几十 GB/中型 genome，且 prototype 可降低 ≥5×；或
- 可删除至少一整次 genome/HSP materialization pass，预计端到端 ≥3×。

**KILL**：

- >80–90% 时间稳定花在不可避免的 RMBlast/RECON 核心计算；
- 工程层优化 ceiling <2×；
- 更现代主流 pipeline 已经完全绕开 RepeatModeler2，并且其使用量明显衰退。

---

# Track C — HLA typing（HLA*LA / HLA-HD）

## C1. 核心问题

目标不是“发明新的 HLA caller”，而是判断：

> 某个今天仍被广泛使用且准确性有竞争力的 HLA typing backend，是否有一个占端到端大头的 graph/DP/read-processing kernel 可以在**不改变 genotype inference 语义**的情况下大幅现代化？

不要投入 OptiType alignment 路线；近期已有 Columba 等明显优化。

## C2. 活跃度与 collision audit

确认 2025–2026：

- HLA-HD、HLA*LA 在 nf-core、临床/WGS/RNA-seq benchmark 中仍有真实使用；
- T1K、SpecHLA、HISAT-genotype、OptiType 等替代工具的准确率/速度现状；
- 是否已有 HLA*LA/HLA-HD 的 fast fork / rewrite；
- 若某旧工具已被同精度更快工具普遍替代，则不要“救古董”。

输出：

| tool | active? | current users/workflows | typing resolution | typical runtime/RAM | obvious faster substitute? | worth saving? |

## C3. 小型 benchmark

必须准备三种输入，但都用小样本数（1–3 sample）：

1. WES；
2. WGS（可用 chr6/HLA-region extracted BAM/FASTQ 做 kernel profiling，但另需至少一次完整 realistic prefilter path）；
3. RNA-seq。

优先选择 Genome in a Bottle / 1000 Genomes /公开细胞系中有高置信 HLA ground truth 的样本。

## C4. Profiling 问题

回答：

- read extraction/prefilter 占比；
- generic alignment 占比；
- graph/DP 占比；
- allele scoring/inference 占比；
- I/O + BAM sort/index 占比；
- peak RAM 来自 reference graph、DP matrix、read cache 还是中间文件；
- threads 是否真正扩展。

## C5. 允许的原型

只对占比 ≥30–40% 的阶段测试：

- banded / sparse / SIMD DP；
- read routing/minimizer prefilter（必须保证超高 recall，不以损失 genotype accuracy 换速度）；
- graph/reference compact representation；
- batch DP；
- zero-copy BAM/FASTQ parsing；
- 避免 BAM→intermediate→reparse 的 pass fusion。

## C6. 结果一致性

主指标是：

- 2-field / 4-field allele call concordance；
- ambiguous calls；
- homozygous/heterozygous；
- difficult alleles / pseudogene-rich loci；
- runtime / RAM。

任何加速如果靠改变 inference heuristic，必须单独标记，不算纯 infra success。

## C7. PASS / KILL

**PASS**：

- 单一关键阶段占 ≥50% 且 prototype 预计/实测 ≥4×；使端到端有 ≥2–3×；或
- peak RAM 降 ≥5× 且维持 allele call；或
- 能去掉一个完整 alignment/intermediate pass。

**KILL**：

- HLA*LA/HLA-HD 已被同准确率更快 caller 实质替代；
- 热点分散，任何单点优化端到端 ceiling <2×；
- 只有牺牲 typing accuracy 才能获得显著加速。

---

# Track D — GATK GenomicsDBImport（低优先级）

## D1. 本路线只有一个继续理由

它必须不是“Java 慢，重写一下快 30%”，而是发现：

> cohort gVCF import 的数据访问模式存在结构性重复，使 sample-major → locus-major 转置可以通过 streaming / partitioning / direct native path 删除大量 random seek、反复 open/close、fragment materialization。

否则直接 KILL，不在此项目投入。

## D2. collision / baseline

必须比较：

- GATK `GenomicsDBImport`；
- GenomicsDB native `vcf2genomicsdb`；
- 若适用，bcftools/GLnexus 等不完全同语义工具，仅作为工程参照，不混淆科学功能。

若 native importer 已经解决大多数开销，GATK wrapper 的优化上限不足，则 KILL。

## D3. 轻量数据设计

构造或下载小型公开 gVCF cohort：

- 50 samples；
- 200 samples；
- 若本地资源允许再做 500–1000 samples；
- 仅 chr20 / chr22 或 exome interval，避免 kill-test 本身变成多小时任务。

同时测试：

1. local NVMe；
2. 若方便，普通共享/NFS 文件系统。

## D4. Profiling

记录：

- open/close/seek 次数；
- bgzip/tabix decode；
- Java vs native 时间；
- interval 数变化对 runtime 的影响；
- batch size 对 RAM/runtime；
- fragment 写入/consolidation；
- sequential read fraction；
- file handle 数。

## D5. 最小原型

若重复 interval I/O 是主因：

做一个**只支持 chr20、小 cohort 的 streaming transposition proof-of-concept**：

```text
sample gVCFs sequentially decoded once
        ↓
partitioned locus buffers
        ↓
direct bulk write / compact intermediate
```

不要求兼容全部 VCF corner cases，只用于测理论 ceiling。

## D6. PASS / KILL

**PASS 仅当**：

- 端到端 ≥3×；或
- RAM / file handles / shared-FS I/O 降 ≥5×；或
- cohort scaling 从明显 superlinear/seek-bound 变为接近 linear，并能形成 population-scale claim。

否则 KILL。此路线论文价值门槛高于其他三条。

---

# 统一输出格式

最终只提交以下文件，不要提交大段无结论日志：

1. `00_EXECUTIVE_DECISION.md`
2. `A_pplacer/REPORT.md`
3. `B_repeatmodeler/REPORT.md`
4. `C_hla/REPORT.md`
5. `D_genomicsdb/REPORT.md`
6. `benchmarks/summary.tsv`
7. `profiles/`（必要的 flamegraph/原始 profile；大文件不要入 git）

`00_EXECUTIVE_DECISION.md` 必须包含：

| Route | Current usage | Collision risk | Critical-path share | Theoretical ceiling | Prototype result | Scientific equivalence risk | Decision |
|---|---|---|---:|---:|---:|---|---|

Decision 只能使用：

- **GO**：值得下一阶段施工；
- **HOLD**：有空间，但证据不足/需特定硬件或数据再测；
- **KILL**：停止投入。

最后强制选出最多 **1 条 GO 主路线**；不要为了让项目“有结果”而把多个弱路线都标 GO。

## 论文级判断标准

GO 路线必须能用一句非技术化语言解释价值，例如：

- “把原本需要 100+ GB RAM 的标准分类步骤降到普通工作站可运行。”
- “把原本以天计的 repeat modeling 中大量中间数据 pass 删除，缩短到小时级。”

如果一句话只能写成：

> “我们优化了实现，使其快 35%。”

则默认不适合作为 011 后继项目。

