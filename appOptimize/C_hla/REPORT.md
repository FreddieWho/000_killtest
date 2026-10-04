# C — HLA typing

2026-10-02：**HOLD，补测完成，仍有范围缺项，无 GO**。完整 WGS/WES prefilter 与 HLA-LA 分型已完成。此前单线程回放因未使用官方 src 工作目录而找不到 hla_nom_g.txt；本次在 `profiles/C_scaling_20261002` 使用正确目录回放，进展与结果见末节，失败证据保留。

## 基线、一次性成本和真值

同一 NA12878 公开细胞系：历史 low-coverage WGS、WES，另已准备 RNA FASTQ。WGS 不是现代 30x WGS；同一 donor 的两个 assay 不是独立外部样本。图准备 7107.22 s；输入 mapped+unmapped BAM 合并另行计时，未把区域 mini 输入当完整 prefilter。

| 完整 HLA-LA 分型 | wall s | backend s | 首次 BWA index s | BWA mem s | backend 其余 s |
|---|---:|---:|---:|---:|---:|
| WGS | 9244.61 | 9064.10 | 6713.81 | 1219.59 | 1130.70 |
| WES | 2002.81 | 1903.13 | 0 | 429.05 | 1474.08 |

首次 WGS 包含共享参考 BWA 建索引，后续线程回放复用该索引；不得用冷启动基线与复用索引的回放计算“新引擎加速”。backend 其余包含 graph、inference、sort 等，不是纯 DP 占比。WES 其余约占 73.6%，需要函数证据进一步区分。实际回放选择较慢 WGS 后端，不能外推 WES 或全人群线程扩展。

公开 20181129 两字段真值覆盖 A/B/C/DQB1/DRB1 共 10 个 allele。严格保留原始 ambiguity sets：WGS 明确匹配 1/10、含糊集合兼容 5/10；WES 明确匹配 1/10、集合兼容 9/10。兼容不是明确正确，两者不能混为准确率；没有唯一四字段真值。该小型输入不能建立现代工具的准确性竞争力，也不能用它宣称临床有效。

证据：`profiles/closure_20261002/C_baseline_evidence.json`，`profiles/C/*/analysis/exec_intervals.json`，原始 bestguess 和 G-group calls 均保留。线程前后调用一致性将另报；不把 G-group 名称截断后当唯一两字段 genotype。

## 活跃度与替代工具


| tool | active / usage evidence | 分辨率与比较边界 | runtime / RAM | worth saving? |
|---|---|---|---|---|
| HLA*LA | nf-core 当前 dev workflow 集成；2025 benchmark 仍纳入 | 主要 G-group，不能无依据扩成唯一 4-field call | 见上表完整基线 | 待 profile；不可认定已被普遍替代 |
| HLA-HD | 2025 consHLA 临床研究采用；nf-core 支持用户提供软件包 | 3-field 等；不能与 class-I-only caller 等价比较 | consHLA 有长耗时，但不是本机基线 | 待合法软件包及 profile |
| T1K | 原作者 RNA/WES/WGS 评估；HLA/KIR | 推断算法不同 | 不引用跨硬件数字作为本机加速比 | 强 comparator，不是同语义 drop-in |
| SpecHLA | nf-core 当前 dev 集成；有全分辨率研究 | 组装/推断路径不同 | 本轮未测 | 必须比较 loci、输入与 ambiguity |
| HISAT-genotype | 现有 graph-based 方法 | 科学推断不同 | 本轮未测 | 未证明普遍取代 HLA*LA/HLA-HD |
| OptiType | 活跃 comparator | 主要 class I | 本任务明确排除对齐优化路线 | 不施工 |

来源：[nf-core usage](https://nf-co.re/hlatyping/dev/docs/usage)、[2025 benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC11839181/)、[T1K primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10519407/)、[SpecHLA primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10545945/)。nf-core dev 中存在模块不等于生产部署数量。

## 已找到的碰撞

[consHLA, 2025](https://link.springer.com/article/10.1186/s12859-025-06223-z)已使用 Bowtie2 HLA read prefilter；文中平均 runtime 降幅为 WGS 31%、RNA 27%，不是数量级加速，也不是我们的测量。该研究证明 HLA-HD 仍有实际使用；不支持“更快工具已普遍淘汰它”。

[HLA-LA 官方 1.0.8-fast](https://github.com/DiltheyLab/HLA-LA/releases/tag/1.0.8-fast)是 2024 预发布。已下载相对 v1.0.4 的四个 commit diff：新增批量 BAM/sampleID、`fastHLAReadExtraction`、按 HLA intervals 选择 reads。它修改 read selection；没有独立结果前不能称纯数据布局、完全等价或已获得固定加速比。证据：`sources/HLA-LA_fast_compare.json` 与 `sources/HLA_fast_diff/`。


## 明确缺项与对 PLAN 的影响

HLA-HD/RNA 分型 **BLOCKED_SOFTWARE**，未用其他算法冒充补齐。2026-10-02 核查的[官方下载请求页](https://www.gckyoto.com/https/wwwgenomemedkyoto-uacjp/hla-hd)仍通过申请提供软件；用户尚未提供安装路径。RNA 原始输入已准备，不等于 RNA benchmark 已运行。

唯一四字段 truth 不可用，且原型 **NOT_RUN**。H-C 的执行成本已有实测，但目前无保持 inference 语义且达到门槛的优化实证；一部分高耗时是已有索引复用能避免的一次性准备，不能作为新贡献。保留 HOLD，并明确双工具/三输入范围尚未全部完成。

## 2026-10-02 完成后核查

1/8 线程回放和 raw allele/G-call 比较均已成功，但 perf_8 退出 134：比对、排序后，`mapper/processBAM.cpp:1214` 触发 `sequencesStream.is_open()` 断言。该段描述旧 perf_8 失败现场，已保留。恢复采样采用 perf_8_retry_20261002，下方结果来自恢复后的完整运行；不再将旧失败片段当完整热点。

<!-- closure-results -->
## 本次补测结果

控制器状态：`COMPLETED_SCHEDULED_COMPUTATION`。

1 线程实际后端回放：wall 5133.00 s，peak RSS 30131964 KiB。

8 线程实际后端回放：wall 1012.01 s，peak RSS 30196536 KiB。

1/8 线程 wall 比 5.072×；同一现有二进制线程对照，不能称新原型加速。

基线与单线程/8 线程 raw allele sets 精确相同：True / True；G calls 精确相同：True / True。

- WGS_baseline：两字段明确匹配 1/10，歧义集合兼容 5/10；五个位点明确 genotype 匹配 0/5。
- WES_baseline：两字段明确匹配 1/10，歧义集合兼容 9/10；五个位点明确 genotype 匹配 0/5。
- WGS_kernel_1：两字段明确匹配 1/10，歧义集合兼容 5/10；五个位点明确 genotype 匹配 0/5。
- WGS_kernel_8：两字段明确匹配 1/10，歧义集合兼容 5/10；五个位点明确 genotype 匹配 0/5。

四字段独立真值仍不可用。完整原始歧义集合与每个位点的杂合/纯合/不确定状态见 call_comparison.json。

恢复后函数采样的调用比较：{"R1_bestguess.txt": {"baseline_rows": 34, "perf_rows": 34, "allele_calls_exact": true, "file_bytes_exact": false}, "R1_bestguess_G.txt": {"baseline_rows": 26, "perf_rows": 26, "allele_calls_exact": true, "file_bytes_exact": false}}

函数采样（49 Hz 用户态 CPU，非端到端 wall 占比）：

```text
19.76%  bwa       bwa                               [.] ksw_u8
15.27%  bwa       bwa                               [.] bwt_occ4
12.46%  bwa       bwa                               [.] ksw_extend2
11.34%  bwa       bwa                               [.] bwt_sa
5.91%  bwa       bwa                               [.] bwt_2occ4
2.72%  bwa       bwa                               [.] bns_get_seq
2.11%  bwa       libc.so.6                         [.] __default_morecore
2.10%  bwa       bwa                               [.] bwt_occ
2.08%  HLA-LA    libz.so.1.2.13                    [.] inflate_fast
1.85%  bwa       bwa                               [.] ks_introsort_mem_ars
1.85%  bwa       bwa                               [.] ksw_global2
1.56%  bwa       bwa                               [.] bwt_extend
1.55%  bwa       bwa                               [.] bwt_smem1a
1.51%  samtools  libdeflate.so.0                   [.] deflate_compress_lazy
1.41%  bwa       bwa                               [.] ks_introsort_mem_ars2
1.08%  bwa       bwa                               [.] mem_chain_flt
0.92%  HLA-LA    HLA-LA                            [.] std::_Rb_tree<Node*, Node*, std::_Identity<Node*>, std::less<Node*>, std::allocator<Node*> >::find
0.76%  bwa       bwa                               [.] mem_chain2aln
```

原型仍为 NOT_RUN；本轮只验证现有程序和证据完整性。新架构收益未得到验证，维持 HOLD，不给 GO。

