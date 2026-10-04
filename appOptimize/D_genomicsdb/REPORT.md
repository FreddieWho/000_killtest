# D — GATK GenomicsDBImport

2026-10-02 收尾结论：**KILL 当前新导入引擎方向**。50/200 个真实样本、完整 chr22 的三组导入均成功，已有 native 路径与 reader 选项未留下已证明达到本任务门槛的独立架构机会。结论限于已测旋转盘，不宣称所有硬件下不存在空间。

## 完整导入基线

GATK 4.6.2.0 内嵌 GenomicsDB 1.5.5-addd795；独立 native 为 1.5.4，因此独立 CLI 与 GATK 不是纯 wrapper 消融。Java/native reader 对照固定相同 GATK 与底层版本。每组 8 CPU、64 GB 容器上限；GATK Xmx12g、batch size 50、reader threads 1，输入缓存未控制。

| 样本数 | 方法 | wall s | peak RSS KiB |
|---|---|---:|---:|
| 50 | native_init | 1.04 | 21072 |
| 50 | native_import | 1334.87 | 510580 |
| 50 | gatk_java | 2814.93 | 1446988 |
| 50 | gatk_native_reader | 2006.95 | 723720 |
| 200 | native_init | 15.56 | 20856 |
| 200 | native_import | 9549.00 | 1184204 |
| 200 | gatk_java | 10101.00 | 2999916 |
| 200 | gatk_native_reader | 8229.00 | 724728 |

独立 native 端到端导入比较需加初始化：50 样本 1335.91 s，200 样本 9564.56 s。相对 Java 的速度比为 2.107×、1.056×；GATK 内部 bypass reader 为 1.403×、1.227×。200 样本 Java/bypass 的 RSS 比为 4.139×，未达到 5×，且是已有选项。并发运行和未控制缓存使这些数字不适合推广为稳定跨机器收益。

样本数增加 4 倍时各方法增长倍数不同，不能仅凭两个规模拟合 population-scale 复杂度。没有新 streaming 原型，也没有用现有 reader 配置充当创新。

## 输出差异已经定位的范围

同一 GATK reader 导出五个 1001 bp 区间：50 样本 3412 条记录、200 样本 4437 条记录。GATK 两种 reader 规范化内容一致；独立 native 分别有 34、53 条差异，所查范围均仅为 INFO/DP，未见 FORMAT/genotype/phase 差异。

独立 native `vidmap.json` 为 DP 配置 `VCF_field_combine_operation=sum`，GATK 没有。在影子工作区只移除这项、底层数组及原始基线保持不变，再用同一 reader 导出，两个规模分别与 GATK **逐记录、逐样本规范化内容完全相同**。200 样本规范化 SHA256：`573d6355fc2e836ea1a736217c220cd1831e8dc12ca97aaf02ee581a3a273965`。

证据：`profiles/D_schema_diagnostic_20260928/comparison.json`、`profiles/D_schema_diagnostic_20261002/comparison.json`，以及两个原始导出比较。该因果诊断只覆盖五个小区间，**不是全染色体/下游 joint calling 等价性证明**，也不是加速原型。

## 强碰撞：已有 native 路径


[GATK performance guidelines，2025 更新](https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines)已覆盖 shared POSIX 优化、batch size、interval merge、contig partitions 和 `--bypass-feature-reader`。该指南报告 bypass reader 的典型收益约 10–15%，不能当作本轮测量或普遍上限。

[GenomicsDB 官方 ETL 文档](https://genomicsdb.readthedocs.io/en/latest/import-etl.html)推荐 native `vcf2genomicsdb`。v1.5.4 源码 `src/main/cpp/src/loader/tiledb_loader.cc` 已使用 circular/ping-pong buffers；不能把“从 sample 到 locus 的 streaming 转置”这一概念本身当全新设计。必须进一步证明存在未被这些实现覆盖的重复访问。

GLnexus 的 joint calling 与 bcftools 的 merge 可作工程参照，不是 GATK joint-genotyping 科学语义的直接替换。


## 数据与未完成项

使用 1000 Genomes DRAGEN 3.5.7b 的 200 个不同样本完整 chr22 gVCF 子集；逐样本源 URL、记录数、字节数及哈希在 `data/D_genomicsdb/*.receipt.json`，已登记数据索引。既有公共 GRCh38 参考只读复用。

原任务的 **NVMe 对照 NOT_RUN（当前主机无对应设备）**；NFS、可选 500–1000 样本未运行。interval 数与 batch-size 扫描、全流程 open/close/seek 次数、decode 和 sequential-read fraction 尚未完成；过程采样记录 FD/IO/RSS，但不等价于这些全量指标。以上缺项限制对普遍 I/O 瓶颈的结论，不隐藏为“全项通过”。

## 对 PLAN 的影响

H-D 在已测范围没有获得支持新架构的证据；现有 native/bypass 与 buffering 构成强碰撞，观察到的收益低于该路线门槛。以 KILL 停止低优先级新引擎投入，保留硬件及 profiling 缺项，不将负面投资判断等同于完整硬件实验验收。
