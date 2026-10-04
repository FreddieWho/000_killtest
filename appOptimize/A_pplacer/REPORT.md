# A — GTDB-Tk / pplacer

2026-10-02：**HOLD，补测完成，仍有范围缺项，无 GO**。原先“0 placement”结论已纠正：GTDB-Tk 默认清理中间文件，旧后处理把文件缺失误计为零；最终表中 46 个 ani_screen、36 个 topology、18 个 RED，实际 54 个经过 bacterial placement。原始记录保留，纠正说明见 `profiles/A_scaling/W2_placement_denominator_correction_20261002.json`。

## 当前完整基线

W1 是 E. coli K-12 GCF_000005845.2；W2 为 100 个真实近期 environmental/uncultured bacterial MAG，输入及来源 receipt 不变。R232 已完整下载、校验并解包。GTDB-Tk 2.7.2、pplacer alpha19、skani 0.3.1，标准 split-tree classify_wf，cpus16 / pplacer8。

W1 wall 82.21 s；W2 wall 9089 s，GNU time peak RSS 138,821,336 KiB（约 132.4 GiB）。该峰值口径不是同时存活进程的独占内存之和。W2 trace 中有 9 次真实 pplacer 调用，区间并集 8189.64 s，占完整分类 **90.10%**。这不是把版本探针算作 placement；也不把 100 输入当作 100 placement。若所有 pplacer 工作真能加速 3 倍，条件端到端约 2.50 倍；尚无原型证实这一局部速度。

新目录 `benchmarks/A_W2_retained_20261002` 使用同一输入与参数，仅加 `--keep_intermediates` 保留证据。补测内容为核对最终 taxonomy、验证 jplace 分母，并对实测最慢的一次 placement 做 1/8 线程回放与函数采样；此单调用覆盖范围不会冒充全体 54 个 query 的线程对照。

碰撞审计见 [A_collision.md](A_collision.md)。现有 mmap、gap masking、split-tree 与 fork 前 reference 共享不能重新包装成创新。运行资源与完整 trace 见 `profiles/A/W1_8`、`profiles/A/W2_8`，汇总见 `profiles/closure_20261002/A_baseline_evidence.json`。

## 解释边界


alpha19 源码已确认：`pplacer_run.ml` 在 fork 前创建、计算 `darr/parr/snodes` 三套 reference likelihood arrays；`--mmap-file` 分支使用 shared mapping。不能先假设每个 worker 都完整复制这三套数组。`gstar_support.ml` 的 float64 likelihood 与 native-int exponent 布局，在本机 64 位下，CAT 模型三数组的名义大小为 `3 × n_glvs × n_sites × 8 × (n_states+1)` bytes；mixture 模型用 `n_states × n_rates` 替代 `n_states`。这只是源码分配量，未含对象、alignment、query scratch、页共享或驻留率。具体名义分配需结合实际 reference 与 masked alignment；未测部分不填推测 RAM 值。内置 `--pretend` 可输出该名义分配估计，后续与 RSS/PSS 区分。

Amdahl 上限按 `1 / ((1-f)+f/s)` 计算。f=40%、局部 3 倍仅得 1.36 倍端到端；不能因超过阶段占比候选值便给 GO。现有 mmap、masking、fork 复用也不能重新包装成创新。

30 行 summary 非逐字一致均只涉及 other_related_references 列的参考条目顺序：按分号拆分、去除分隔空格后，完整条目集合相同，未更改数值或参考条目。证据：`profiles/A_scaling_20261002/summary_field_diagnosis.json`。

## 对 PLAN 的影响

H-A 的“placement 仍是关键路径”获得支持；但 H-A 的可压缩资源/执行空间尚未证实。无改写原型，无新引擎速度或内存缩减声明。线程回放与 taxonomy 复现是证据补充，不自动构成 GO。

<!-- closure-results -->
## 本次补测结果

控制器状态：`COMPLETED_SCHEDULED_COMPUTATION`。

保留的 jplace 分母：54 个不同 bacterial query，9 次调用。

线程回放选择最慢实际调用：基线 1141.46 s，覆盖 6 个 query；不能外推为所有 query 的线程等价性。

100 样本重跑的 taxonomy 完全一致：True；全部 summary 字段完全一致：False。这是同一工作流保留中间文件的复现。

回放 query 集合相同：True；基线与 1/8 线程的全部 jplace 数值精确相同：True / True。每条 query 的 edge、LWR、数值差异及 top-two gap 已写入 jplace_comparison.json；保留多候选 edge 的 query 数 0。

1 线程实际后端回放：wall 1240.91 s，peak RSS 138821336 KiB。

8 线程实际后端回放：wall 1124.55 s，peak RSS 138821336 KiB。

1/8 线程 wall 比 1.103×；同一现有二进制线程对照，不能称新原型加速。

函数采样（49 Hz 用户态 CPU，非端到端 wall 占比）：

```text
24.99%  pplacer  pplacer        [.] cblas_dgemv
12.03%  pplacer  pplacer        [.] mark_slice
9.82%  pplacer  pplacer        [.] caml_page_table_lookup
5.73%  pplacer  pplacer        [.] sweep_slice
4.81%  pplacer  pplacer        [.] __ieee754_log_avx
4.17%  pplacer  pplacer        [.] mat_masked_logdot_c
3.29%  pplacer  pplacer        [.] gsl_vector_max
2.72%  pplacer  pplacer        [.] caml_adjust_gc_speed
2.49%  pplacer  pplacer        [.] caml_ba_get_N
2.49%  pplacer  pplacer        [.] caml_ba_slice
2.03%  pplacer  pplacer        [.] caml_ba_offset
1.69%  pplacer  pplacer        [.] dediagonalize
1.65%  pplacer  pplacer        [.] caml_ba_alloc
1.46%  pplacer  pplacer        [.] caml_ba_finalize
1.23%  pplacer  pplacer        [.] caml_c_call
1.23%  pplacer  pplacer        [.] camlGstar_support__masked_total_twoexp_10340
1.16%  pplacer  pplacer        [.] caml_fl_merge_block
1.12%  pplacer  pplacer        [.] caml_alloc_shr
```

原型仍为 NOT_RUN；本轮只验证现有程序和证据完整性。新架构收益未得到验证，维持 HOLD，不给 GO。

