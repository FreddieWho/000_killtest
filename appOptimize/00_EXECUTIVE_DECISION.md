# 四路线 Kill-Test 决策

本机补测与报告收尾完成；原任务全范围仍有明确缺项

**当前 0 条 GO；A/C HOLD，B/D KILL。** HOLD 不是建议进入大重写阶段，KILL 是停止当前范围的新项目投入，不是证明任何算法或硬件下都无优化空间。

| Route | Current usage | Collision risk | Critical-path share | Theoretical ceiling | Prototype result | Scientific equivalence risk | Decision |
|---|---|---|---|---|---|---|---|
| A GTDB-Tk / pplacer | 2.7.2/R232 标准工作流 | 已有 mmap、gap masking、split-tree/reference 共享 | W2 pplacer 90.10%，54 个实际 placement | 假设全部 pplacer 3×，端到端约 2.50×；未证实 | 无新架构原型；线程/热点补测见 A 报告 | 单调用 edge/LWR 与全流程 taxonomy 比较范围分开 | HOLD |
| B RepeatModeler2 | EBP/Earl Grey 的现有使用证据 | 近期 refinement 已优化；不允许重写 aligner | 完整 small comparison 85.64% | comparison 不变时删除其余工作仅 1.168× | 未做；两批函数热点支持对齐主导 | 无修改版，无 paired remask 等价性结论 | KILL |
| C HLA-LA / HLA-HD | 现有 workflow/benchmark 使用证据 | 官方 fast 分支与既有 prefilter/索引复用 | WGS 首次 index 6713.81 s；WES backend 其余 73.6%，非纯 DP | 一次性建索引不能作新收益；单阶段可加速空间未证实 | 无新原型；1/8 线程与调用比较见 C 报告 | 低覆盖 WGS 真值表现弱；ambiguity 与唯一 allele 不等价 | HOLD |
| D GenomicsDB | GATK/native 标准路径 | 高：native buffering 与 bypass reader 已有 | 已完成 50/200 完整 chr22 三组导入 | 已有选项最多观察到 2.107× 时间、4.139× RSS 比；非理论普遍上界 | 无新原型，DP 元数据诊断仅解释输出 | 两规模五区间差异由 DP sum 配置解释；非全染色体等价 | KILL |

## 已闭合的证据与仍未完成的范围

B 完整果蝇基线 4:05:01、461 families，水稻按原任务许可仅前三轮。D 50/200 样本完整导入与五区间差异诊断完成；200 样本 53 条 INFO/DP 差异在仅改合并元数据后全部消失。

A 旧零计数已纠正：GTDB-Tk 清理中间文件导致后处理误判，不是 100 个 MAG 都被 ANI 分流。C 旧回放缺失 hla_nom_g.txt 已定位为工作目录问题。A/C 本次状态：`COMPLETED_SCHEDULED_COMPUTATION` / `COMPLETED_SCHEDULED_COMPUTATION`；实际测量与一致性细节见两份报告及 profiles/closure_20261002/closure_assessment.json。

HLA-HD 软件/RNA 分型仍 BLOCKED_SOFTWARE；当前主机无 NVMe，对照 NOT_RUN。唯一四字段 HLA 真值不可用。D interval/batch 扫描、全量 syscall/decode/sequential-read 分解及全染色体等价性未完成；B parse/IO/scheduling 精确分离未完成。上述缺项均未用 proxy 或文献数字替代，不宣称原任务所有必做项已验收。

A/C 没有得到新架构达到门槛的原型实证，暂不升级 GO。B/D 在已有测量、碰撞与任务允许的修改范围内停止投入；不为获得正面结果补写替代引擎。

## 对 PLAN 的影响

H-A 的关键路径前提获得支持，但可压缩空间未证实；H-B 被完整基线削弱；H-C 尚受科学准确性/输入范围和原型证据限制；H-D 未支持独立新架构。PLAN 保持原文。最终最多一条 GO 的要求满足，本轮实际为零。
