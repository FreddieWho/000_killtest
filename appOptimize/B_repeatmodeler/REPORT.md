# B — RepeatModeler2

2026-10-02 收尾结论：**KILL 当前搬运/调度重写方向**。完整 small 基线已经完成；不再投入 RMBlast 替代引擎。此结论针对当前输入与允许的工程方向，不是所有 repeat 算法都无改进空间。

## 完整基线与加速空间

软件 2.0.9、16 threads、seed 20260928，未启用可选 LTRStruct；Dfam 4.0 curated consensus。果蝇 GCF_000001215.4，143,726,002 bp，5 轮完成并生成 461 条最终 consensus families，51 条标为 Unknown。原始证据在 `profiles/B_dmel_full_20260928/`，最终库与 seed alignment 在 `benchmarks/B_dmel_full_20260928/dmel-families.fa` / `.stk`。

- 完整 wall 14,701 s（4:05:01），user+system CPU 186,853.32 s，平均约 12.71 cores。
- comparison 合计 12,590 s，占全流程 **85.64%**。在 comparison 不变的条件下，删除其余全部工作最多 **1.168×**。
- GNU time peak RSS 969,260 KiB 是最大进程/后代口径；30 秒采样的同时进程 RSS 总和最大约 4.69 GiB，含共享页重复计数，不能称准确总峰值。
- 5 分钟目录采样最大 664,244,750 bytes、10,436 files，包含输入和保留输出且可能漏掉短峰；未观察到支持“几十 GB 临时文件、可降 5 倍”的证据，不把采样最大值当绝对上界。

| round | wall s | comparison s | comparison 占本轮 |
|---|---:|---:|---:|
| 1 | 387 | NA | RepeatScout 路径 |
| 2 | 331 | 269 | 81.27% |
| 3 | 1549 | 1366 | 88.19% |
| 4 | 11587 | 10638 | 91.81% |
| 5 | 418 | 317 | 75.84% |

## 阶段分解与边界

| 不重叠日志分组 | 秒 |
|---|---:|
| sampling 标签（包含 TRF/TE masking） | 835 |
| 第一轮 extraction | 22 |
| comparison | 12,590 |
| RepeatScout | 219 |
| RECON cluster | 20 |
| instance gathering | 2 |
| refinement | 439 |
| 未归因残余（含最终分类、解析、I/O、调度等） | 574 |

合计 14,701 s。不能把 sampling 或残余整体称为纯 I/O。parse、I/O、scheduling 三者没有各自独立的 wall-time 计量；这里明确保留不可辨识，而不人为分摊。每个已完成 external exec 的区间在 `analysis/exec_intervals.json`，父子区间可能重叠，不可直接相加为流程 wall。

两个真实 RMBlast 批次的用户态 CPU 采样分别有 94.13%（果蝇）和 94.37%（水稻）集中于对齐扩展。果蝇批次文件系统调用计账总计 0.002959 s；未覆盖 mmap 页错误与全流程。支持 comparison 计算主导的解释，但这些批次不能证明所有 comparison 时间均不可优化，1.168× 仍是条件上限。

## medium 与科学输出

水稻 GCF_001433935.1，374,422,835 bp，按原任务许可完成前三轮后有界停止，第四轮 partial 排除。10 Mb 到 30 Mb 的 comparison 从 224 s 增至 1425 s，推断完整流程超过 3 h；这不是实测全程。第一轮受 trace 干扰，后续轮次已证明位于 detach 之后，第三轮 comparison 占 73.8%。

无架构原型，未改变 sampling 或 family discovery。因此 paired family identity、分类一致性及同一 RepeatMasker 的前后 remask **NOT_RUN / 无修改版可比较**，不能标为等价性通过。最终 small 库的条数、类别及 SHA256 见 `profiles/closure_20261002/B_evidence.json`；medium 没有完整最终库。完整数据搬运分解与 medium 全流程均未冒充已完成。

## 当前使用与碰撞


RepeatModeler2 没有被主流流程完全淘汰：[EBP 2026 年 5 月推荐工具](https://www.earthbiogenome.org/report-on-annotation-recommended-tools)仍列出它；[Earl Grey 当前实现](https://github.com/TobyBaril/EarlGrey)仍调用它。Earl Grey 已处理部分小基因组反复失败重跑问题，不能把同一修复当新贡献。

[2.0.9 release](https://github.com/Dfam-consortium/RepeatModeler/releases/tag/2.0.9)支持 RepeatMasker 4.2.4 / Dfam 4.0，修复旧 RECON 兼容；2.0.8 修复 consensus extension 正反链混接，不能用有该问题的 2.0.7 做科学等价基线。2025 起 Refiner 已引入 RepeatAfterMe，旧版本时间结构不能直接沿用。

EDTA 含结构化 TE 发现路线且部分使用 RepeatModeler；RepeatExplorer 面向 reads 聚类；RED 的 masking 结果不是 de novo consensus library。[BRAKER4](https://github.com/Gaius-Augustus/BRAKER4)提供 RED 路线只说明部分 annotation 场景可替代，不证明同科学输出的替代。有限检索未确认一个已独立验证、同语义的快速完整重写。


## 对 PLAN 的影响

完整 small 分母及跨基因组批次热点共同削弱 H-B。依原任务“不重新发明 aligner”的范围，现有证据不足以支持 ≥2–3× 的工程重写；停止该方向。无需为获得正结果追加更大 genome 或重写原型。
