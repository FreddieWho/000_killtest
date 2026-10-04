# A：碰撞审计（2026-09-28）

本表是官方发布、源码及论文的有限范围审计；不是本机性能测量。检索窗口为 2025-01-01 至今，较早方法作为现有基线补充。

| 方法 | 年份 | GTDB 实际可用性 | 是否改变 placement algorithm | RAM / runtime 证据 | 冲突判断 |
|---|---:|---|---|---|---|
| GTDB-Tk 2.7.2 / R232 | 2026 | 标准路径 | ANI 先筛；需要时仍调用 pplacer | 官方要求约 140 GB bacterial RAM，约 100 GB 参考存储；不是本机实测 | prescreen、预建 skani sketch、split tree 已存在 |
| pplacer alpha20–22 | 2025–2026 | GTDB-Tk 固定 alpha19，不能随意升级 | 发布说明主要为编译兼容、并行死锁和数值边界修复 | 未发布可据此认定数量级收益的测量 | “现代编译环境重写”本身不足以构成新贡献 |
| EPA-ng | 既有基线 | 可做 placement；GTDB 完整接入与结果等价未证明 | 搜索和优化实现不等同于 pplacer | 未取得同一 R232 / MAG 条件的可比测量 | 相同 jplace 格式不能证明 likelihood 或 taxonomy 等价 |
| SCAMPP / BSCAMPP | 2022–2025 | 可包装 placement；并非已证实的 GTDB 原位替换 | 使用子树，批量摊销搜索 | 2025 BSCAMPP 论文讨论大型树的扩展；不可直接外推 R232 | 子树与 batching 方向已有实质先例 |
| APPLES / APPLES-2 | 既有基线 | 未证明 GTDB workflow 等价 | 距离法，不是相同 likelihood inference | 未取得同输入的可比数值 | 属算法替代，不可冒充相同语义执行引擎 |
| atlas-place | 2026 | PyPI 自称替代；未独立核实完整可用性 | embedding、最近原型、conformal decision | 没有本轮独立数值 | 输出 jplace 不代表 likelihood-equivalent；不能据营销描述直接 KILL |
| GTDB-Tk scratch / pplacer mmap | 既有源码 | 已内置 | 改存储方式 | RSS 与页错误权衡待测 | 不能将“增加 mmap”本身当成新方案 |

源码已核实：alpha19 的 `pplacer_src/glvm.ml` 包含 `mask_into`、`masked_logdot` 与 `mmap_glv_arrays`；`pplacer_run.ml:306` 已分配三个 mmap GLV arrays，并在后续 fork。故 gap masking、mmap、进程间 reference 复用不能未经审计便假定完全不存在。是否有进一步布局空间仍须测量。

来源：
- [GTDB-Tk releases](https://github.com/Ecogenomics/GTDBTk/releases)，API 快照见 `sources/GTDBTk_releases.json`。
- [GTDB-Tk installation](https://ecogenomics.github.io/GTDBTk/installing/index.html)。官方文档中部分旧页面数字不一致，本轮固定 2.7.2/R232。
- [pplacer releases](https://github.com/matsen/pplacer/releases)，快照见 `sources/pplacer_releases.json`。
- [pplacer alpha19 source](https://github.com/matsen/pplacer/tree/v1.1.alpha19/pplacer_src)。
- [EPA-ng official implementation](https://github.com/Pbdas/epa-ng)。
- [BSCAMPP 2025](https://pubmed.ncbi.nlm.nih.gov/40811324/)。
- [APPLES official implementation](https://github.com/balabanmetin/apples)。
- [atlas-place publisher description](https://pypi.org/project/atlas-place/)。
- [GTDB maintainer response to 32 GB workstation limitation, 2026](https://forum.gtdb.ecogenomic.org/t/gtdb-error-55-gb/823)。该案例用旧 R214，不能拿其 55 GB 当 R232 测量。
