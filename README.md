# 000_killtest

收录各类 kill-test 项目的集合。每个子目录是一个独立项目，各自带 `PLAN.md` / `ROADMAP.md` / `DECISIONS.md` / `TODO.md` / `STATUS.md` / `LEADS.md`，当前进展从各自的 `STATUS.md` 读起。

| 子目录 | 内容 |
|---|---|
| [`cci/`](cci/) | CCI（细胞通讯推断）kill-test：碰撞矩阵、统计检验、IC 删除、工程画像、空间同源性，六份科学报告 + `CCI_killtest_final_report.md`。入口见 [cci/README.md](cci/README.md)。 |
| [`appOptimize/`](appOptimize/) | 四条生信应用路线（pplacer 放置、RepeatModeler、HLA-LA、GenomicsDB）的 kill-test 与工程对照，路线报告在 `A_pplacer/`～`D_genomicsdb/`，总览见 `011_killtest_4routes.md`。 |

## 数据与产物的存放约定

仓库只入库文档、脚本与数据索引。体积大的原始数据、运行产物与环境目录留在本机磁盘，不进 git：

- 排除：`data/`、`benchmarks/`、`profiles/`、`recovery/`、`appOptimize/sources/`、`cci/results/`、`.cache/`、`.envs/`、`.deps/`、`.pi/`
- 入库：`*.md`、`scripts/`、`infra/bioinf-data-index/`（含 `BIOINF_DATA_INDEX.tsv` 与 `BIOINF_DATA_INVENTORY.md`）

各项目的 `scripts/reproduce.sh` 与 `results/selected_runs.json` 依赖本机已有数据；换机器复算前先读对应 `README.md`。