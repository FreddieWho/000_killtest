# 当前待办（仓库级）

本文件是仓库层面的共读入口，只放跨项目或仓库自身的事项。各子项目的技术待办仍在 `cci/TODO.md` 与 `appOptimize/TODO.md`，不搬到此处。

- [x] 建立仓库并首次推送到 GitHub
    - [x] 写根 `.gitignore`，排除 341G 本地数据与环境目录
    - [x] 写根 `README.md` 说明两个子项目与数据存放约定
    - [x] `git init` + 首次提交 + 推送 `FreddieWho/000_killtest`

- [ ] 确认两个子项目在 GitHub 上的可读性
    - [ ] 抽查 `README.md` 中的相对链接与表格在网页端是否正常
    - [ ] 确认 `cci/collision_matrix.md` 符号链接在网页端的行为（必要时改为普通文件）

## 分支记录

- 2026-10-04：入库范围选择「仅文档+脚本+索引」；`data/`、`benchmarks/`、`profiles/`、`recovery/`、`appOptimize/sources/`、`cci/results/` 全部留在本地磁盘。子项目的停止/收尾分支维持原状，本仓库不重启任何科学分支。

## 变更记录

- 2026-10-04：建立根 `README.md` 与本文件；新增符号链接检查待办（起因：符号链接在 GitHub 网页端不会跟随）。