# 当前待办（仓库级）

本文件是仓库层面的共读入口，只放跨项目或仓库自身的事项。各子项目的技术待办仍在 `cci/TODO.md` 与 `appOptimize/TODO.md`，不搬到此处。

- [x] 建立仓库并推送到 GitHub（完成 2026-10-04）
    - [x] 写根 `.gitignore`，排除 341G 本地数据与环境目录
    - [x] 写根 `README.md` 说明两个子项目与数据存放约定
    - [x] 放行 `appOptimize/benchmarks/*.py` 与 `profiles/README.md`
    - [x] 改用 SSH 推送（`git@github.com:FreddieWho/000_killtest.git`）

- [ ] 确认两个子项目在 GitHub 上的可读性
    - [ ] 抽查 `README.md` 中的相对链接与表格在网页端是否正常
    - [ ] 确认 `cci/collision_matrix.md` 符号链接在网页端的行为（必要时改为普通文件）

## 分支记录

- 2026-10-04：入库范围选择「仅文档+脚本+索引」；`data/`、`benchmarks/` 运行产物、`profiles/`、`recovery/`、`appOptimize/sources/`、`cci/results/` 全部留在本地磁盘。子项目的停止/收尾分支维持原状，本仓库不重启任何科学分支。

## 变更记录

- 2026-10-04：建立根 `README.md` 与本文件；新增符号链接检查待办（起因：符号链接在 GitHub 网页端不会跟随）。
- 2026-10-04：推送改走 SSH。本机原全局配置 `url.https://github.com/.insteadof=ssh://git@github.com/` 会把 SSH 地址改写成 HTTPS，导致认证失败，已用 `git config --global --unset-all url.https://github.com/.insteadof` 删除；恢复用 `git config --global url."https://github.com/".insteadOf "ssh://git@github.com/"`。凭据存储中只有 hf-mirror，未新增任何 token。
- 2026-10-04：`benchmarks/*.py` 改为入库，原因是 `appOptimize/STATUS.md` 的「复现汇总」直接引用 `benchmarks/finalize_closure.py`，排除后报告不可复算；其下运行产物仍排除。