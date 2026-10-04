# 探索线索
## L-001 2026-09-28 HLA*LA 的线程参数分阶段生效
稳定版入口给 read alignment 传入 `threads=1`，但 inference 另用 `threads_for_HLAtyping=maxThreads`；不能笼统称整个 caller 单线程。
若真实 DP 阶段占比足够高，可审计其共享状态与并行安全性；最低代价是现有 8 线程基线后的同输入阶段计时与 1 线程对照。
状态：待挖掘；当前不改实现，等待实际阶段占比。
