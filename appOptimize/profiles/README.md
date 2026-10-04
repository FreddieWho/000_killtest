# 原始证据与测量口径

这里保存运行命令、软件版本、资源采样、子进程 trace、退出回执及准备过程日志。四路线尚在执行，文件存在不等于完整基线已完成。

- `environment_versions.json`：实际安装的软件版本。
- A/C/D 的 `*_controller.log`：自动衔接输入完成与计算；等待依据进程退出或文件完成事件，不使用分钟级状态轮询。
- `B_dmel.log` / `B_dmel.time.txt`：small 基线；完整 GNU time 在进程结束后才写完。
- `B_rice/`、`A/`、`C/`：`invocation.json`、`samples.jsonl`、`process.trace`、`time.txt`、`status.json`。资源每 30 秒采样，临时文件 metadata 每 300 秒采样。进程 trace 用于外部子命令边界，不是每个内部函数的 CPU profile。
- `D_20260928/`：计时工具在容器内执行；launcher wall 包含启动开销，单独保存，不能拿 Docker 客户端 RSS 当 importer RSS。

`benchmarks/summary.tsv` 由 `python3 benchmarks/summarize_profiles.py` 从回执生成。`EXIT_0` 只表示进程正常退出，科学等价性另行判断；运行中值为 NA，不能填写预期结果。

GNU time 的 filesystem input/output 按 Linux 512-byte block 单位转换，是内核记账的存储 I/O，不等同于逻辑 read()/write() bytes。page cache 命中可使物理读取很低。单进程/后代最大 RSS 不是进程树同时驻留内存总和；共享页不能简单重复相加。

本机存储为 sda 旋转盘，未测 NVMe/NFS。普通用户 perf 权限受 `perf_event_paranoid=4` 限制；内部函数热点和 cache misses 不能假装已采集。

`D_download.log` 保留初次 CA 路径失败记录；`D_download_retry.log` 是修复后的正式准备。`A_controller_startup_failed.log` 保留 Python pidfd API 兼容问题；控制器已使用 Linux syscall 事件等待修复。准备失败不计入基线速度。

原型仅在实际 profile 支持时启动；目前没有原型，也没有等价性通过结论。

2026-09-28 测量修正：水稻初始 strace 耗用显著 CPU，已仅停止观察进程、保留 RepeatModeler 计算，并用 `B_rice_after_detach/` 继续低开销采样。`B_rice/tracer_detach.json` 保存边界；该次全程 wall time 不作干净基线。后续 process trace 启用 `--seccomp-bpf`。

perf 能力检查：PERFMON-only 容器仍被内核拒绝；无网络、只读、无 host PID 挂载、仅 SYS_ADMIN 的容器可以采样其自身子进程，见 `perf_capability_confined_admin.txt`。主机 sysctl 未改变。这只是工具能力检查，不是性能实验证据。

后续 `perf record` 初始发生工具自身崩溃，首次故障定位于 BPF 元数据查询，见 `perf_record_first_fault.txt`。仅关闭 BPF 事件或去掉调用栈均未解决。容器 seccomp 禁止 BPF 系统调用后录制成功，见 `perf_record_capability_bpf_denied.txt`；科学任务的采样仍须另看对应回执。`benchmarks/perf_replay.py` 使用新输出目录、49 Hz 用户态平面采样；采样回放不是干净计时基线。

`benchmarks/B_round_metrics.tsv` 从日志提取已完成轮次与原始阶段标签；未完成轮次不估计 ceiling。`benchmarks/analyze_process_trace.py` 提取 exec 到 exit 的跨度，包含子进程等待，禁止把重叠跨度直接相加为流程 wall time。

后续自动衔接：`A_post_controller.log` / `C_post_controller.log` 等待完整基线后执行真实后端线程回放、函数采样与输出比较；当前日志为空表示尚在等待，不能据文件存在宣称已完成。`D_export_comparison_controller.log` 对每个成功导入执行同读取器五区间比较。控制器一次性快照在 `controller_manifest.json`，不是实时完成状态。

medium 按任务许可有界结束：`B_rice/bounded_stop.json` 是范围回执，前三轮保留、第四轮 partial 排除；底层程序即使返回 0，也不代表完整建库。`B_rice/post_detach_round_boundary.json` 由启动时刻与 round 1 时长下界证明 rounds 2–3 在观察器退出后，允许使用其阶段日志，不恢复初始全程计时的有效性。

新增数据在项目与中央索引同步；`index_completion_events.log` 在现有下载/分类控制器退出事件后登记新增完成资产，不把稀疏下载文件算完整数据。

恢复版本由 `appopt-20260928-a.service`、`appopt-20260928-b.service`、`appopt-20260928-c.service`、`appopt-20260928-d.service` 承载；命令和最终状态在 `persistent_20260928/`，不依赖 API 工具会话存活。`B_dmel_full_20260928/` 是新连续基线；原 `B_dmel.log` 是中断历史。D 的首次 200 初始化因用户管理器缺失 docker 附加组而未启动容器，见 `N200_native_init_failed_docker_group/`；通过已有组成员资格恢复后，正式初始化已成功。

`D_schema_diagnostic_20260928/` 是只改变复制元数据的 DP 合并规则诊断，源数组与基线不变；五区间输出与 GATK 匹配。该结果不等于全染色体等价证明。多个路线共享同一旋转盘与主机，缓存、并发负载未隔离；时间用于探索定位，不直接外推独占硬件的精确速度比。
