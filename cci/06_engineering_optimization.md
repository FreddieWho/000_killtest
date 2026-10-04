# 工程优化结果 · 2026-09-30

结论：相同两库分数与配对模型的输入→聚合→推断流程，内存显著降低，单次实测耗时也下降。统一新方法的KILL结论不变。

| 路径 | 计算耗时（秒） | 计算阶段进程峰值RSS（MiB） | 保留基因 |
|---|---:|---:|---:|
| 原版，全矩阵 | 11.598 | 995.508 | 15706 |
| 分块1024行 | 13.434 | 390.574 | 575 |
| 分块4096行（推荐） | 8.810 | 429.383 | 575 |

推荐路径较原版内存减少56.9%，计算耗时减少24.0%。三者都从同一Kang原始H5AD开始，计算CellChat与LIANA/CellPhoneDB两套资源的全部分数、支持集合及相同配对模型；未缩减供体、细胞类型、LR或被检验特征。计时包括读取和聚合，排除Python导入、结果验证和压缩落盘；RSS在计算结束采样，是进程累计峰值，包含导入。宿主机非独占、仅单次测量，部分运行与独立的完整Kang检查重叠，速度百分比不是稳定性能保证。这一比较不包含原生工具的specificity输出，也不是对四个原生工具的全面胜出。

## 改动与数值证据

- 新增 `scripts/streaming_input.py`：backed H5AD分块读取，先从所有15706基因计算每个细胞library总量，再保留所需基因。默认4096行；较小块更省内存，但读取更慢。
- `cci_core.aggregate`允许传入全库library总量；原API默认行为不变。复合体仍在细胞层取子单元最小值，然后按样本汇总。
- `scripts/run_kang.py --stream`接入实际完整流程。LR资源从项目本地来源文件读取，移除对个人目录下LIANA安装文件的依赖，并停止每次运行重写来源文件。
- 分块只限制输入阶段临时矩阵，仍保留全部细胞×所需基因稀疏矩阵与结果；并非对任意细胞数恒定内存。
- 数值检查覆盖所有分数、支持集合、实际被检验特征、beta/SE/p值、保留基因均值/方差/检出率/原始计数、全库library量及细胞数（容差1e-12）。完整流程的IC、配对检验、256条件随机化、因子模型及合成功效与旧运行一致。不同等价最小覆盖可选不同靶点，核对的是相同最小数目和程序，而非强制同一最优解代表。

证据：`results/optimization_20260930/{full,stream,stream4096}/receipt.json`、`checks.json`、`kang/receipt.json`。完整Kang验证使用1024块，推荐4096块已独立核对全部分数及模型输出；分块大小不改变推断。历史六份报告与原结果保持原状；旧artifact manifest是历史快照，不代表修改后的当前源码，当前变更见新目录manifest。

## 复现

先按README设置本机动态库与Python路径。必须使用新输出目录：

```bash
python scripts/run_kang.py results/new_kang_stream --stream
python scripts/benchmark_memory.py full results/new_full
python scripts/benchmark_memory.py stream results/new_stream --chunk-size 4096
```

`benchmark_memory.py`依赖原冻结结果作对照，`check_optimization.py`复核本轮固定交付目录，二者不是可替代科学验证的通用测试框架。默认 `run_kang.py` 保留全基因统计；`--stream` 的sufficient_statistics仅含所需基因，基因名随文件保存，不能将其用于全转录组分析。

## 剩余阻塞及对PLAN的影响

H4获得更强的同口径工程证据。H2的一般设计创新不足、H3的count→IC留一供体R²约0.918不会因省内存而改变。真实阻断验证仍NOT_VALIDATED，空间仍只有一切片，MOSANIC全文缺口仍在；本轮未新增外部数据或重启这些分支。跨越这些科学阻塞需要匹配的独立扰动/多样本数据或有实证支持的新科学问题，不应修改判据包装为通过。
