# CCI kill-test：最小复现入口

控制要求：`011_cci_killtest.md`。`collision_matrix.md` 是第2节名称的兼容链接，指向第12节规范名称 `01_collision_matrix.md`，不是另一份报告。结论与范围从 `00_DECISION.md` 开始；六份科学报告、`results/benchmark.tsv` 与本目录脚本为交付主体。项目管理文件按 AGENTS.md 保留。

## 当前冻结结果

`results/selected_runs.json` 指定被报告采用的运行。其他早期目录保留失败/修订线索，不混入最终数字。

- `kang_run_03`：全部原型统计、IC、100个complex ILP和256次条件交换。
- `scaling_run_01`：2k/5k/10k/full、共同450特征的规模曲线。
- `pooled_full_run_01`：六主要细胞类型、4951特征的100×199实际条件置换。
- `ic_sanity_01`：物理gene去重、图删除合成证据、IFN方向检查。
- `spatial_run_04`：全部2695spots、官方100µm中心距标定、双核、typed结果、全LR dense对照。
- `native_*_01`：四个scRNA原生包及Spatial CellChat的输出、参数和receipt。

数据来源、bytes、SHA256和列数见 `infra/bioinf-data-index/BIOINF_DATA_INDEX.tsv`。raw数据约几十MB量级；未下载FASTQ/BAM/百万细胞atlas。10x提供的spatial归档包含图像，只提取坐标和scale factors用于计算。作者注释对象和当前filtered counts有1条barcode差异，见空间报告。

## 环境

当前主机 Python：`/opt/anaconda3/bin/python`；R：`/usr/bin/R`。准确包版本与平台见 `results/environment.json`、`results/R_environment.txt`。R端需要已安装的CellChat、Matrix、future、jsonlite；Python基础环境需要NumPy/SciPy/pandas/anndata/statsmodels/sklearn/networkx/patsy/h5py/matplotlib。

本轮只将额外Python包放入项目`.deps`；清单为 `sources/extra_requirements.txt`。当前主机需要指定conda编译库和可写缓存路径：

```bash
export LD_LIBRARY_PATH=/opt/anaconda3/lib
export PYTHONPATH="$PWD/.deps"
export NUMBA_CACHE_DIR="$PWD/.cache/numba"
export MPLCONFIGDIR="$PWD/.cache/mpl"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
```

## 复算

输入已在磁盘。完整数值重算使用**新的**结果目录，旧目录存在即停止；原生CPDB/CellChat会花费分钟级时间，不需高频轮询：

```bash
bash scripts/reproduce.sh results/reproduce_new
```

脚本重建派生输入并产生新数值结果，不自动替换当前六份报告。数据缺失时可先运行 `python scripts/fetch_inputs.py --download`；无该参数时仅复用和检查已有原始数据，下载后仍需更新数据索引。源hash改变会停止。

单独复算核心路径（仍应使用新输出目录）：

```bash
python scripts/run_kang.py results/new_kang
python scripts/run_ic_sanity.py results/new_kang results/new_ic_sanity
python scripts/run_spatial.py results/new_spatial
python scripts/finish_spatial.py results/new_spatial
python scripts/run_scaling.py results/new_full --full-only
```

当前交付报告从已选择、已完成的receipts生成：

```bash
python scripts/make_reports.py
```

该命令拒绝未完成运行，写入六份报告中的五份（碰撞矩阵保留）、benchmark和三张诊断图。`scripts/check_algebra.py` 是本任务的全输出一致性、原始计数守恒和power真值检查；没有另建通用测试框架。

## 解释边界

CPDB/CellChat原生specificity p不是condition-difference p。原生包的样本分数接同一配对OLS进行比较；真实供体重复数为8。Visium只有1张切片。IC100与hitting set是固定支持图上的理想删除数量，不是药物效应预测；真实blockade验证未运行。可选IC90、skin、LMM和50k/100k复制规模未做，没有用它们支撑结论。

## 2026-09-30 内存优化

推荐在新输出目录运行 `python scripts/run_kang.py results/new_kang_stream --stream`，环境设置同上。全基因library分母不变，输出统计仅覆盖资源所需基因；原默认全基因路径保留。相同两库分数与配对模型实测峰值996→429MiB，11.60→8.81秒（单次非独占测量）。完整结果与边界见[优化报告](06_engineering_optimization.md)。历史manifest不表示当前修改后的源码，新manifest位于 `results/optimization_20260930/artifact_manifest.tsv`。
