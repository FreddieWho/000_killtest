# Spatial homology

## 数据与单位

10x Mouse Brain Serial Section 1, Sagittal Anterior标准数据，官方filtered H5与spatial包；CellChat作者Figshare条目23621151/file41446839提供既有注释。保留全部**2695 spots**；原始32285 gene IDs按重复symbol求和后32245 symbols。作者对象有1073 spots，其中1072条barcode匹配；1条不在此次filtered counts。剩余1623标作UNANNOTATED，未重新聚类或做deconvolution。

只有**1个生物样本**，不做spatial replicate-level显著性结论，也不把复制坐标块作为重复。作者对象仅用于注释，表达使用官方全filtered矩阵，不用其小型示例表达子集替换全数据。

官方名义v1 spot直径55µm、中心距100µm：[10x术语表](https://www.10xgenomics.com/support/cn/software/space-ranger/3.1/getting-started/space-ranger-glossary)。此数据JSON直径配55µm却得到近邻84.216µm；CellChat教程用65µm约99.528µm。因此用实际近邻像素间距137.000与官方100µm pitch标定：0.729927µm/pixel。保留冲突及来源，不盲目套直径换算。

## 同一接口和代数

两模态调用 `cci_core.communication`，相同resource loader、sample/sender/receiver/lr_index schema、相同active规则、同一IC和下游样本线性模型。scRNA仅population mixing：Ω=1/(nA nB)，C=mean(L)mean(R)，不补造坐标。空间Ω=K(d)/(nA nB)，归一化分母固定为所有可能A/B配对数，避免稀疏后改分母掩盖机会减少。

Kang为human、Visium为mouse，使用同一CellChatDB家族的物种资源；不冒称两物种是相同基因图。真正同资源的对照是在同一Visium矩阵上分别应用mixing/contact/secreted Ω。其all-to-all dense机会算子恢复mixing，包括complex组件。

- contact：110µm radius、二值、无self loops，表示spot邻接proxy，不能解释成单细胞物理接触。
- secreted：250µm截断exponential，scale100µm、无self loops。
- 正式typed结果按资源annotation：Cell-Cell Contact用contact，Secreted Signaling用secreted。分别把两种核应用全部LR的结果仅是同资源计算对照。

## 结果与数值证据

| kernel                   |   max_absolute_difference | all_LR_and_all_group_pairs   | edges   | dense_pairs   |
|:-------------------------|--------------------------:|:-----------------------------|:--------|:--------------|
| contact                  |                 6.939e-18 | True                         | 15656.0 | 7263025.0     |
| exponential              |                 6.939e-18 | True                         | 46202.0 | 7263025.0     |
| all_to_all_equals_mixing |                 6.217e-15 | True                         | —       | —             |

全部可映射LR及全部group pairs均与dense comparator对照，而非只抽几条边。Kernel替换后支持通道从25,350降至16,162，移除9,188、新增0；有支持的group pairs从81降至61。这支持空间约束确实移除了部分不邻近/不共现的机会；不等于这些被保留通道已在实验中通信。

图的距离矩阵对照、建图和aggregation的分开runtime/存储见 `results/spatial_run_04/profile.tsv`。1/2/4个分离坐标块测到2695/5390/10780 spots，仅检验complexity。正式路径不构造dense N²距离；full dense只用于基准和数值对照。

## Spatial CellChat 语义对照

作者教程采用空间距离与接触范围约束并计算群体概率。本轮实际运行安装版本2.2.0，全部2695spots和既有/未知标签，保留全部已测LR与cofactor，20次bootstrap仅作profiling，耗时85.26s；native net保存在 `results/native_SpatialCellChat_01/net.rds`。这是CellChat 2.2.0的spatial mode，未实跑v3；v3纳入了碰撞审计。本原型的逐边核和传统product与其群体距离/Hill/cofactor模型不同，不要求数值相同，不称为CellChat复现准确度。

KT6：接口、资源表示、IC和数学退化关系通过；空间稀疏同口径优势见工程报告。**Spatial replicate calibration NOT_TESTABLE（只有一切片）**。对PLAN的影响：H4工程同构成立且无需虚构scRNA空间；不支持跨模态生物效应等价或泛化声明。
