# CCI 外部生信数据索引

2026-09-28；项目范围，不代表全主机资产。

已登记 23 个原始输入或导出文件，合计 86,220,503 bytes（包含派生文件，不是独立数据量）。

Kang：24,673 cells，8 donors，16 sample units；Visium：2695 spots、一张切片；1072条作者注释匹配，1623未注释，另1条作者barcode不在filtered counts。不能用于供体层推断。作者注释对象与 10x counts 按 barcode 对齐；不得将副本或空间块当成独立重复。

各文件来源、SHA256、字节数及用途见 BIOINF_DATA_INDEX.tsv。源归档中的其他示例不会自动纳入分析；只提取需要的数据库/坐标。
