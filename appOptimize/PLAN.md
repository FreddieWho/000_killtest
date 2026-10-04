# 四路线 Kill-Test
依据：011_killtest_4routes.md。目标是判断保持科学语义的执行层改进是否具备论文级空间。

- H-A：现代 GTDB-Tk 的 placement 仍存在足够大的时间或内存可压缩空间。
- H-B：RepeatModeler2 的数据搬运、解析、采样与调度占足够大的成本。
- H-C：活跃 HLA backend 的关键阶段可在保持推断语义下大幅提速。
- H-D：GenomicsDB import 存在尚未被 native importer 覆盖的结构性重复。

判据沿用任务文件，不新增门槛。最多一条 GO；无实测支持时不得把预期收益写成实测或完成。
