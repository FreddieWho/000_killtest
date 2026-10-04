"""Evidence-driven report refresh; never promote missing/failed measurements to completion."""
import csv,datetime,hashlib,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];E=ROOT/'profiles/closure_20261002'
def read(path):
 p=ROOT/path
 if not p.exists():return None
 try:return json.loads(p.read_text())
 except json.JSONDecodeError:return {'state':'INVALID_RECEIPT'}
def timing(path):
 p=ROOT/path/'time.txt'
 if not p.exists():return {}
 d={}
 for line in p.read_text().splitlines():
  if ': ' in line:
   k,v=line.strip().split(': ',1);d[k]=v
 wall=next((v for k,v in d.items() if k.startswith('Elapsed (wall clock)')),None)
 def num(k):
  try:return float(d[k])
  except (KeyError,ValueError):return None
 return {'wall_s':sum(float(v)*60**i for i,v in enumerate(reversed(wall.split(':')))) if wall else None,'user_s':num('User time (seconds)'),'system_s':num('System time (seconds)'),'peak_rss_kib':num('Maximum resident set size (kbytes)')}
def supplemental(path,text):
 p=ROOT/path;s=p.read_text();s=s.split('\n<!-- closure-results -->')[0]
 p.write_text(s+'\n<!-- closure-results -->\n'+text+'\n')
def stage_status(r):
 recovery=read('profiles/closure_20261002/C_recovery/status.json') if r=='C' else None
 return (recovery or read(f'profiles/closure_20261002/{r}/status.json') or {}).get('state','NOT_RUN')
a_state=stage_status('A');c_state=stage_status('C');done=all(s=='COMPLETED_SCHEDULED_COMPUTATION' for s in [a_state,c_state]);failed=any(s=='FAILED' for s in [a_state,c_state])
aq=read('profiles/A_scaling_20261002/W2_placement_denominator.json');ac=read('profiles/A_scaling_20261002/jplace_comparison.json');at=read('profiles/A_scaling_20261002/taxonomy_comparison.json');asel=read('profiles/A_scaling_20261002/case_selection.json');cc=read('profiles/C_scaling_20261002/call_comparison.json')
required_present=all(x is not None for x in [aq,ac,at,asel,cc])
if done and not required_present:raise RuntimeError('Completed controller but required comparison artifacts absent')
checks={'A_qualified':bool(aq and aq.get('task_requires_at_least_50')),'A_same_query_sets':bool(ac and ac.get('query_sets_equal')),'A_same_1thread_jplace':bool(ac and ac.get('baseline_vs_1_exact')),'A_same_8thread_jplace':bool(ac and ac.get('baseline_vs_8_exact')),'A_same_taxonomy':bool(at and at.get('taxonomy_exact')),'C_same_1thread_alleles':bool(cc and cc.get('baseline_kernel1_exact_allele_sets')),'C_same_8thread_alleles':bool(cc and cc.get('baseline_kernel8_exact_allele_sets'))}
for key in list(checks):
 source=aq if key=='A_qualified' else at if key=='A_same_taxonomy' else ac if key.startswith('A_') else cc
 if source is None:checks[key]=None  # Not run is not a failed scientific comparison.
# Presence/completion and scientific equivalence are separate; a negative comparison is retained.
measurements={}
for r,base in [('A','profiles/A_scaling_20261002'),('C','profiles/C_scaling_20261002')]:
 measurements[r]={str(t):timing(f'{base}/threads_{t}') for t in [1,8]}
 perf_dir='perf_8_retry_20261002' if r=='C' and (ROOT/'profiles/closure_20261002/C_recovery/status.json').exists() else 'perf_8'
 measurements[r]['perf_path']=f'{base}/{perf_dir}'
 measurements[r]['perf_status']=read(f'{base}/{perf_dir}/status.json')
 v=measurements[r]
 if v['1'].get('wall_s') and v['8'].get('wall_s'):v['one_over_eight_wall_ratio']=v['1']['wall_s']/v['8']['wall_s']
 # Keep function samples separate from wall-time fractions.
 p=ROOT/base/perf_dir/'perf_report.txt'
 if p.exists():v['flat_profile_lines']=[l.strip() for l in p.read_text(errors='replace').splitlines() if '%' in l and not l.lstrip().startswith('#')][:18]
a_text='## 本次补测结果\n\n控制器状态：`'+a_state+'`。\n\n'
if aq:a_text+=f"保留的 jplace 分母：{aq.get('distinct_bacterial_placed_queries')} 个不同 bacterial query，{aq.get('pplacer_calls')} 次调用。\n\n"
if asel:a_text+=f"线程回放选择最慢实际调用：基线 {asel['wall_s']:.2f} s，覆盖 {asel['queries']} 个 query；不能外推为所有 query 的线程等价性。\n\n"
if at:a_text+=f"100 样本重跑的 taxonomy 完全一致：{at['taxonomy_exact']}；全部 summary 字段完全一致：{at['all_summary_fields_exact']}。这是同一工作流保留中间文件的复现。\n\n"
if ac:
 gaps=[x for x in ac['per_query'] if x['top_two_weight_gap'] is not None];gaps.sort(key=lambda x:x['top_two_weight_gap'])
 a_text+=f"回放 query 集合相同：{ac['query_sets_equal']}；基线与 1/8 线程的全部 jplace 数值精确相同：{ac['baseline_vs_1_exact']} / {ac['baseline_vs_8_exact']}。每条 query 的 edge、LWR、数值差异及 top-two gap 已写入 jplace_comparison.json；保留多候选 edge 的 query 数 {len(gaps)}。\n\n"
 if gaps:a_text+='最接近的候选 edge（按已有输出的 top-two LWR gap 排序，无额外阈值）：\n\n'+''.join(f"- {x['query']}：gap={x['top_two_weight_gap']:.6g}。\n" for x in gaps[:5])+'\n'
for r,text,report in [('A',a_text,'A_pplacer/REPORT.md'),('C','## 本次补测结果\n\n控制器状态：`'+c_state+'`。\n\n','C_hla/REPORT.md')]:
 m=measurements[r]
 for t in ['1','8']:
  v=m[t]
  if v.get('wall_s') is not None:text+=f"{t} 线程实际后端回放：wall {v['wall_s']:.2f} s，peak RSS {v['peak_rss_kib']:.0f} KiB。\n\n"
 if 'one_over_eight_wall_ratio' in m:text+=f"1/8 线程 wall 比 {m['one_over_eight_wall_ratio']:.3f}×；同一现有二进制线程对照，不能称新原型加速。\n\n"
 if r=='C' and cc:
  text+=f"基线与单线程/8 线程 raw allele sets 精确相同：{cc['baseline_kernel1_exact_allele_sets']} / {cc['baseline_kernel8_exact_allele_sets']}；G calls 精确相同：{cc['baseline_kernel1_exact_G_calls']} / {cc['baseline_kernel8_exact_G_calls']}。\n\n"
  for name,s in cc['truth_summary'].items():text+=f"- {name}：两字段明确匹配 {s['confirmed_two_field_alleles']}/10，歧义集合兼容 {s['compatible_two_field_alleles']}/10；五个位点明确 genotype 匹配 {s['unambiguous_matching_genotypes']}/5。\n"
  text+='\n四字段独立真值仍不可用。完整原始歧义集合与每个位点的杂合/纯合/不确定状态见 call_comparison.json。\n\n'
 if r=='C':
  recovery_comparison=read('profiles/closure_20261002/C_recovery/call_comparison.json')
  if recovery_comparison:text+='恢复后函数采样的调用比较：'+json.dumps(recovery_comparison,ensure_ascii=False)+'\n\n'
 if (m.get('perf_status') or {}).get('exit_code') not in [None,0]:text+='注意：本次 perf 失败，以下仅为失败前片段，不是完整后端热点。\n\n'
 if m.get('flat_profile_lines'):text+='函数采样（49 Hz 用户态 CPU，非端到端 wall 占比）：\n\n```text\n'+'\n'.join(m['flat_profile_lines'])+'\n```\n'
 text+='\n原型仍为 NOT_RUN；本轮只验证现有程序和证据完整性。新架构收益未得到验证，维持 HOLD，不给 GO。\n'
 supplemental(report,text)
 p=ROOT/report;s=p.read_text();s=s.replace('补测进行中','补测未完成').replace('补测完成，仍有范围缺项','补测未完成')
 if done:
  s=s.replace('补测未完成','补测完成，仍有范围缺项')
  if r=='C':s=s.replace('**下方函数采样表仅覆盖失败前执行片段，不能当完整 graph/inference 热点分布。** 当前无计算在运行，本次仅检查，没有重启。完整函数诊断仍待单独恢复。','该段描述旧 perf_8 失败现场，已保留。恢复采样采用 perf_8_retry_20261002，下方结果来自恢复后的完整运行；不再将旧失败片段当完整热点。')
 p.write_text(s)
report={'updated_at':datetime.datetime.now().astimezone().isoformat(),'local_measurement_jobs_complete':done,'any_job_failed':failed,'route_controller_states':{'A':a_state,'C':c_state},'scientific_comparisons':checks,'measurements':measurements,'decisions':{'A':'HOLD','B':'KILL','C':'HOLD','D':'KILL'},'GO_count':0,'external_gaps':['HLA-HD software / RNA typing','NVMe benchmark'],'scope_gaps':['No unique four-field HLA truth','D full-chromosome equivalence and interval/batch/IO sweeps NOT_RUN','B fine-grained parse/IO/scheduling wall split unresolved; no modified prototype for remask comparison','A/C no architecture prototype; thread replay is not an innovation']}
(E/'closure_assessment.json').write_text(json.dumps(report,indent=2))
subtitle='本机补测与报告收尾完成；原任务全范围仍有明确缺项' if done else '补充收尾进行中；A/C 尚未全部完成'
(ROOT/'00_EXECUTIVE_DECISION.md').write_text('''# 四路线 Kill-Test 决策

'''+subtitle+'''

**当前 0 条 GO；A/C HOLD，B/D KILL。** HOLD 不是建议进入大重写阶段，KILL 是停止当前范围的新项目投入，不是证明任何算法或硬件下都无优化空间。

| Route | Current usage | Collision risk | Critical-path share | Theoretical ceiling | Prototype result | Scientific equivalence risk | Decision |
|---|---|---|---|---|---|---|---|
| A GTDB-Tk / pplacer | 2.7.2/R232 标准工作流 | 已有 mmap、gap masking、split-tree/reference 共享 | W2 pplacer 90.10%，54 个实际 placement | 假设全部 pplacer 3×，端到端约 2.50×；未证实 | 无新架构原型；线程/热点补测见 A 报告 | 单调用 edge/LWR 与全流程 taxonomy 比较范围分开 | HOLD |
| B RepeatModeler2 | EBP/Earl Grey 的现有使用证据 | 近期 refinement 已优化；不允许重写 aligner | 完整 small comparison 85.64% | comparison 不变时删除其余工作仅 1.168× | 未做；两批函数热点支持对齐主导 | 无修改版，无 paired remask 等价性结论 | KILL |
| C HLA-LA / HLA-HD | 现有 workflow/benchmark 使用证据 | 官方 fast 分支与既有 prefilter/索引复用 | WGS 首次 index 6713.81 s；WES backend 其余 73.6%，非纯 DP | 一次性建索引不能作新收益；单阶段可加速空间未证实 | 无新原型；1/8 线程与调用比较见 C 报告 | 低覆盖 WGS 真值表现弱；ambiguity 与唯一 allele 不等价 | HOLD |
| D GenomicsDB | GATK/native 标准路径 | 高：native buffering 与 bypass reader 已有 | 已完成 50/200 完整 chr22 三组导入 | 已有选项最多观察到 2.107× 时间、4.139× RSS 比；非理论普遍上界 | 无新原型，DP 元数据诊断仅解释输出 | 两规模五区间差异由 DP sum 配置解释；非全染色体等价 | KILL |

## 已闭合的证据与仍未完成的范围

B 完整果蝇基线 4:05:01、461 families，水稻按原任务许可仅前三轮。D 50/200 样本完整导入与五区间差异诊断完成；200 样本 53 条 INFO/DP 差异在仅改合并元数据后全部消失。

A 旧零计数已纠正：GTDB-Tk 清理中间文件导致后处理误判，不是 100 个 MAG 都被 ANI 分流。C 旧回放缺失 hla_nom_g.txt 已定位为工作目录问题。A/C 本次状态：`'''+a_state+'` / `'+c_state+'''`；实际测量与一致性细节见两份报告及 profiles/closure_20261002/closure_assessment.json。

HLA-HD 软件/RNA 分型仍 BLOCKED_SOFTWARE；当前主机无 NVMe，对照 NOT_RUN。唯一四字段 HLA 真值不可用。D interval/batch 扫描、全量 syscall/decode/sequential-read 分解及全染色体等价性未完成；B parse/IO/scheduling 精确分离未完成。上述缺项均未用 proxy 或文献数字替代，不宣称原任务所有必做项已验收。

A/C 没有得到新架构达到门槛的原型实证，暂不升级 GO。B/D 在已有测量、碰撞与任务允许的修改范围内停止投入；不为获得正面结果补写替代引擎。

## 对 PLAN 的影响

H-A 的关键路径前提获得支持，但可压缩空间未证实；H-B 被完整基线削弱；H-C 尚受科学准确性/输入范围和原型证据限制；H-D 未支持独立新架构。PLAN 保持原文。最终最多一条 GO 的要求满足，本轮实际为零。
''')
subprocess.run([sys.executable,str(ROOT/'benchmarks/summarize_profiles.py')],check=True)
p=ROOT/'benchmarks/summary.tsv'
with p.open() as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;rows=list(r)
for route,workload,path,n,unit in [('A','W2_retained_20261002','profiles/A/W2_retained_20261002',100,'input_genomes')]+[(r,f'closure_threads_{t}',f'profiles/{r}_scaling_20261002/threads_{t}',(asel or {}).get('queries','NA') if r=='A' else 1,'placed_queries' if r=='A' else 'sample_assay') for r in ['A','C'] for t in [1,8]]:
 status=read(path+'/status.json');row={'route':route,'workload':workload,'n':n,'unit':unit,'run_state':status.get('status','UNKNOWN') if status else ('FAILED_OR_INCOMPLETE' if stage_status(route)=='FAILED' else 'PENDING'),'profile_path':path,'storage':'local_rotational_sda','measurement_scope':'same_binary_thread_replay_not_new_prototype' if 'threads' in workload else 'repeat_classify_keep_intermediates','equivalence':'SEE_CLOSURE_ASSESSMENT',**timing(path)};rows.append(row)
with p.open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows({k:r.get(k,'NA') for k in fields} for r in rows)
if done:
 status='''# 当前状态
本机可执行的 A/C 补测与 B/D 报告收尾已经完成，四份报告、决策表和汇总均已更新。当前 A/C 保留观察，B/D 停止当前新引擎方向，没有推荐施工的主路线。

当前补测服务已完成，未安排新的计算。A 的实际 placement 分母为 54，旧零计数已纠正；B 完整果蝇基线与 D 两种规模的导入、差异诊断已有完整证据，A/C 线程与输出比较详见报告。

原任务全范围仍缺 HLA-HD/RNA 分型与 NVMe 条件，唯一四字段 HLA 真值也不可用。报告列明尚未完成的细分 profiling、全染色体等价性和原型，不把本机收尾等同于原任务所有项目通过。

目前以明确缺项的阶段报告交付，不启动大重写。后续只有取得缺失软件/硬件或新的关键证据后才恢复对应补测；H-A 的关键路径得到支持，H-B 被削弱，H-C/H-D 没有达到推荐施工的证据。

'''
else:
 status='''# 当前状态
B 完整果蝇基线分析及 D 200 样本差异诊断已完成，两条路线停止当前新引擎投入。A 旧零计数已纠正，实际有 54 个 placement；C 回放工作目录问题已修正。

A/C 补测仍未全部完成，具体状态见下方控制器信息。新结果写入独立目录，成功旧基线与失败证据均保留。

HLA-HD/RNA 分型缺软件，NVMe 对照缺硬件；这些条件尚未补齐。C 一次性 BWA 建索引成本与重复推断分开报告，严格两字段真值并未建立高准确率。

下一步依据完成事件汇总 A/C 线程、输出一致性与函数热点，自动更新报告并保留未完成项。H-A 的关键路径得到支持，H-B 被完整基线削弱；仍无 GO，不启动大重写。

'''
 status+='2026-10-02\nROADMAP [#--] 1/3 节点；R2 补测尚未完成。\n'
if done:status+='2026-10-02\nROADMAP [##-] 2/3 节点：R1、R3 完成，R2 有外部条件和明确范围缺项。\n'
status+='本周投入：未统计工时；以基线解释、科学输出比较及报告为主。\n偏离程度：中\n偏离位置：R2 缺软件/硬件及部分细分测量，不能宣称全范围验收。\n建议：不扩展替代引擎，按报告保留缺项边界。\n\n## 接手信息\n活跃节点：'+('无新计算；R2 缺项、R3 阶段报告已交付。' if done else 'R2 补测，R3 报告随证据更新。')+'\n核心文件：benchmarks/finalize_closure.py、benchmarks/post_a.py、benchmarks/post_c.py。\n复现汇总：python3 benchmarks/finalize_closure.py。\n控制器：A '+a_state+'；C '+c_state+'。\n证据入口：profiles/closure_20261002/closure_assessment.json。\n最近决策：D017。\n下一步：'+('需软件/硬件或新证据才恢复对应 HOLD 分支。' if done else '等待完成事件；若 FAILED 先读对应 controller.log。')+'\n'
(ROOT/'STATUS.md').write_text(status)
if done:
 p=ROOT/'TODO.md';s=p.read_text();s=s.replace('- [ ] 修正 A 中间文件清理导致的零计数，保留 54 个实际 placement 的回放证据。','- [x] 修正 A 分母并完成保留中间文件、线程及输出比较。（2026-10-02）').replace('- [ ] 修复 C 单线程回放查找 hla_nom_g.txt 的工作目录，继续后续比较。','- [x] 修复 C 回放路径，完成线程、调用及真值比较。（2026-10-02）');s=s.replace('- [ ] 完成 A/C 补测后的报告、决策表与汇总收尾。','- [x] 完成 A/C 补测后的报告、决策表与汇总收尾。（2026-10-02；外部缺项保留）')
 s=s.replace('- [ ] 单独重跑 C 完整函数采样，比较输出并更新最终报告。','- [x] C 完整函数采样恢复完成，调用比较与最终报告已更新。（2026-10-02）')
 if '2026-10-02 完成事件：' not in s:s+='\n2026-10-02 完成事件：A/C 计划补测退出成功，报告与 summary 更新；无 GO，原任务 HLA-HD/RNA、NVMe 等缺项仍未完成。\n';p.write_text(s)
 p=ROOT/'ROADMAP.md';p.write_text('''# 路线
1. [infra] R1 本地资产盘点与四路线近期碰撞审计；服务 H-A 至 H-D。已完成首轮有限审计。
2. R2 现实工作负载 profiling 与必要的最小原型；直接测量 H-A 至 H-D。本机补测完成，但 HLA-HD/RNA、NVMe 及报告明确列出的细分测量仍有缺项；未做新架构原型，不算全范围验收。
3. R3 汇总科学语义风险、条件 Amdahl 上限与 GO/HOLD/KILL；判断 H-A 至 H-D。阶段报告完成：A/C HOLD、B/D KILL、零 GO，所有缺项随结论保留。
''')
paths=['00_EXECUTIVE_DECISION.md','A_pplacer/REPORT.md','B_repeatmodeler/REPORT.md','C_hla/REPORT.md','D_genomicsdb/REPORT.md','benchmarks/summary.tsv','profiles/closure_20261002/closure_assessment.json','profiles/closure_20261002/B_evidence.json','profiles/closure_20261002/D_evidence.json','profiles/closure_20261002/C_baseline_evidence.json','profiles/D_schema_diagnostic_20261002/comparison.json']
for candidate in ['profiles/closure_20261002/C_recovery/permission_probe.json','profiles/closure_20261002/C_recovery/call_comparison.json','profiles/closure_20261002/C_recovery/status.json','profiles/C_scaling_20261002/perf_8_retry_20261002/status.json','profiles/C_scaling_20261002/perf_retry_spec.json','profiles/A_scaling_20261002/summary_field_diagnosis.json','profiles/A_scaling_20261002/jplace_comparison.json','profiles/C_scaling_20261002/call_comparison.json']:
 if (ROOT/candidate).exists():paths.append(candidate)
(E/'SHA256SUMS').write_text(''.join(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()+'  '+p+'\n' for p in paths))
print(json.dumps({'jobs_complete':done,'failed':failed,'decisions':report['decisions'],'reports_updated':True}))
