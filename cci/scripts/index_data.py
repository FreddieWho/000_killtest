"""Register downloaded sources and local exports; hashes are actual file content."""
from pathlib import Path
import hashlib,csv,json
root=Path(__file__).resolve().parents[1]
base='https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/'
entries=[
 ('data/raw/kang_2018.h5ad','https://ndownloader.figshare.com/files/34464122','Kang GSE96583 processed counts',24673,'8 donors;16 paired sample units','statistics;IC;profiling'),
 ('data/raw/visium_filtered_feature_bc_matrix.h5',base+'V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5','10x anterior1 filtered counts',2695,'one section;not biological replication','spatial operator;profiling'),
 ('data/raw/visium_spatial.tar.gz',base+'V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz','10x anterior1 coordinates/scales','NA','same one section','spatial coordinates and physical scale'),
 ('data/raw/visium_mouse_cortex_annotated.RData','https://ndownloader.figshare.com/files/41446839','CellChat author processed annotation',1073,'same one section','existing spot labels only;full expression from 10x'),
 ('data/raw/cellphonedb_v5_source.zip','https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v5.0.0.zip','CellPhoneDB v5 official archive','NA','not observations','native baselines'),
 ('data/resources/cellchat_human.tsv','R installed CellChat; scripts/export_cellchat_resources.R','CellChat human export','NA','LR prior;not data replication','IC and reference'),
 ('data/resources/cellchat_mouse.tsv','R installed CellChat; scripts/export_cellchat_resources.R','CellChat mouse export','NA','LR prior;not data replication','same operator spatial comparison'),
 ('data/resources/liana_cellphonedb.tsv','installed liana 1.6.1 resource/omni_resource.csv; cellphonedb subset','LIANA packaged CPDB resource','NA','distinct packaged catalogue from native v5','database sensitivity'),
 ('data/resources/visium_annotation.tsv','data/raw/visium_mouse_cortex_annotated.RData','exported annotations (1072 match official filtered counts)','NA','same one section','barcode-matched labels')]
entries += [
 ('data/resources/cellphonedb_v5.zip','data/raw/cellphonedb_v5_source.zip','official nested database archive','NA','LR prior','native CPDB'),
 ('data/resources/tissue_positions_list.csv','data/raw/visium_spatial.tar.gz','coordinates including off-tissue rows','NA','one section','barcode-matched geometry'),
 ('data/resources/scalefactors_json.json','data/raw/visium_spatial.tar.gz','image scale metadata','NA','one section','scale audit'),
 ('sources/GSE96583.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583&targ=self&form=text&view=brief','official GEO study metadata','NA','source accession','data source verification')]
for p in sorted((root/'data/resources/cpdb_v5').glob('*.csv')):
 entries.append((str(p.relative_to(root)),'data/resources/cellphonedb_v5.zip','extracted official CPDB table','NA','LR prior','native FastCCC'))
for p in sorted((root/'data/resources').glob('*cofactor_genes.txt')):
 entries.append((str(p.relative_to(root)),'R installed CellChat DB cofactor table','native cofactor gene export','NA','LR prior','complete native baseline input'))
rows=[]
for rel,url,title,cols,units,use in entries:
 p=root/rel
 if not p.is_file():continue
 rows.append(dict(path=rel,source=url,description=title,bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),expression_columns=cols,biological_units=units,intended_use=use,registered_date='2026-09-28'))
out=root/'infra/bioinf-data-index';out.mkdir(parents=True,exist_ok=True)
with (out/'BIOINF_DATA_INDEX.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(out/'BIOINF_DATA_INVENTORY.md').write_text('# CCI 外部生信数据索引\n\n2026-09-28；项目范围，不代表全主机资产。\n\n'+f'已登记 {len(rows)} 个原始输入或导出文件，合计 {sum(r["bytes"] for r in rows):,} bytes（包含派生文件，不是独立数据量）。\n\n'+'Kang：24,673 cells，8 donors，16 sample units；Visium：2695 spots、一张切片；1072条作者注释匹配，1623未注释，另1条作者barcode不在filtered counts。不能用于供体层推断。作者注释对象与 10x counts 按 barcode 对齐；不得将副本或空间块当成独立重复。\n\n各文件来源、SHA256、字节数及用途见 BIOINF_DATA_INDEX.tsv。源归档中的其他示例不会自动纳入分析；只提取需要的数据库/坐标。\n')
print('registered',len(rows),'files')
