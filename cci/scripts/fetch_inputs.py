"""Fetch only named small task inputs; reuse existing files; atomic download."""
from pathlib import Path
import urllib.request,zipfile,tarfile,hashlib,csv,sys
root=Path(__file__).resolve().parents[1]
base='https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/'
files={
 'data/raw/kang_2018.h5ad':'https://ndownloader.figshare.com/files/34464122',
 'data/raw/visium_filtered_feature_bc_matrix.h5':base+'V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5',
 'data/raw/visium_spatial.tar.gz':base+'V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz',
 'data/raw/visium_mouse_cortex_annotated.RData':'https://ndownloader.figshare.com/files/41446839',
 'data/raw/cellphonedb_v5_source.zip':'https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v5.0.0.zip'}
index=root/'infra/bioinf-data-index/BIOINF_DATA_INDEX.tsv'
expected={r['path']:r['sha256'] for r in csv.DictReader(index.open(),delimiter='\t')} if index.exists() else {}
for rel,url in files.items():
 p=root/rel;p.parent.mkdir(parents=True,exist_ok=True)
 if not p.exists():
  if '--download' not in sys.argv:raise FileNotFoundError(str(p)+' (use --download)')
  part=p.with_suffix(p.suffix+'.part');urllib.request.urlretrieve(url,part);part.rename(p)
 actual=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
 if rel in expected and actual!=expected[rel]:raise ValueError('Source hash changed: '+rel)
 print(rel,p.stat().st_size,actual)
res=root/'data/resources';res.mkdir(exist_ok=True)
with zipfile.ZipFile(root/'data/raw/cellphonedb_v5_source.zip') as z:
 n=next(n for n in z.namelist() if n.endswith('/cellphonedb.zip'));(res/'cellphonedb_v5.zip').write_bytes(z.read(n))
(res/'cpdb_v5').mkdir(exist_ok=True)
with zipfile.ZipFile(res/'cellphonedb_v5.zip') as z:
 for n in z.namelist():
  if '/' not in n and n.endswith('.csv'):(res/'cpdb_v5'/n).write_bytes(z.read(n))
with tarfile.open(root/'data/raw/visium_spatial.tar.gz') as z:
 for n in z.getmembers():
  if n.isfile() and Path(n.name).name in ['tissue_positions_list.csv','scalefactors_json.json']:(res/Path(n.name).name).write_bytes(z.extractfile(n).read())
