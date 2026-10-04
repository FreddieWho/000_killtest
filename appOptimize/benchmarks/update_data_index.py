"""Register completed download receipts; partial downloads are excluded."""
import argparse, csv, datetime, fcntl, io, json, pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--sync-central',action='store_true');a=p.parse_args()
dest=ROOT/'infra/bioinf-data-index';dest.mkdir(parents=True,exist_ok=True)
lock=(dest/'.index.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
index=dest/'BIOINF_DATA_INDEX.tsv'
rows={}
if index.exists():
    with index.open() as f:
        rows={r['path']:r for r in csv.DictReader(f,delimiter='\t')}
for receipt in (ROOT/'data').rglob('*.receipt.json'):
    r=json.loads(receipt.read_text())
    if not r.get('path') or not r.get('sha256'):
        continue
    file=ROOT/r['path']
    if not file.exists() or file.stat().st_size!=int(r['bytes']) or pathlib.Path(str(file)+'.aria2').exists():
        continue
    rows[r['path']]=dict(r,origin='download_or_remote_region_extract',status='COMPLETE')
fields=['name','path','origin','status','url','bytes','sha256','md5','sequences','bases','records','region','downloaded_utc','purpose']
temporary=index.with_suffix('.tsv.tmp')
with temporary.open('w') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows[k] for k in sorted(rows))
temporary.replace(index)
total=sum(int(r['bytes']) for r in rows.values())
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
summary=f'已完成外部资产 {len(rows)} 个，按已注册文件计 {total} bytes；下载中的稀疏文件不计入。'
(dest/'BIOINF_DATA_INVENTORY.md').write_text(f'# 本轮新增外部数据\n\n更新时间：{now}\n\n{summary}\n\n来源、哈希与用途见 BIOINF_DATA_INDEX.tsv；派生解压副本不重复计为独立数据集。参考图和索引也不是独立样本。现有公共 CTAT GRCh38 参考只读复用，见 data/D_genomicsdb/reference.reuse.json。\n')
if a.sync_central:
    central=pathlib.Path('/home/huyudi/Infra/bioinf-data-index')
    def append_new(file,key,records):
        with file.open('r+') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            reader=csv.DictReader(f,delimiter='\t');header=reader.fieldnames;seen={r[key] for r in reader}
            f.seek(0,2);w=csv.DictWriter(f,fieldnames=header,delimiter='\t')
            for r in records:
                if r[key] not in seen:w.writerow(r);seen.add(r[key])
            f.flush()
    ext=[];inventory=[]
    for r in rows.values():
        file=ROOT/r['path'];relative=file.relative_to('/home/huyudi')
        ext.append(dict(dataset_id='appOptimize_'+r['name'],accession=r['name'],source_url=r['url'],
            license_or_consent='SEE_PUBLIC_SOURCE_TERMS',role=r.get('purpose','benchmark_asset'),local_path=str(file),
            sha256='sha256:'+r['sha256'],bytes=r['bytes'],registered_utc=now))
        inventory.append(dict(project_root='000_killtest',relative_path=str(relative),format=''.join(file.suffixes).lstrip('.'),
            category='public_bioinformatics_benchmark_asset',bytes=r['bytes'],
            modified_utc=datetime.datetime.fromtimestamp(file.stat().st_mtime,datetime.timezone.utc).isoformat(),same_size_group=''))
    append_new(central/'BIOINF_DATA_EXTERNAL.tsv','local_path',ext)
    append_new(central/'BIOINF_DATA_INDEX.tsv','relative_path',inventory)
    with (central/'BIOINF_DATA_INVENTORY.md').open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        f.write(f'\n### appOptimize 登记更新 {now}\n\n{summary} 当前细目：`/home/huyudi/000_killtest/appOptimize/infra/bioinf-data-index/BIOINF_DATA_INDEX.tsv`。\n')
print(json.dumps({'registered_assets':len(rows),'bytes':total,'central_synced':a.sync_central}))
