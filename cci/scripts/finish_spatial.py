import sys,json
from pathlib import Path
import numpy as np,pandas as pd
from cci_core import exact_cover
p=Path(sys.argv[1]);r=pd.read_csv(p/'mapped_resource.tsv',sep='\t');f=pd.read_csv(p/'features.tsv',sep='\t');labels=f.sender.unique().tolist();k=len(r)
c=np.load(p/'contact_communication.npz');s=np.load(p/'secreted_communication.npz');mix=np.load(p/'mixing_communication.npz')
is_contact=np.tile((r.annotation=='Cell-Cell Contact').to_numpy(),len(labels)**2)[None,:]
score=np.where(is_contact,c['scores'],s['scores']);support=np.where(is_contact,c['support'],s['support'])
np.savez_compressed(p/'typed_communication.npz',scores=score,support=support)
rows=[];simple=~(r.ligand_subunits.str.contains('|',regex=False)|r.receptor_subunits.str.contains('|',regex=False))
for a in range(len(labels)):
 for b in range(len(labels)):
  sl=slice((a*len(labels)+b)*k,(a*len(labels)+b+1)*k);active=support[0,sl]&simple.to_numpy();z=r[active];ic,targets=exact_cover(z.ligand_subunits,z.receptor_subunits)
  rows.append(dict(mode='typed_kernel',sample='visium_anterior1',sender=labels[a],receiver=labels[b],count=len(z),strength=float(score[0,sl][active].sum()),ic100=ic,targets=';'.join(targets)))
pd.DataFrame(rows).to_csv(p/'typed_programs.tsv',sep='\t',index=False)
summary={'mixing_supported_channels':int(mix['support'].sum()),'typed_supported_channels':int(support.sum()),'removed_by_geometry':int((mix['support']&~support).sum()),'added_by_geometry':int((~mix['support']&support).sum()),'kernel_rule':'Cell-Cell Contact:110um binary radius; Secreted Signaling:250um truncated exponential scale100um','unannotated_group':'retained, not a biological cell-type claim','active_program_pairs_mixing':int(mix['support'].reshape(-1,k).any(axis=1).sum()),'active_program_pairs_typed':int(support.reshape(-1,k).any(axis=1).sum())}
(p/'typed_summary.json').write_text(json.dumps(summary,indent=2))
print(summary)
