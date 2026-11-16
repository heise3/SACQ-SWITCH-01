"""SACQ-SWITCH-01 Phase 1. Read-only Phase-0 adapters, fixed gene modules."""
from pathlib import Path
import ast,hashlib,json,sys,os
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
SRC=Path(os.environ.get('SACQ_INPUTS',str(ROOT/'inputs')))
CACHE=SRC/'cache'
SEED=20261116

def save(obj,name):
 def clean(x):
  if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
  if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
  if isinstance(x,(np.integer,)):return int(x)
  if isinstance(x,(np.floating,float)):return float(x) if np.isfinite(x) else None
  if isinstance(x,np.bool_):return bool(x)
  return x
 (ROOT/name).write_text(json.dumps(clean(obj),indent=2,allow_nan=False))

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()

# Original pure parsers, separated from the historical execution driver.
import phase0_parsers
ns=vars(phase0_parsers)
ns['CACHE_DIR']=str(CACHE)
MODULES=ns['CANONICAL_MODULES'];score=ns['compute_module_scores'];parse=ns['parse_series_matrix'];normalize=ns['normalize_gene_symbol']

def load65():
 fn=ROOT/'work/gse65391.pkl'
 if fn.exists():return pd.read_pickle(fn)
 meta,gene,info=parse(str(CACHE/'GSE65391_series_matrix.txt.gz'),ns['load_gpl10558_map']())
 raw=len(meta);keep=~meta.batch_replicate.str.lower().isin(['yes','true','1']);meta=meta.loc[keep].copy();gene=gene.loc[keep]
 meta['visit_num']=pd.to_numeric(meta['visit'],errors='coerce');meta['sledai_num']=pd.to_numeric(meta.sledai,errors='coerce')
 hc=meta['disease state'].str.lower().eq('healthy');sle=meta['disease state'].eq('SLE')&meta.visit_num.notna()&meta.sledai_num.notna()
 meta=meta.loc[hc|sle].copy();gene=gene.loc[meta.index]
 meta['_dedup_key']=np.where(meta['disease state'].eq('Healthy'),'HC_'+meta.subject,'SLE_'+meta.subject+'_V'+meta.visit_num.astype(str))
 meta=meta.loc[~meta._dedup_key.duplicated()].copy();gene=gene.loc[meta.index]
 hc=meta['disease state'].eq('Healthy').values;scores,cov=score(gene,hc);d=pd.concat([meta,scores],axis=1)
 for c in ['c3','c4','ds_dna','age','renal','cumulative_time','days_since_last_visit','neutrophil_count','oral_steroids_category','steroid_iv_category','mycophenolate_category']:
  d[c+'_num']=pd.to_numeric(d[c],errors='coerce')
 d['Sero_SLEDAI']=2*d.increased_dna_binding.isin(['2','2.0','Yes','yes','1'])+2*d.low_complement.isin(['2','2.0','Yes','yes','1'])
 d['cSLEDAI']=(d.sledai_num-d.Sero_SLEDAI).clip(lower=0)
 d['state']=np.select([(d.cSLEDAI==0)&(d.Sero_SLEDAI==0),(d.cSLEDAI==0)&(d.Sero_SLEDAI>0)],['CSQ','SACQ'],default='ClinActive')
 d.to_pickle(fn);save({'raw_samples':raw,'retained':len(d),'sle_visits':int((~hc).sum()),'sle_patients':d.loc[~hc,'subject'].nunique(),'unique_hc':int(hc.sum()),'coverage':cov,'processing':info},'results/gse65391_processing.json')
 return d

def pair_summary(a,b):
 from scipy import stats
 a,b=np.asarray(a,float),np.asarray(b,float);keep=np.isfinite(a)&np.isfinite(b);a,b=a[keep],b[keep];z=b-a;n=len(z)
 if n<2:return {'n':n}
 se=z.std(ddof=1)/np.sqrt(n);ci=stats.t.interval(.95,n-1,loc=z.mean(),scale=se)
 return {'n':n,'mean_base':a.mean(),'mean_followup':b.mean(),'delta':z.mean(),'ci95':list(ci),'paired_dz':z.mean()/z.std(ddof=1) if z.std(ddof=1)>0 else None,'p':stats.ttest_rel(b,a).pvalue,'wilcoxon_p':stats.wilcoxon(z).pvalue if np.any(z) else 1}
