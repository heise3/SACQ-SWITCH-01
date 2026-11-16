"""Task 3 adult open external therapy validation; preserve original sensitivities."""
from common import *
import re,gzip
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
assert (ROOT/'results/task2_results.json').exists()
# Verify complete gzip before analysis, fetch task must complete.
with gzip.open(ROOT/'sources/GSE72747_series_matrix.txt.gz','rb') as f:
 while f.read(1024*1024):pass
map570={}
with gzip.open(ROOT/'sources/GPL570.annot.gz','rt') as f:
 active=False;header=None
 for line in f:
  if line.startswith('!platform_table_begin'):active=True;continue
  if line.startswith('!platform_table_end'):break
  if not active:continue
  a=line.rstrip('\n').split('\t')
  if header is None:header=a;gi=header.index('Gene symbol');continue
  gsym=normalize(a[gi]) if len(a)>gi else None
  if gsym:map570[a[0].strip()]=gsym
m,g,inf=parse(str(ROOT/'sources/GSE72747_series_matrix.txt.gz'),map570)
assert g.shape[1]>40, 'insufficient probe matching';scores,cov=score(g,np.zeros(len(m),bool));d=pd.concat([m,scores],axis=1)
d['subject']=d.patient;d['month']=d.time.map({'prior to administration of immunosuppressive therapy':0,'3 months':3,'6 months':6})
for c,k in [('sledai','SLEDAI'),('glomerular filtration rate','GFR'),('dsdna ab titers (u/l by farr assay)','DNA'),('age','age_num')]:d[k]=pd.to_numeric(d[c],errors='coerce')
assert d.subject.nunique()==10 and len(d)==30;assert d.age_num.min()>=18
result={'adult_LN_blood':{'n_patients':d.subject.nunique(),'n_visits':len(d),'age_range':[d.age_num.min(),d.age_num.max()],'coverage':cov,'scale':inf,'anchor':'whole cohort mean/SD; no healthy controls, same fixed gene modules','available_labels':['SLEDAI total','quantitative Farr anti-dsDNA','GFR','treatment','dates'],'unavailable_labels':['C3','C4','full SLEDAI individual items'],'SACQ_replication':False},'paired':[],'FE':[],'sensitivities':{},'partial_kidney':{}}
for tp in [3,6]:
 for k in ['M_IFN','M_PB','M_LDG','ESI','NC_HK','SLEDAI','GFR','DNA']:
  p=d.pivot(index='subject',columns='month',values=k)[[0,tp]].dropna();r=pair_summary(p[0],p[tp]);r.update({'timepoint_month':tp,'measure':k});result['paired'].append(r)
# Small-cluster t inference, participant FE. Total disease activity, not clinical-only SLEDAI.
d['logDNA']=np.log1p(d.DNA)
for k in ['M_IFN','M_PB','M_LDG']:
 for outcome in ['SLEDAI','GFR']:
  ss=d.dropna(subset=[k,outcome,'logDNA']);fit=smf.ols(f'{outcome} ~ {k} + logDNA + C(month) + C(subject)',data=ss).fit(cov_type='cluster',cov_kwds={'groups':ss.subject},use_t=True)
  result['FE'].append({'outcome':outcome,'module':k,'n':len(ss),'patients':ss.subject.nunique(),'beta':fit.params[k],'ci95':fit.conf_int().loc[k].tolist(),'p':fit.pvalues[k]})
d.to_csv(ROOT/'results/task3_GSE72747_source.csv',index_label='gsm')
# GSE49454: reproduce exact frozen operational state, patient-level paired mean contrasts.
m,g,info=parse(str(CACHE/'GSE49454_series_matrix.txt.gz'),ns['load_gpl10558_map']());hc=m.group.str.contains('Healthy',case=False).values;sc,coverage=score(g,hc);q=pd.concat([m,sc],axis=1);q=q.loc[~hc].copy()
q['subject']=q._title.str.extract(r'^(C\d*PL\d+|PL\d+)',expand=False).fillna(q._title);q['sledai_num']=pd.to_numeric(q.sledai,errors='coerce');q['DNA']=pd.to_numeric(q['a-dsd'],errors='coerce')
pos=lambda x:x.astype(str).str.lower().isin(['yes','y','1','true','pos'])
q['sero']=(q.DNA.fillna(0)>=20)|pos(q['low c3'])|pos(q['low c4'])|pos(q['a-dsd'])
q['active']=pos(q.renal)|pos(q.joints)|pos(q.skin)|pos(q.hemato)|(q.sledai_num>4)
q['state']=np.where(~q.active&q.sero,'SACQ',np.where(q.active,'ClinActive','CSQ'));pp=q.groupby(['subject','state'])[['M_IFN','M_LDG','M_PB','ESI','NC_HK']].mean().unstack('state');rs={}
for k in ['M_IFN','M_LDG','M_PB','ESI','NC_HK']:
 z=pp[k][['SACQ','ClinActive']].dropna();rs[k]=pair_summary(z.SACQ,z.ClinActive)
result['sensitivities']['GSE49454']={'n_patients':q.subject.nunique(),'n_visits':len(q),'paired':rs,'definition':'Frozen operational SACQ rule; organ flags plus total SLEDAI>4; not a complete validated cSLEDAI reconstruction'}
q.to_csv(ROOT/'results/task3_GSE49454_source.csv',index_label='gsm')
# GSE121239: total-SLEDAI operational low/high; not SACQ and not causal counterfactual.
m,g,info=parse(str(CACHE/'GSE121239_series_matrix.txt.gz'),ns['load_gpl13158_map']());hc=m['disease state'].eq('Healthy').values;sc,cov121=score(g,hc);q=pd.concat([m,sc],axis=1);q=q.loc[~hc].copy();q['subject']=q['patient id'];q['sledai_num']=pd.to_numeric(q.sledai,errors='coerce');q['state']=np.select([q.sledai_num<=4,q.sledai_num>=6],['Low','High'],default='Intermediate');pp=q.groupby(['subject','state'])[['M_IFN','M_LDG','M_PB','ESI','NC_HK']].mean().unstack('state');rs={}
for k in ['M_IFN','M_LDG','M_PB','ESI','NC_HK']:
 z=pp[k][['Low','High']].dropna();rs[k]=pair_summary(z.Low,z.High)
result['sensitivities']['GSE121239']={'n_patients':q.subject.nunique(),'n_visits':len(q),'paired':rs,'definition':'total SLEDAI<=4 vs >=6; no item decomposition; PBMC; descriptive sensitivity, not identified counterfactual'}
q.to_csv(ROOT/'results/task3_GSE121239_source.csv',index_label='gsm')
# Adult serial kidney cohort: use sample-table expression and metadata, no ambiguous supplementary column mapping.
p=ROOT/'sources/GSE200306_series_matrix.txt.gz'
with gzip.open(p,'rt') as f:
 for ln in f:
  if ln.startswith('"ID_REF"'):break
 for _ in range(3):print('200306 matrix row',next(f)[:150],flush=True)
# Gene IDs in this matrix are symbols. Retain only fixed target genes.
probe={gene:gene for gs in MODULES.values() for gene in gs};m,g,info=parse(str(p),probe)
if g.shape[1]>0:
 hc=m._source.eq('Healthy Control').values;sc,cov200=score(g,hc);q=pd.concat([m,sc],axis=1)
 q['subject']=q._title.str.extract(r'(?i)patient\s+(\d+)',expand=False);q['biopsy']=q._title.str.lower().map(lambda x:'First' if 'first' in x else 'Second' if 'second' in x else 'Control');q['compartment']=q._title.str.lower().map(lambda x:'Glomeruli' if 'glom' in x else 'Tubulointerstitium')
 q['proteinuria']=pd.to_numeric(q['proteinuria (g/24hr)'],errors='coerce');q['creatinine']=pd.to_numeric(q['serum creatinine (mg/dl)'],errors='coerce');q['age_num']=pd.to_numeric(q['age (yrs)'],errors='coerce')
 # Collapse any technical duplicates at donor/biopsy/compartment level before inference.
 adult=q.loc[q.subject.notna()&(q.age_num>=18)].copy();qq=adult.groupby(['subject','biopsy','compartment'],as_index=False).agg({**{k:'mean' for k in ['M_IFN','M_LDG','M_PB','ESI','proteinuria','creatinine']},'clinical response group':'first'})
 prs=[]
 for (comp,response),ss in qq.groupby(['compartment','clinical response group']):
  for k in ['M_IFN','M_LDG','M_PB']:
   u=ss.pivot(index='subject',columns='biopsy',values=k)
   if {'First','Second'}.issubset(u.columns):
    u=u[['First','Second']].dropna();r=pair_summary(u.First,u.Second);r.update({'compartment':comp,'response':response,'module':k});prs.append(r)
 result['partial_kidney']={'raw_n':len(q),'adult_donors':adult.subject.nunique(),'retained_donor_biopsy_compartment_rows':len(qq),'coverage':cov200,'label':'partial gene-panel sensitivity; does not validate full frozen modules','paired':prs,'age_min':adult.age_num.min()}
 qq.to_csv(ROOT/'results/task3_GSE200306_partial_source.csv',index=False)
for family in ['paired','FE']:
 qvals=multipletests([r['p'] for r in result[family]],method='fdr_bh')[1]
 for r,v in zip(result[family],qvals):r['BH_q']=v
save(result,'results/task3_results.json')
pd.DataFrame(result['paired']).to_csv(ROOT/'results/task3_adult_paired.csv',index=False)
print('TASK3 COMPLETE',result['adult_LN_blood'],result['partial_kidney'].get('coverage'),flush=True)
