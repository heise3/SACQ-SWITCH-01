"""Task 2: adult trial quadrants and continuous serology gradients."""
from common import *
import re,subprocess
from collections import defaultdict
from scipy import stats
import statsmodels.api as sm,statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from patsy import dmatrix

def load88():
 fn=ROOT/'work/gse88884.pkl'
 if fn.exists():return pd.read_pickle(fn)
 # Literal, read-only Phase-0 expression loading block preserves score definition.
 text=original[original.index('meta_88884, _, _ ='):original.index('# --- Load GSE65391')]
 env=dict(ns);env.update({'gpl17586_map':ns['load_gpl17586_map'](),'CACHE_DIR':str(CACHE),'subprocess':subprocess})
 exec(text,env)
 d=env['df88'];d.to_pickle(fn)
 save({'n_samples':len(d),'coverage':env['cov_88884'],'q99':env['q99_88'],'sle':int(env['sle_mask_88884'].sum()),'healthy':int(env['hc_mask_88884'].sum())},'results/gse88884_processing.json')
 return d
if __name__=='__main__':
 assert (ROOT/'results/task1_results.json').exists(),'Task 1 must finish before Task 2'
 d=load88();s=d.loc[d.group=='SLE'].copy()
 for k,c in [('sledai_num','sledai_at_baseline'),('age_num','age_at_baseline'),('dsdna_iu','antidsdna_at_baseline (iu)'),('c3_val','c3_at_baseline (g/l)'),('c4_val','c4_at_baseline (g/l)')]:s[k]=pd.to_numeric(s[c],errors='coerce')
 adna=next(c for c in s.columns if c.startswith('adna_at_baseline'));c3=next(c for c in s.columns if c.startswith('c3base'));c4=next(c for c in s.columns if c.startswith('c4base'))
 dna=s[adna].str.lower().str.endswith('positive')|(s.dsdna_iu.fillna(0)>=30)
 low=s[c3].str.upper().isin(['L','LOW'])|s[c4].str.upper().isin(['L','LOW'])|(s.c3_val<.9)|(s.c4_val<.1)
 s['Sero_SLEDAI']=2*dna.astype(float)+2*low.astype(float);s['Est_cSLEDAI']=(s.sledai_num-s.Sero_SLEDAI).clip(lower=0)
 s['sero']=s.Sero_SLEDAI>0;s['quadrant']=np.select([s.sero&(s.Est_cSLEDAI<=4),s.sero&(s.Est_cSLEDAI>=6),~s.sero&(s.Est_cSLEDAI<=4),~s.sero&(s.Est_cSLEDAI>=6)],['Sero+ / low clinical','Sero+ / high clinical','Sero- / low clinical','Sero- / high clinical'],default='Intermediate clinical (5) / missing')
 s['serology_complete']=s[['dsdna_iu','c3_val','c4_val']].notna().all(axis=1)
 s.loc[~s.serology_complete,'quadrant']='Incomplete serology'
 s['exact_estimated_SACQ']=(s.Est_cSLEDAI==0)&s.sero
 s.to_csv(ROOT/'results/task2_trial_source.csv',index_label='gsm')
 result={'trials':{},'quadrants':[],'slopes':[],'phase0_reproduction':[],'contrasts':[]}
 for trial,t in s.groupby('batch'):
  sub=t.dropna(subset=['sledai_num','age_num','M_IFN','M_PB','M_LDG']).copy()
  result['trials'][trial]={'raw':len(t),'patients':t.subject_id.nunique(),'complete_model':len(sub),'missing_age':int(t.age_num.isna().sum()),'estimated_SACQ':int(t.exact_estimated_SACQ.sum()),'quadrants':t.quadrant.value_counts().to_dict(),'missing_serology':t[['dsdna_iu','c3_val','c4_val']].isna().sum().to_dict()}
  for outcome in ['Sero_SLEDAI','Est_cSLEDAI','sledai_num']:
   fit=smf.ols(f'{outcome} ~ M_IFN + M_LDG + M_PB + age_num + C(Sex) + C(race)',data=sub).fit(cov_type='HC3')
   for k in ['M_IFN','M_PB','M_LDG']:
    result['phase0_reproduction'].append({'trial':trial,'outcome':outcome,'module':k,'beta':fit.params[k],'se':fit.bse[k],'ci95':fit.conf_int().loc[k].tolist(),'p':fit.pvalues[k],'n':len(sub)})
  for q,ss in t.groupby('quadrant'):
   for m in ['M_IFN','M_PB','M_LDG']:
    y=ss[m].dropna();se=stats.sem(y);ci=stats.t.interval(.95,len(y)-1,loc=y.mean(),scale=se) if len(y)>1 else [None,None]
    result['quadrants'].append({'trial':trial,'quadrant':q,'module':m,'n':len(y),'mean':y.mean(),'sd':y.std(),'ci95':ci})
  # Fixed reference Sero+ low; robust adjusted subgroup effects with CI.
  ss=sub.loc[sub.quadrant.isin(['Sero+ / low clinical','Sero+ / high clinical','Sero- / low clinical','Sero- / high clinical'])].copy()
  ss['quadrant']=pd.Categorical(ss.quadrant,categories=['Sero+ / low clinical','Sero+ / high clinical','Sero- / low clinical','Sero- / high clinical'])
  for m in ['M_IFN','M_PB','M_LDG']:
   fit=smf.ols(f'{m} ~ C(quadrant) + age_num + C(Sex) + C(race)',data=ss).fit(cov_type='HC3')
   for k in fit.params.index:
    if k.startswith('C(quadrant)'):result['contrasts'].append({'trial':trial,'module':m,'contrast':k,'beta':fit.params[k],'ci95':fit.conf_int().loc[k].tolist(),'p':fit.pvalues[k],'n':len(ss)})
  for marker in ['logDNA','negative_C3','negative_C4']:
   v=np.log1p(sub.dsdna_iu) if marker=='logDNA' else -sub.c3_val if marker=='negative_C3' else -sub.c4_val
   ss=sub.copy();ss['x']=(v-v.mean())/v.std(ddof=1);ss=ss.dropna(subset=['x']);B=dmatrix('bs(x,df=4,degree=3,include_intercept=False)-1',ss,return_type='dataframe').values
   A=np.column_stack([np.ones(len(ss)),ss.x]);res=B-A@np.linalg.lstsq(A,B,rcond=None)[0];U,S,V=np.linalg.svd(res,full_matrices=False);rank=int((S>1e-8).sum())
   for j in range(rank):ss[f'nl{j}']=U[:,j]*np.sqrt(len(ss))
   for m in ['M_IFN','M_PB','M_LDG']:
    ss['y']=(ss[m]-ss[m].mean())/ss[m].std(ddof=1)
    lin=smf.ols('y ~ x + age_num + C(Sex) + C(race)',data=ss).fit(cov_type='HC3')
    non=smf.ols('y ~ x + '+' + '.join(f'nl{j}' for j in range(rank))+' + age_num + C(Sex) + C(race)',data=ss).fit(cov_type='HC3')
    R=np.zeros((rank,len(non.params)))
    for j in range(rank):R[j,list(non.params.index).index(f'nl{j}')]=1
    pn=float(non.wald_test(R,scalar=True).pvalue)
    result['slopes'].append({'trial':trial,'marker':marker,'module':m,'n':len(ss),'beta_SD_per_SD':lin.params.x,'ci95':lin.conf_int().loc['x'].tolist(),'p':lin.pvalues.x,'nonlinear_p':pn,'nonlinear_df':rank})
 for family in ['slopes','contrasts','phase0_reproduction']:
  arr=result[family];q=multipletests([r['p'] for r in arr],method='fdr_bh')[1]
  for r,v in zip(arr,q):r['BH_q']=v
 if result['slopes']:
  q=multipletests([r['nonlinear_p'] for r in result['slopes']],method='fdr_bh')[1]
  for r,v in zip(result['slopes'],q):r['nonlinear_BH_q']=v
 save(result,'results/task2_results.json')
 for k in ['quadrants','slopes','contrasts','phase0_reproduction']:pd.DataFrame(result[k]).to_csv(ROOT/f'results/task2_{k}.csv',index=False)
 print('TASK2 COMPLETE',result['trials'],flush=True)
