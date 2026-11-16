"""Task 1: patient-cluster mediation and landmark/discrete survival."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from common import *
from scipy import stats,linalg,special
import statsmodels.api as sm
from lifelines import CoxPHFitter,KaplanMeierFitter
from lifelines.statistics import logrank_test,proportional_hazard_test
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore',category=RuntimeWarning)
rng=np.random.default_rng(SEED)
d=load65();s=d.loc[d['disease state']=='SLE'].sort_values(['subject','visit_num']).copy()
records=[];intervals=[];landmarks=[];timecheck=[]
for sid,g in s.groupby('subject',sort=True):
 rows=g.to_dict('records')
 for j in range(len(rows)-1):
  a,b=rows[j:j+2]
  if a['state']!='SACQ':continue
  r={'subject':sid,'gsm_t':g.index[j],'visit_t':a['visit_num'],'visit_t1':b['visit_num'],'days':b['days_since_last_visit_num'],'time_diff':b['cumulative_time_num']-a['cumulative_time_num'],'any':int(b['cSLEDAI']>=2),'renal':int(b['renal_num']>0),'cSLEDAI_next':b['cSLEDAI']}
  for k,c in [('PB','M_PB'),('IFN','M_IFN'),('LDG','M_LDG'),('ESI','ESI'),('C3','c3_num'),('C4','c4_num'),('DNA','ds_dna_num'),('age','age_num'),('oral','oral_steroids_category_num'),('iv','steroid_iv_category_num'),('mmf','mycophenolate_category_num')]:r[k]=a[c]
  r['male']=int(a['gender']=='M');records.append(r)
 eligible=[i for i,a in enumerate(rows) if a['state']=='SACQ' and i<len(rows)-1]
 if not eligible:continue
 j=eligible[0];a=rows[j];later=rows[j+1:];base=a['cumulative_time_num']
 r={'subject':sid,'ESI':a['ESI'],'IFN':a['M_IFN'],'PB':a['M_PB'],'LDG':a['M_LDG'],'age':a['age_num'],'male':int(a['gender']=='M'),'oral':a['oral_steroids_category_num'],'iv':a['steroid_iv_category_num'],'mmf':a['mycophenolate_category_num'],'landmark_time':base,'C3':a['c3_num'],'C4':a['c4_num'],'DNA':a['ds_dna_num']}
 for outcome in ['any','renal']:
  found=next((b for b in later if b['cSLEDAI']>=2 if outcome=='any'),None) if outcome=='any' else next((b for b in later if b['renal_num']>0),None)
  end=found or later[-1];r[outcome+'_event']=int(found is not None);r[outcome+'_time']=end['cumulative_time_num']-base
 landmarks.append(r)
 for outcome in ['any','renal']:
  for k in range(j,len(rows)-1):
   b,c=rows[k:k+2];dt=c['cumulative_time_num']-b['cumulative_time_num'];ev=int(c['cSLEDAI']>=2) if outcome=='any' else int(c['renal_num']>0)
   ir={'subject':sid,'outcome':outcome,'days':dt,'event':ev,'start':b['cumulative_time_num']-base,'ESI':b['ESI'],'IFN':b['M_IFN'],'C3':b['c3_num'],'C4':b['c4_num'],'DNA':b['ds_dna_num'],'age':b['age_num'],'male':int(b['gender']=='M'),'oral':b['oral_steroids_category_num'],'iv':b['steroid_iv_category_num'],'mmf':b['mycophenolate_category_num']}
   intervals.append(ir)
   if ev:break
tr=pd.DataFrame(records);lm=pd.DataFrame(landmarks);it=pd.DataFrame(intervals)
assert len(tr)==274 and tr.subject.nunique()==78,(len(tr),tr.subject.nunique())
time_mismatch=int((~np.isclose(tr.days,tr.time_diff,equal_nan=True)).sum())
tr['days_since_last_visit']=tr.days
tr['days']=tr.time_diff
assert (tr.days>0).all()
assert (it.days>0).all()
tr.to_csv(ROOT/'results/task1_transitions.csv',index=False);lm.to_csv(ROOT/'results/task1_landmarks.csv',index=False);it.to_csv(ROOT/'results/task1_intervals.csv',index=False)
print('TASK1 transitions',len(tr),'patients',tr.subject.nunique(),'events',tr['any'].sum(),flush=True)
result={'transitions':{'n':len(tr),'patients':tr.subject.nunique(),'any_events':tr['any'].sum(),'renal_events':tr.renal.sum(),'metadata_interval_mismatch_n':time_mismatch,'elapsed_time_source':'difference in cumulative_time; original days_since_last_visit retained'},'mediation':[],'cox':[],'discrete':[],'next_visit_models':[],'landmark':{},'missingness':[]}
# Gaussian linear mediator, natural effects integrated by 20-point Gauss-Hermite quadrature.
# Cluster bootstrap samples patient labels with replacement; patient multiplicities weight all their rows.
conf=['age','male','oral','iv','mmf'];vars=['PB','IFN','LDG','C3','days']+conf
med=tr.dropna(subset=vars).copy();loss=1-len(med)/len(tr)
result['missingness'].append({'analysis':'mediation','eligible':len(tr),'complete':len(med),'loss':loss,'patients':med.subject.nunique()})
if loss>.15:raise RuntimeError('Mediation loss >15%; multiple-imputation branch required')
for col in ['PB','IFN','LDG','C3','age','days']:
 med[col]=(med[col]-med[col].mean())/med[col].std(ddof=1)
med['lowC3']=-med.C3
C=np.column_stack([np.ones(len(med)),med[conf].values,med.days.values]);yM=med.lowC3.values
ids=np.unique(med.subject);grp=pd.Categorical(med.subject,categories=ids).codes
nodes,weights=np.polynomial.hermite.hermgauss(20);weights=weights/np.sqrt(np.pi)
def fit_lm(X,y,w):return np.linalg.lstsq(X*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=None)[0]
def fit_bin(X,y,w):
 b=np.zeros(X.shape[1])
 for _ in range(60):
  eta=X@b;p=special.expit(eta);ww=np.maximum(p*(1-p),1e-9)*w
  H=X.T@(ww[:,None]*X);grad=X.T@(w*(y-p));step=np.linalg.lstsq(H,grad,rcond=None)[0]
  b+=step
  if np.max(np.abs(step))<1e-7:break
 if not np.all(np.isfinite(b)) or np.max(np.abs(b))>40:raise ValueError('unstable logit')
 return b
# Return order: ACME0, ACME1, mean ACME, ADE0, ADE1, mean ADE, total, mediated fraction.
def decompose(exposure,outcome,w):
 x=med[exposure].values;other=med['IFN' if exposure=='PB' else 'PB'].values
 X=np.column_stack([C,other,x]);am=fit_lm(X,yM,w);M0=np.column_stack([C,other,np.zeros(len(x))])@am;M1=M0+am[-1]
 Yo=med[outcome].values;Xy=np.column_stack([C,other,x,yM]);by=fit_lm(Xy,Yo,w) if outcome=='LDG' else fit_bin(Xy,Yo,w)
 if outcome=='LDG':
  ac=am[-1]*by[-1];ade=by[-2];total=ac+ade
  return np.array([ac,ac,ac,ade,ade,ade,total,ac/total if abs(total)>1e-9 else np.nan])
 sigma=np.sqrt(np.average((yM-X@am)**2,weights=w));fixed=np.column_stack([C,other])@by[:-2]
 def risk(x,m):
  p=special.expit(fixed[:,None]+by[-2]*x+by[-1]*(m[:,None]+np.sqrt(2)*sigma*nodes[None,:]))@weights
  return np.average(p,weights=w)
 r00,r01,r10,r11=risk(0,M0),risk(0,M1),risk(1,M0),risk(1,M1)
 ac0,ac1=r01-r00,r11-r10;ad0,ad1=r10-r00,r11-r01;total=r11-r00;ac=(ac0+ac1)/2
 return np.array([ac0,ac1,ac,ad0,ad1,(ad0+ad1)/2,total,ac/total if abs(total)>1e-9 else np.nan])
keys=['ACME0','ACME1','ACME','ADE0','ADE1','ADE','total','proportion']
bootweights=[np.bincount(rng.integers(len(ids),size=len(ids)),minlength=len(ids))[grp].astype(float) for _ in range(2000)]
if (ROOT/'results/task1_checkpoint.json').exists():
 result=json.loads((ROOT/'results/task1_checkpoint.json').read_text())
else:
 for x in ['PB','IFN']:
  for y in ['LDG','any','renal']:
   point=decompose(x,y,np.ones(len(med)));boots=[];failed=0
   for w in bootweights:
    try:boots.append(decompose(x,y,w))
    except (ValueError,np.linalg.LinAlgError):failed+=1
   B=np.array(boots);r={'exposure':x,'outcome':y,'n':len(med),'patients':len(ids),'bootstrap_requested':2000,'bootstrap_success':len(B),'bootstrap_failed':failed,'causal_identification':'NOT_IDENTIFIED: concurrent exposure/mediator and observational assumptions','unit':'SD LDG' if y=='LDG' else 'risk difference'}
   for i,k in enumerate(keys):
    v=B[:,i];v=v[np.isfinite(v)];r[k]={'estimate':point[i],'ci95':np.quantile(v,[.025,.975]),'bootstrap_sign_p':min(1,2*min((1+(v<=0).sum())/(len(v)+1),(1+(v>=0).sum())/(len(v)+1)))}
   r['proportion_interpretable']=bool(r['total']['ci95'][0]*r['total']['ci95'][1]>0)
   result['mediation'].append(r);pd.DataFrame(B,columns=keys).to_csv(ROOT/f'results/task1_boot_{x}_{y}.csv',index=False)
   print('mediation',x,y,'ACME',r['ACME'],'success',len(B),flush=True)
 save(result,'results/task1_checkpoint.json')
# Reproduce Phase-0 univariate prospective association and add full adjustment.
for c in ['ESI','IFN','PB','LDG','C3','C4','age']:
 tr[c+'z']=(tr[c]-tr[c].mean())/tr[c].std(ddof=1)
tr['logDNA']=np.log1p(tr.DNA);tr['logdays']=np.log(tr.days)
for outcome in ['any','renal']:
 for predictor in ['ESI','IFN']:
  for model,covs in [('unadjusted',[predictor+'z']),('demography_treatment_interval',[predictor+'z','agez','male','oral','iv','mmf','logdays']),('serology_adjusted',[predictor+'z']+(['IFNz'] if predictor=='ESI' else [])+['C3z','C4z','logDNA','agez','male','oral','iv','mmf','logdays'])]:
   sub=tr.dropna(subset=covs+[outcome]);loss=1-len(sub)/len(tr)
   fit=sm.GLM(sub[outcome],sm.add_constant(sub[covs]),family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':sub.subject})
   k=predictor+'z';ci=fit.conf_int().loc[k]
   result['next_visit_models'].append({'outcome':outcome,'predictor':predictor,'model':model,'n':len(sub),'patients':sub.subject.nunique(),'OR':np.exp(fit.params[k]),'ci95':np.exp(ci).tolist(),'p':fit.pvalues[k]})
   result['missingness'].append({'analysis':'next_visit_'+outcome+'_'+model,'eligible':len(tr),'complete':len(sub),'loss':loss})
   if loss>.15:print('MI sensitivity required',outcome,model,loss,flush=True)
# One record per patient for landmark Cox and KM.
for c in ['ESI','IFN','PB','LDG','age']:
 lm[c+'z']=(lm[c]-lm[c].mean())/lm[c].std(ddof=1)
median=lm.ESI.median();lm['high_ESI']=(lm.ESI>median).astype(int)
for outcome in ['any','renal']:
 valid=lm.loc[np.isfinite(lm[outcome+'_time'])&(lm[outcome+'_time']>0)].copy()
 result['landmark'][outcome]={'eligible':len(lm),'retained':len(valid),'events':int(valid[outcome+'_event'].sum()),'median_followup_days':valid[outcome+'_time'].median(),'median_ESI_cutoff':median,'n_low':int((valid.high_ESI==0).sum()),'n_high':int((valid.high_ESI==1).sum())}
 lr=logrank_test(valid.loc[valid.high_ESI==0,outcome+'_time'],valid.loc[valid.high_ESI==1,outcome+'_time'],valid.loc[valid.high_ESI==0,outcome+'_event'],valid.loc[valid.high_ESI==1,outcome+'_event'])
 result['landmark'][outcome]['logrank_p']=lr.p_value
 for predictor in ['ESI','IFN']:
  for model,covs in [('unadjusted',[predictor+'z']),('adjusted',[predictor+'z','agez','male','oral','iv','mmf'])]:
   sub=valid[['subject',outcome+'_time',outcome+'_event']+covs].dropna().copy();loss=1-len(sub)/len(lm)
   result['missingness'].append({'analysis':'Cox_'+outcome+'_'+model,'eligible':len(lm),'complete':len(sub),'loss':loss})
   if loss>.15:raise RuntimeError('Cox loss >15%; MI required')
   # Patient is both independent record and robust cluster.
   cph=CoxPHFitter().fit(sub,duration_col=outcome+'_time',event_col=outcome+'_event',cluster_col='subject',robust=True)
   row=cph.summary.loc[predictor+'z'];ph=proportional_hazard_test(cph,sub,time_transform='rank')
   result['cox'].append({'outcome':outcome,'predictor':predictor,'model':model,'n':len(sub),'events':int(sub[outcome+'_event'].sum()),'HR':row['exp(coef)'],'ci95':[row['exp(coef) lower 95%'],row['exp(coef) upper 95%']],'p':row['p'],'PH_test_p':ph.summary.loc[predictor+'z','p']})
# Discrete-time piecewise exponential approximation: cloglog plus log(interval length).
for outcome in ['any','renal']:
 sub=it.loc[(it.outcome==outcome)&(it.days>0)].copy()
 for c in ['ESI','IFN','age','C3','C4']:
  sub[c+'z']=(sub[c]-sub[c].mean())/sub[c].std(ddof=1)
 sub['logDNA']=np.log1p(sub.DNA);sub['logstart']=np.log1p(sub.start/90)
 for model,covs in [('adjusted',['ESIz','IFNz','agez','male','oral','iv','mmf','logstart']),('serology_adjusted',['ESIz','IFNz','agez','male','oral','iv','mmf','logstart','C3z','C4z','logDNA'])]:
  ss=sub.dropna(subset=covs);loss=1-len(ss)/len(sub)
  result['missingness'].append({'analysis':'discrete_'+outcome+'_'+model,'eligible':len(sub),'complete':len(ss),'loss':loss})
  if loss>.15:print('MI sensitivity required discrete',outcome,model,loss,flush=True)
  fit=sm.GLM(ss.event,sm.add_constant(ss[covs]),family=sm.families.Binomial(sm.families.links.CLogLog()),offset=np.log(ss.days)).fit(cov_type='cluster',cov_kwds={'groups':ss.subject})
  for k in ['ESIz','IFNz']:
   result['discrete'].append({'outcome':outcome,'predictor':k,'model':model,'n_intervals':len(ss),'patients':ss.subject.nunique(),'events':int(ss.event.sum()),'rate_ratio':np.exp(fit.params[k]),'ci95':np.exp(fit.conf_int().loc[k]).tolist(),'p':fit.pvalues[k]})
# Missing covariate sensitivity: 20 stochastic chained-equation imputations.
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
result['multiple_imputation']=[]
def pooled(estimates,variances):
 Q=np.asarray(estimates);U=np.asarray(variances);m=len(Q);q=Q.mean();u=U.mean();b=Q.var(ddof=1);v=u+(1+1/m)*b
 df=(m-1)*(1+u/((1+1/m)*b))**2 if b>0 else np.inf
 ci=q+np.array([-1,1])*stats.t.ppf(.975,df)*np.sqrt(v)
 return {'estimate':np.exp(q),'ci95':np.exp(ci).tolist(),'p':2*stats.t.sf(abs(q)/np.sqrt(v),df),'df_Rubin':df,'m':m}
# Include outcomes and cluster-level means as auxiliary predictors; never impute outcomes.
mi_columns=['ESIz','IFNz','C3z','C4z','logDNA','agez','male','oral','iv','mmf','logdays']
aux=tr[mi_columns+['any','renal']].copy()
for c in ['C3z','C4z','logDNA']:aux['patient_mean_'+c]=tr.groupby('subject')[c].transform('mean')
mi_arrays=[]
for j in range(20):
 imp=IterativeImputer(random_state=SEED+j,sample_posterior=True,max_iter=20)
 mi_arrays.append(imp.fit_transform(aux))
for outcome in ['any','renal']:
 for predictor in ['ESI','IFN']:
  covs=[predictor+'z']+(['IFNz'] if predictor=='ESI' else [])+['C3z','C4z','logDNA','agez','male','oral','iv','mmf','logdays']
  qs=[];us=[]
  for ar in mi_arrays:
   ss=pd.DataFrame(ar[:,:len(mi_columns)],columns=mi_columns,index=tr.index)
   fit=sm.GLM(tr[outcome],sm.add_constant(ss[covs]),family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':tr.subject})
   qs.append(fit.params[predictor+'z']);us.append(fit.bse[predictor+'z']**2)
  result['multiple_imputation'].append({'analysis':'next_visit','outcome':outcome,'predictor':predictor,'n':len(tr),'patients':tr.subject.nunique(),'scale':'OR',**pooled(qs,us)})
for outcome in ['any','renal']:
 ss=it.loc[(it.outcome==outcome)&(it.days>0)].copy()
 for c in ['ESI','IFN','age','C3','C4']:ss[c+'z']=(ss[c]-ss[c].mean())/ss[c].std(ddof=1)
 ss['logDNA']=np.log1p(ss.DNA);ss['logstart']=np.log1p(ss.start/90);ss['logdays']=np.log(ss.days)
 cols=['ESIz','IFNz','agez','male','oral','iv','mmf','logstart','C3z','C4z','logDNA']
 aux=ss[cols+['event','logdays']].copy()
 for c in ['C3z','C4z','logDNA']:aux['patient_mean_'+c]=ss.groupby('subject')[c].transform('mean')
 out={k:([],[]) for k in ['ESIz','IFNz']}
 for j in range(20):
  ar=IterativeImputer(random_state=SEED+j,sample_posterior=True,max_iter=20).fit_transform(aux)
  xx=pd.DataFrame(ar[:,:len(cols)],columns=cols,index=ss.index)
  fit=sm.GLM(ss.event,sm.add_constant(xx),family=sm.families.Binomial(sm.families.links.CLogLog()),offset=ss.logdays).fit(cov_type='cluster',cov_kwds={'groups':ss.subject})
  for k in out:out[k][0].append(fit.params[k]);out[k][1].append(fit.bse[k]**2)
 for k,(qs,us) in out.items():result['multiple_imputation'].append({'analysis':'discrete','outcome':outcome,'predictor':k,'n':len(ss),'patients':ss.subject.nunique(),'scale':'rate ratio',**pooled(qs,us)})
result['consecutive_visit_sensitivity']=[]
ss=tr.loc[(tr.visit_t1-tr.visit_t)==1].copy()
for outcome in ['any','renal']:
 for predictor in ['ESIz','IFNz']:
  covs=[predictor,'agez','male','oral','iv','mmf','logdays'];q=ss.dropna(subset=covs)
  fit=sm.GLM(q[outcome],sm.add_constant(q[covs]),family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':q.subject})
  result['consecutive_visit_sensitivity'].append({'outcome':outcome,'predictor':predictor,'n':len(q),'OR':np.exp(fit.params[predictor]),'ci95':np.exp(fit.conf_int().loc[predictor]).tolist(),'p':fit.pvalues[predictor]})


lm.to_csv(ROOT/'results/task1_landmarks.csv',index=False)
# Paired Phase-0 means with inferential CIs for writing.
p=s.groupby(['subject','state'])[['M_IFN','M_LDG','M_PB','ESI','NC_HK']].mean().unstack('state');paired={}
for m in ['M_IFN','M_LDG','M_PB','ESI','NC_HK']:
 q=p[m][['SACQ','ClinActive']].dropna();paired[m]=pair_summary(q.SACQ,q.ClinActive)
save(paired,'results/task1_state_paired.json')
save(result,'results/task1_results.json')
pd.DataFrame(result['cox']).to_csv(ROOT/'results/task1_cox.csv',index=False);pd.DataFrame(result['discrete']).to_csv(ROOT/'results/task1_discrete.csv',index=False);pd.DataFrame(result['next_visit_models']).to_csv(ROOT/'results/task1_next_visit.csv',index=False)
print('TASK1 COMPLETE',json.dumps(result['landmark']),flush=True)
