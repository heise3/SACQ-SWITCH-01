"""Task 4: donor-level granulocyte subtype trajectories and renal compartment tests."""
from common import *
import h5py,itertools
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
assert (ROOT/'results/task3_results.json').exists(),'Complete Task 3 first'
rng=np.random.default_rng(SEED)

def hread(grp,k):
 q=grp[k]
 if isinstance(q,h5py.Group):
  cats=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in q['categories'][:]])
  codes=q['codes'][:];return np.where(codes>=0,cats[np.maximum(codes,0)],'NA')
 a=q[:];return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])

def extract(kind):
 fn=ROOT/f'work/sc_{kind}_aggregates.pkl'
 if fn.exists():return pd.read_pickle(fn)
 with h5py.File(CACHE/f'natmed2026_{kind}.h5ad') as f:
  obs=pd.DataFrame({k:hread(f['obs'],k) for k in ['donor_id','timepoint','cell_type','sample']});obs.index.name='cell'
  # All cells are used, no subsampling; retain all target genes and full-row library sizes.
  raw=f['raw/X'];vg=f['raw/var'];key='feature_name' if 'feature_name' in vg else 'gene_symbols';genes=hread(vg,key)
  targets={i:normalize(g) for i,g in enumerate(genes) if normalize(g)};idx=np.array(list(targets));names=[targets[i] for i in idx]
  def read_matrix(group,normalize_counts):
   shape=tuple(group.attrs['shape']);ptr=group['indptr'][:];out=np.zeros((shape[0],len(idx)),dtype=np.float32)
   for lo in range(0,shape[0],4096):
    hi=min(shape[0],lo+4096);a,b=int(ptr[lo]),int(ptr[hi]);z=sparse.csr_matrix((group['data'][a:b],group['indices'][a:b],ptr[lo:hi+1]-a),shape=(hi-lo,shape[1]));v=z[:,idx].toarray()
    if normalize_counts:
     lib=np.asarray(z.sum(axis=1)).ravel();assert (lib>0).all();v=np.log1p(v/lib[:,None]*1e4)
    out[lo:hi]=v
   e=pd.DataFrame(out,columns=names).T.groupby(level=0).mean().T
   return pd.concat([obs,e],axis=1)
  xgenes=hread(f['var'],'feature_name');assert np.array_equal(genes,xgenes), 'raw and X feature order mismatch'
  cell=read_matrix(raw,True);xx=read_matrix(f['X'],False)
  assert len(cell)==(77112 if kind=='ldg' else 4829)
  group=['donor_id','timepoint','cell_type'];gn=sorted(set(names));pb=cell.groupby(group)[gn].mean();pbX=xx.groupby(group)[gn].mean()
  cnt=obs.groupby(group).size().rename('cells');dt=cnt.groupby(level=[0,1]).sum().rename('compartment_cells');frac=cnt.to_frame().join(dt);frac['fraction']=frac.cells/frac.compartment_cells
  # Add absent subtypes as zero fractions only in observed donor/timepoint samples.
  grid=pd.MultiIndex.from_tuples([(d,t,c) for d,t in dt.index for c in sorted(obs.cell_type.unique())],names=group)
  frac=frac.reindex(grid);frac['compartment_cells']=frac['compartment_cells'].fillna(pd.Series([dt.loc[(d,t)] for d,t,c in grid],index=grid));frac['cells']=frac.cells.fillna(0).astype(int);frac['fraction']=frac.cells/frac.compartment_cells
  score_tables={}
  for label,expr in [('log1p_CPTT',pb),('original_X',pbX)]:
   baseline=expr.loc[expr.index.get_level_values('timepoint')=='Base'];mu=baseline.mean();sd=baseline.std(ddof=1).replace(0,np.nan);z=(expr-mu)/sd
   sc=pd.DataFrame(index=expr.index);coverage={}
   for m,gs in MODULES.items():
    avail=[g for g in gs if g in z and z[g].notna().any()];coverage[m]=avail;sc[m]=z[avail].mean(axis=1)
   sc['ESI']=.5*(sc.M_PB+sc.M_LDG);score_tables[label]=sc
  # Primary score uses common Base donor/subtype gene reference, not subtype-specific z-scores.
  out=frac.join(score_tables['log1p_CPTT']);out=out.join(score_tables['original_X'].add_suffix('_X'))
  out=out.reset_index();out.to_pickle(fn);out.to_csv(ROOT/f'results/task4_{kind}_donor_subtype_source.csv',index=False)
  pb.reset_index().to_csv(ROOT/f'results/task4_{kind}_pseudobulk_genes.csv',index=False)
  save({'cells':len(obs),'donors':obs.donor_id.nunique(),'sample_counts':obs.groupby(['donor_id','timepoint']).size().to_dict(),'coverage':coverage,'normalization':'raw per-cell counts / library size * 10000 then log1p, mean per donor/timepoint/subtype; common Base reference','citation':f['uns/citation'][()].decode()},f'results/task4_{kind}_processing.json')
  print('EXTRACTED',kind,len(obs),flush=True)
  return out
ld=extract('ldg');bc=extract('bcells');res={'paired':[],'baseline_subtype':[],'change_contrasts':[],'correlations':[],'kidney':{},'claim_limits':['cell proportions are compositional, not absolute depletion','uncontrolled pre/post therapy does not identify a treatment effect','frozen LDG granule program is not a functional NETosis assay','baseline subtype and mediator are contemporaneous where measured']}
for kind,tab in [('granulocytes',ld),('B_cells',bc)]:
 for ct,ss in tab.groupby('cell_type'):
  for measure in ['fraction','M_LDG','M_IFN','M_PB','M_LDG_X','M_IFN_X']:
   u=ss.pivot(index='donor_id',columns='timepoint',values=measure)
   for tp in ['D28','M3']:
    if {'Base',tp}.issubset(u.columns):
     pp=u[['Base',tp]].dropna();r=pair_summary(pp.Base,pp[tp]);r.update({'compartment':kind,'subtype':ct,'measure':measure,'timepoint':tp});res['paired'].append(r)
# Baseline donor-paired subtype contrasts across common-reference scores.
base=ld.loc[ld.timepoint=='Base']
for m in ['M_LDG','M_IFN','M_LDG_X','M_IFN_X']:
 u=base.pivot(index='donor_id',columns='cell_type',values=m)
 for other in ['mature neutrophil','neutrophil']:
  p=u[[other,'immature neutrophil']].dropna();r=pair_summary(p[other],p['immature neutrophil']);r.update({'measure':m,'reference':other,'comparison':'immature neutrophil'});res['baseline_subtype'].append(r)
for measure in ['fraction','M_LDG','M_IFN']:
 u=ld.pivot(index=['donor_id','cell_type'],columns='timepoint',values=measure)
 delta=(u.M3-u.Base).unstack('cell_type')
 for other in ['mature neutrophil','neutrophil']:
  p=delta[[other,'immature neutrophil']].dropna();r=pair_summary(p[other],p['immature neutrophil']);r.update({'measure':measure,'reference_change':other,'comparison_change':'immature neutrophil'});res['change_contrasts'].append(r)
# Cross-compartment proportions; absent B compartment timepoint stays missing.
u=ld.loc[ld.cell_type=='immature neutrophil'].pivot(index='donor_id',columns='timepoint',values='fraction');contraction=(u.Base-u.M3).rename('immature_contraction')
corr_source=[]
for ct in ['plasmablast','class switched memory B cell']:
 for tp in ['D28','M3']:
  a=bc.loc[(bc.cell_type==ct)&(bc.timepoint==tp)].set_index('donor_id').fraction.rename('B_fraction');q=pd.concat([a,contraction],axis=1).dropna();n=len(q)
  if n<4:continue
  x,y=q.B_fraction.values,q.immature_contraction.values
  if np.std(x)==0 or np.std(y)==0:
   res['correlations'].append({'B_subtype':ct,'timepoint':tp,'n_donors':n,'status':'NOT_ESTIMABLE_CONSTANT_FRACTION','spearman_rho':None,'spearman_permutation_p':None})
   for donor,row in q.iterrows():corr_source.append({'donor':donor,'B_subtype':ct,'timepoint':tp,**row.to_dict()})
   continue
  rho=stats.spearmanr(x,y).statistic;pear=stats.pearsonr(x,y);xr=stats.rankdata(x);yr=stats.rankdata(y);xr=(xr-xr.mean())/np.linalg.norm(xr-xr.mean());yr=(yr-yr.mean())/np.linalg.norm(yr-yr.mean())
  if n<=7:perm=np.array(list(itertools.permutations(range(n))));pperm=float(np.mean(np.abs(yr[perm]@xr)>=abs(rho)-1e-12));method='all donor permutations'
  else:
   perm=np.array([rng.permutation(n) for _ in range(9999)]);pperm=float((1+(np.abs(yr[perm]@xr)>=abs(rho)-1e-12).sum())/10000);method='9999 donor permutations + observed'
  ii=rng.integers(0,n,size=(5000,n));rx=stats.rankdata(x[ii],axis=1);ry=stats.rankdata(y[ii],axis=1);rx-=rx.mean(axis=1)[:,None];ry-=ry.mean(axis=1)[:,None];rr=(rx*ry).sum(axis=1)/np.sqrt((rx*rx).sum(axis=1)*(ry*ry).sum(axis=1));rr=rr[np.isfinite(rr)]
  r={'B_subtype':ct,'timepoint':tp,'n_donors':n,'spearman_rho':rho,'spearman_ci95':np.quantile(rr,[.025,.975]),'spearman_permutation_p':pperm,'permutation_method':method,'pearson_r':pear.statistic,'pearson_ci95':list(pear.confidence_interval()),'pearson_p':pear.pvalue};res['correlations'].append(r)
  for donor,row in q.iterrows():corr_source.append({'donor':donor,'B_subtype':ct,'timepoint':tp,**row.to_dict()})
pd.DataFrame(corr_source).to_csv(ROOT/'results/task4_B_granulocyte_correlations_source.csv',index=False)
# GSE32591 matched renal compartments: two observations per donor; compare LN-HC differences in donor contrast.
m,g,info=parse(str(CACHE/'GSE32591_series_matrix.txt.gz'),ns['load_gpl14663_map']());hc=m['disease status'].eq('control').values;sc,cov=score(g,hc);k=pd.concat([m,sc],axis=1);k['is_LN']=(~hc).astype(int);k['donor']=k['is_LN'].astype(str)+'_'+k['patient number'].astype(str);k['compartment']=np.where(k.tissue.str.startswith('Glomeruli'),'Glomeruli','Tubulointerstitium')
k.to_csv(ROOT/'results/task4_kidney_source.csv',index_label='gsm')
assert not k.duplicated(['donor','compartment']).any()
res['kidney']['clinical_labels_available']=False;res['kidney']['coverage']=cov;res['kidney']['compartment_results']=[];res['kidney']['interaction']=[]
for mod in ['M_IFN','M_LDG','M_PB','ESI','NC_HK']:
 for comp,ss in k.groupby('compartment'):
  a=ss.loc[ss.is_LN==1,mod].values;b=ss.loc[ss.is_LN==0,mod].values;delta=a.mean()-b.mean();se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b));df=(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))**2/((a.var(ddof=1)/len(a))**2/(len(a)-1)+(b.var(ddof=1)/len(b))**2/(len(b)-1));ci=delta+stats.t.ppf(.975,df)*np.array([-1,1])*se
  # Independent donor bootstrap Cohen d CI, not cell/visit bootstrap.
  ds=[]
  for _ in range(2000):ds.append(ns['cohens_d'](rng.choice(a,len(a)),rng.choice(b,len(b))))
  res['kidney']['compartment_results'].append({'module':mod,'compartment':comp,'n_LN':len(a),'n_controls':len(b),'difference':delta,'ci95':ci,'cohens_d':ns['cohens_d'](a,b),'cohens_d_ci95':np.quantile(ds,[.025,.975]),'p':stats.ttest_ind(a,b,equal_var=False).pvalue})
 u=k.pivot(index=['donor','is_LN'],columns='compartment',values=mod).dropna();u['delta']=u.Glomeruli-u.Tubulointerstitium
 a=u.loc[u.index.get_level_values('is_LN')==1,'delta'];b=u.loc[u.index.get_level_values('is_LN')==0,'delta'];fit=stats.ttest_ind(a,b,equal_var=False);se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b));df=(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))**2/((a.var(ddof=1)/len(a))**2/(len(a)-1)+(b.var(ddof=1)/len(b))**2/(len(b)-1));ci=(a.mean()-b.mean())+np.array([-1,1])*stats.t.ppf(.975,df)*se
 res['kidney']['interaction'].append({'module':mod,'n_LN_pairs':len(a),'n_HC_pairs':len(b),'difference_in_compartment_differences':a.mean()-b.mean(),'ci95':ci,'p':fit.pvalue})
for family in ['paired','baseline_subtype','change_contrasts']:
 arr=[r for r in res[family] if 'p' in r];q=multipletests([r['p'] for r in arr],method='fdr_bh')[1]
 for r,v in zip(arr,q):r['BH_q']=v
for family in ['compartment_results','interaction']:
 arr=res['kidney'][family];q=multipletests([r['p'] for r in arr],method='fdr_bh')[1]
 for r,v in zip(arr,q):r['BH_q']=v
estimable=[r for r in res['correlations'] if r.get('spearman_permutation_p') is not None]
q=multipletests([r['spearman_permutation_p'] for r in estimable],method='fdr_bh')[1] if estimable else []
for r,v in zip(estimable,q):r['BH_q']=v
save(res,'results/task4_results.json')
print('TASK4 COMPLETE',len(res['paired']),'donor tests',res['correlations'],flush=True)
