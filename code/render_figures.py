# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# Fig1b-e, Fig2c-d, Fig4a/d/e, S2c-d: Forest via cross-type Bar -> parameter inherit.
# Fig1a: original vector schematic -> vector primitives at the user-specified panel map.
# Fig2a-b, Fig3: PairedBoxScatter/LineTrend -> parameter inherit, all observed donor/patient points.
# S1: GroupedViolin -> parameter inherit; all four-quadrant observations and box summaries.
# S2a-b: scatter -> parameter inherit; all ten donor points with donor labels.
# Fig4b-c: frozen vector KM geometry copied from the original SVG; NO estimator is refitted.
# Existing templates embed synthetic or incompatible data; native PNG composition would violate editable SVG.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})
# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch,FancyBboxPatch
from matplotlib.lines import Line2D
from pathlib import Path
import numpy as np,pandas as pd,json,re,hashlib,xml.etree.ElementTree as ET
from decimal import Decimal, ROUND_HALF_UP
from scipy import stats
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'rendered_figures';OUT.mkdir(exist_ok=True)
SUP=ROOT/'documentation';SUP.mkdir(exist_ok=True)
R={i:json.loads((ROOT/f'results/task{i}_results.json').read_text()) for i in [1,2,3,4]}
P=json.loads((ROOT/'results/task1_state_paired.json').read_text())
COL={'M_IFN':'#2166AC','IFN':'#2166AC','M_PB':'#762A83','PB':'#762A83','M_LDG':'#B2182B','LDG':'#B2182B','ESI':'#1B7837','NC_HK':'#777777','SLEDAI':'#222222'}
LAB={'M_IFN':'IFN','M_PB':'PB','M_LDG':'LDG','ESI':'ESI','NC_HK':'HK'}
CTS=['immature neutrophil','neutrophil','mature neutrophil'];CTLAB=['Immature','Intermediate','Mature'];CTCOL=['#B2182B','#2166AC','#1B7837'];CTMARK=['o','^','s']
W=170/25.4;QA=[];DATA={};SOURCE={}
ld=pd.read_csv(ROOT/'results/task4_ldg_donor_subtype_source.csv');ad=pd.read_csv(ROOT/'results/task3_GSE72747_source.csv')
def pstr(x):
 d=Decimal(str(float(x)))
 if not d:return '0'
 q=d.quantize(Decimal(1).scaleb(d.adjusted()-2),rounding=ROUND_HALF_UP)
 return (format(q,'.2e') if d<Decimal('.001') else format(q,'f')).replace('e-0','e−').replace('e-','e−')
def num(x,n=3):
 x=float(x)
 return (f'{x:.2g}' if 0<abs(x)<.0005 else f'{x:.{n}f}').replace('e-0','e−').replace('-','−')
def val(r,est='delta',ci='ci95',p='p',label='',module=None,scale=1,q='BH_q'):
 return {'label':label,'est':float(r[est])*scale,'ci':[float(v)*scale for v in r[ci]],'p':float(r[p]),'q':r.get(q),'color':COL.get(module,'#34495E')}
def addtitle(fig,rect,letter,title,dy=.018):
 x,y,w,h=rect;fig.text(x-.025,y+h+dy,letter,fontsize=10,fontweight='bold',va='bottom');fig.text(x+.014,y+h+dy,title,fontsize=8,fontweight='bold',va='bottom')
def forest(fig,rect,letter,title,rows,xlabel,zero=0,log=False,showq=False,header='Estimate [95% CI]',label_frac=.29,table_frac=.41,fontsize=6.7,dy=.018,add_title=True,left_header=''):
 x,y,w,h=rect;fw=1-label_frac-table_frac
 l=fig.add_axes([x,y,w*label_frac,h]);ax=fig.add_axes([x+w*label_frac,y,w*fw,h]);t=fig.add_axes([x+w*(1-table_frac)+.006,y,w*table_frac-.006,h]);n=len(rows)
 for a in [l,ax,t]:a.set_ylim(n-.45,-1.05)
 l.set_axis_off();t.set_axis_off();
 if left_header:l.text(.01,-.72,left_header,fontsize=fontsize,fontweight='bold',transform=l.get_yaxis_transform())
 ax.spines[['left','right','top']].set_visible(False);ax.set_yticks([])
 lo=min([r['ci'][0] for r in rows]+[zero]);hi=max([r['ci'][1] for r in rows]+[zero]);pad=(hi-lo)*.12
 if log:
  ax.set_xscale('log');ax.set_xlim(max(.05,lo/1.15),hi*1.15)
 else:ax.set_xlim(lo-pad,hi+pad)
 for i,r in enumerate(rows):
  if i%2==0:
   for a in [l,ax,t]:a.axhspan(i-.42,i+.42,color='#F4F6F8',zorder=-5,lw=0)
  l.text(.01,i,r['label'],va='center',fontsize=fontsize,transform=l.get_yaxis_transform())
  ax.plot(r['ci'],[i,i],color=r['color'],lw=1.05,zorder=2);ax.plot(r['est'],i,'o',ms=3.5,color=r['color'],zorder=3)
  t.text(.60 if showq else .73,i,f"{num(r['est'])} [{num(r['ci'][0])}, {num(r['ci'][1])}]",ha='right',va='center',fontsize=fontsize,transform=t.get_yaxis_transform())
  t.text(.81 if showq else .99,i,pstr(r['p']),ha='right',va='center',fontsize=fontsize,transform=t.get_yaxis_transform())
  if showq:t.text(.99,i,pstr(r['q']) if r['q'] is not None else '—',ha='right',va='center',fontsize=fontsize,transform=t.get_yaxis_transform())
 t.text(.60 if showq else .73,-.72,header,ha='right',fontsize=fontsize,fontweight='bold',transform=t.get_yaxis_transform())
 t.text(.81 if showq else .99,-.72,'p',ha='right',fontsize=fontsize,fontweight='bold',transform=t.get_yaxis_transform())
 if showq:t.text(.99,-.72,'q',ha='right',fontsize=fontsize,fontweight='bold',transform=t.get_yaxis_transform())
 ax.axvline(zero,color='#999999',ls=(0,(3,3)),lw=.65,zorder=1);ax.set_xlabel(xlabel,fontsize=7,labelpad=2);ax.tick_params(axis='x',labelsize=6.5,pad=2)
 if log:
  ax.xaxis.set_minor_locator(mpl.ticker.NullLocator());ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
 if not log:ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
 else:ax.set_xticks([v for v in [.5,1,1.5,2,3] if ax.get_xlim()[0]<=v<=ax.get_xlim()[1]],labels=[str(v) for v in [.5,1,1.5,2,3] if ax.get_xlim()[0]<=v<=ax.get_xlim()[1]])
 if add_title:addtitle(fig,rect,letter,title,dy=dy)
 if letter:DATA[CURRENT+letter]=rows
 return ax

def save(fig,name):
 fig.canvas.draw();ren=fig.canvas.get_renderer();outside=[]
 for t in fig.findobj(mpl.text.Text):
  if not t.get_visible() or not t.get_text():continue
  # Major ticks outside the actual view are not drawn by Matplotlib.
  bb=t.get_window_extent(ren)
  if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.x1+.5 or bb.y1>fig.bbox.y1+.5:outside.append({'text':t.get_text(),'bounds':[bb.x0,bb.y0,bb.x1,bb.y1]})
 with mpl.rc_context({'savefig.bbox':None}):
  for ext in ['svg','pdf','png']:fig.savefig(OUT/f'{name}.{ext}',dpi=300,bbox_inches=None,facecolor='white')
 QA.append({'figure':name,'width_mm':fig.get_figwidth()*25.4,'height_mm':fig.get_figheight()*25.4,'outside_canvas_text':outside})
 plt.close(fig);print('RENDERED',name,flush=True)
def trajectory(ax,frame,idcol,timecol,measure,times,color,marker='o',scale=1,alpha=.32):
 u=frame.pivot(index=idcol,columns=timecol,values=measure).reindex(columns=times)*scale
 for _,row in u.iterrows():ax.plot(range(len(times)),row.values,color=color,alpha=alpha,lw=.5,marker=marker,ms=2.0,zorder=1)
 means=u.mean();n=u.notna().sum();ci=stats.t.ppf(.975,np.maximum(n-1,1))*u.sem()
 ax.errorbar(range(len(times)),means,yerr=ci,color=color,marker=marker,ms=3.5,lw=1.2,capsize=2,zorder=4)
 ax.set_xticks(range(len(times)),[str(v) for v in times]);ax.set_xlim(-.22,len(times)-.78);ax.tick_params(labelsize=6.5,pad=2)
 return {'means':means.tolist(),'n':n.tolist(),'pointwise_t_CI_halfwidth':ci.tolist(),'observed_records':int(u.notna().sum().sum()),'n_unique_units':len(u)}
def find(rows,**criteria):return next(r for r in rows if all(r.get(k)==v for k,v in criteria.items()))

# FIGURE 1. All evidence belongs to the first Results section.
CURRENT='Figure1';fig=plt.figure(figsize=(W,190/25.4))
ax=fig.add_axes([.04,.825,.935,.145]);ax.set_axis_off();ax.text(-.026,1.03,'a',fontweight='bold',fontsize=10,transform=ax.transAxes);ax.text(.015,1.04,'Cross-scale study design | eight public resources',fontsize=8.5,fontweight='bold',transform=ax.transAxes)
cards=[('GSE65391','157 patients · 899 visits'),('GSE88884','1,760 SLE · complete n=876/880'),('GSE49454','46 patients · 157 visits'),('GSE121239','65 patients · 292 PBMC visits'),('GSE72747','10 adults · 30 treatment visits'),('GSE32591','32 LN · 14/15 controls'),('GSE200306','53 adults · 162 biopsy rows'),('CD3×CD19 trial','11 donors · 81,941 cells')]
for i,(name,detail) in enumerate(cards):
 x=(i%4)*.25;y=.60-(i//4)*.32
 ax.add_patch(FancyBboxPatch((x,y),.236,.25,boxstyle='round,pad=0.008,rounding_size=.02',fc='#F0F5F8',ec='#D1DEE7',lw=.6,transform=ax.transAxes))
 ax.text(x+.012,y+.158,name,fontsize=6.6,fontweight='bold',transform=ax.transAxes);ax.text(x+.012,y+.061,detail,fontsize=5.9,transform=ax.transAxes)
ax.text(.01,.14,'SLEDAI = serological DNA/complement points + clinical items',fontsize=7,fontweight='bold',color='#24384C',transform=ax.transAxes)
ax.text(.01,.015,'Fixed IFN / PB / LDG / ESI / HK → state • cell/tissue • treatment • prospective risk',fontsize=6.8,transform=ax.transAxes)
rows=[]
for out,short in [('Sero_SLEDAI','Serology'),('Est_cSLEDAI','Clinical')]:
 for trial in ['ILLUMINATE-1','ILLUMINATE-2']:
  for m in ['M_IFN','M_PB','M_LDG']:
   r=find(R[2]['phase0_reproduction'],trial=trial,outcome=out,module=m);rows.append(val(r,'beta',label=f'{short} · T{trial[-1]} · {LAB[m]}',module=m))
forest(fig,[.04,.589,.935,.187],'b','ILLUMINATE: mutually adjusted module loadings (T1 n=876; T2 n=880)',rows,'',showq=True,header='β [95% CI]',fontsize=6.6)
rows=[]
for marker,short in [('logDNA','log DNA'),('negative_C3','−C3'),('negative_C4','−C4')]:
 for trial in ['ILLUMINATE-1','ILLUMINATE-2']:
  for m in ['M_PB','M_LDG']:
   r=find(R[2]['slopes'],trial=trial,marker=marker,module=m);rows.append(val(r,'beta_SD_per_SD',label=f'{short} · T{trial[-1]} · {LAB[m]}',module=m))
forest(fig,[.04,.365,.935,.166],'c','Continuous serology: PB / LDG slopes (module SD per serology SD)',rows,'',showq=True,header='β [95% CI]',fontsize=6.6)
rows=[val(P[m],label=LAB[m],module=m) for m in ['M_IFN','M_PB','M_LDG','ESI','NC_HK']]
forest(fig,[.04,.215,.935,.090],'d','Clinically active − SACQ: 75 patient-paired state contrasts',rows,'',header='Δ [95% CI]')
f=pd.read_csv(ROOT/'results/S0_frozen_FE_with_CIs.csv');rect=[.04,.027,.935,.128];addtitle(fig,rect,'e','Patient fixed effects: treatment and neutrophil-adjusted coefficients',dy=.022)
for j,(family,short,y) in enumerate([('GSE65391_patient_FE_seropositive_SACQ_vs_ClinActive','Seropositive: clinical activity indicator',.104),('GSE65391_patient_FE_adjusted_steroids_MMF_neutrophils_all_visits','All visits: clinical SLEDAI',.024)]):
 rows=[]
 for m in ['M_IFN','M_PB','M_LDG']:
  r=f.loc[(f.family==family)&(f.term==m)].iloc[0];rows.append({'label':LAB[m],'est':float(r.beta),'ci':[float(r.CI_low),float(r.CI_high)],'p':float(r.p),'q':None,'color':COL[m]})
 forest(fig,[.04,y,.935,.051],'e'+str(j),short,rows,'',header='β [95% CI]',fontsize=6.3,add_title=False,left_header='Sero+ (activity indicator)' if j==0 else 'All visits (clinical SLEDAI)')
save(fig,'Figure_1_score_decomposition_and_states')

# FIGURE 2. Subtype localization and the distinct renal compartment test.
CURRENT='Figure2';fig=plt.figure(figsize=(W,180/25.4))
for j,m in enumerate(['M_LDG','M_IFN']):
 rect=[.095+j*.49,.603,.365,.275];a=fig.add_axes(rect);base=ld.loc[ld.timepoint=='Base'];u=base.pivot(index='donor_id',columns='cell_type',values=m).reindex(columns=CTS)
 for donor,row in u.iterrows():a.plot(range(3),row.values,color='#B4BDC5',lw=.6,alpha=.7,zorder=1)
 for i,ct in enumerate(CTS):
  vals=u[ct].dropna();offset=np.linspace(-.045,.045,len(vals));a.scatter(i+offset,vals,color=CTCOL[i],marker=CTMARK[i],s=16,edgecolor='white',lw=.3,zorder=3)
 a.set_xticks(range(3),[f'{s}\nn={u[c].notna().sum()}' for s,c in zip(CTLAB,CTS)]);a.set_xlim(-.25,2.25);a.set_ylabel(LAB[m]+' module score',fontsize=7);a.tick_params(labelsize=6.5)
 r=find(R[4]['baseline_subtype'],measure=m,reference='mature neutrophil',comparison='immature neutrophil')
 addtitle(fig,rect,chr(97+j),'Baseline '+LAB[m]+' across paired donor states',dy=.073)
 a.text(0,1.075,f"Immature − mature: Δ={num(r['delta'])}; p={pstr(r['p'])}; n={r['n']} pairs",fontsize=6.4,transform=a.transAxes)
 DATA[CURRENT+chr(97+j)]={'contrast':r,'n_per_state':{ct:int(u[ct].notna().sum()) for ct in CTS}}
rows=[]
for m in ['M_IFN','M_PB','M_LDG']:
 for comp,short in [('Glomeruli','Glom · 32/14'),('Tubulointerstitium','TI · 32/15')]:
  r=find(R[4]['kidney']['compartment_results'],module=m,compartment=comp);rows.append(val(r,'cohens_d','cohens_d_ci95',label=f'{LAB[m]} · {short}',module=m))
forest(fig,[.04,.277,.935,.175],'c','Renal microdissection: LN versus controls (2,000 donor bootstraps)',rows,"Cohen d",showq=True,header='d [95% CI]')
rows=[]
for m in ['M_IFN','M_PB','M_LDG','ESI','NC_HK']:
 r=find(R[4]['kidney']['interaction'],module=m);rows.append(val(r,'difference_in_compartment_differences',label=LAB[m],module=m))
forest(fig,[.04,.045,.935,.13],'d','Direct renal compartment interaction (32 LN + 14 control pairs)',rows,'Glom − TI: LN minus control (module units)',showq=True,header='Δ [95% CI]')
save(fig,'Figure_2_cell_states_and_renal_compartments')

# FIGURE 3. Induction paired changes and CD3×CD19 cell-state/composition trajectories.
CURRENT='Figure3';fig=plt.figure(figsize=(W,190/25.4))
for i,m in enumerate(['M_IFN','M_PB','M_LDG','SLEDAI']):
 row=i//2;col=i%2;rect=[.095+col*.49,.754-row*.284,.365,.16 if row==0 else .17];a=fig.add_axes(rect)
 DATA[CURRENT+chr(97+i)]=trajectory(a,ad,'subject','month',m,[0,3,6],COL[m]);a.set_xlabel('Months after induction',fontsize=6.6,labelpad=2);a.set_ylabel('Total SLEDAI' if m=='SLEDAI' else LAB[m]+' module score',fontsize=7)
 addtitle(fig,rect,chr(97+i),('Total SLEDAI' if m=='SLEDAI' else LAB[m]+' program')+' | 10 adults',dy=.059)
 for k,t in enumerate([3,6]):
  rr=find(R[3]['paired'],measure=m,timepoint_month=t);note=f"M{t}: Δ={num(rr['delta'])}; p={pstr(rr['p'])}"
  if m=='M_LDG' and t==3:note+='; q='+pstr(rr['BH_q'])
  a.text(0,1.165-k*.082,note,fontsize=6.2,transform=a.transAxes,color=COL[m] if m in ['M_PB','M_LDG'] else '#222222',fontweight='bold' if m=='M_LDG' and t==3 else 'normal')
# e: all subtype fractions on one compact axis; equal-donor summaries, all observed points.
rect=[.095,.265,.865,.125];a=fig.add_axes(rect)
for i,ct in enumerate(CTS):
 DATA[CURRENT+'e_'+ct]=trajectory(a,ld.loc[ld.cell_type==ct],'donor_id','timepoint','fraction',['Base','D28','M3'],CTCOL[i],CTMARK[i],scale=100,alpha=.22)
a.set_ylabel('Granulocytes\n(% sampled)',fontsize=7);a.legend([Line2D([],[],color=c,marker=m,lw=1,ms=3) for c,m in zip(CTCOL,CTMARK)],CTLAB,ncol=3,loc='upper right',bbox_to_anchor=(1,-.125),fontsize=6,borderaxespad=0,frameon=False)
addtitle(fig,rect,'e','CD3×CD19 therapy: composition trajectories (10 Base–M3 donor pairs)',dy=.027)
imm=find(R[4]['paired'],compartment='granulocytes',subtype=CTS[0],measure='fraction',timepoint='M3');mat=find(R[4]['paired'],compartment='granulocytes',subtype=CTS[2],measure='fraction',timepoint='M3');inter=find(R[4]['change_contrasts'],measure='fraction',comparison_change=CTS[0],reference_change=CTS[2])
a.text(0,1.12,f"Immature {imm['mean_base']*100:.1f}% → {imm['mean_followup']*100:.1f}%: p={pstr(imm['p'])}, q={pstr(imm['BH_q'])} | Mature Δ={mat['delta']*100:+.1f} pp: p={pstr(mat['p'])}",transform=a.transAxes,fontsize=6.4)
a.text(0,1.035,f"Immature − mature change: {inter['delta']*100:.1f} pp; p={pstr(inter['p'])}, q={pstr(inter['BH_q'])}",transform=a.transAxes,fontsize=6.4)
for j,m in enumerate(['M_LDG','M_IFN']):
 rect=[.095+j*.49,.052,.365,.125];a=fig.add_axes(rect)
 for i,ct in enumerate(CTS):DATA[CURRENT+chr(102+j)+'_'+ct]=trajectory(a,ld.loc[ld.cell_type==ct],'donor_id','timepoint',m,['Base','D28','M3'],CTCOL[i],CTMARK[i],alpha=.22)
 a.set_ylabel(LAB[m]+' score',fontsize=7);addtitle(fig,rect,chr(102+j),'Within-state '+LAB[m]+' trajectories',dy=.04)
 a.legend([Line2D([],[],color=c,marker=mark,lw=1,ms=3) for c,mark in zip(CTCOL,CTMARK)],CTLAB,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.22),fontsize=5.7,columnspacing=.6,handlelength=1,frameon=False)
 rr=find(R[4]['paired'],compartment='granulocytes',subtype='immature neutrophil',measure=m,timepoint='M3');a.text(0,-.30,f"Immature Base–M3: Δ={num(rr['delta'])}; p={pstr(rr['p'])}",transform=a.transAxes,fontsize=6.2)
save(fig,'Figure_3_treatment_and_cellular_redistribution')

# Frozen KM geometry: convert old SVG paths to data coordinates; NO fitting or re-estimation.
def svg_vertices(d):
 toks=re.findall(r'[MLZmlz]|[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?',d);verts=[];codes=[];mode=None;i=0
 while i<len(toks):
  tok=toks[i]
  if tok in ['M','L','m','l']:mode=tok.upper();i+=1;continue
  if tok in ['Z','z']:
   verts.append(verts[0]);codes.append(MPath.CLOSEPOLY);i+=1;continue
  x=float(tok);y=float(toks[i+1]);verts.append((x,y));codes.append(MPath.MOVETO if mode=='M' else MPath.LINETO);mode='L';i+=2
 return np.array(verts),codes
ns={'s':'http://www.w3.org/2000/svg'};svgpath=ROOT/'figures/Figure_1_longitudinal_and_mediation.svg';svgroot=ET.parse(svgpath).getroot() if svgpath.exists() else None;idmap={el.attrib['id']:el for el in svgroot.iter() if 'id' in el.attrib} if svgroot is not None else {}
KMFILE=SUP/'KM_GEOMETRY.json';KMCACHE=json.loads(KMFILE.read_text()) if KMFILE.exists() else {}
def style(s):return {k.strip():v.strip() for k,v in (z.split(':',1) for z in s.split(';') if ':' in z)}
def km_frozen(a,axisid,out):
 if axisid in KMCACHE:
  geo=KMCACHE[axisid]
  for r in geo['bands']:a.add_patch(PathPatch(MPath(r['vertices'],r['codes']),fc=r['color'],ec='none',alpha=r['alpha'],zorder=1))
  for r in geo['lines']:a.add_patch(PathPatch(MPath(r['vertices'],r['codes']),fill=False,edgecolor=r['color'],lw=1.05,zorder=3))
  for r in geo['censors']:a.plot(r['xy'][0],r['xy'][1],'+',color=r['color'],ms=3,mew=.7,zorder=4)
  a.set_xlim(*geo['xlim']);a.set_ylim(0,1.02);a.set_ylabel('Probability without activity',fontsize=7);a.set_xlabel('Days after first SACQ',fontsize=7,labelpad=2);a.tick_params(labelsize=6.5,pad=2)
  a.text(.025,.045,'Log-rank p='+pstr(R[1]['landmark'][out]['logrank_p']),transform=a.transAxes,fontsize=6.5)
  a.legend([Line2D([],[],color=c,lw=1.1) for c in ['#2166AC','#B2182B']],['Lower ESI (n=39)','Higher ESI (n=39)'],loc='upper right',bbox_to_anchor=(1.01,1.0),fontsize=6,frameon=False)
  return geo['provenance']
 geo={'bands':[],'lines':[],'censors':[]}
 g=idmap[axisid];rectpath=g.find('s:g/s:path',ns);v,_=svg_vertices(rectpath.attrib['d']);left,right=v[:,0].min(),v[:,0].max();top,bottom=v[:,1].min(),v[:,1].max()
 axis=next(c for c in g if c.attrib.get('id','').startswith('matplotlib.axis_'));ticks=[]
 for tick in axis:
  if not tick.attrib.get('id','').startswith('xtick'):continue
  use=tick.find('.//s:use',ns);text=tick.find('.//s:text',ns)
  if use is not None and text is not None:
   try:ticks.append((float(use.attrib['x']),float(''.join(text.itertext()).replace('−','-'))))
   except ValueError:pass
 assert len(ticks)>=2
 scale=(ticks[1][1]-ticks[0][1])/(ticks[1][0]-ticks[0][0]);offset=ticks[0][1]-ticks[0][0]*scale
 def data(verts,dy=0):return np.column_stack((verts[:,0]*scale+offset,(bottom-(verts[:,1]+dy))/(bottom-top)*1.02))
 counts={'CI_paths':0,'survival_paths':0,'censor_marks':0}
 for child in g:
  cid=child.attrib.get('id','')
  if cid.startswith('FillBetweenPolyCollection'):
   use=child.find('.//s:use',ns);ref=use.attrib['{http://www.w3.org/1999/xlink}href'][1:];path=idmap[ref];v,c=svg_vertices(path.attrib['d']);v=data(v,float(use.attrib.get('y',0)));st=style(use.attrib.get('style',''))
   a.add_patch(PathPatch(MPath(v,c),fc=st.get('fill','#999999'),ec='none',alpha=float(st.get('fill-opacity',.25)),zorder=1));counts['CI_paths']+=1;geo['bands'].append({'vertices':v.tolist(),'codes':[int(z) for z in c],'color':st.get('fill','#999999'),'alpha':float(st.get('fill-opacity',.25))})
  elif cid.startswith('line2d_'):
   for path in child.findall('s:path',ns):
    if 'clip-path' not in path.attrib:continue
    v,c=svg_vertices(path.attrib['d']);st=style(path.attrib.get('style',''));a.add_patch(PathPatch(MPath(data(v),c),fill=False,edgecolor=st.get('stroke','#222222'),lw=1.05,zorder=3));counts['survival_paths']+=1;geo['lines'].append({'vertices':data(v).tolist(),'codes':[int(z) for z in c],'color':st.get('stroke','#222222')})
   for use in child.findall('.//s:use',ns):
    if 'x' not in use.attrib or 'y' not in use.attrib:continue
    xy=data(np.array([[float(use.attrib['x']),float(use.attrib['y'])]]))[0];st=style(use.attrib.get('style',''));a.plot(xy[0],xy[1],'+',color=st.get('stroke','#222222'),ms=3,mew=.7,zorder=4);counts['censor_marks']+=1;geo['censors'].append({'xy':xy.tolist(),'color':st.get('stroke','#222222')})
 a.set_xlim(left*scale+offset,right*scale+offset);a.set_ylim(0,1.02);a.set_ylabel('Probability without activity',fontsize=7);a.set_xlabel('Days after first SACQ',fontsize=7,labelpad=2);a.tick_params(labelsize=6.5,pad=2)
 a.text(.025,.045,'Log-rank p='+pstr(R[1]['landmark'][out]['logrank_p']),transform=a.transAxes,fontsize=6.5)
 a.legend([Line2D([],[],color=c,lw=1.1) for c in ['#2166AC','#B2182B']],['Lower ESI (n=39)','Higher ESI (n=39)'],loc='upper right',bbox_to_anchor=(1.01,1.0),fontsize=6,frameon=False)
 assert counts['CI_paths']==2 and counts['survival_paths']==2 and counts['censor_marks']>0,counts
 geo['xlim']=[left*scale+offset,right*scale+offset]
 geo['provenance']={'method':'frozen SVG path reuse; no Kaplan–Meier fit','source':str(svgpath.relative_to(ROOT)),'source_sha256':hashlib.sha256(svgpath.read_bytes()).hexdigest(),'axis_id':axisid,**counts}
 KMCACHE[axisid]=geo;KMFILE.write_text(json.dumps(KMCACHE,indent=2));return geo['provenance']
CURRENT='Figure4';fig=plt.figure(figsize=(W,190/25.4));rows=[]
for mod,label in [('unadjusted','Unadjusted · n=274'),('demography_treatment_interval','Treatment / interval · n=274'),('serology_adjusted','Complete serology · n=208')]:
 r=find(R[1]['next_visit_models'],outcome='any',predictor='ESI',model=mod);rows.append(val(r,'OR',label=label,module='ESI'))
r=find(R[1]['multiple_imputation'],analysis='next_visit',outcome='any',predictor='ESI');rows.append(val(r,'estimate',label='20 imputations · n=274',module='ESI'))
forest(fig,[.04,.775,.935,.115],'a','Next-visit clinical activity: sequential adjustment (78 patients)',rows,'OR per SD of ESI',zero=1,log=True,header='OR [95% CI]',label_frac=.32,table_frac=.40)
for j,(axisid,out) in enumerate([('axes_2','any'),('axes_3','renal')]):
 rect=[.10+j*.49,.478,.365,.205];a=fig.add_axes(rect);DATA[CURRENT+chr(98+j)]=km_frozen(a,axisid,out);addtitle(fig,rect,chr(98+j),'Clinical activity' if out=='any' else 'Renal activity',dy=.025)
rows=[]
for out,short in [('any','Clinical'),('renal','Renal')]:
 for model,lab in [('unadjusted','unadjusted'),('adjusted','adjusted')]:
  r=find(R[1]['cox'],outcome=out,predictor='ESI',model=model);rows.append(val(r,'HR',label=f'{short} · {lab}',module='ESI'))
forest(fig,[.04,.254,.935,.112],'d','First-SACQ landmark ESI and time to activity (n=78)',rows,'HR per SD of ESI',zero=1,log=True,header='HR [95% CI]')
rows=[]
for exp in ['PB','IFN']:
 for out,short in [('any','clinical'),('renal','renal')]:
  r=find(R[1]['mediation'],exposure=exp,outcome=out);rows.append(val(r['ACME'],'estimate',p='bootstrap_sign_p',label=f'{exp} → {short}',module=exp,scale=100))
forest(fig,[.04,.05,.935,.112],'e','Lower-baseline-C3 decomposition (267 transitions; 77 patient clusters)',rows,'Model indirect risk difference (percentage points)',header='pp [95% CI]')
save(fig,'Figure_4_prospective_models_and_C3_decomposition')

# S1: all observations from the original quadrant displays.
CURRENT='S1';s=pd.read_csv(ROOT/'results/task2_trial_source.csv');quads=['Sero+ / low clinical','Sero+ / high clinical','Sero- / low clinical','Sero- / high clinical'];qcol=['#2166AC','#B2182B','#777777','#F1A340'];rng=np.random.default_rng(20261116)
fig,axs=plt.subplots(3,2,figsize=(W,187/25.4),layout='constrained');fig.set_constrained_layout_pads(w_pad=.05,h_pad=.11,wspace=.12,hspace=.08)
for row,m in enumerate(['M_IFN','M_PB','M_LDG']):
 for col,trial in enumerate(['ILLUMINATE-1','ILLUMINATE-2']):
  a=axs[row,col];ss=s.loc[s.batch==trial];ns_=[]
  for i,q in enumerate(quads):
   vals=ss.loc[ss.quadrant==q,m].dropna().values;ns_.append(len(vals))
   if len(vals):
    a.scatter(i+rng.normal(0,.055,len(vals)),vals,s=3.2,color=qcol[i],alpha=.38,edgecolor='none');a.boxplot([vals],positions=[i],widths=.48,showfliers=False,patch_artist=True,boxprops={'facecolor':'none','edgecolor':qcol[i],'linewidth':.8},medianprops={'color':'#222222','linewidth':1},whiskerprops={'color':qcol[i],'linewidth':.6},capprops={'color':qcol[i],'linewidth':.6})
  a.set_xticks(range(4),['S+ L','S+ H','S− L','S− H']);a.set_xlim(-.5,3.5);a.set_ylabel(LAB[m]+' module score');a.set_title(trial+' · '+LAB[m]+'\nn='+', '.join(map(str,ns_)),loc='left',fontsize=8,pad=8);a.text(-.19,1.13,chr(97+row*2+col),fontweight='bold',fontsize=10,transform=a.transAxes);DATA[CURRENT+chr(97+row*2+col)]={'n_by_quadrant':dict(zip(quads,ns_)),'all_points':sum(ns_)}
save(fig,'Supplementary_Figure_S1_trial_quadrants')
# S2: observed donor correlations and the complete original adult sensitivities.
CURRENT='S2';fig=plt.figure(figsize=(W,193/25.4));cs=pd.read_csv(ROOT/'results/task4_B_granulocyte_correlations_source.csv')
for j,ct in enumerate(['plasmablast','class switched memory B cell']):
 rect=[.105+j*.49,.64,.355,.25];a=fig.add_axes(rect);sub=cs.loc[(cs.B_subtype==ct)&(cs.timepoint=='M3')];a.scatter(sub.B_fraction*100,sub.immature_contraction*100,color='#762A83',s=20,edgecolors='white',lw=.3,zorder=3)
 for _,r in sub.iterrows():
  off=(-5,-8) if j==1 and r.donor=='P7' else (3,3)
  a.annotate(r.donor,(r.B_fraction*100,r.immature_contraction*100),xytext=off,textcoords='offset points',fontsize=6,color='#222222',ha='right' if off[0]<0 else 'left')
 rr=find(R[4]['correlations'],B_subtype=ct,timepoint='M3');a.set_xlabel('M3 '+('plasmablast' if j==0 else 'memory B')+' fraction (%)',fontsize=7);a.set_ylabel('Immature contraction (pp)',fontsize=7);a.tick_params(labelsize=6.5);a.margins(.15)
 addtitle(fig,rect,chr(97+j),'B-lineage association | 10 donors',dy=.052);a.text(0,1.075,'Spearman ρ='+num(rr['spearman_rho'])+'; p='+pstr(rr['spearman_permutation_p']),transform=a.transAxes,fontsize=6.5);DATA[CURRENT+chr(97+j)]=rr
for j,(cohort,title) in enumerate([('GSE49454','Adult whole blood: operational state contrast (10 paired adults)'),('GSE121239','Adult PBMC: high − low total SLEDAI (37 paired adults)')]):
 rows=[val(R[3]['sensitivities'][cohort]['paired'][m],label=LAB[m],module=m) for m in ['M_IFN','M_PB','M_LDG','ESI','NC_HK']]
 forest(fig,[.04,.345-j*.262,.935,.14],chr(99+j),title,rows,'Paired module difference',header='Δ [95% CI]')
save(fig,'Supplementary_Figure_S2_B_lineage_and_adult_sensitivities')
(SUP/'FIGURE_VALUES_v1.3.json').write_text(json.dumps(DATA,indent=2,allow_nan=False))
(SUP/'FIGURE_LAYOUT_v1.3.json').write_text(json.dumps(QA,indent=2))
print('Six figures exported; numerical drawing payload and layout records saved.')
