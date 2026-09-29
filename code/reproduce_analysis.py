from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DATA=ROOT/'data'; DATA.mkdir(exist_ok=True)
sys.path.insert(0,str(HERE))
from cardioreg_bounds import (
    lp_bounds, lp_bounds_interval_signatures,
    lp_bounds_terminal_prevalence,
    shared_rate_contrast_exact, arm_specific_rate_contrast_exact,
    shared_rate_contrast_from_difference_intervals,
    constrained_least_squares,
)
from simulation_utils import (
    S_REF,X_CENTER,TARGET_PRODUCTIVE,TARGET_ONE,TARGET_REPEAT,TARGET_DEATH,
    simultaneous_t_intervals,simulate_paired_panel,STATES,ASSAYS,
)

SEED=20260929
rng=np.random.default_rng(SEED)

# 1) Published Leone/Musa/Engel summary and contrast sensitivities ----------------
leone=pd.DataFrame([
    ['FBS','two_sided',.25,.13,5],['FBS','one_sided',.33,.14,5],['FBS','no_sided',.42,.09,5],
    ['FGF1_p38i','two_sided',.43,.10,5],['FGF1_p38i','one_sided',.26,.13,5],['FGF1_p38i','no_sided',.31,.05,5],
],columns=['condition','furrow_category','mean_fraction','sd_fraction','n_independent_experiments'])
leone.to_csv(DATA/'leone_published_summary.csv',index=False)
cats=['two_sided','one_sided','no_sided']
C=leone[leone.condition.eq('FBS')].set_index('furrow_category').loc[cats]
T=leone[leone.condition.eq('FGF1_p38i')].set_index('furrow_category').loc[cats]
meanC=C.mean_fraction.to_numpy(); meanT=T.mean_fraction.to_numpy()
rate_low=np.array([.90,0,0.]); rate_high=np.array([1,.15,.15])
shared=shared_rate_contrast_exact(meanC,meanT,rate_low,rate_high)
arm_specific=arm_specific_rate_contrast_exact(meanC,meanT,rate_low,rate_high)

# Simultaneous 95% intervals for the three *differences*, not four marginal means.
# Independent-arm Welch standard errors are used because raw pairing/covariance is unavailable.
alpha=.05; rows=[]
for k,cat in enumerate(cats):
    mc,mt=float(C.loc[cat,'mean_fraction']),float(T.loc[cat,'mean_fraction'])
    sc,st=float(C.loc[cat,'sd_fraction']),float(T.loc[cat,'sd_fraction'])
    nc,nt=int(C.loc[cat,'n_independent_experiments']),int(T.loc[cat,'n_independent_experiments'])
    diff=mt-mc
    se=np.sqrt(sc**2/nc+st**2/nt)
    df=(sc**2/nc+st**2/nt)**2/((sc**2/nc)**2/(nc-1)+(st**2/nt)**2/(nt-1))
    crit=stats.t.ppf(1-alpha/(2*3),df)
    rows.append([cat,diff,se,df,crit,diff-crit*se,diff+crit*se])
diff_df=pd.DataFrame(rows,columns=['furrow_category','difference_T_minus_C','welch_se','welch_df','bonferroni_critical_t','simultaneous_lower','simultaneous_upper'])
diff_df.to_csv(DATA/'leone_difference_intervals.csv',index=False)
diff_intervals=diff_df[['simultaneous_lower','simultaneous_upper']].to_numpy()
summary_shared=shared_rate_contrast_from_difference_intervals(diff_intervals,rate_low,rate_high)

contrast=pd.DataFrame([
    ['published_means_shared_rates',shared['lower'],shared['upper'],'Sensitivity set conditional on one shared morphology-to-outcome rate vector; not a confidence interval.'],
    ['published_means_arm_specific_rates',arm_specific['lower'],arm_specific['upper'],'Same rate boxes allowed to differ between treatment arms.'],
    ['summary_difference_intervals_shared_rates',summary_shared['lower'],summary_shared['upper'],'Uses all three published composition SDs and simultaneous Welch intervals for T-C differences; assumes independent arms because raw covariance/pairing is unavailable.'],
],columns=['analysis','lower','upper','interpretation'])
contrast.to_csv(DATA/'leone_contrast_sensitivity.csv',index=False)

# 2) Endpoint denominator and net-gain identification ----------------------------
x=np.array([.60,.25,.15])
d=np.array([1.,2.,0.])
S_terminal=np.array([[.10,.80,0.],[.60,.90,0.]])
p=(S_terminal*d[None,:])@x/(d@x)
net_c=d-1
no_count=lp_bounds_terminal_prevalence(S_terminal,d,np.c_[p,p],net_c)
with_count=lp_bounds_terminal_prevalence(S_terminal,d,np.c_[p,p],net_c,(1.05,1.15))
exact_count=lp_bounds_terminal_prevalence(S_terminal,d,np.c_[p,p],net_c,(1.10,1.10))
naive=float(S_terminal[0]@x); correct=float(p[0]); truth=float(net_c@x)
endpoint=pd.DataFrame([
    ['naive_baseline_weighted_prevalence',naive,np.nan,np.nan],
    ['correct_endpoint_cell_prevalence',correct,np.nan,np.nan],
    ['net_gain_truth',truth,np.nan,np.nan],
    ['net_gain_from_endpoint_prevalences_only',np.nan,no_count['lower'],no_count['upper']],
    ['net_gain_plus_cell_count_ratio_1.05_to_1.15',np.nan,with_count['lower'],with_count['upper']],
    ['net_gain_plus_exact_cell_count_ratio_1.10',np.nan,exact_count['lower'],exact_count['upper']],
],columns=['quantity','point_value','lower','upper'])
endpoint.to_csv(DATA/'endpoint_denominator_and_net_gain.csv',index=False)

# 3) Structural partial identification ------------------------------------------
obs=S_REF@X_CENTER
struct=[]
for name,c in [('productive_any',TARGET_PRODUCTIVE),('one_division',TARGET_ONE),('repeated_division',TARGET_REPEAT),('death',TARGET_DEATH)]:
    ans=lp_bounds(S_REF,np.c_[obs,obs],c)
    struct.append([name,float(c@X_CENTER),ans['lower'],ans['upper'],ans['upper']-ans['lower']])
# Public-audit consequence: if every signature entry is unsupported [0,1], the target is vacuous.
uninformative=lp_bounds_interval_signatures(np.zeros_like(S_REF),np.ones_like(S_REF),np.tile([[0,1]],(S_REF.shape[0],1)),TARGET_PRODUCTIVE)
struct.append(['productive_any_all_signatures_0_to_1',float(TARGET_PRODUCTIVE@X_CENTER),uninformative['lower'],uninformative['upper'],uninformative['upper']-uninformative['lower']])
struct_df=pd.DataFrame(struct,columns=['target','truth','lower','upper','width'])
struct_df.to_csv(DATA/'structural_target_bounds.csv',index=False)

# 4) Cluster-aware synthetic design ---------------------------------------------
# The signatures are hypothetical and are not presented as biological calibration.
# The misspecified scenario lowers nonproductive-state signatures by 0.25 and then
# analyzes with a narrow +/-0.01 box around that wrong matrix; this is deliberately
# severe so coverage failure becomes observable as biological-unit precision rises.
B=100
simrows=[]
miss_center=S_REF.copy(); miss_center[:,[0,1,2,3,6]]=np.clip(miss_center[:,[0,1,2,3,6]]-.25,0,1)
scenarios={
    'exact':(S_REF,S_REF),
    'hypothetical_box_0.05':(np.clip(S_REF-.05,0,1),np.clip(S_REF+.05,0,1)),
    'hypothetical_box_0.15':(np.clip(S_REF-.15,0,1),np.clip(S_REF+.15,0,1)),
    'deliberately_misspecified':(np.clip(miss_center-.01,0,1),np.clip(miss_center+.01,0,1)),
}
true=float(TARGET_PRODUCTIVE@X_CENTER)
for n in [3,5,8,12,20,30]:
    for scen,(Sl,Sh) in scenarios.items():
        covered=[]; infeas=[]; lows=[]; highs=[]; nnls_err=[]; clipped=[]
        for _ in range(B):
            props,counts=simulate_paired_panel(rng,n,S_REF)
            ints,raw=simultaneous_t_intervals(props)
            clipped.append(np.mean((raw[:,0]<0)|(raw[:,1]>1)))
            ans=lp_bounds_interval_signatures(Sl,Sh,ints,TARGET_PRODUCTIVE)
            infeas.append(ans['status']!='OK')
            if ans['status']=='OK':
                lows.append(ans['lower']); highs.append(ans['upper']); covered.append(ans['lower']<=true<=ans['upper'])
            else:
                covered.append(False)
            xhat=constrained_least_squares((Sl+Sh)/2,props.mean(axis=1))
            nnls_err.append(abs(float(TARGET_PRODUCTIVE@xhat)-true))
        simrows.append([n,scen,B,float(np.mean(covered)),float(np.mean(infeas)),float(np.median(lows)) if lows else np.nan,float(np.median(highs)) if highs else np.nan,float(np.median(np.array(highs)-np.array(lows))) if lows else np.nan,float(np.median(nnls_err)),float(np.mean(clipped))])
sim_df=pd.DataFrame(simrows,columns=['n_biological_units','scenario','n_sim','coverage_true_target','infeasible_fraction','median_lower_feasible','median_upper_feasible','median_width_feasible','median_constrained_LS_target_abs_error','mean_fraction_assay_intervals_clipped'])
sim_df.to_csv(DATA/'cluster_simulation.csv',index=False)

# 5) Signature invariance: population-level and finite-sample null contrast -------
# Both groups have the same latent x. Treatment changes transient visibility only.
# This asks how much apparent contrast can arise from using S_control for treatment.
inv_rows=[]
for shift in [0,.05,.10,.15,.20]:
    St=S_REF.copy(); St[:3,4:6]=np.clip(St[:3,4:6]+shift,0,1)
    pC=S_REF@X_CENTER; pT=St@X_CENTER
    bC=lp_bounds(S_REF,np.c_[pC,pC],TARGET_PRODUCTIVE)
    bT_wrong=lp_bounds(S_REF,np.c_[pT,pT],TARGET_PRODUCTIVE)
    pop_lo=bT_wrong['lower']-bC['upper']; pop_hi=bT_wrong['upper']-bC['lower']
    for n in [5,8,12,30]:
        false_effect=0; widths=[]; valid=0
        for _ in range(60):
            pc,_=simulate_paired_panel(rng,n,S_REF)
            pt,_=simulate_paired_panel(rng,n,St)
            ic,_=simultaneous_t_intervals(pc); it,_=simultaneous_t_intervals(pt)
            c=lp_bounds(S_REF,ic,TARGET_PRODUCTIVE); t=lp_bounds(S_REF,it,TARGET_PRODUCTIVE)
            if c['status']!='OK' or t['status']!='OK': continue
            lo=t['lower']-c['upper']; hi=t['upper']-c['lower']; valid+=1
            false_effect += int(lo>0 or hi<0); widths.append(hi-lo)
        inv_rows.append([shift,n,pop_lo,pop_hi,0.0,false_effect/max(valid,1),float(np.median(widths)) if widths else np.nan])
inv_df=pd.DataFrame(inv_rows,columns=['treatment_visibility_shift','n_biological_units_per_arm','population_shared_S_contrast_lower','population_shared_S_contrast_upper','true_contrast','finite_sample_false_effect_fraction','median_finite_sample_contrast_width'])
inv_df.to_csv(DATA/'signature_invariance_null_sensitivity.csv',index=False)

summary={
    'seed':SEED,
    'leone_shared_rate_means':[shared['lower'],shared['upper']],
    'leone_arm_specific_rate_means':[arm_specific['lower'],arm_specific['upper']],
    'leone_summary_shared_rate_interval':[summary_shared['lower'],summary_shared['upper']],
    'endpoint_prevalence_example':{'naive':naive,'correct':correct},
    'net_gain_bounds_without_count':[no_count['lower'],no_count['upper']],
    'net_gain_bounds_with_count_1.05_1.15':[with_count['lower'],with_count['upper']],
    'productive_structural_bounds':struct_df.iloc[0][['lower','upper']].tolist(),
    'public_uninformative_signature_bounds':[uninformative['lower'],uninformative['upper']],
    'scope':'AJP version is a literature-based Perspective with no unpublished simulation results; full manuscript includes computational design analyses but claims no biological validation.'
}
(DATA/'analysis_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
