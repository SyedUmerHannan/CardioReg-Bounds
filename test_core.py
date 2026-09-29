import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from cardioreg_bounds import (
    lp_bounds, lp_bounds_interval_signatures,
    lp_bounds_terminal_prevalence, lp_bounds_terminal_interval_signatures,
    shared_rate_contrast_exact, arm_specific_rate_contrast_exact,
    shared_rate_contrast_from_difference_intervals,
    target_globally_identified,
)
from simulation_utils import S_REF,X_CENTER,TARGET_PRODUCTIVE,TARGET_ONE,TARGET_REPEAT,simulate_paired_panel,simultaneous_t_intervals


def test_interval_signature_point_equivalence():
    S=np.array([[.1,.8],[.2,.6]]); obs=np.array([[.38,.42],[.34,.38]]); c=np.array([0.,1.])
    a=lp_bounds(S,obs,c); b=lp_bounds_interval_signatures(S,S,obs,c)
    assert a['status']==b['status']=='OK'
    assert abs(a['lower']-b['lower'])<1e-10 and abs(a['upper']-b['upper'])<1e-10


def test_signature_widening_cannot_narrow():
    S=np.array([[.05,.8,.9],[.01,.25,.7]]); obs=np.array([[.30,.34],[.16,.20]]); c=np.array([0,0,1.])
    t=lp_bounds_interval_signatures(S,S,obs,c)
    w=lp_bounds_interval_signatures(np.clip(S-.05,0,1),np.clip(S+.05,0,1),obs,c)
    assert w['lower']<=t['lower']+1e-10 and w['upper']>=t['upper']-1e-10


def test_terminal_reduces_when_descendants_one():
    S=np.array([[.1,.8,.2],[.3,.4,.9]]); obs=np.array([[.30,.35],[.45,.50]]); c=np.array([0,1,0.])
    a=lp_bounds(S,obs,c); b=lp_bounds_terminal_prevalence(S,np.ones(3),obs,c)
    assert a['status']==b['status']=='OK'
    assert abs(a['lower']-b['lower'])<1e-9 and abs(a['upper']-b['upper'])<1e-9


def test_endpoint_denominator_truth_and_count_constraint():
    x=np.array([.60,.25,.15]); d=np.array([1.,2.,0.]); S=np.array([[.10,.80,0.],[.60,.90,0.]])
    p=(S*d[None,:])@x/(d@x); c=d-1; truth=float(c@x)
    no=lp_bounds_terminal_prevalence(S,d,np.c_[p,p],c)
    exact=lp_bounds_terminal_prevalence(S,d,np.c_[p,p],c,(1.10,1.10))
    assert no['lower']<=truth<=no['upper']
    assert abs(exact['lower']-truth)<1e-9 and abs(exact['upper']-truth)<1e-9


def test_all_death_is_closure_not_artificial_epsilon():
    d=np.array([1.,2.,0.]); S=np.array([[.1,.8,0.]])
    # Endpoint prevalence alone cannot exclude the all-death limit; net gain infimum is -1.
    ans=lp_bounds_terminal_prevalence(S,d,np.array([[.4,.5]]),d-1)
    assert ans['status']=='OK' and abs(ans['lower']+1)<1e-9


def test_terminal_interval_point_equivalence():
    S=np.array([[.1,.8,0.],[.6,.9,0.]]); d=np.array([1.,2.,0.]); obs=np.array([[.4,.45],[.7,.8]]); c=np.array([0,1,-1.])
    a=lp_bounds_terminal_prevalence(S,d,obs,c,(.5,2.0))
    b=lp_bounds_terminal_interval_signatures(S,S,d,obs,c,(.5,2.0))
    assert a['status']==b['status']=='OK'
    assert abs(a['lower']-b['lower'])<1e-9 and abs(a['upper']-b['upper'])<1e-9


def test_leone_shared_and_arm_specific_means():
    C=np.array([.25,.33,.42]); T=np.array([.43,.26,.31]); lo=np.array([.9,0,0]); hi=np.array([1,.15,.15])
    s=shared_rate_contrast_exact(C,T,lo,hi); a=arm_specific_rate_contrast_exact(C,T,lo,hi)
    assert abs(s['lower']-.135)<1e-12 and abs(s['upper']-.18)<1e-12
    assert abs(a['lower']-.0245)<1e-12 and abs(a['upper']-.2905)<1e-12


def test_difference_interval_contrast_uses_sum_zero():
    D=np.array([[-.04491701,.40491701],[-.32802221,.18802221],[-.25920904,.03920904]])
    ans=shared_rate_contrast_from_difference_intervals(D,np.array([.9,0,0]),np.array([1,.15,.15]))
    assert ans['status']=='OK'
    assert ans['lower']<0<ans['upper']


def test_underdetermined_targets_differ_in_resolution():
    p=S_REF@X_CENTER
    prod=lp_bounds(S_REF,np.c_[p,p],TARGET_PRODUCTIVE)
    one=lp_bounds(S_REF,np.c_[p,p],TARGET_ONE)
    rep=lp_bounds(S_REF,np.c_[p,p],TARGET_REPEAT)
    assert prod['upper']-prod['lower'] < .02
    assert one['upper']-one['lower'] > .08
    assert rep['upper']-rep['lower'] > .07


def test_uninformative_signatures_are_vacuous():
    obs=np.tile([[0,1]],(S_REF.shape[0],1))
    ans=lp_bounds_interval_signatures(np.zeros_like(S_REF),np.ones_like(S_REF),obs,TARGET_PRODUCTIVE)
    assert abs(ans['lower'])<1e-12 and abs(ans['upper']-1)<1e-12


def test_global_rowspace_boundary_counterexample():
    A=np.array([[1,1,1],[2,3,5]],float)
    assert target_globally_identified(A,2*A[0]-.5*A[1])
    assert not target_globally_identified(A,np.array([1,0,0.]))
    B=np.array([[0,1]],float)
    ans=lp_bounds(B,np.array([[0,0]],float),np.array([1,0.]))
    assert ans['status']=='OK' and abs(ans['lower']-1)<1e-12 and abs(ans['upper']-1)<1e-12


def test_cluster_simulator_shapes_and_clipping():
    rng=np.random.default_rng(99)
    props,counts=simulate_paired_panel(rng,5,S_REF)
    ints,raw=simultaneous_t_intervals(props)
    assert props.shape==(4,5) and counts.shape==(4,5) and ints.shape==(4,2)
    assert np.all((props>=0)&(props<=1)) and counts.min()>=120 and counts.max()<=800
    assert np.all(ints[:,0]<=ints[:,1])
