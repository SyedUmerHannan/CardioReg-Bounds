"""Core LP utilities for CardioReg-Bounds.

These functions implement baseline-cohort prevalence models, descendant-weighted
endpoint prevalence models, interval-valued signatures, and contrast bounds.
They deliberately separate biological assumptions from optimization.
"""
from __future__ import annotations
import itertools
import numpy as np
from scipy.optimize import linprog, minimize


def _solve(Aub, bub, Aeq, beq, bounds, c):
    rlo = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq,
                  bounds=bounds, method="highs")
    if not rlo.success:
        return {"status":"INFEASIBLE","lower":np.nan,"upper":np.nan}
    rhi = linprog(-np.asarray(c), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq,
                  bounds=bounds, method="highs")
    if not rhi.success:
        return {"status":"OPTIMIZATION_FAILURE","lower":np.nan,"upper":np.nan}
    return {"status":"OK","lower":float(rlo.fun),"upper":float(-rhi.fun)}


def lp_bounds(S, intervals, c):
    """Target bounds for baseline-unit prevalences p = Sx."""
    S=np.asarray(S,float); intervals=np.asarray(intervals,float); c=np.asarray(c,float)
    A,K=S.shape
    if intervals.shape!=(A,2) or c.shape!=(K,): raise ValueError("shape mismatch")
    if np.any(intervals[:,0]>intervals[:,1]): raise ValueError("invalid intervals")
    Aub=np.vstack([S,-S]); bub=np.r_[intervals[:,1],-intervals[:,0]]
    return _solve(Aub,bub,np.ones((1,K)),np.array([1.0]),[(0,1)]*K,c)


def lp_bounds_interval_signatures(S_low,S_high,intervals,c):
    """Baseline-unit bounds with row-wise independent signature boxes."""
    L=np.asarray(S_low,float); U=np.asarray(S_high,float)
    if L.shape!=U.shape or np.any(L>U) or np.any(L<0) or np.any(U>1):
        raise ValueError("invalid signature box")
    intervals=np.asarray(intervals,float); c=np.asarray(c,float)
    A,K=L.shape
    if intervals.shape!=(A,2) or c.shape!=(K,): raise ValueError("shape mismatch")
    Aub=np.vstack([L,-U]); bub=np.r_[intervals[:,1],-intervals[:,0]]
    return _solve(Aub,bub,np.ones((1,K)),np.array([1.0]),[(0,1)]*K,c)


def lp_bounds_terminal_prevalence(S, descendants, intervals, c, descendant_ratio_bounds=None):
    """Bounds when assays are prevalences among descendants present at H.

    p_a = sum_j x_j d_j s_aj / sum_j x_j d_j.

    The LP uses the closure of the positive-denominator feasible set. Thus if no
    cell-count/survival lower bound is supplied, an all-death limit can appear as
    an infimum. This is mathematically informative: endpoint prevalences alone do
    not identify net gain. Supply descendant_ratio_bounds=(L,U) when N(H)/N0 is
    observed or otherwise bounded.
    """
    S=np.asarray(S,float); d=np.asarray(descendants,float); intervals=np.asarray(intervals,float); c=np.asarray(c,float)
    A,K=S.shape
    if d.shape!=(K,) or c.shape!=(K,) or intervals.shape!=(A,2): raise ValueError("shape mismatch")
    if np.any(d<0): raise ValueError("descendant counts must be nonnegative")
    Aub=[]; bub=[]
    for a in range(A):
        lo,hi=intervals[a]
        Aub.append(d*(S[a]-hi)); bub.append(0.0)
        Aub.append(d*(lo-S[a])); bub.append(0.0)
    if descendant_ratio_bounds is not None:
        L,U=map(float,descendant_ratio_bounds)
        if L<0 or L>U: raise ValueError("invalid descendant ratio bounds")
        Aub.append(d); bub.append(U)
        Aub.append(-d); bub.append(-L)
    return _solve(np.asarray(Aub),np.asarray(bub),np.ones((1,K)),np.array([1.0]),[(0,1)]*K,c)


def lp_bounds_terminal_interval_signatures(S_low,S_high,descendants,intervals,c,descendant_ratio_bounds=None):
    """Endpoint-cell prevalence model with row-wise interval signatures."""
    L=np.asarray(S_low,float); U=np.asarray(S_high,float); d=np.asarray(descendants,float)
    if L.shape!=U.shape or np.any(L>U) or np.any(L<0) or np.any(U>1):
        raise ValueError("invalid signature box")
    intervals=np.asarray(intervals,float); c=np.asarray(c,float)
    A,K=L.shape
    if d.shape!=(K,) or c.shape!=(K,) or intervals.shape!=(A,2): raise ValueError("shape mismatch")
    Aub=[]; bub=[]
    for a in range(A):
        lo,hi=intervals[a]
        Aub.append(d*(L[a]-hi)); bub.append(0.0)
        Aub.append(d*(lo-U[a])); bub.append(0.0)
    if descendant_ratio_bounds is not None:
        loD,hiD=map(float,descendant_ratio_bounds)
        Aub.append(d); bub.append(hiD)
        Aub.append(-d); bub.append(-loD)
    return _solve(np.asarray(Aub),np.asarray(bub),np.ones((1,K)),np.array([1.0]),[(0,1)]*K,c)


def shared_rate_contrast_from_difference_intervals(delta_intervals, rate_low, rate_high):
    """Bounds Delta = r' delta_y with shared rates and sum(delta_y)=0.

    delta_intervals are simultaneous intervals for category differences T-C.
    Extrema are obtained by enumerating rate-box corners and solving a linear LP
    in delta_y for each corner.
    """
    D=np.asarray(delta_intervals,float); rl=np.asarray(rate_low,float); rh=np.asarray(rate_high,float)
    K=len(rl)
    if D.shape!=(K,2) or rh.shape!=(K,): raise ValueError("shape mismatch")
    lows=[]; highs=[]
    for bits in itertools.product([0,1], repeat=K):
        r=np.array([rh[k] if bits[k] else rl[k] for k in range(K)])
        bnds=[tuple(D[k]) for k in range(K)]
        lo=linprog(r,A_eq=np.ones((1,K)),b_eq=[0.0],bounds=bnds,method="highs")
        hi=linprog(-r,A_eq=np.ones((1,K)),b_eq=[0.0],bounds=bnds,method="highs")
        if lo.success and hi.success:
            lows.append(lo.fun); highs.append(-hi.fun)
    if not lows: return {"status":"INFEASIBLE","lower":np.nan,"upper":np.nan}
    return {"status":"OK","lower":float(min(lows)),"upper":float(max(highs))}


def shared_rate_contrast_exact(control, treatment, rate_low, rate_high):
    """Exact-composition contrast with rate box shared across arms."""
    delta=np.asarray(treatment,float)-np.asarray(control,float)
    return shared_rate_contrast_from_difference_intervals(np.c_[delta,delta],rate_low,rate_high)


def arm_specific_rate_contrast_exact(control,treatment,rate_low,rate_high):
    """Contrast when each arm may take a different rate vector in the same box."""
    C=np.asarray(control,float); T=np.asarray(treatment,float); rl=np.asarray(rate_low,float); rh=np.asarray(rate_high,float)
    valsC=[]; valsT=[]
    for bits in itertools.product([0,1], repeat=len(rl)):
        r=np.array([rh[k] if bits[k] else rl[k] for k in range(len(rl))])
        valsC.append(float(r@C)); valsT.append(float(r@T))
    return {"status":"OK","lower":min(valsT)-max(valsC),"upper":max(valsT)-min(valsC)}


def constrained_least_squares(S, p):
    """Simplex-constrained least-squares point fit using the same S and p."""
    S=np.asarray(S,float); p=np.asarray(p,float); K=S.shape[1]
    fun=lambda x: float(np.sum((S@x-p)**2))
    res=minimize(fun,np.ones(K)/K,bounds=[(0,1)]*K,constraints=[{"type":"eq","fun":lambda x: np.sum(x)-1}],method="SLSQP")
    if not res.success: raise RuntimeError(res.message)
    return res.x


def target_globally_identified(A,c,tol=1e-10):
    A=np.asarray(A,float); c=np.asarray(c,float)
    return np.linalg.matrix_rank(A,tol)==np.linalg.matrix_rank(np.vstack([A,c]),tol)
