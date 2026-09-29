from __future__ import annotations
import numpy as np
from scipy import stats

STATES=['no_cycle','dna_only','endoreplication','multinucleation','one_division','repeated_division','death']
ASSAYS=['DNA_synthesis','mitosis','cytokinesis_associated','persistent_event_label']
S_REF=np.array([
    [0.01,0.90,0.85,0.80,0.88,0.92,0.02],
    [0.001,0.02,0.08,0.65,0.70,0.75,0.01],
    [0.001,0.01,0.03,0.30,0.75,0.82,0.01],
    [0.02,0.15,0.25,0.55,0.85,0.95,0.05],
],float)
X_CENTER=np.array([0.60,0.08,0.07,0.06,0.07,0.02,0.10],float)
TARGET_PRODUCTIVE=np.array([0,0,0,0,1,1,0.],float)
TARGET_ONE=np.array([0,0,0,0,1,0,0.],float)
TARGET_REPEAT=np.array([0,0,0,0,0,1,0.],float)
TARGET_DEATH=np.array([0,0,0,0,0,0,1.],float)


def simultaneous_t_intervals(unit_props,alpha=0.05):
    unit_props=np.asarray(unit_props,float); A,n=unit_props.shape
    if n<2: raise ValueError('at least two biological units required')
    crit=stats.t.ppf(1-alpha/(2*A),df=n-1)
    mean=unit_props.mean(axis=1); se=unit_props.std(axis=1,ddof=1)/np.sqrt(n)
    raw=np.c_[mean-crit*se,mean+crit*se]
    clipped=np.c_[np.clip(raw[:,0],0,1),np.clip(raw[:,1],0,1)]
    return clipped,raw


def simulate_paired_panel(rng,n_units,S_true,x_center=X_CENTER,kappa=80,cell_low=120,cell_high=800):
    """Same biological units contribute all assays; cells scored can differ by assay."""
    S_true=np.asarray(S_true,float); x_center=np.asarray(x_center,float)
    A,K=S_true.shape
    xs=rng.dirichlet(kappa*x_center,size=n_units)
    props=np.zeros((A,n_units)); counts=np.zeros((A,n_units),int)
    for i in range(n_units):
        for a in range(A):
            ncell=int(rng.integers(cell_low,cell_high+1)); counts[a,i]=ncell
            p=float(S_true[a]@xs[i])
            props[a,i]=rng.binomial(ncell,p)/ncell
    return props,counts
