"""
Physical constants used in the model.
Also a tiny value and NOFLOW value which are not really constants but here
anyway.
"""
import numpy as np

TINY = np.finfo(float).eps
RHO_ICE = 910
RHO_WATER = 1000
LATENT_HEAT = 334000
GRAVITY = 9.81
SPY = 3.154e7
Ss = 0.000982977696
NOFLOW_VALUE = 1e-14 # value used to mimic noflow condition

NOFLOW_FLAG = 2
COMPUTE_FLAG = 0
DIRICHLET_FLAG = 1      # connection to the ocean
DIRICHLET_LAKE_FLAG = 3 # connected to periglacial lake, not tidal influenced
INVALID_FLAG = 5        # from gen_bnd_mask.py