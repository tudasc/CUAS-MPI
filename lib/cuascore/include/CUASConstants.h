#ifndef CUAS_PHYSICALCONSTANTS_H
#define CUAS_PHYSICALCONSTANTS_H

#include <cfloat>

// Also a tiny value and NOFLOW value which are not really constants but here anyway
#define TINY FLT_EPSILON
#define NOFLOW_VALUE 1e-14  // value used to mimic noflow condition

#define RHO_ICE 910
#define RHO_WATER 1000
#define LATENT_HEAT 334000
#define GRAVITY 9.81
#define SPY 3.154e7
#define Ss 0.000982977696

#define COMPUTE_FLAG 0          // active cuas grid point
#define DIRICHLET_FLAG 1        // lateral boundary
#define NOFLOW_FLAG 2           // lateral boundary
#define DIRICHLET_OCEAN_FLAG 3  // connected to ocean, tidal influenced
#define DIRICHLET_LAKE_FLAG 4   // connected to periglacial lake, not tidal influenced
#define INVALID_FLAG 5          // from gen_bnd_mask.py

#endif
