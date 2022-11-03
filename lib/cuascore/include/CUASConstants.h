#ifndef CUAS_PHYSICALCONSTANTS_H
#define CUAS_PHYSICALCONSTANTS_H

// Also a tiny value and NOFLOW value which are not really constants but here anyway
#ifndef NOFLOW_VALUE
#define NOFLOW_VALUE 1e-14  // value used to mimic noflow condition
#endif
#ifndef TINY
#define TINY 1.0e-20  // much smaller than NOFLOW_VALUE and consistent with CUAS-python version
#endif
#ifndef RHO_ICE
#define RHO_ICE 910.0
#endif
#ifndef SPY
#define SPY 31536000.0  //  365 day calendar used in lib/cuascore/include/timeparse.h
#endif
#ifndef RHO_WATER
#define RHO_WATER 1000.0
#endif
#ifndef LATENT_HEAT
#define LATENT_HEAT 334000.0
#endif
#ifndef GRAVITY
#define GRAVITY 9.81
#endif
#ifndef PSI_BELOW_ZERO_REPLACE_VALUE
#define PSI_BELOW_ZERO_REPLACE_VALUE 0.01
#endif

#define COMPUTE_FLAG 0          // active cuas grid point
#define DIRICHLET_FLAG 1        // lateral boundary
#define NOFLOW_FLAG 2           // lateral boundary
#define DIRICHLET_OCEAN_FLAG 3  // connected to ocean, tidal influenced
#define DIRICHLET_LAKE_FLAG 4   // connected to periglacial lake, not tidal influenced
#define INVALID_FLAG 5          // from gen_bnd_mask.py

#endif
