#include <float.h>
// Physical constants
// Also a tiny value and NOFLOW value which are not really constants but here anyway
#define TINY FLT_EPSILON
#define RHO_ICE 910
#define RHOW_WATER 1000
#define LATENT_HEAT 334000
#define GRAVITY 9.81
#define SPY 3.154e7
#define Ss 0.000982977696
#define NOFLOW_VALUE 1e-14  // value used to mimic noflow condition
#define NOFLOW_FLAG 2
#define COMPUTE_FLAG 0
#define DIRICHLET_FLAG 1       // connection to the ocean
#define DIRICHLET_LAKE_FLAG 3  // connected to periglacial lake, not tidal influenced
#define INVALID_FLAG 5         // from gen_bnd_mask.py
