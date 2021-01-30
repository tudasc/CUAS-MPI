#ifndef CUAS_HELPER_H
#define CUAS_HELPER_H

#include "physicalConstants.h"

#include "PetscGrid.h"

namespace CUAS {

inline void pressure2head(PetscGrid &result, PetscGrid &pressure, PetscGrid &bed_elevation,
                          PetscScalar const sea_level) {
  auto result2d = result.getAsGlobal2dArr();
  auto pressure2d = pressure.getAsGlobal2dArr();
  auto bed_elevation2d = bed_elevation.getAsGlobal2dArr();
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d[j][i] - sea_level;
      result2d[j][i] = pressure2d[j][i] / (RHO_WATER * GRAVITY) + effective_bed_elevation;
    }
  }
  result.setAsGlobal2dArr(result2d);
  pressure.restoreGlobal2dArr(pressure2d);
  bed_elevation.restoreGlobal2dArr(bed_elevation2d);
}

inline void head2pressure(PetscGrid &result, PetscGrid &head, PetscGrid &bed_elevation, PetscScalar const sea_level) {
  auto result2d = result.getAsGlobal2dArr();
  auto head2d = head.getAsGlobal2dArr();
  auto bed_elevation2d = bed_elevation.getAsGlobal2dArr();
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d[j][i] - sea_level;
      result2d[j][i] = RHO_WATER * GRAVITY * (head2d[j][i] - effective_bed_elevation);
    }
  }
  result.setAsGlobal2dArr(result2d);
  head.restoreGlobal2dArr(head2d);
  bed_elevation.restoreGlobal2dArr(bed_elevation2d);
}

inline void overburdenPressure(PetscGrid &result, PetscGrid &thk) {
  auto result2d = result.getAsGlobal2dArr();
  auto thk2d = thk.getAsGlobal2dArr();
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      result2d[j][i] = thk2d[j][i] * RHO_ICE * GRAVITY;
    }
  }
  result.setAsGlobal2dArr(result2d);
  thk.restoreGlobal2dArr(thk2d);
}

inline void cavityOpenB(PetscGrid &result, PetscScalar const beta, PetscScalar const v_b, PetscGrid &K) {
  auto resultGlobal = result.getAsGlobal2dArr();
  auto KGlobal = K.getAsGlobal2dArr();
  PetscScalar betaVb = beta * v_b;
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      resultGlobal[j][i] = betaVb * KGlobal[j][i];
    }
  }
  result.setAsGlobal2dArr(resultGlobal);
  K.restoreGlobal2dArr(KGlobal);
}

inline void compute_melt(PetscGrid &result, PetscScalar const r, PetscScalar const g, PetscScalar const rho_w,
                         PetscGrid &T, PetscGrid &K, PetscGrid &gradh2, PetscScalar const rho_i, PetscScalar const L,
                         PetscScalar const bt) {
  auto resultGlobal = result.getAsGlobal2dArr();
  auto KGlobal = K.getAsGlobal2dArr();
  auto TGlobal = T.getAsGlobal2dArr();
  auto gradh2Global = gradh2.getAsGlobal2dArr();
  PetscScalar r_g_rhow = r * g * rho_w;
  PetscScalar rhoi_L = rho_i * L;
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      resultGlobal[j][i] = r_g_rhow * TGlobal[j][i] * KGlobal[j][i] * gradh2Global[j][i] / rhoi_L;
    }
  }
  result.setAsGlobal2dArr(resultGlobal);
  K.restoreGlobal2dArr(KGlobal);
  T.restoreGlobal2dArr(TGlobal);
  gradh2.restoreGlobal2dArr(gradh2Global);
}

}  // namespace CUAS

#endif
