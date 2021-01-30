#ifndef CUAS_HELPER_H
#define CUAS_HELPER_H

#include "physicalConstants.h"

#include "PetscGrid.h"

namespace CUAS {

inline void pressure2head(PetscGrid &result, PetscGrid const &pressure, PetscGrid const &bed_elevation,
                          PetscScalar const sea_level) {
  auto result2d = result.getWriteHandle();
  auto &pressure2d = pressure.getReadHandle();
  auto &bed_elevation2d = bed_elevation.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d(j, i) - sea_level;
      result2d(j, i) = pressure2d(j, i) / (RHO_WATER * GRAVITY) + effective_bed_elevation;
    }
  }
}

inline void head2pressure(PetscGrid &result, PetscGrid const &head, PetscGrid const &bed_elevation,
                          PetscScalar const sea_level) {
  auto result2d = result.getWriteHandle();
  auto &head2d = head.getReadHandle();
  auto &bed_elevation2d = bed_elevation.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d(j, i) - sea_level;
      result2d(j, i) = RHO_WATER * GRAVITY * (head2d(j, i) - effective_bed_elevation);
    }
  }
}

inline void overburdenPressure(PetscGrid &result, PetscGrid const &thk) {
  auto result2d = result.getWriteHandle();
  auto &thk2d = thk.getReadHandle();
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      result2d(j, i) = thk2d(j, i) * RHO_ICE * GRAVITY;
    }
  }
}

inline void cavityOpenB(PetscGrid &result, PetscScalar const beta, PetscScalar const v_b, PetscGrid const &K) {
  auto resultGlobal = result.getWriteHandle();
  auto &KGlobal = K.getReadHandle();
  PetscScalar betaVb = beta * v_b;
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      resultGlobal(j, i) = betaVb * KGlobal(j, i);
    }
  }
}

inline void compute_melt(PetscGrid &result, PetscScalar const r, PetscScalar const g, PetscScalar const rho_w,
                         PetscGrid const &T, PetscGrid const &K, PetscGrid const &gradh2, PetscScalar const rho_i,
                         PetscScalar const L, PetscScalar const bt) {
  auto resultGlobal = result.getWriteHandle();
  auto &KGlobal = K.getReadHandle();
  auto &TGlobal = T.getReadHandle();
  auto &gradh2Global = gradh2.getReadHandle();
  PetscScalar r_g_rhow = r * g * rho_w;
  PetscScalar rhoi_L = rho_i * L;
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      resultGlobal(j, i) = r_g_rhow * TGlobal(j, i) * KGlobal(j, i) * gradh2Global(j, i) / rhoi_L;
    }
  }
}

}  // namespace CUAS

#endif
