#ifndef CUAS_KERNELS_H
#define CUAS_KERNELS_H

#include "CUASConstants.h"
#include "specialgradient.h"

#include "Logger.h"
#include "PETScGrid.h"

#include <algorithm>

namespace CUAS {

/*
 * Converts hydraulic head to water pressure
 * bedElevation and seaLevel are relative to the mean sea level = 0.0
 */
inline void head2pressure(PETScGrid &result, PETScGrid const &head, PETScGrid const &bedElevation,
                          PetscScalar const seaLevel = 0.0) {
  if (!result.isCompatible(head) || !result.isCompatible(bedElevation)) {
    Logger::instance().error("CUASKernels.h: head2pressure was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto result2d = result.getWriteHandleGhost();
  auto &head2d = head.getReadHandle();
  auto &bedElevation2d = bedElevation.getReadHandle();

  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      double effective_bed_elevation = bedElevation2d(j, i, true) - seaLevel;
      result2d(j, i) = RHO_WATER * GRAVITY * (head2d(j, i, true) - effective_bed_elevation);
    }
  }
}

/*
 * Convert water pressure to hydraulic head
 * bed_elevation and sea_level are relative to the mean sea level = 0.0
 * pressure is ice pressure
 */
inline void pressure2head(PETScGrid &result, PETScGrid const &pressure, PETScGrid const &bedElevation,
                          PetscScalar const seaLevel = 0.0) {
  if (!result.isCompatible(pressure) || !result.isCompatible(bedElevation)) {
    Logger::instance().error("CUASKernels.h: pressure2head was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto result2d = result.getWriteHandleGhost();
  auto &pressure2d = pressure.getReadHandle();
  auto &bedElevation2d = bedElevation.getReadHandle();

  constexpr double rgMultiplicator = 1.0 / (RHO_WATER * GRAVITY);

  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      double effective_bed_elevation = bedElevation2d(j, i, true) - seaLevel;
      result2d(j, i) = pressure2d(j, i, true) * rgMultiplicator + effective_bed_elevation;
    }
  }
}

/*
 * Compute ice overburden pressure
 */
inline void overburdenPressure(PETScGrid &result, PETScGrid const &thk) {
  if (!result.isCompatible(thk)) {
    Logger::instance().error("CUASKernels.h: overburdenPressure was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto result2d = result.getWriteHandleGhost();
  auto &thk2d = thk.getReadHandle();
  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      result2d(j, i) = thk2d(j, i, true) * RHO_ICE * GRAVITY;
    }
  }
}

/*
 * From Werder2013 / summers2018:
 * beta = (b r − b)/l r for b < b r , beta = 0 for b ≥ b r
 */
inline void cavityOpenB(PETScGrid &result, PetscScalar const beta, PetscScalar const v_b, PETScGrid const &K) {
  if (!result.isCompatible(K)) {
    Logger::instance().error("CUASKernels.h: cavityOpenB was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto resultGlobal = result.getWriteHandleGhost();
  auto &KGlobal = K.getReadHandle();
  PetscScalar betaVb = beta * v_b;
  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      resultGlobal(j, i) = betaVb * KGlobal(j, i, true);
    }
  }
}

/*
 * like deFleurian2016
 */
inline void computeMelt(PETScGrid &melt, PetscScalar const r, PetscScalar const g, PetscScalar const rho_w,
                        PETScGrid const &T, PETScGrid const &K, PETScGrid const &gradh2, PetscScalar const rho_i,
                        PetscScalar const L, PetscScalar const bt) {
  if (!melt.isCompatible(T) || !melt.isCompatible(T) || !melt.isCompatible(K) || !melt.isCompatible(gradh2)) {
    Logger::instance().error("CUASKernels.h: computeMelt was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto resultGlobal = melt.getWriteHandle();
  auto &KGlobal = K.getReadHandle();
  auto &TGlobal = T.getReadHandle();
  auto &gradh2Global = gradh2.getReadHandle();
  const PetscScalar r_g_rhow = r * g * rho_w;
  const PetscScalar rhoi_L_inv = 1.0 / (rho_i * L);
  for (int j = 0; j < melt.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < melt.getLocalNumOfCols(); ++i) {
      resultGlobal(j, i) = r_g_rhow * TGlobal(j, i) * KGlobal(j, i) * gradh2Global(j, i) * rhoi_L_inv;
    }
  }
}
/** morphology operation for expanding the shapes in an image.
 *
 * @param output
 * @param input
 */
inline void binaryDilation(PETScGrid &output, PETScGrid const &input) {
  if (!output.isCompatible(input)) {
    Logger::instance().error("CUASKernels.h: binaryDilation was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto &in = input.getReadHandle();
  auto out = output.getWriteHandle();

  int index_cols = 1;
  int index_rows = 1;

  auto iter_rows = output.getLocalNumOfRows();
  auto iter_cols = output.getLocalNumOfCols();

  for (int i = 0; i < iter_rows; ++i) {
    for (int j = 0; j < iter_cols; ++j) {
      if (in(index_rows, index_cols, GHOSTED))
        out(i, j) = true;
      if (in(index_rows + 1, index_cols, GHOSTED))
        out(i, j) = true;
      if (in(index_rows - 1, index_cols, GHOSTED))
        out(i, j) = true;
      if (in(index_rows, index_cols + 1, GHOSTED))
        out(i, j) = true;
      if (in(index_rows, index_cols - 1, GHOSTED))
        out(i, j) = true;
      ++index_cols;
    }
    index_cols = 1;
    ++index_rows;
  }
}

inline void enableUnconfined(PETScGrid &Teff, PETScGrid &TeffPowTexp, PETScGrid &Sp, PETScGrid const &T_n,
                             PETScGrid const &K, PETScGrid const &topg, PETScGrid const &u_n, PetscScalar const Texp,
                             PetscScalar const unconfSmooth, PetscScalar const bt) {
  if (!Teff.isCompatible(TeffPowTexp) || !Teff.isCompatible(u_n) || !Teff.isCompatible(T_n) ||
      !Teff.isCompatible(topg) || !Teff.isCompatible(K) || !Teff.isCompatible(Sp)) {
    Logger::instance().error("CUASKernels.h: enableUnconfined was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto &u_nGlobal = u_n.getReadHandle();
  auto &T_nGlobal = T_n.getReadHandle();
  auto &topgGlobal = topg.getReadHandle();
  auto &KGlobal = K.getReadHandle();
  auto TeffGlobal = Teff.getWriteHandle();
  auto TeffPowTexpGlobal = TeffPowTexp.getWriteHandle();
  auto SpGlobalWrite = Sp.getWriteHandle();

  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      double psi = u_nGlobal(j, i) - topgGlobal(j, i);
      if (psi < 0) {
        psi = 0.01;
      }
      if (psi < bt) {
        TeffGlobal(j, i) = KGlobal(j, i) * psi;
      } else {
        TeffGlobal(j, i) = T_nGlobal(j, i);
      }
      TeffPowTexpGlobal(j, i) = pow(TeffGlobal(j, i), Texp);

      if (psi < (bt - unconfSmooth)) {
        SpGlobalWrite(j, i) = 0.4;
      } else {
        SpGlobalWrite(j, i) = 0;
      }
    }
  }
}

inline void calculateTeffPowTexp(PETScGrid &Teff, PETScGrid &TeffPowTexp, PETScGrid const &T, PetscScalar const Texp) {
  if (!Teff.isCompatible(TeffPowTexp) || !Teff.isCompatible(T)) {
    Logger::instance().error("CUASKernels.h: calculateTeffPowTexp was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  // copy T to Teff
  auto &TGlobal = T.getReadHandle();
  auto TeffGlobal = Teff.getWriteHandle();
  auto TeffPowTexpGlobal = TeffPowTexp.getWriteHandle();
  for (int j = 0; j < Teff.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Teff.getLocalNumOfCols(); ++i) {
      TeffGlobal(j, i) = TGlobal(j, i);
      TeffPowTexpGlobal(j, i) = pow(TeffGlobal(j, i), Texp);
    }
  }
}

inline void calculateSeValues(PETScGrid &Se, PETScGrid const &Sp, PETScGrid const &S) {  // extra scope for handles
  if (!Se.isCompatible(Sp) || !Se.isCompatible(S)) {
    Logger::instance().error("CUASKernels.h: calculateSeValues was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto SeGlobal = Se.getWriteHandle();
  auto &SpGlobalRead = Sp.getReadHandle();
  auto &SGlobal = S.getReadHandle();
  for (int j = 0; j < Se.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Se.getLocalNumOfCols(); ++i) {
      SeGlobal(j, i) = SGlobal(j, i) + SpGlobalRead(j, i);  // effective storage
    }
  }
}

inline void convolveStar11411(PETScGrid &melt, PETScGrid &result) {
  if (!melt.isCompatible(result)) {
    Logger::instance().error("CUASKernels.h: convolveStar11411 was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  // the stencil is: stencil = np.array([[0, 1, 0],
  //                                    [1, 4, 1],
  //                                    [0, 1, 0]]) / 8

  auto melt2d = melt.getReadHandle();
  auto result2d = result.getWriteHandle();
  int ghostIndexCols = 1;
  int ghostIndexRows = 1;

  for (int i = 0; i < result.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < result.getLocalNumOfCols(); ++j) {
      result2d(i, j) =
          (4.0 * melt2d(ghostIndexRows, ghostIndexCols, true) + melt2d(ghostIndexRows, ghostIndexCols - 1, true) +
           melt2d(ghostIndexRows, ghostIndexCols + 1, true) + melt2d(ghostIndexRows - 1, ghostIndexCols, true) +
           melt2d(ghostIndexRows + 1, ghostIndexCols, true)) /
          8.0;
      ++ghostIndexCols;
    }
    ghostIndexCols = 1;
    ++ghostIndexRows;
  }
}

inline void clamp(PETScGrid &input, PetscScalar const minimum, PetscScalar const maximum) {
  auto inputGlobal = input.getWriteHandle();
  for (int j = 0; j < input.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < input.getLocalNumOfCols(); ++i) {
      inputGlobal(j, i) = std::clamp(inputGlobal(j, i), minimum, maximum);
    }
  }
}

inline void doChannels(PETScGrid &melt, PETScGrid &creep, PETScGrid const &u_n, PETScGrid const &gradMask, PETScGrid &T,
                       PETScGrid const &T_n, PETScGrid const &pIce, PETScGrid const &topg, PETScGrid const &K,
                       PETScGrid const &bndMask, PETScGrid &cavityOpening, PetscScalar const flowConstant,
                       PetscScalar const Texp, PetscScalar const roughnessFactor, bool const noSmoothMelt,
                       PetscScalar const cavityBeta, PetscScalar const basalVelocityIce, PetscScalar const tMin,
                       PetscScalar const tMax, PetscScalar const bt, PetscScalar const dx, PetscScalar const dt_secs,
                       bool doMelt = true, bool doCreep = true, bool doCavity = true) {
  if (!melt.isCompatible(creep) || !melt.isCompatible(u_n)) {
    Logger::instance().error("CUASKernels.h: doChannels was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  PETScGrid maggrad2(u_n.getTotalNumOfCols(), u_n.getTotalNumOfRows());
  gradient2(maggrad2, u_n, dx);
  auto maggrad2Global = maggrad2.getWriteHandle();
  auto &grad_maskGlobal = gradMask.getReadHandle();
  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      // if (maggrad2Global(j, i) == grad_maskGlobal(j, i)) {
      if (grad_maskGlobal(j, i) == 1) {
        maggrad2Global(j, i) = 0;
      }
    }
  }

  auto creepGlobal = creep.getWriteHandle();
  auto &TGlobal = T.getReadHandle();
  auto &p_iceGlobal = pIce.getReadHandle();
  auto &u_nGlobal = u_n.getReadHandle();
  auto &topgGlobal = topg.getReadHandle();
  PetscScalar multValue = 2 * flowConstant;
  PetscScalar onethird = 1.0 / 3.0;
  PetscScalar rhowater_gravity = RHO_WATER * GRAVITY;
  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      double p_w = (u_nGlobal(j, i) - topgGlobal(j, i) - bt) * rhowater_gravity;
      double N = p_iceGlobal(j, i) - p_w;
      creepGlobal(j, i) = multValue * pow(N * onethird, 3) * TGlobal(j, i);
    }
  }

  if (Texp != 1) {
    PETScGrid T_nExp(u_n.getTotalNumOfCols(), u_n.getTotalNumOfRows());
    auto T_nExpGlobal = T_nExp.getWriteHandle();
    auto &T_nGlobal = T_n.getReadHandle();
    for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
        T_nExpGlobal(j, i) = pow(T_nGlobal(j, i), Texp);
      }
    }
    computeMelt(melt, roughnessFactor, GRAVITY, RHO_WATER, T_nExp, K, maggrad2, RHO_ICE, LATENT_HEAT, bt);
  } else {
    computeMelt(melt, roughnessFactor, GRAVITY, RHO_WATER, T_n, K, maggrad2, RHO_ICE, LATENT_HEAT, bt);
  }

  if (!noSmoothMelt) {
    PETScGrid result(melt.getTotalNumOfCols(), melt.getTotalNumOfRows());
    convolveStar11411(melt, result);
    melt.copy(result);
  }
  cavityOpenB(cavityOpening, cavityBeta, basalVelocityIce, K);

  // Update T with melt cavityOpening and creep
  {
    auto TGlobal = T.getWriteHandle();
    auto T_nGlobal = T_n.getReadHandle();
    auto meltGlobal = melt.getReadHandle();
    auto creepGlobal = creep.getReadHandle();
    auto cavityGlobal = cavityOpening.getReadHandle();
    if (doMelt && doCavity && doCreep) {
      for (int j = 0; j < T.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < T.getLocalNumOfCols(); ++i) {
          TGlobal(j, i) = T_nGlobal(j, i) + (meltGlobal(j, i) + cavityGlobal(j, i) - creepGlobal(j, i)) * dt_secs;
        }
      }
    } else {
      for (int j = 0; j < T.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < T.getLocalNumOfCols(); ++i) {
          TGlobal(j, i) = TGlobal(j, i) =
              T_nGlobal(j, i) + ((doMelt ? meltGlobal(j, i) : 0.0) + (doCavity ? cavityGlobal(j, i) : 0.0) -
                                 (doCreep ? creepGlobal(j, i) : 0.0)) *
                                    dt_secs;
        }
      }
    }
  }

  clamp(T, tMin, tMax);

  {
    auto TGlobal = T.getWriteHandle();
    auto mask = bndMask.getReadHandle();
    for (int j = 0; j < T.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < T.getLocalNumOfCols(); ++i) {
        if (mask(j, i) == (PetscScalar)NOFLOW_FLAG) {
          TGlobal(j, i) = NOFLOW_VALUE;
        }
      }
    }
  }
}

inline void noChannels(PETScGrid &melt, PETScGrid &creep, PETScGrid &cavityOpening) {
  melt.setZero();
  creep.setZero();
  cavityOpening.setZero();

  // this is probably faster than setting zero
  /*auto meltGlobal = melt.getWriteHandle();
  auto creepGlobal = creep.getWriteHandle();
  for (int j = 0; j < melt.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < melt.getLocalNumOfCols(); ++i) {
      meltGlobal(j, i) = 0;
      creepGlobal(j, i) = 0;
      cavityOpening(j, i) = 0;
    }
  }*/
}

}  // namespace CUAS

#endif
