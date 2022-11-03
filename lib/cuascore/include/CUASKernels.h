#ifndef CUAS_KERNELS_H
#define CUAS_KERNELS_H

#include "CUASConstants.h"
#include "specialgradient.h"

#include "Logger.h"
#include "PETScGrid.h"

#include <algorithm>

namespace CUAS {

/** Converts hydraulic head to water pressure
 *
 *    zw = 0.0 -> water pressure at the base of the aquifer
 *    zw = layer thickness -> ... at top of aquifer
 *
 * @param result water pressure (Pa)
 * @param hydraulicHead (m)
 * @param bedElevation  (m)
 * @param zw \f$0 \le z_w \le layerThickness\f$
 */
inline void headToPressure(PETScGrid &result, PETScGrid const &hydraulicHead, PETScGrid const &bedElevation,
                           PetscScalar const zw = 0.0) {
  if (!result.isCompatible(hydraulicHead) || !result.isCompatible(bedElevation)) {
    CUAS_ERROR("CUASKernels.h: headToPressure was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandleGhost();
  auto &head = hydraulicHead.getReadHandle();
  auto &topg = bedElevation.getReadHandle();

  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      res(j, i) = (head(j, i, GHOSTED) - topg(j, i, GHOSTED) - zw) * RHO_WATER * GRAVITY;
    }
  }
}

/** Converts hydraulic head to effective pressure
 *
 * This uses water pressure at the aquifer ice interface and ice pressure.
 *
 * \f$ N = p_i - p_w\f$
 *
 * @param result effective pressure (Pa)
 * @param hydraulicHead  (m)
 * @param bedElevation   (m)
 * @param icePressure    (m)
 * @param layerThickness (m)
 */
inline void headToEffectivePressure(PETScGrid &result, PETScGrid const &hydraulicHead, PETScGrid const &bedElevation,
                                    PETScGrid const &icePressure, PetscScalar const layerThickness = 0.0) {
  if (!result.isCompatible(hydraulicHead) || !result.isCompatible(bedElevation)) {
    CUAS_ERROR("CUASKernels.h: headToPressure was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandle();  // model->pIce is without ghosts
  auto &head = hydraulicHead.getReadHandle();
  auto &topg = bedElevation.getReadHandle();
  auto &pi = icePressure.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      // see headToPressure()
      double pw = (head(j, i) - topg(j, i) - layerThickness) * RHO_WATER * GRAVITY;
      res(j, i) = pi(j, i) - pw;
    }
  }
}

/** Convert water pressure to hydraulic head
 *
 * @param result
 * @param waterPressure (Pa)
 * @param bedElevation (m)
 */
inline void pressureToHead(PETScGrid &result, PETScGrid const &waterPressure, PETScGrid const &bedElevation) {
  if (!result.isCompatible(waterPressure) || !result.isCompatible(bedElevation)) {
    CUAS_ERROR("CUASKernels.h: pressureToHead was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandleGhost();
  auto &pw = waterPressure.getReadHandle();
  auto &topg = bedElevation.getReadHandle();

  constexpr double rgMultiplicator = 1.0 / (RHO_WATER * GRAVITY);

  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      res(j, i) = pw(j, i, GHOSTED) * rgMultiplicator + topg(j, i, GHOSTED);
    }
  }
}

/** Compute ice overburden pressure
 *
 * @param result ice pressure (Pa)
 * @param iceThickness (m)
 */
inline void overburdenPressure(PETScGrid &result, PETScGrid const &iceThickness) {
  if (!result.isCompatible(iceThickness)) {
    CUAS_ERROR("CUASKernels.h: overburdenPressure was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandleGhost();
  auto &thk = iceThickness.getReadHandle();
  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      res(j, i) = thk(j, i, GHOSTED) * RHO_ICE * GRAVITY;
    }
  }
}

/** Implements dT/dt_cavity as in Beyer et al., 2018 Eq 5.
 *
 * From Werder2013 / summers2018:
 * beta = (b r − b)/l r for b < b r , beta = 0 for b ≥ b r
 *
 * @param result cavity opening (m^2/s^2)
 * @param beta cavity opening parameter (1)
 * @param hydraulicConductivity (m/s)
 * @param basalVelocity (m/s)
 */
inline void computeCavityOpening(PETScGrid &result, PetscScalar const beta, PetscScalar const hydraulicConductivity,
                                 PETScGrid const &basalVelocity) {
  if (!result.isCompatible(basalVelocity)) {
    CUAS_ERROR("CUASKernels.h: computeCavityOpening was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }
  auto res = result.getWriteHandle();
  auto &vb = basalVelocity.getReadHandle();
  PetscScalar betaK = beta * hydraulicConductivity;
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      res(j, i) = betaK * vb(j, i);
    }
  }
}

/** Implements dT/dt_creep as in Beyer et al., 2018 Eq 5.
 *
 * Note, here we use "creep opening" instead of creep closure by using -1 * creep_closure.
 * Thus, all contributions to dT/dt have the same sign
 *
 * @param result creep opening (m^2/s^2)
 * @param rateFactor (Pa^-3 / s)
 * @param effectivePressure (Pa)
 * @param hydraulicTransmissivity (m^2/s)
 */
inline void computeCreepOpening(PETScGrid &result, PETScGrid const &rateFactor, PETScGrid const &effectivePressure,
                                PETScGrid const &hydraulicTransmissivity) {
  if (!result.isCompatible(rateFactor) || !result.isCompatible(effectivePressure) ||
      !result.isCompatible(hydraulicTransmissivity)) {
    CUAS_ERROR("CUASKernels.h: computeCreepOpening was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandle();
  auto &A = rateFactor.getReadHandle();
  auto &N = effectivePressure.getReadHandle();
  auto &T = hydraulicTransmissivity.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      res(j, i) = -2.0 * A(j, i) * pow(N(j, i) / 3.0, 3) * T(j, i);
    }
  }
}

/** Implements dT/dt_melt as in Beyer et al., 2018 Eq 5.
 *
 * @param result
 * @param roughnessFactor
 * @param hydraulicConductivity
 * @param hydraulicTransmissivity
 * @param gradientHeadSquared
 */
inline void computeMeltOpening(PETScGrid &result, PetscScalar const roughnessFactor,
                               PetscScalar const hydraulicConductivity, PETScGrid const &hydraulicTransmissivity,
                               PETScGrid const &gradientHeadSquared) {
  if (!result.isCompatible(hydraulicTransmissivity) || !result.isCompatible(gradientHeadSquared)) {
    CUAS_ERROR("CUASKernels.h: computeMeltOpening was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandle();
  auto &T = hydraulicTransmissivity.getReadHandle();
  auto &gradh2 = gradientHeadSquared.getReadHandle();
  const PetscScalar r_g_rhow_K = roughnessFactor * GRAVITY * RHO_WATER * hydraulicConductivity;
  const PetscScalar rhoi_L_inv = 1.0 / (RHO_ICE * LATENT_HEAT);
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      res(j, i) = r_g_rhow_K * T(j, i) * gradh2(j, i) * rhoi_L_inv;
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
    CUAS_ERROR("CUASKernels.h: binaryDilation was called with incompatible PETScGrids. Exiting.");
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

/** implements Ehlig & Halepaska eq. 5a,b for effective transmissivity and eq. 7 for effective storativity
 *
 * @param effectiveStorativity (1)
 * @param effectiveTransmissivity (m^2/s)
 * @param hydraulicTransmissivity (m^2/s)
 * @param hydraulicHead (m)
 * @param bedElevation (m)
 * @param bndMask (1)
 * @param layerThickness (m)
 * @param specificStorage (m^-1)
 * @param specificYield (1)
 * @param unconfinedSmooth confined - unconfined transition length, \f$ 0 \le d \le b \f$ (m)
 */
inline void getEffectiveAquiferProperties(PETScGrid &effectiveStorativity, PETScGrid &effectiveTransmissivity,
                                          PETScGrid const &hydraulicTransmissivity, PETScGrid const &hydraulicHead,
                                          PETScGrid const &bedElevation, PETScGrid const &bndMask,
                                          PetscScalar const layerThickness, PetscScalar const specificStorage,
                                          PetscScalar const specificYield, PetscScalar const unconfinedSmooth) {
  if (!effectiveStorativity.isCompatible(bndMask) || !effectiveTransmissivity.isCompatible(bndMask) ||
      !hydraulicTransmissivity.isCompatible(bndMask) || !hydraulicHead.isCompatible(bndMask) ||
      !bedElevation.isCompatible(bndMask)) {
    CUAS_ERROR("CUASKernels.h: getEffectiveAquiferProperties was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  if (unconfinedSmooth < 0.0 || unconfinedSmooth > layerThickness) {
    CUAS_ERROR(
        "CUASKernels.h: getEffectiveAquiferProperties was called with incompatible parameter unconfinedSmooth={}. "
        "Exiting.",
        unconfinedSmooth);
    exit(1);
  }

  auto Seff = effectiveStorativity.getWriteHandle();
  auto Teff = effectiveTransmissivity.getWriteHandle();
  auto &T = hydraulicTransmissivity.getReadHandle();
  auto &head = hydraulicHead.getReadHandle();
  auto &topg = bedElevation.getReadHandle();
  auto &mask = bndMask.getReadHandle();
  for (int j = 0; j < effectiveStorativity.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < effectiveStorativity.getLocalNumOfCols(); ++i) {
      // Seff = S = Ss * b in the confined case
      Seff(j, i) = specificStorage * layerThickness;

      if (mask(j, i) == (PetscScalar)COMPUTE_FLAG) {
        auto psi_org = head(j, i) - topg(j, i);
        // fixme: CUAS is using PSI_BELOW_ZERO_REPLACE_VALUE without justification.
        auto psi = (psi_org < 0.0) ? PSI_BELOW_ZERO_REPLACE_VALUE : psi_org;  // this is needed for Seff, why?
        // implements Ehlig & Halepaska Eq. 5a,b
        Teff(j, i) = (psi < layerThickness) ? T(j, i) / layerThickness * psi : T(j, i);
        // implements Ehlig & Halepaska Eq. 7
        PetscScalar Sprime;
        if (unconfinedSmooth > 0.0) {
          // psi could be negative and needs to be handled next
          Sprime = (psi < layerThickness) ? specificYield / unconfinedSmooth * (layerThickness - psi) : 0.0;
          Sprime = (psi < (layerThickness - unconfinedSmooth)) ? specificYield : Sprime;
        } else {
          Sprime = (psi < layerThickness) ? specificYield : 0.0;
        }
        Seff(j, i) += Sprime;
      } else {
        Teff(j, i) = T(j, i);
      }
    }
  }
}

/** Calls getEffectiveAquiferProperties() if needed
 *
 * @param effectiveStorativity (1)
 * @param effectiveTransmissivity (m^2/s)
 * @param hydraulicTransmissivity (m^2/s)
 * @param hydraulicHead (m)
 * @param bedElevation (m)
 * @param bndMask (1)
 * @param layerThickness (m)
 * @param specificStorage (m^-1)
 * @param specificYield (1)
 * @param unconfinedSmooth confined - unconfined transition length, \f$ 0 \le d \le b \f$ (m)
 * @param disableUnconfined (boolean flag) See also CUASArgs::parseArgs()
 * @param doAnyChannel (boolean flag) See also CUASArgs::parseArgs()
 *
 */
inline void updateEffectiveAquiferProperties(PETScGrid &effectiveStorativity, PETScGrid &effectiveTransmissivity,
                                             PETScGrid const &hydraulicTransmissivity, PETScGrid const &hydraulicHead,
                                             PETScGrid const &bedElevation, PETScGrid const &bndMask,
                                             PetscScalar const layerThickness, PetscScalar const specificStorage,
                                             PetscScalar const specificYield, PetscScalar const unconfinedSmooth,
                                             bool const disableUnconfined, bool const doAnyChannel) {
  /*
   * Do not check input field compatibility here
   */

  if (disableUnconfined) {
    // We don't need to update Seff until we have some evolution in this aquifer property.
    if (doAnyChannel) {
      effectiveTransmissivity.copy(hydraulicTransmissivity);
    }
  } else {
    getEffectiveAquiferProperties(effectiveStorativity, effectiveTransmissivity, hydraulicTransmissivity, hydraulicHead,
                                  bedElevation, bndMask, layerThickness, specificStorage, specificYield,
                                  unconfinedSmooth);
  }
}

/** applies 3x3 smoothing kernel
 *
 * @param input
 * @param result
 */
inline void convolveStar11411(PETScGrid &input, PETScGrid &result) {
  if (!input.isCompatible(result)) {
    CUAS_ERROR("CUASKernels.h: convolveStar11411 was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  // the stencil is: stencil = np.array([[0, 1, 0],
  //                                    [1, 4, 1],
  //                                    [0, 1, 0]]) / 8

  auto &inp = input.getReadHandle();
  auto res = result.getWriteHandle();
  int ghostIndexCols = 1;
  int ghostIndexRows = 1;

  for (int i = 0; i < result.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < result.getLocalNumOfCols(); ++j) {
      res(i, j) = (4.0 * inp(ghostIndexRows, ghostIndexCols, true) + inp(ghostIndexRows, ghostIndexCols - 1, true) +
                   inp(ghostIndexRows, ghostIndexCols + 1, true) + inp(ghostIndexRows - 1, ghostIndexCols, true) +
                   inp(ghostIndexRows + 1, ghostIndexCols, true)) /
                  8.0;
      ++ghostIndexCols;
    }
    ghostIndexCols = 1;
    ++ghostIndexRows;
  }
}

/**
 *
 * @param inout
 * @param minimum
 * @param maximum
 */
inline void clamp(PETScGrid &inout, PetscScalar const minimum, PetscScalar const maximum) {
  auto res = inout.getWriteHandle();
  for (int j = 0; j < inout.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < inout.getLocalNumOfCols(); ++i) {
      res(j, i) = std::clamp(res(j, i), minimum, maximum);
    }
  }
}

/** explicit transmissivity update step within Tmin, Tmax limits
 *
 * Note, creep is usually considered as creep closure but here we use creep opening (opposite sign) to
 * consider all terms leading to transmissivity change with the same sign.
 *
 * T = T_n + (creep + melt + cavity) * dt_secs
 *
 * @param newHydraulicTransmissivity (m^2/s)
 * @param hydraulicTransmissivity (m^2/s)
 * @param aCreep transmissivity change due to creep opening (m^2/s^2)
 * @param aMelt transmissivity change due melt opening (m^2/s^2)
 * @param aCavity transmissivity change due cavity opening (m^2/s^2)
 * @param bndMask (1)
 * @param Tmin (m^2/s)
 * @param Tmax (m^2/s)
 * @param dt_secs time step length (s)
 */
inline void doChannels(PETScGrid &newHydraulicTransmissivity, PETScGrid const &hydraulicTransmissivity,
                       PETScGrid const &aCreep, PETScGrid const &aMelt, PETScGrid const &aCavity,
                       PETScGrid const &bndMask, PetscScalar const Tmin, PetscScalar const Tmax,
                       PetscScalar const dt_secs) {
  if (!aMelt.isCompatible(hydraulicTransmissivity) || !aCreep.isCompatible(hydraulicTransmissivity) ||
      !aCavity.isCompatible(hydraulicTransmissivity) ||
      !newHydraulicTransmissivity.isCompatible(hydraulicTransmissivity)) {
    CUAS_ERROR("CUASKernels.h: doChannels was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  // Update T with melt, creep and cavity
  {
    auto T = newHydraulicTransmissivity.getWriteHandle();
    auto &T_n = hydraulicTransmissivity.getReadHandle();
    auto &melt = aMelt.getReadHandle();      // dT/dt_melt
    auto &creep = aCreep.getReadHandle();    // dT/dt_creep
    auto &cavity = aCavity.getReadHandle();  // dT/dt_cavity
    auto &mask = bndMask.getReadHandle();
    for (int j = 0; j < newHydraulicTransmissivity.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < newHydraulicTransmissivity.getLocalNumOfCols(); ++i) {
        if (mask(j, i) == (PetscScalar)COMPUTE_FLAG) {
          // update with clamp
          auto newValue = T_n(j, i) + (creep(j, i) + melt(j, i) + cavity(j, i)) * dt_secs;
          T(j, i) = std::clamp(newValue, Tmin, Tmax);
        } else if (mask(j, i) == (PetscScalar)NOFLOW_FLAG) {
          T(j, i) = NOFLOW_VALUE;
        } else if (mask(j, i) == (PetscScalar)DIRICHLET_OCEAN_FLAG || mask(j, i) == (PetscScalar)DIRICHLET_LAKE_FLAG) {
          T(j, i) = Tmax;  // todo: better use Tocean here
        } else if (mask(j, i) == (PetscScalar)DIRICHLET_FLAG) {
          // Do nothing. If the head is as Dirichlet BC assume the transmissivity is also given.
        } else {
          CUAS_ERROR("CUASKernels.h: doChannels was called with unknown bndMask value {}. Exiting.", mask(j, i));
          exit(1);
        }
      }
    }
  }
}

/** Compute the hydraulic flux.
 *
 * It is assumed that the gradient mask has already be applied to the gradient!
 *
 * @param result (m^2/s)
 * @param gradientHeadSquared (1)
 * @param hydraulicTransmissivity (m^2/s)
 */
inline void getFluxMagnitude(PETScGrid &result, PETScGrid const &gradientHeadSquared,
                             PETScGrid const &hydraulicTransmissivity) {
  if (!result.isCompatible(gradientHeadSquared) || !result.isCompatible(hydraulicTransmissivity)) {
    CUAS_ERROR("CUASKernels.h: getFluxMagnitude() was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandle();
  auto &T = hydraulicTransmissivity.getReadHandle();
  auto &gradh2 = gradientHeadSquared.getReadHandle();
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      res(j, i) = T(j, i) * PetscSqrtScalar(gradh2(j, i));
    }
  }
}

/**
 * Simple driver for gradient2(...) that also aplies the gradient mask
 * @param result (1)
 * @param hydraulicHead (m)
 * @param dx (m)
 * @param gradMask (1)
 */
inline void getGradHeadSQR(PETScGrid &result, PETScGrid const &hydraulicHead, PetscScalar const dx,
                           PETScGrid const &gradMask) {
  if (!result.isCompatible(hydraulicHead) || !result.isCompatible(gradMask)) {
    CUAS_ERROR("CUASKernels.h: getGradHeadSQR was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  gradient2(result, hydraulicHead, dx);

  // apply mask
  auto res = result.getWriteHandle();
  auto &gmask = gradMask.getReadHandle();
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      if (gmask(j, i) == 1.0) {  // gradMask has inverted meaning
        res(j, i) = 0.0;
      }
    }
  }
}

}  // namespace CUAS

#endif
