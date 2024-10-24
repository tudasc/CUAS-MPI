/**
 * File: CUASKernels.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_KERNELS_H
#define CUAS_KERNELS_H

#include "CUASArgs.h"
#include "CUASConstants.h"
#include "Logger.h"
#include "PETScGrid.h"

#include <algorithm>

namespace CUAS {

namespace  // unnamed namespace
{
// fixme: code duplication from  lib/cuascore/src/systemmatrix.cpp
inline PetscScalar harmonicmean(PetscScalar x1, PetscScalar x2) {  // can only be accessed in this file
  return 2.0 * x1 * x2 / (x1 + x2 + TINY);
}
}  // namespace

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
    CUAS_ERROR("CUASKernels.h: headToPressure was called with incompatible PETScGrids. Exiting.")
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
 * This uses water pressure at the ice interface (top of aquifer) and ice pressure.
 *
 * \f$ N = p_i - p_w(z_w = layerThickness) \f$
 *
 * @param result effective pressure (Pa)
 * @param hydraulicHead  (m)
 * @param bedElevation   (m)
 * @param icePressure    (Pa)
 * @param layerThickness (m)
 */
inline void headToEffectivePressure(PETScGrid &result, PETScGrid const &hydraulicHead, PETScGrid const &bedElevation,
                                    PETScGrid const &icePressure, PetscScalar const layerThickness = 0.0) {
  if (!result.isCompatible(hydraulicHead) || !result.isCompatible(bedElevation)) {
    CUAS_ERROR("CUASKernels.h: headToPressure was called with incompatible PETScGrids. Exiting.")
    exit(1);
  }

  auto res = result.getWriteHandle();  // model->pIce is without ghosts
  auto &head = hydraulicHead.getReadHandle();
  auto &topg = bedElevation.getReadHandle();
  auto &pi = icePressure.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      // see headToPressure()
      auto pw = (head(j, i) - topg(j, i) - layerThickness) * RHO_WATER * GRAVITY;
      res(j, i) = pi(j, i) - pw;
    }
  }
}

/** Convert water pressure to hydraulic head
 *
 *    zw = 0.0 -> water pressure at the base of the aquifer
 *    zw = layer thickness -> ... at top of aquifer
 *
 * @param result
 * @param waterPressure (Pa)
 * @param bedElevation (m)
 * @param zw \f$0 \le z_w \le layerThickness\f$
 *
 */
inline void pressureToHead(PETScGrid &result, PETScGrid const &waterPressure, PETScGrid const &bedElevation,
                           PetscScalar const zw = 0.0) {
  if (!result.isCompatible(waterPressure) || !result.isCompatible(bedElevation)) {
    CUAS_ERROR("CUASKernels.h: pressureToHead was called with incompatible PETScGrids. Exiting.")
    exit(1);
  }

  auto res = result.getWriteHandleGhost();
  auto &pw = waterPressure.getReadHandle();
  auto &topg = bedElevation.getReadHandle();

  constexpr double rgMultiplicator = 1.0 / (RHO_WATER * GRAVITY);

  for (int j = 0; j < result.getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalGhostNumOfCols(); ++i) {
      res(j, i) = pw(j, i, GHOSTED) * rgMultiplicator + topg(j, i, GHOSTED) + zw;
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
    CUAS_ERROR("CUASKernels.h: overburdenPressure was called with incompatible PETScGrids. Exiting.")
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
                                 PETScGrid const &basalVelocity, PETScGrid const &bndMask) {
  if (!result.isCompatible(basalVelocity) || !result.isCompatible(bndMask)) {
    CUAS_ERROR("CUASKernels.h: computeCavityOpening was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }
  auto res = result.getWriteHandle();
  auto &vb = basalVelocity.getReadHandle();
  auto &mask = bndMask.getReadHandle();

  const auto betaK = beta * hydraulicConductivity;
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      // only on active points, todo: Do we need the trivial else case ?
      if (mask(j, i) == COMPUTE_FLAG) {
        res(j, i) = betaK * vb(j, i);
      } else {
        res(j, i) = 0.0;
      }
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
                                PETScGrid const &hydraulicTransmissivity, PETScGrid const &bndMask) {
  if (!result.isCompatible(rateFactor) || !result.isCompatible(effectivePressure) ||
      !result.isCompatible(hydraulicTransmissivity) || !result.isCompatible(bndMask)) {
    CUAS_ERROR("CUASKernels.h: computeCreepOpening was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto res = result.getWriteHandle();
  auto &A = rateFactor.getReadHandle();
  auto &N = effectivePressure.getReadHandle();
  auto &T = hydraulicTransmissivity.getReadHandle();
  auto &mask = bndMask.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      // only on active points, todo: Do we need the trivial else case ?
      if (mask(j, i) == COMPUTE_FLAG) {
        res(j, i) = -2.0 * A(j, i) * pow(N(j, i) / 3.0, 3) * T(j, i);
      } else {
        res(j, i) = 0.0;
      }
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
    CUAS_ERROR("CUASKernels.h: binaryDilation was called with incompatible PETScGrids. Exiting.")
    exit(1);
  }

  // otherwise the output is not correctly updated,
  // when called a second time with smaller masked area
  output.setZero();

  auto &in = input.getReadHandle();
  auto out = output.getWriteHandle();

  int index_cols = 1;
  int index_rows = 1;

  auto iter_rows = output.getLocalNumOfRows();
  auto iter_cols = output.getLocalNumOfCols();

  for (int i = 0; i < iter_rows; ++i) {
    for (int j = 0; j < iter_cols; ++j) {
      // TODO we only have to evaluate this until we found one true
      if (0.0 != in(index_rows, index_cols, GHOSTED)) {
        out(i, j) = true;
      }
      if (0.0 != in(index_rows + 1, index_cols, GHOSTED)) {
        out(i, j) = true;
      }
      if (0.0 != in(index_rows - 1, index_cols, GHOSTED)) {
        out(i, j) = true;
      }
      if (0.0 != in(index_rows, index_cols + 1, GHOSTED)) {
        out(i, j) = true;
      }
      if (0.0 != in(index_rows, index_cols - 1, GHOSTED)) {
        out(i, j) = true;
      }
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
    CUAS_ERROR("CUASKernels.h: getEffectiveAquiferProperties was called with incompatible PETScGrids. Exiting.")
    exit(1);
  }

  if (unconfinedSmooth < 0.0 || unconfinedSmooth > layerThickness) {
    CUAS_ERROR(
        "CUASKernels.h: getEffectiveAquiferProperties was called with incompatible parameter unconfinedSmooth={}. "
        "Exiting.",
        unconfinedSmooth)
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
inline void convolveStar11411(PETScGrid const &input, PETScGrid &result) {
  if (!input.isCompatible(result)) {
    CUAS_ERROR("CUASKernels.h: convolveStar11411 was called with incompatible PETScGrids. Exiting.")
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
    CUAS_ERROR("CUASKernels.h: doChannels was called with incompatible PETScGrids. Exiting.")
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
          // T(j, i) = T_n(j, i); // fixme: here or in CUASSolver.cpp ? actual done in CUASSolver
        } else {
          CUAS_ERROR("CUASKernels.h: doChannels was called with unknown bndMask value {}. Exiting.", mask(j, i))
          exit(1);
        }
      }
    }
  }
}

// args.initialHead == 'Nopc':
//  See issm v4.23: trunk/src/c/classes/Loads/Friction.cpp:
//                      p_water = max(0.,rho_water*gravity*(sealevel-base));
//                      Neff = p_ice - p_water;
//
//  Note, ISSM uses the base = iceSurface - iceThickness and base = bed only for grounded ice
//
inline void setInitialHeadFromArgs(PETScGrid &result, PETScGrid const &bndMask, PETScGrid const &bedElevation,
                                   PETScGrid const &icePressure, CUASArgs const &args, PETScGrid &workSpace) {
  if (!result.isCompatible(bedElevation) || !result.isCompatible(bndMask) || !result.isCompatible(icePressure) ||
      !result.isCompatible(workSpace)) {
    CUAS_ERROR("{}:was called with incompatible PETScGrids. Exiting.", __PRETTY_FUNCTION__)
    exit(1);
  }

  // Use this to keep the implementation consistent with all the post-EGU runs for Greenland
  // This might be changed in the future
  bool positivePreserving = false;
  constexpr PetscScalar seaLevel = 0.0;  // sea level with respect to datum. This may change in future applications.

  //
  // fixme: I would like to write pressureToHead(result, 0.1 * icePressure, ...)
  //
  if (args.initialHead == "topg") {
    result.copy(bedElevation);
    // positivePreserving = true;  // if true -> this would be the same as Nopc
  } else if (args.initialHead == "low") {
    workSpace.copy(icePressure);  // store waterPressure = icePressure in workSpace
    workSpace.applyMultiplier(0.1);
    pressureToHead(result, workSpace, bedElevation, args.layerThickness);  // 10% of ice overburden pressure
    positivePreserving = true;
  } else if (args.initialHead == "mid") {
    workSpace.copy(icePressure);
    workSpace.applyMultiplier(0.5);
    pressureToHead(result, workSpace, bedElevation, args.layerThickness);  // 50% of ice overburden pressure
    positivePreserving = true;
  } else if (args.initialHead == "high") {
    workSpace.copy(icePressure);
    workSpace.applyMultiplier(0.9);
    pressureToHead(result, workSpace, bedElevation, args.layerThickness);  // 90% of ice overburden pressure
    positivePreserving = true;
  } else if (args.initialHead == "Nzero") {
    pressureToHead(result, icePressure, bedElevation, args.layerThickness);  // 100% of ice overburden pressure
    positivePreserving = true;
  } else if (args.initialHead == "Nopc") {
    // implements Wolovick et al., 2023, eq.6 (doi: 10.5194/tc-17-5027-2023)
    {
      auto waterPressure = workSpace.getWriteHandle();
      auto &topg = bedElevation.getReadHandle();
      for (int row = 0; row < workSpace.getLocalNumOfRows(); ++row) {
        for (int col = 0; col < workSpace.getLocalNumOfCols(); ++col) {
          waterPressure(row, col) = std::max(0.0, RHO_WATER * GRAVITY * (seaLevel - topg(row, col)));
        }
      }
    }
    pressureToHead(result, workSpace, bedElevation, args.layerThickness);  // Nopc
    positivePreserving = true;
  } else {
    // set from numeric value
    try {
      auto initialHeadValue = (PetscScalar)std::stod(args.initialHead);
      result.setConst(initialHeadValue);
    } catch (std::invalid_argument const &ex) {
      CUAS_ERROR("{}: args->initialHead invalid_argument: {}. Exiting.", __PRETTY_FUNCTION__, args.initialHead)
      exit(1);
    } catch (std::out_of_range const &ex) {
      CUAS_ERROR("{}: args->initialHead out_of_range: {}. Exiting.", __PRETTY_FUNCTION__, args.initialHead)
      exit(1);
    } catch (...) {
      CUAS_ERROR(
          "{}: args->initialHead needs to be 'Nzero', 'Nopc', 'low', 'mid', 'high', 'topg' or valid number. Exiting.",
          __PRETTY_FUNCTION__)
      exit(1);
    }
  }

  // positive preserving (head => 0, psi => 0)
  if (positivePreserving) {
    auto head = result.getWriteHandle();
    auto &topg = bedElevation.getReadHandle();
    for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
        head(row, col) = std::max({head(row, col), topg(row, col), 0.0});  // positive preserving head
      }
    }
  } else {
    CUAS_WARN_RANK0("{}: positivePreserving is false (off)", __PRETTY_FUNCTION__)
  }

  // duplicate some code from prepare() to make sure the ocean head is at sea level
  {
    auto head = result.getWriteHandle();
    auto &mask = bndMask.getReadHandle();
    for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
        if (mask(row, col) == (PetscScalar)DIRICHLET_OCEAN_FLAG) {
          // ensure proper ocean bc's after restart
          head(row, col) = seaLevel;
        }
      }
    }
  }
}

/** Compute the hydraulic flux using upwind difference scheme (UDS)
 *
 * We use compass notation to simplify indexing, thus "E" means to the East i+1 while "e" means
 * at the eastern interface at i+1/2 and in a similar fashion for N:North, S:South, W:West.
 * East = col + 1, West = col - 1
 * North = row + 1, South = row - 1
 *
 *
 * @param fluxMagnitude (m^2/s)
 * @param hydraulicHead (1)
 * @param hydraulicTransmissivity (m^2/s)
 */
inline void getFlux(PETScGrid &fluxMagnitude, PETScGrid &fluxXDir, PETScGrid &fluxYDir, PETScGrid const &bndMask,
                    PETScGrid const &hydraulicHead, PETScGrid const &effTransEast, PETScGrid const &effTransWest,
                    PETScGrid const &effTransNorth, PETScGrid const &effTransSouth, PetscScalar dx) {
  if (!fluxMagnitude.isCompatible(hydraulicHead) || !fluxMagnitude.isCompatible(fluxXDir) ||
      !fluxMagnitude.isCompatible(fluxYDir) || !fluxMagnitude.isCompatible(bndMask) ||
      !fluxMagnitude.isCompatible(effTransEast) || !fluxMagnitude.isCompatible(effTransWest) ||
      !fluxMagnitude.isCompatible(effTransNorth) || !fluxMagnitude.isCompatible(effTransSouth)) {
    CUAS_ERROR("CUASKernels.h: getFlux() was called with incompatible PETScGrids. Exiting.");
    exit(1);
  }

  auto fl = fluxMagnitude.getWriteHandle();
  auto flx = fluxXDir.getWriteHandle();
  auto fly = fluxYDir.getWriteHandle();

  auto &T_up_x = effTransEast.getReadHandle();   // interface (eff.) transmissivity in x-dir
  auto &T_dn_x = effTransWest.getReadHandle();   // interface (eff.) transmissivity in x-dir
  auto &T_up_y = effTransNorth.getReadHandle();  // interface (eff.) transmissivity in y-dir
  auto &T_dn_y = effTransSouth.getReadHandle();  // interface (eff.) transmissivity in y-dir

  auto &head = hydraulicHead.getReadHandle();
  auto &mask = bndMask.getReadHandle();

  auto const cornerX = hydraulicHead.getCornerX();
  auto const cornerY = hydraulicHead.getCornerY();
  auto cornerXGhost = hydraulicHead.getCornerXGhost();
  auto cornerYGhost = hydraulicHead.getCornerYGhost();

  const auto dx_inv = 1.0 / dx;

  for (int row = 0; row < hydraulicHead.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < hydraulicHead.getLocalNumOfCols(); ++col) {
      // as in lib/cuascore/src/systemmatrix.cpp
      auto j = row + (cornerY - cornerYGhost);
      auto i = col + (cornerX - cornerXGhost);

      if (mask(row, col) == (PetscScalar)COMPUTE_FLAG) {
        // some abbreviations
        auto head_E = head(j, i + 1, GHOSTED);
        auto head_W = head(j, i - 1, GHOSTED);
        auto head_N = head(j + 1, i, GHOSTED);
        auto head_S = head(j - 1, i, GHOSTED);
        auto head_P = head(j, i, GHOSTED);

        auto flux_e = -T_up_x(row, col) * (head_E - head_P) * dx_inv;
        auto flux_w = -T_dn_x(row, col) * (head_P - head_W) * dx_inv;
        auto flux_n = -T_up_y(row, col) * (head_N - head_P) * dx_inv;
        auto flux_s = -T_dn_y(row, col) * (head_P - head_S) * dx_inv;

        // Results:
        flx(row, col) = 0.5 * (flux_e + flux_w);
        fly(row, col) = 0.5 * (flux_n + flux_s);
        fl(row, col) = PetscSqrtScalar((flx(row, col) * flx(row, col)) + (fly(row, col) * fly(row, col)));
      }
    }
  }
}

/** Implements dT/dt_melt
 *
 */
inline void computeMeltOpening(PETScGrid &result, PetscScalar const roughnessFactor,
                               PetscScalar const hydraulicConductivity, PetscScalar const dx,
                               PETScGrid const &effTransEast, PETScGrid const &effTransWest,
                               PETScGrid const &effTransNorth, PETScGrid const &effTransSouth,
                               PETScGrid const &hydraulicHead, PETScGrid const &bedElevation, PETScGrid const &bndMask,
                               const PetscScalar alpha = 1.0, const PetscScalar beta = 2.0) {
  if (!result.isCompatible(effTransEast) || !result.isCompatible(hydraulicHead) || !result.isCompatible(bedElevation)) {
    CUAS_ERROR("CUASKernels.h: computeMeltOpening() was called with incompatible PETScGrids. Exiting.")
    exit(1);
  }
  if (dx <= 0.0) {
    CUAS_ERROR("CUASKernels.h: computeMeltOpening() was called with invalid grid spacing dx = {}. Exiting.", dx)
    exit(1);
  }
  if (alpha < 1.0) {
    CUAS_ERROR("CUASKernels.h: computeMeltOpening() was called with invalid alpha = {} < 1. Exiting.", alpha)
    exit(1);
  }

  // todo: We need a proof that this implementation is actually correct or find a reference
  if (beta != 2.0) {
    CUAS_WARN_RANK0("CUASKernels.h: computeMeltOpening() was called with beta = {} != 2.", beta)
  }

  const PetscScalar r_g_rhow_K = roughnessFactor * GRAVITY * RHO_WATER * hydraulicConductivity;
  const PetscScalar rhoi_L_inv = 1.0 / (RHO_ICE * LATENT_HEAT);

  auto res = result.getWriteHandle();
  auto const cornerX = result.getCornerX();
  auto const cornerY = result.getCornerY();
  auto cornerXGhost = result.getCornerXGhost();
  auto cornerYGhost = result.getCornerYGhost();

  auto &T_up_x = effTransEast.getReadHandle();   // interface (eff.) transmissivity in x-dir
  auto &T_dn_x = effTransWest.getReadHandle();   // interface (eff.) transmissivity in x-dir
  auto &T_up_y = effTransNorth.getReadHandle();  // interface (eff.) transmissivity in y-dir
  auto &T_dn_y = effTransSouth.getReadHandle();  // interface (eff.) transmissivity in y-dir

  auto &head = hydraulicHead.getReadHandle();
  auto &topg = bedElevation.getReadHandle();
  auto &mask = bndMask.getReadHandle();

  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      // We only want opening at active cuas locations.
      // We don't need to check if result.isOnRealBoundary(row, col), because we don't
      // have active points on the real boundary
      if (mask(row, col) != (PetscScalar)COMPUTE_FLAG || head(row, col) - topg(row, col) <= 0.0) {
        res(row, col) = 0.0;
      } else {
        // as in lib/cuascore/src/systemmatrix.cpp
        auto j = row + (cornerY - cornerYGhost);
        auto i = col + (cornerX - cornerXGhost);

        // some abbreviations
        auto head_E = head(j, i + 1, GHOSTED);
        auto head_W = head(j, i - 1, GHOSTED);
        auto head_N = head(j + 1, i, GHOSTED);
        auto head_S = head(j - 1, i, GHOSTED);
        auto head_P = head(j, i, GHOSTED);

        auto dhdx_up_x = (head_E - head_P) / dx;
        auto dhdx_dn_x = (head_P - head_W) / dx;
        auto dhdy_up_y = (head_N - head_P) / dx;
        auto dhdy_dn_y = (head_P - head_S) / dx;

        // extended version of gradient2() from specialgradient.cpp. See Beyer et al., 2018 Eq. B9
        // if beta is 2 we get (dhdx*dhdx)^(beta/2) = (dhdx^2)^(2/2) = dhdx^2 as expected
        res(row, col) = r_g_rhow_K * rhoi_L_inv * 0.5 *
                        (pow(T_up_x(row, col), alpha) * pow(dhdx_up_x * dhdx_up_x, beta / 2.0) +
                         pow(T_dn_x(row, col), alpha) * pow(dhdx_dn_x * dhdx_dn_x, beta / 2.0) +
                         pow(T_up_y(row, col), alpha) * pow(dhdy_up_y * dhdy_up_y, beta / 2.0) +
                         pow(T_dn_y(row, col), alpha) * pow(dhdy_dn_y * dhdy_dn_y, beta / 2.0));
      }  // end mask
    }
  }
}

inline void updateInterfaceTransmissivityCDS(PETScGrid &effTransEast, PETScGrid &effTransWest, PETScGrid &effTransNorth,
                                             PETScGrid &effTransSouth, const PETScGrid &bndMask,
                                             const PETScGrid &bedElevation, const PETScGrid &hydraulicHead,
                                             PETScGrid const &effectiveTransmissivity) {
  if (!bndMask.isCompatible(hydraulicHead) || !bndMask.isCompatible(bedElevation) ||
      !bndMask.isCompatible(effectiveTransmissivity) || !bndMask.isCompatible(effTransEast) ||
      !bndMask.isCompatible(effTransWest) || !bndMask.isCompatible(effTransNorth) ||
      !bndMask.isCompatible(effTransSouth)) {
    CUAS_ERROR("CUASKernels.h: updateInterfaceTransmissivityCDS was called with incompatible PETScGrids. Exiting.")
    exit(1);
  }

  // read from
  auto &mask = bndMask.getReadHandle();
  auto &transmissivity = effectiveTransmissivity.getReadHandle();
  // auto &head = hydraulicHead.getReadHandle();
  // auto &topg = bedElevation.getReadHandle();

  // write to
  auto T_e = effTransEast.getWriteHandle();   // interface (eff.) transmissivity in pos. x-dir
  auto T_w = effTransWest.getWriteHandle();   // interface (eff.) transmissivity in neg. x-dir
  auto T_n = effTransNorth.getWriteHandle();  // interface (eff.) transmissivity in pos. y-dir
  auto T_s = effTransSouth.getWriteHandle();  // interface (eff.) transmissivity in neg. y-dir

  auto cornerX = bndMask.getCornerX();
  auto cornerY = bndMask.getCornerY();
  auto cornerXGhost = bndMask.getCornerXGhost();
  auto cornerYGhost = bndMask.getCornerYGhost();

  // we take localnumOfCols, so that we can iterate over the normal boundaries
  for (int row = 0; row < bndMask.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < bndMask.getLocalNumOfCols(); ++col) {
      if (mask(row, col) == (PetscScalar)COMPUTE_FLAG) {
        // needed for transmissivity and head in order to use LocalArray rather than GlobalArray
        auto j = row + (cornerY - cornerYGhost);
        auto i = col + (cornerX - cornerXGhost);

        // same for transmissivity aka effective transmissivity
        auto T_E = transmissivity(j, i + 1, GHOSTED);
        auto T_W = transmissivity(j, i - 1, GHOSTED);
        auto T_N = transmissivity(j + 1, i, GHOSTED);
        auto T_S = transmissivity(j - 1, i, GHOSTED);
        auto T_P = transmissivity(j, i, GHOSTED);

        T_e(row, col) = harmonicmean(T_P, T_E);
        T_w(row, col) = harmonicmean(T_P, T_W);
        T_s(row, col) = harmonicmean(T_P, T_S);
        T_n(row, col) = harmonicmean(T_P, T_N);
      }
    }
  }
}

inline void updateInterfaceTransmissivityUDS(PETScGrid &effTransEast, PETScGrid &effTransWest, PETScGrid &effTransNorth,
                                             PETScGrid &effTransSouth, const PETScGrid &bndMask,
                                             const PETScGrid &bedElevation, const PETScGrid &hydraulicHead,
                                             PETScGrid const &effectiveTransmissivity, PetscScalar threshold) {
  if (!bndMask.isCompatible(hydraulicHead) || !bndMask.isCompatible(bedElevation) ||
      !bndMask.isCompatible(effectiveTransmissivity) || !bndMask.isCompatible(effTransEast) ||
      !bndMask.isCompatible(effTransWest) || !bndMask.isCompatible(effTransNorth) ||
      !bndMask.isCompatible(effTransSouth)) {
    CUAS_ERROR("CUASKernels.h: updateInterfaceTransmissivityUDS was called with incompatible PETScGrids. Exiting.")
    exit(1);
  }

  if (threshold < 0.0) {
    CUAS_ERROR("CUASSolver::updateInterfaceTransmissivityUDS() called with threshold = {} < 0. Exiting.", threshold);
    exit(1);
  }

  // read from
  auto &mask = bndMask.getReadHandle();
  auto &transmissivity = effectiveTransmissivity.getReadHandle();
  auto &head = hydraulicHead.getReadHandle();
  auto &topg = bedElevation.getReadHandle();

  // write to
  auto T_up_x = effTransEast.getWriteHandle();   // interface (eff.) transmissivity in pos. x-dir
  auto T_dn_x = effTransWest.getWriteHandle();   // interface (eff.) transmissivity in neg. x-dir
  auto T_up_y = effTransNorth.getWriteHandle();  // interface (eff.) transmissivity in pos. y-dir
  auto T_dn_y = effTransSouth.getWriteHandle();  // interface (eff.) transmissivity in neg. y-dir

  auto cornerX = bndMask.getCornerX();
  auto cornerY = bndMask.getCornerY();
  auto cornerXGhost = bndMask.getCornerXGhost();
  auto cornerYGhost = bndMask.getCornerYGhost();

  // PetscScalar T_up_x, T_dn_x;  // interface (eff.) transmissivity in x-dir
  // PetscScalar T_up_y, T_dn_y;  // interface (eff.) transmissivity in y-dir

  PetscScalar T_e, T_w, T_s, T_n;

  // we take localnumOfCols, so that we can iterate over the normal boundaries
  for (int row = 0; row < bndMask.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < bndMask.getLocalNumOfCols(); ++col) {
      if (mask(row, col) == (PetscScalar)COMPUTE_FLAG) {
        // needed for transmissivity and head in order to use LocalArray rather than GlobalArray
        auto j = row + (cornerY - cornerYGhost);
        auto i = col + (cornerX - cornerXGhost);

        // some abbreviations
        auto head_E = head(j, i + 1, GHOSTED);
        auto head_W = head(j, i - 1, GHOSTED);
        auto head_N = head(j + 1, i, GHOSTED);
        auto head_S = head(j - 1, i, GHOSTED);
        auto head_P = head(j, i, GHOSTED);

        // we also need psi = head - topg for all
        auto psi_E = head_E - topg(j, i + 1, GHOSTED);
        auto psi_W = head_W - topg(j, i - 1, GHOSTED);
        auto psi_N = head_N - topg(j + 1, i, GHOSTED);
        auto psi_S = head_S - topg(j - 1, i, GHOSTED);
        auto psi_P = head_P - topg(j, i, GHOSTED);

        // same for transmissivity aka effective transmissivity
        auto T_E = transmissivity(j, i + 1, GHOSTED);
        auto T_W = transmissivity(j, i - 1, GHOSTED);
        auto T_N = transmissivity(j + 1, i, GHOSTED);
        auto T_S = transmissivity(j - 1, i, GHOSTED);
        auto T_P = transmissivity(j, i, GHOSTED);

        //
        // eastern eff. boundary transmissivity
        //
        if (mask(j, i + 1, GHOSTED) == NOFLOW_FLAG) {
          // This was also harmonic mean in the CDS scheme
          T_e = NOFLOW_VALUE;  // fixme: should we use zero instead?
        } else {
          if (psi_E > threshold && psi_P > threshold) {
            T_e = harmonicmean(T_P, T_E);
          } else {
            // arithmetic mean or better harmonic mean?
            T_e = 0.5 * (T_P + T_E);  // Teff(psi_e) at (i+1/2, j)
            T_e = 1e-1 * T_P;         // 1% of the center
          }
        }
        //
        // western eff. boundary transmissivity
        //
        if (mask(j, i - 1, GHOSTED) == NOFLOW_FLAG) {
          // This was also harmonic mean in the CDS scheme
          T_w = NOFLOW_VALUE;  // fixme: should we use zero instead?
        } else {
          if (psi_W > threshold && psi_P > threshold) {
            T_w = harmonicmean(T_P, T_W);
          } else {
            // arithmetic mean or better harmonic mean
            T_w = 0.5 * (T_P + T_W);  // Teff(psi_w) at (i-1/2, j)
            T_w = 1e-1 * T_P;
          }
        }
        //
        // northern eff. boundary transmissivity
        //
        if (mask(j + 1, i, GHOSTED) == NOFLOW_FLAG) {
          // This was also harmonic mean in the CDS scheme
          T_n = NOFLOW_VALUE;  // fixme: should we use zero instead?
        } else {
          if (psi_N > threshold && psi_P > threshold) {
            T_n = harmonicmean(T_P, T_N);
          } else {
            // arithmetic mean or better harmonic mean?
            T_n = 0.5 * (T_P + T_N);  // Teff(psi_w) at (i, j+1/2)
            T_n = 1e-1 * T_P;
          }
        }
        //
        // southern eff. boundary transmissivity
        //
        if (mask(j - 1, i, GHOSTED) == NOFLOW_FLAG) {
          // This was also harmonic mean in the CDS scheme
          T_s = NOFLOW_VALUE;  // fixme: should we use zero instead?
        } else {
          if (psi_S > threshold && psi_P > threshold) {
            T_s = harmonicmean(T_P, T_S);
          } else {
            // arithmetic mean or better harmonic mean
            T_s = 0.5 * (T_P + T_S);  // Teff(psi_s) at (i, j-1/2)
            T_s = 1e-1 * T_P;
          }
        }

        // no flux out of dry cells
        T_up_x(row, col) = (psi_P <= 0.0 && head_P > head_E && psi_P <= psi_E) ? NOFLOW_VALUE : T_e;
        T_dn_x(row, col) = (psi_P <= 0.0 && head_P > head_W && psi_P <= psi_W) ? NOFLOW_VALUE : T_w;
        T_up_y(row, col) = (psi_P <= 0.0 && head_P > head_N && psi_P <= psi_N) ? NOFLOW_VALUE : T_n;
        T_dn_y(row, col) = (psi_P <= 0.0 && head_P > head_S && psi_P <= psi_S) ? NOFLOW_VALUE : T_s;
      }
    }
  }
}

}  // namespace CUAS

#endif
