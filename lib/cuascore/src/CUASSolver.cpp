/**
 * File: CUASSolver.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASSolver.h"

#include "CUASConstants.h"
#include "CUASKernels.h"
#include "specialgradient.h"
#include "systemmatrix.h"
#include "timing.h"

#include "Logger.h"
#include "PETScSolver.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace CUAS {

void CUASSolver::setup() {
  melt->setZero();
  creep->setZero();
  cavity->setZero();
  pEffective->setZero();
  gradHeadSquared->setZero();
  fluxMagnitude->setZero();

  rateFactorIce->setConst(args->flowConstant);         // todo: read 2d field from file (optional)
  basalVelocityIce->setConst(args->basalVelocityIce);  // todo: read 2d field from file (optional)

  // as in CUAS-python
  if (args->doAnyChannel) {
    currTransmissivity->setConst(args->Tinit);
  } else {
    currTransmissivity->setConst(args->layerThickness * args->conductivity);
  }

  {
    auto T = currTransmissivity->getWriteHandle();
    auto &mask = model->bndMask->getReadHandle();
    for (int i = 0; i < currTransmissivity->getLocalNumOfRows(); ++i) {
      for (int j = 0; j < currTransmissivity->getLocalNumOfCols(); ++j) {
        if (mask(i, j) == (PetscScalar)NOFLOW_FLAG) {
          T(i, j) = NOFLOW_VALUE;
        }
      }
    }
  }

  //
  // initialize the head
  //
  if (args->initialHead == "zero") {
    currHead->setZero();
  } else if (args->initialHead == "Nzero") {
    // np.maximum(np.maximum(helpers.pressure_to_head(pwater=pice, topg=topg), topg), 0.0)
    pressureToHead(*currHead, *model->pIce, *model->topg, args->layerThickness);
    {
      auto head = currHead->getWriteHandle();
      auto &topg = model->topg->getReadHandle();
      for (int j = 0; j < currHead->getLocalNumOfRows(); ++j) {
        for (int i = 0; i < currHead->getLocalNumOfCols(); ++i) {
          head(j, i) = std::max({head(j, i), topg(j, i), 0.0});
        }
      }
    }
  } else if (args->initialHead == "topg") {
    currHead->copy(*model->topg);
  } else {
    try {
      auto initialHeadValue = (PetscScalar)std::stod(args->initialHead);
      currHead->setConst(initialHeadValue);
    } catch (std::invalid_argument const &ex) {
      CUAS_ERROR("CUASSolver.cpp: solve(): args->initialHead invalid_argument: '" + args->initialHead + "'. Exiting.")
      exit(1);
    } catch (std::out_of_range const &ex) {
      CUAS_ERROR("CUASSolver.cpp: solve(): args->initialHead out_of_range: '" + args->initialHead + "'. Exiting.")
      exit(1);
    } catch (...) {
      CUAS_ERROR("CUASSolver.cpp: solve(): args->initialHead needs to be zero, Nzero, topg or valid number. Exiting.")
      exit(1);
    }
  }

  // local noFlowMask
  {
    auto cols = model->bndMask->getTotalNumOfCols();
    auto rows = model->bndMask->getTotalNumOfRows();
    auto grid = std::make_unique<PETScGrid>(cols, rows);
    auto glob = grid->getWriteHandle();
    auto &mask = model->bndMask->getReadHandle();

    for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
      for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
        // should we use true = 1.0 and false = 0.0 to match type of PetscScalar?
        glob(i, j) = (mask(i, j) == (PetscScalar)NOFLOW_FLAG);
      }
    }
    glob.setValues();
    binaryDilation(*gradMask, *grid);  // fixme: check gradMask for all Dirichlet BCs in Poisson Test
    gradMask->setRealBoundary(1.0);    // do not allow the gradient at the outermost grid points
  }

  // TODO: Time dependent forcing (Interpolierung)
  //
  // TODO: Time dependent tidal forcing

  // SET PETSc Options for direct solver MUMPS+PARDISO
  if (args->directSolver) {
#if (PETSC_HAVE_MUMPS == 1) && (PETSC_HAVE_PARMETIS == 1)
    CUAS_INFO("CUASSolver::setup(): Setup MUMPS+PARMETIS")
    // See: https://petsc.org/release/docs/manualpages/Mat/MATSOLVERMUMPS.html
    PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_type", "preonly");
    PetscOptionsSetValue(PETSC_NULLPTR, "-pc_type", "lu");
    PetscOptionsSetValue(PETSC_NULLPTR, "-pc_factor_mat_solver_type", "mumps");
    PetscOptionsSetValue(PETSC_NULLPTR, "-mat_mumps_icntl_14", "120");  // needed?
    PetscOptionsSetValue(PETSC_NULLPTR, "-mat_mumps_icntl_28", "2");    // parallel analysis
    PetscOptionsSetValue(PETSC_NULLPTR, "-mat_mumps_icntl_29", "2");    // parallel ordering: 2 = parmetis
    // if MUMPS factorization fails inside a KSP solve, give information about the failure
    PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_error_if_not_converged", nullptr);
#else
    CUAS_ERROR("CUASSolver.cpp: setup(): args->directSolver requires PETSc compiled with MUMPS and PARMETIS. Exiting.");
    exit(1);
#endif
  }
}

void CUASSolver::solve(std::vector<CUAS::timeSecs> &timeSteps) {
  prepare();

  beginSolverTiming();

  for (int timeStepIndex = 0; timeStepIndex < timeSteps.size(); ++timeStepIndex) {
    auto [currTime, dt] = getTimeStepInformation(timeSteps, timeStepIndex);

    // time dependent forcing
    auto &currentQ = model->Q->getCurrentQ(currTime);
    preComputation();

    storeData(currentQ, dt, timeSteps, timeStepIndex);

    updateHeadAndTransmissivity(dt, currentQ, timeSteps, timeStepIndex);
  }

  endSolverTiming();
}

std::pair<timeSecs, timeSecs> CUASSolver::getTimeStepInformation(std::vector<CUAS::timeSecs> const &timeSteps,
                                                                 int timeStepIndex) {
  timeSecs dt = 0;
  auto currTime = timeSteps[timeStepIndex];
  if (timeSteps.size() > 1) {
    if (timeStepIndex < timeSteps.size() - 1) {
      auto nextTime = timeSteps[timeStepIndex + 1];
      dt = nextTime - currTime;
    } else {
      // This is the last iteration only used to recompute the diagnostic fields for saving
      dt = -9999;  // Something stupid to let the solver crash if used
    }
  }

  return std::make_pair(currTime, dt);
}

void CUASSolver::prepare() {
  //
  // SOLVER PREPARATION
  //

  // only if we really want to solve something, the choice of theta is relevant
  if (args->timeSteppingTheta < 0.0 || args->timeSteppingTheta > 1.0) {
    CUAS_ERROR("CUASSolver::prepare(): timeSteppingTheta = {} is invalid for calling systemmatrix() . Exiting.",
               args->timeSteppingTheta)
    exit(1);
  }

  // after restart apply checks and set values consistent to cuas mask
  {
    auto trans = currTransmissivity->getWriteHandle();
    auto head = currHead->getWriteHandle();
    auto &mask = model->bndMask->getReadHandle();

    auto rows = currTransmissivity->getLocalNumOfRows();
    auto cols = currTransmissivity->getLocalNumOfCols();

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        if (mask(i, j) == (PetscScalar)NOFLOW_FLAG) {
          // no-flow must be obtained after restart
          trans(i, j) = NOFLOW_VALUE;
        } else if (mask(i, j) == (PetscScalar)DIRICHLET_OCEAN_FLAG) {
          // ensure proper ocean bc's after restart
          trans(i, j) = args->Tmax;
          head(i, j) = 0.0;
        }
      }
    }
  }
  dirichletValues->copy(*currHead);               // store initial values as dirichlet values for this run
  nextTransmissivity->copy(*currTransmissivity);  // otherwise invalid for bnd:mask == DIRICHLET_FLAG

  // initialized as confined only and thus Seff = S = Ss * b, and Teff = T
  Seff->setConst(args->layerThickness * args->specificStorage);
  Teff->copy(*currTransmissivity);

  if (args->verbose) {
    CUAS_INFO_RANK0("  S = Ss * b = {}", args->specificStorage * args->layerThickness)
    CUAS_INFO_RANK0("          Sy = {}", args->specificYield)
    CUAS_INFO_RANK0("        TINY = {}", TINY)
    CUAS_INFO_RANK0("NOFLOW_VALUE = {}", NOFLOW_VALUE)
    CUAS_INFO_RANK0("     RHO_ICE = {}", RHO_ICE)
    CUAS_INFO_RANK0("         SPY = {}", SPY)
    CUAS_INFO_RANK0("       input = {}", args->input)
    CUAS_INFO_RANK0("      output = {}", args->output)
    CUAS_INFO_RANK0("Starting CUASSolver")
  }
}

void CUASSolver::preComputation() {
  // TODO: basal velocity and rate factor fields are probably also time dependent
  // auto &currentRateFactor = XXX->getCurrent(currTime);
  // auto &currentBasalVelocity = XXX->getCurrent(currTime);

  // if (args->seaLevelForcing) {
  // TODO
  //}

  // Update Seff, Teff or both if needed.
  updateEffectiveAquiferProperties(*Seff, *Teff, *currTransmissivity, *currHead, *model->topg, *model->bndMask,
                                   args->layerThickness, args->specificStorage, args->specificYield, args->unconfSmooth,
                                   args->disableUnconfined, args->doAnyChannel);

  // calculate effective pressure for diagnostic output and probably for creep opening/closure
  headToEffectivePressure(*pEffective, *currHead, *model->topg, *model->pIce, args->layerThickness);

  // we need gradient of head for the flux and probably for melt opening
  getGradHeadSQR(*gradHeadSquared, *currHead, model->dx, *gradMask);

  // TODO extract?
  // compute channel opening terms to be ready for netcdf output
  if (args->doAnyChannel) {
    // creep
    if (args->doCreep) {
      computeCreepOpening(*creep, *rateFactorIce, *pEffective, *currTransmissivity);
    }
    // melt
    if (args->doMelt) {
      computeMeltOpening(*melt, args->roughnessFactor, args->conductivity, *currTransmissivity, *gradHeadSquared);
      if (!args->noSmoothMelt) {
        convolveStar11411(*melt, *tmpMelt);
        melt->copy(*tmpMelt);
      }
    }
    // cavity
    if (args->doCavity) {
      computeCavityOpening(*cavity, args->cavityBeta, args->conductivity, *basalVelocityIce);
    }
  }
}

void CUASSolver::storeData(PETScGrid const &currentQ, timeSecs dt, std::vector<CUAS::timeSecs> const &timeSteps,
                           int timeStepIndex) {
  //
  // STORE DATA, IF NEEDED
  //
  if (solutionHandler != nullptr) {
    OutputReason reason = solutionHandler->getOutputReason(timeStepIndex, (int)timeSteps.size(), args->saveEvery);

    auto currTime = timeSteps[timeStepIndex];

    if (reason != OutputReason::NONE) {
      if (!args->verboseSolver && args->verbose) {
        // show only if verboseSolver is off
        CUAS_INFO_RANK0("time({}/{}) = {} s, dt = {} s", timeStepIndex, timeSteps.size() - 1, currTime, dt)
      }
      // Process diagnostic variables for output only. We don't need them in every time step
      getFluxMagnitude(*fluxMagnitude, *gradHeadSquared, *Teff);  // was currTransmissivity for a very long time
      // ... more

      if (reason == OutputReason::INITIAL) {
        // storeInitialSetup() calls storeSolution() to store initial values for time dependent fields
        solutionHandler->storeInitialSetup(*currHead, *currTransmissivity, *model, *fluxMagnitude, *melt, *creep,
                                           *cavity, *pEffective, *Seff, *Teff, currentQ, *args);
      } else {
        solutionHandler->storeSolution(currTime, *currHead, *currTransmissivity, *model, *fluxMagnitude, *melt, *creep,
                                       *cavity, *pEffective, *Seff, *Teff, currentQ, eps, Teps);
      }
    }
  }
}

void CUASSolver::updateHeadAndTransmissivity(timeSecs dt, PETScGrid const &currentQ,
                                             std::vector<CUAS::timeSecs> const &timeSteps, int timeStepIndex) {
  //
  // WE NEED TO SOLVE AGAIN
  //
  if (timeStepIndex < timeSteps.size() - 1) {
    // Sanity check
    if (dt <= 0) {
      CUAS_ERROR("CUASSolver.cpp: solve(): dt={}s is invalid for calling systemmatrix() . Exiting.", dt)
      exit(1);
    }

    //
    // UPDATE HEAD
    //
    systemmatrix(*matA, *bGrid, *Seff, *Teff, model->dx, static_cast<double>(dt), args->timeSteppingTheta, *currHead,
                 currentQ, *dirichletValues, *model->bndMask, *globalIndicesBlocked);

    // solve the equation A*sol = b,
    auto res = PETScSolver::solve(*matA, *bGrid, *solGrid);
    nextHead->copyGlobal(*solGrid);
    // eps_head = np.max(np.abs(nextHead - currHead))
    eps = nextHead->getMaxAbsDiff(*currHead) / static_cast<double>(dt);
    // switch pointers
    currHead.swap(nextHead);

    //
    // UPDATE TRANSMISSIVITY, IF NEEDED
    //
    if (args->doAnyChannel) {
      doChannels(*nextTransmissivity, *currTransmissivity, *creep, *melt, *cavity, *model->bndMask, args->Tmin,
                 args->Tmax, static_cast<double>(dt));
      // eps_T = np.max(np.abs(T - currTransmissivity))
      Teps = nextTransmissivity->getMaxAbsDiff(*currTransmissivity) / static_cast<double>(dt);
      // switch pointers
      currTransmissivity.swap(nextTransmissivity);
    }

    if (args->verboseSolver) {
      // compute max digits needed for verboseSolver output format
      auto maxDigitsIndex = static_cast<int>(std::ceil(std::log10(timeSteps.size())));
      auto maxDigitsTime = static_cast<int>(std::ceil(std::log10(timeSteps.back())));

      auto currTime = timeSteps[timeStepIndex];
      // no fmt given for res.numberOfIterations, because this is only available from PETSc
      CUAS_INFO_RANK0("SOLVER: {:{}d} {:{}d} {:.5e} {:.5e} {:.5e} {} {}", timeStepIndex, maxDigitsIndex, currTime,
                      maxDigitsTime, eps, Teps, res.residualNorm, res.numberOfIterations, res.reason)
    }
  } else {
    // NO NEED TO UPDATE nextHead or T again
  }
}

}  // namespace CUAS
