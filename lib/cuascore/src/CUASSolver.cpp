/**
 * File: CUASSolver.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASSolver.h"

#include "CUASConstants.h"
#include "CUASKernels.h"
#include "CUASTimeIntegrator.h"
#include "systemmatrix.h"
#include "timing.h"

#include "Logger.h"
#include "PETScSolver.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace CUAS {

#define CUAS_MAX_TIMESTEP std::numeric_limits<typeof(timeSecs)>::max()

CUASSolver::CUASSolver(CUASModel *model, CUASArgs const *args, SolutionHandler *solutionHandler)
    : model(model), args(args), solutionHandler(solutionHandler), numOfCols(model->Ncols), numOfRows(model->Nrows) {
  /***** setup grids to work with *****/
  nextHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  currHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  nextTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  currTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);

  dirichletValues = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  // rate factor from flow law (ice rheology)
  rateFactorIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  // basal velocity of ice
  basalVelocityIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  melt = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  tmpMelt = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  creep = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  cavity = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  pEffective = std::make_unique<PETScGrid>(numOfCols, numOfRows);  // same as model->pIce
  fluxXDir = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  fluxYDir = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  fluxMagnitude = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  workSpace = std::make_unique<PETScGrid>(numOfCols, numOfRows);

  Seff = std::make_unique<PETScGrid>(numOfCols, numOfRows);  // effective Storativity
  Teff = std::make_unique<PETScGrid>(numOfCols, numOfRows);  // effective Transmissivity

  effTransEast = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  effTransWest = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  effTransNorth = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  effTransSouth = std::make_unique<PETScGrid>(numOfCols, numOfRows);

  /***** setup water source *****/
  waterSource = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  addWaterSource(dynamic_cast<WaterSource *>(model));

  /***** setup equation system *****/
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, numOfCols, numOfRows,
               PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr, &dm);
  DMSetFromOptions(dm);
  DMSetUp(dm);
  Mat petMat;
  DMCreateMatrix(dm, &petMat);
  matA = std::make_unique<PETScMatrix>(petMat);
  bGrid = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  solGrid = std::make_unique<PETScGrid>(numOfCols, numOfRows);

  /***** setup helper global indices *****/
  globalIndicesBlocked = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  fillGlobalIndicesBlocked(*globalIndicesBlocked);
}

void CUASSolver::setup() {
  melt->setZero();
  creep->setZero();
  cavity->setZero();
  pEffective->setZero();
  fluxMagnitude->setZero();
  waterSource->setZero();
  workSpace->setZero();

  rateFactorIce->setConst(args->flowConstant);         // todo: read 2d field from file (optional)
  basalVelocityIce->setConst(args->basalVelocityIce);  // todo: read 2d field from file (optional)

  currTransmissivity->setConst(args->Tinit);
  {
    // See also CUASSolver::solve() for args->applyRestartChecks
    auto trans = currTransmissivity->getWriteHandle();
    auto &mask = model->bndMask->getReadHandle();
    auto &topg = model->topg->getReadHandle();
    for (int row = 0; row < currHead->getLocalNumOfRows(); ++row) {
      for (int col = 0; col < currHead->getLocalNumOfCols(); ++col) {
        if (mask(row, col) == (PetscScalar)NOFLOW_FLAG) {
          // no-flow must be obtained after restart
          trans(row, col) = NOFLOW_VALUE;
        } else if (mask(row, col) == (PetscScalar)DIRICHLET_OCEAN_FLAG) {
          // ensure proper ocean bc's after restart from e.g. interpolated fields
          trans(row, col) = args->Tmax;
        } else if (mask(row, col) == (PetscScalar)DIRICHLET_LAKE_FLAG) {
          trans(row, col) = args->Tmax;
        }
      }
    }
  }

  //
  // initialize the head
  setInitialHeadFromArgs(*currHead, *model->bndMask, *model->topg, *model->pIce, *args, *workSpace);

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

void CUASSolver::solve(std::vector<timeSecs> &timeSteps) {
  prepare();

  CUASTimeIntegrator timeIntegrator(timeSteps);

  beginSolverTiming();

  bool continueSolverLoop = true;

  while (continueSolverLoop) {
    auto [currTime, dt] = timeIntegrator.getTimestepInformation(CUAS_MAX_TIMESTEP);

    updateWaterSource(currTime);

    preComputation();

    storeData(timeIntegrator);

    continueSolverLoop = updateHeadAndTransmissivity(timeIntegrator);

    timeIntegrator.finalizeTimestep(dt);
  }

  endSolverTiming();
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

void CUASSolver::updateWaterSource(timeSecs currTime) {
  if (waterSources.empty()) {
    // if we do not have a water source in water sources, water source is set to 0
    waterSource->setZero();
  } else if (waterSources.size() == 1) {
    // if we have one water source, we only check this water source
    // if it provides a water source, we copy the current field in water source
    // if it does not provide a water source, we set wat source to 0
    if (waterSources[0]->providesWaterSource()) {
      waterSource->copy(waterSources[0]->getCurrentWaterSource(currTime));
    } else {
      waterSource->setZero();
    }
  } else {
    // if multiple water sources are provided we iterate all of them
    // all water sources are checked if they provide a water source
    // if they do not provide a water source they are skipped
    waterSource->setZero();
    auto writeHandle = waterSource->getWriteHandle();
    for (auto w : waterSources) {
      if (!w->providesWaterSource()) {
        continue;
      }
      auto &curr = w->getCurrentWaterSource(currTime);
      auto &readHandle = curr.getReadHandle();
      for (int i = 0; i < waterSource->getLocalNumOfRows(); ++i) {
        for (int j = 0; j < waterSource->getLocalNumOfCols(); ++j) {
          writeHandle(i, j) = writeHandle(i, j) + readHandle(i, j);
        }
      }
    }
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
  if (args->enableUDS) {
    // experimental for unconfined flow mainly
    updateInterfaceTransmissivityUDS(*effTransEast, *effTransWest, *effTransNorth, *effTransSouth, *model->bndMask,
                                     *model->topg, *currHead, *Teff, args->thresholdThicknessUDS);
  } else {
    // this is the unmodified old scheme
    updateInterfaceTransmissivityCDS(*effTransEast, *effTransWest, *effTransNorth, *effTransSouth, *model->bndMask,
                                     *model->topg, *currHead, *Teff);
  }

  // Use unified version for the flux (CDS and UDS). So far only needed for output.
  getFlux(*fluxMagnitude, *fluxXDir, *fluxYDir, *model->bndMask, *currHead, *effTransEast, *effTransWest,
          *effTransNorth, *effTransSouth, model->dx);
  // ... more diagnostic variables ?

  // calculate effective pressure for diagnostic output and probably for creep opening/closure
  headToEffectivePressure(*pEffective, *currHead, *model->topg, *model->pIce, args->layerThickness);

  // compute channel opening terms to be ready for netcdf output
  if (args->doAnyChannel) {
    // creep
    if (args->doCreep) {
      computeCreepOpening(*creep, *rateFactorIce, *pEffective, *currTransmissivity, *model->bndMask);
    }
    // melt
    if (args->doMelt) {
      // use unified version of computeMeltOpening
      computeMeltOpening(*melt, args->roughnessFactor, args->conductivity, model->dx, *effTransEast, *effTransWest,
                         *effTransNorth, *effTransSouth, *currHead, *model->topg, *model->bndMask);
    }
    // cavity
    if (args->doCavity) {
      computeCavityOpening(*cavity, args->cavityBeta, args->conductivity, *basalVelocityIce, *model->bndMask);
    }
  }
}

void CUASSolver::storeData(CUASTimeIntegrator const &timeIntegrator) {
  //
  // STORE DATA, IF NEEDED
  //
  if (solutionHandler) {
    solutionHandler->storeData(*this, *model, *args, *waterSource, timeIntegrator);
  }
}

bool CUASSolver::updateHeadAndTransmissivity(CUASTimeIntegrator const &timeIntegrator) {
  auto dt = timeIntegrator.getCurrentDt();
  auto currTime = timeIntegrator.getCurrentTime();
  auto timeSteps = timeIntegrator.getTimesteps();
  auto timeStepIndex = timeIntegrator.getTimestepIndex();

  // Sanity check
  if (dt < 0) {
    CUAS_ERROR("CUASSolver.cpp: solve(): dt={}s is invalid for calling systemmatrix() . Exiting.", dt)
    exit(1);
  }

  // do we need to solve again?
  bool solve = dt != 0;

  if (solve) {
    //
    // UPDATE HEAD
    //
    //    systemmatrixDeprecated(*matA, *bGrid, *Seff, *Teff, model->dx, static_cast<double>(dt),
    //    args->timeSteppingTheta,
    //                           *currHead, *waterSource, *dirichletValues, *model->bndMask, *globalIndicesBlocked);
    systemmatrix(*matA, *bGrid, *Seff, *effTransEast, *effTransWest, *effTransNorth, *effTransSouth, model->dx,
                 static_cast<double>(dt), args->timeSteppingTheta, *currHead, *waterSource, *dirichletValues,
                 *model->bndMask, *globalIndicesBlocked);

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

      // no fmt given for res.numberOfIterations, because this is only available from PETSc
      CUAS_INFO_RANK0("SOLVER: {:{}d} {:{}d} {:.5e} {:.5e} {:.5e} {} {}", timeStepIndex, maxDigitsIndex, currTime,
                      maxDigitsTime, eps, Teps, res.residualNorm, res.numberOfIterations, res.reason)
    }
  } else {
    // NO NEED TO UPDATE nextHead or T again
  }

  return solve;
}

void CUASSolver::addWaterSource(WaterSource *newWaterSource) {
  if (std::find(waterSources.begin(), waterSources.end(), newWaterSource) != waterSources.end()) {
    CUAS_WARN("New water source was already added to the CUAS Solver before and is refused to be added again.")
    return;
  }
  waterSources.emplace_back(newWaterSource);
}

}  // namespace CUAS
