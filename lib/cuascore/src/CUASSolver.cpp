#include "CUASSolver.h"

#include "CUASConstants.h"
#include "CUASKernels.h"
#include "specialgradient.h"
#include "systemmatrix.h"

#include "Logger.h"
#include "PETScSolver.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace CUAS {

void CUASSolver::setup() {
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
    pressureToHead(*currHead, *model->pIce, *model->topg);
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
      CUAS_ERROR("CUASSolver.cpp: solve(): args->initialHead invalid_argument: '" + args->initialHead + "'. Exiting.");
      exit(1);
    } catch (std::out_of_range const &ex) {
      CUAS_ERROR("CUASSolver.cpp: solve(): args->initialHead out_of_range: '" + args->initialHead + "'. Exiting.");
      exit(1);
    } catch (...) {
      CUAS_ERROR("CUASSolver.cpp: solve(): args->initialHead needs to be zero, Nzero, topg or valid number. Exiting.");
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
        glob(i, j) = (mask(i, j) == (PetscScalar)NOFLOW_FLAG) ? true : false;
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
    CUAS_INFO("CUASSolver::setup(): Setup MUMPS+PARMETIS");
    // See: https://petsc.org/release/docs/manualpages/Mat/MATSOLVERMUMPS.html
    PetscOptionsSetValue(PETSC_NULL, "-ksp_type", "preonly");
    PetscOptionsSetValue(PETSC_NULL, "-pc_type", "lu");
    PetscOptionsSetValue(PETSC_NULL, "-pc_factor_mat_solver_type", "mumps");
    PetscOptionsSetValue(PETSC_NULL, "-mat_mumps_icntl_14", "120");  // needed?
    PetscOptionsSetValue(PETSC_NULL, "-mat_mumps_icntl_28", "2");    // parallel analysis
    PetscOptionsSetValue(PETSC_NULL, "-mat_mumps_icntl_29", "2");    // parallel ordering: 2 = parmetis
    // if MUMPS factorization fails inside a KSP solve, give information about the failure
    PetscOptionsSetValue(PETSC_NULL, "-ksp_error_if_not_converged", nullptr);
#else
    CUAS_ERROR("CUASSolver.cpp: setup(): args->directSolver requires PETSc compiled with MUMPS and PARMETIS. Exiting.");
    exit(1);
#endif
  }
}

void CUASSolver::solve(std::vector<CUAS::timeSecs> &timeSteps) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //
  // SOLVER PREPARATION
  //

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
  dirichletValues->copy(*currHead);  // store initial values as dirichlet values for this run

  PETScGrid melt(model->Ncols, model->Nrows);
  PETScGrid creep(model->Ncols, model->Nrows);
  PETScGrid cavity(model->Ncols, model->Nrows);
  PETScGrid pEffective(model->Ncols, model->Nrows);  // same as model->pIce
  PETScGrid gradHeadSquared(model->Ncols, model->Nrows);
  PETScGrid fluxMagnitude(model->Ncols, model->Nrows);

  PETScGrid Seff(model->Ncols, model->Nrows);  // effective Storativity
  PETScGrid Teff(model->Ncols, model->Nrows);  // effectife Transmissivity

  melt.setZero();
  creep.setZero();
  cavity.setZero();
  pEffective.setZero();
  gradHeadSquared.setZero();
  fluxMagnitude.setZero();

  // initialized as confined only and thus Seff = S = Ss * b, and Teff = T
  Seff.setConst(args->layerThickness * args->specificStorage);
  Teff.copy(*currTransmissivity);

  clock_t t;
  if (rank == 0) {
    t = clock();
  }

  // start
  // creating grids outside of loop to save time
  int size = model->Ncols * model->Nrows;
  PetscScalar eps = 0.0;
  PetscScalar Teps = 0.0;

  // Why is this not an arg?
  // tkleiner: Because other values are not working! This indicates that the choice of initial
  // conditions is not appropriate.
  const int theta = 1;  // 1 means fully implicit, 0 means fully explicit, 0.5 is Crank-Nicholson

  if (args->verbose) {
    CUAS_INFO_RANK0("  S = Ss * b = {}", args->specificStorage * args->layerThickness);
    CUAS_INFO_RANK0("          Sy = {}", args->specificYield);
    CUAS_INFO_RANK0("        TINY = {}", TINY);
    CUAS_INFO_RANK0("NOFLOW_VALUE = {}", NOFLOW_VALUE);
    CUAS_INFO_RANK0("     RHO_ICE = {}", RHO_ICE);
    CUAS_INFO_RANK0("         SPY = {}", SPY);
    CUAS_INFO_RANK0("       input = {}", args->input);
    CUAS_INFO_RANK0("      output = {}", args->output);
    CUAS_INFO_RANK0("Starting CUASSolver");
  }

  timeSecs dt = 0;
  for (int timeStepIndex = 0; timeStepIndex < timeSteps.size(); ++timeStepIndex) {
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

    // time dependent forcing
    auto &currentQ = model->Q->getCurrentQ(currTime);

    // TODO: basal velocity and rate factor fields are probably also time dependent
    // auto &currentRateFactor = XXX->getCurrent(currTime);
    // auto &currentBasalVelocity = XXX->getCurrent(currTime);

    // if (args->seaLevelForcing) {
    // TODO
    //}

    if (args->doAnyChannel) {
      if (args->disableUnconfined) {
        // We don't need to update Seff until we have also some kind of evolution in this parameter.
        Teff.copy(*currTransmissivity);
      } else {
        getEffectiveAquiferProperties(Seff, Teff, *currTransmissivity, *currHead, *model->topg, *model->bndMask,
                                      args->layerThickness, args->specificStorage, args->specificYield,
                                      args->unconfSmooth);
      }
    }

    // calculate effective pressure for diagnostic output and probably for creep opening/closure
    headToEffectivePressure(pEffective, *currHead, *model->topg, *model->pIce, args->layerThickness);

    // we need gradient of head for the flux and probably for melt opening
    getGradHeadSQR(gradHeadSquared, *currHead, model->dx, *gradMask);

    // compute channel opening terms to be ready for netcdf output
    if (args->doAnyChannel) {
      // creep
      if (args->doCreep) {
        computeCreepOpening(creep, *rateFactorIce, pEffective, *currTransmissivity);
      }
      // melt
      if (args->doMelt) {
        computeMeltOpening(melt, args->roughnessFactor, args->conductivity, *currTransmissivity, gradHeadSquared);
        if (!args->noSmoothMelt) {
          PETScGrid tmp(melt.getTotalNumOfCols(), melt.getTotalNumOfRows());
          convolveStar11411(melt, tmp);
          melt.copy(tmp);
        }
      }
      // cavity
      if (args->doCavity) {
        computeCavityOpening(cavity, args->cavityBeta, args->conductivity, *basalVelocityIce);
      }
    }

    //
    // STORE DATA, IF NEEDED
    //
    if (solutionHandler != nullptr) {
      OutputReason reason = solutionHandler->getOutputReason(timeStepIndex, (int)timeSteps.size(), args->saveEvery);
      if (reason != OutputReason::NONE) {
        if (args->verbose) {
          CUAS_INFO_RANK0("time({}/{}) = {} s, dt = {} s", timeStepIndex, timeSteps.size() - 1, currTime, dt);
        }
        // Process diagnostic variables for output only. We don't need them in every time step
        getFluxMagnitude(fluxMagnitude, gradHeadSquared, Teff);  // was currTransmissivity for a very long time
        // ... more

        if (reason == OutputReason::INITIAL) {
          // storeInitialSetup() calls storeSolution() to store initial values for time dependent fields
          solutionHandler->storeInitialSetup(*currHead, *currTransmissivity, *model, fluxMagnitude, melt, creep, cavity,
                                             pEffective, Seff, Teff, currentQ, *args);
        } else {
          solutionHandler->storeSolution(currTime, *currHead, *currTransmissivity, *model, fluxMagnitude, melt, creep,
                                         cavity, pEffective, Seff, Teff, currentQ, eps, Teps);
        }
      }
    }  // solutionHandler

    //
    // WE NEED TO SOLVE AGAIN
    //
    if (timeStepIndex < timeSteps.size() - 1) {
      // Sanity check
      if (dt <= 0) {
        CUAS_ERROR("CUASSolver.cpp: solve(): dt={}s is invalid for calling systemmatrix() . Exiting.", dt);
        exit(1);
      }

      //
      // UPDATE HEAD
      //
      systemmatrix(*matA, *bGrid, Seff, Teff, model->dx, dt, theta, *currHead, currentQ, *dirichletValues,
                   *model->bndMask, *globalIndicesBlocked);

      // solve the equation A*sol = b,
      // todo: - return number of iterations and rnorm for logging
      PETScSolver::solve(*matA, *bGrid, *solGrid, args->verboseSolver && !args->directSolver);
      nextHead->copyGlobal(*solGrid);
      // todo: eps_head = np.max(np.abs(nextHead - currHead))
      eps = nextHead->getMaxAbsDiff(*currHead) / dt;
      // switch pointers
      currHead.swap(nextHead);

      //
      // UPDATE TRANSMISSIVITY, IF NEEDED
      //
      if (args->doAnyChannel) {
        doChannels(*nextTransmissivity, *currTransmissivity, creep, melt, cavity, *model->bndMask, args->Tmin,
                   args->Tmax, dt);
        // todo: eps_T = np.max(np.abs(T - currTransmissivity))
        Teps = nextTransmissivity->getMaxAbsDiff(*currTransmissivity) / dt;
        // switch pointers
        currTransmissivity.swap(nextTransmissivity);
      }

    } else {
      // NO NEED TO UPDATE nextHead or T again
    }
  }

  // end
  if (rank == 0) {
    t = clock() - t;
    CUAS_INFO_RANK0("CUASSolver.cpp: solve(): computation took: {} seconds.", ((float)t) / CLOCKS_PER_SEC);
  }
}

}  // namespace CUAS
