#include "CUASSolver.h"

#include "CUASConstants.h"
#include "CUASKernels.h"
#include "specialgradient.h"
#include "systemmatrix.h"

#include "Logger.h"
#include "PETScMatrix.h"
#include "PETScSolver.h"
#include "PETScVector.h"

#include <cmath>
#include <memory>

namespace CUAS {

void CUASSolver::setup() {
  Sp->setZero();
  K->setConst(args->conductivity);

  auto bt = args->layerThickness;

  nextTransmissivity->setConst(0.2);

  S->setConst(Ss * bt * args->Ssmulti);

  auto K_arr = K->getWriteHandle();
  auto T_arr = nextTransmissivity->getWriteHandle();
  auto &mask = model->bndMask->getReadHandle();

  auto rows = nextTransmissivity->getLocalNumOfRows();
  auto cols = nextTransmissivity->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)NOFLOW_FLAG) {
        K_arr(i, j) = NOFLOW_VALUE;
        T_arr(i, j) = NOFLOW_VALUE;
      }
    }
  }

  // needed for function call later on
  T_arr.setValues();

  currTransmissivity->copy(*nextTransmissivity);

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
    binaryDilation(*gradMask, *grid);
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

  // why is this not an arg?
  const int theta = 1;  // 1 means fully implicit, 0 means fully explicit, 0.5 is Crank-Nicholson

  PETScVector b(model->Nrows * model->Ncols);

  if (args->initialHead == "zero") {
    currHead->setZero();
  } else if (args->initialHead == "Nzero") {
    pressure2head(*currHead, *model->pIce, *model->topg, 0.0);
  } else if (args->initialHead == "topg") {
    currHead->copy(*model->topg);
  } else {
    CUAS_ERROR("CUASSolver.cpp: solve(): args->initialHead needs to be zero, Nzero or topg. Exiting.");
    exit(1);
  }

  // TODO: restart needs to read currTransmissivity (T_n), currHead (u_n)
  //  if(args->restart){
  //    if(args->verbose){
  //      std::cout << "Read restart from file " << args->restart << std::endl;
  //    }
  // TODO
  //  }

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

  // TODO!! solution init (part of saving to netcdf, see original-python main: 272-274)
  // melt, creep and Q are supposed to be part of the solution class.
  PETScGrid melt(model->Ncols, model->Nrows);
  PETScGrid creep(model->Ncols, model->Nrows);
  PETScGrid cavity(model->Ncols, model->Nrows);

  melt.setZero();
  creep.setZero();
  cavity.setZero();

  PetscScalar currTime = 0.0;
  clock_t t;
  if (rank == 0) {
    t = clock();
  }

  // start
  // creating grids outside of loop to save time
  PETScGrid Se(Sp->getTotalNumOfCols(), Sp->getTotalNumOfRows());
  int size = model->Ncols * model->Nrows;
  PETScMatrix A(size, size);
  PETScGrid Teff(nextTransmissivity->getTotalNumOfCols(), nextTransmissivity->getTotalNumOfRows());
  PETScGrid TeffPowTexp(nextTransmissivity->getTotalNumOfCols(), nextTransmissivity->getTotalNumOfRows());

  // initial condition
  if (solutionHandler != nullptr) {
    solutionHandler->storeInitialSetup(*nextHead, *nextTransmissivity, *model, melt, creep, cavity, *args);
  }

  if (args->verbose) {
    CUAS_INFO_RANK0("Starting CUASSolver");
  }

  for (int timeStepIndex = 1; timeStepIndex < timeSteps.size(); ++timeStepIndex) {
    auto currTime = timeSteps[timeStepIndex];
    auto prevTime = timeSteps[timeStepIndex - 1];
    auto dt = currTime - prevTime;

    // time dependent forcing
    auto &currentQ = model->Q->getCurrentQ(currTime);

    // if (args->seaLevelForcing) {
    // TODO
    //}

    // copying T to Teff is done inside if and else
    // enable_unconfied
    if (!args->disableUnconfined) {
      enableUnconfined(Teff, TeffPowTexp, *Sp, *currTransmissivity, *K, *model->topg, *currHead, args->Texp,
                       args->unconfSmooth, args->layerThickness);
    } else {
      Sp->setZero();
      calculateTeffPowTexp(Teff, TeffPowTexp, *nextTransmissivity, args->Texp);
    }

    calculateSeValues(Se, *Sp, *S);

    systemmatrix(A, b, model->Nrows, model->Ncols, Se, TeffPowTexp, model->dx, dt, theta, *currHead, currentQ,
                 *dirichletValues, *model->bndMask);
    // solve the equation A*sol = b
    PETScSolver::solve(A, b, *sol, args->verboseSolver && !args->directSolver);

    nextHead->setGlobalVecColMajor(*sol);

    if (args->doAllChannels) {
      doChannels(melt, creep, *currHead, *gradMask, *nextTransmissivity, *currTransmissivity, *model->pIce,
                 *model->topg, *K, *model->bndMask, cavity, args->flowConstant, args->Texp, args->roughnessFactor,
                 args->noSmoothMelt, args->cavityBeta, args->basalVelocityIce, args->Tmin, args->Tmax,
                 args->layerThickness, model->dx, dt);
    } else if (args->doAnyChannel) {
      doChannels(melt, creep, *currHead, *gradMask, *nextTransmissivity, *currTransmissivity, *model->pIce,
                 *model->topg, *K, *model->bndMask, cavity, args->flowConstant, args->Texp, args->roughnessFactor,
                 args->noSmoothMelt, args->cavityBeta, args->basalVelocityIce, args->Tmin, args->Tmax,
                 args->layerThickness, model->dx, dt, args->doMelt, args->doCreep, args->doCavity);
    } else {
      noChannels(melt, creep, cavity);
    }

    // switch pointers
    currHead.swap(nextHead);
    currTransmissivity.swap(nextTransmissivity);

    if (solutionHandler != nullptr &&
        ((timeStepIndex % args->saveEvery == 0) || (timeStepIndex == timeSteps.size() - 1))) {
      solutionHandler->storeSolution(currTime, *nextHead, *nextTransmissivity, *model, melt, creep, cavity);
    }
  }

  // end
  if (rank == 0) {
    t = clock() - t;
    CUAS_INFO_RANK0("CUASSolver.cpp: solve(): computation took: {} seconds.", ((float)t) / CLOCKS_PER_SEC);
  }
}

}  // namespace CUAS
