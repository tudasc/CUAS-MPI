#include "CUASSolver.h"

#include "CUASConstants.h"
#include "CUASKernels.h"
#include "savetimestep.h"
#include "specialgradient.h"
#include "systemmatrix.h"

#include "Logger.h"
#include "PETScMat.h"
#include "PETScSolver.h"
#include "PETScVec.h"

#include <cmath>
#include <memory>

namespace CUAS {

void CUASSolver::setup() {
  Sp->setZero();
  K->setConst(args->conductivity);

  auto bt = args->layerThickness;

  T->setConst(0.2);

  S->setConst(Ss * bt * args->Ssmulti);

  auto global_mask = noFlowMask->getWriteHandle();
  auto K_arr = K->getWriteHandle();
  auto T_arr = T->getWriteHandle();
  auto &mask = model->bndMask->getReadHandle();

  auto rows = noFlowMask->getLocalNumOfRows();
  auto cols = noFlowMask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)NOFLOW_FLAG) {
        global_mask(i, j) = true;
        K_arr(i, j) = NOFLOW_VALUE;
        T_arr(i, j) = NOFLOW_VALUE;
      } else {
        global_mask(i, j) = false;
      }
    }
  }

  // needed for function call later on
  global_mask.setValues();
  T_arr.setValues();

  T_n->copy(*T);

  auto global_dir_mask = dirichletMask->getWriteHandle();

  rows = dirichletMask->getLocalNumOfRows();
  cols = dirichletMask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)DIRICHLET_FLAG) {
        global_dir_mask(i, j) = true;
      } else {
        global_dir_mask(i, j) = false;
      }
    }
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)DIRICHLET_LAKE_FLAG) {
        global_dir_mask(i, j) = true;
      }
    }
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      global_dir_mask(i, j) = global_dir_mask(i, j) || global_mask(i, j);
    }
  }

  pressure2head(*dirichletValues, *model->pIce, *model->topg, 0.0);

  auto global_sea_mask = seaLevelForcingMask->getWriteHandle();

  rows = seaLevelForcingMask->getLocalNumOfRows();
  cols = seaLevelForcingMask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)DIRICHLET_FLAG) {
        global_sea_mask(i, j) = true;
      } else {
        global_sea_mask(i, j) = false;
      }
    }
  }

  binaryDialation(*gradMask, *noFlowMask);

  // TODO: Time dependent forcing (Interpolierung)
  //
  // TODO: Time dependent tidal forcing
}

void CUASSolver::solve(std::unique_ptr<PETScGrid> &u, std::unique_ptr<PETScGrid> &u_n, int const Nt,
                       PetscScalar const totaltime_secs, PetscScalar const dt_secs) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //
  // SOLVER PREPARATION
  //
  PETScVec sol(model->Nrows * model->Ncols);

  PetscScalar It[Nt + 1];
  for (int i = 0; i < Nt + 1; ++i) {
    It[i] = i;
  }

  // auto u = std::make_unique<PETScGrid>(model->Ncols, model->Nrows);    // unknown u at new time level
  // auto u_n = std::make_unique<PETScGrid>(model->Ncols, model->Nrows);  // u at the previous time level

  // why is this not an arg?
  const int theta = 1;  // 1 means fully implicit, 0 means fully explicit, 0.5 is Crank-Nicholson

  PETScVec b(model->Nrows * model->Ncols);

  if (args->initialHead == "zero") {
    u_n->setZero();
  } else if (args->initialHead == "Nzero") {
    pressure2head(*u_n, *model->pIce, *model->topg, 0.0);
  } else if (args->initialHead == "topg") {
    u_n->copy(*model->topg);
  } else {
    Logger::instance().error("CUASSolver.cpp: solve(): args->initialHead needs to be zero, Nzero or topg. Exiting.");
    exit(1);
  }

  //  if(args->restart){
  //    if(args->verbose){
  //      std::cout << "Read restart from file " << args->restart << std::endl;
  //    }
  // TODO
  //  }

  // save timedependent values
  PetscScalar Ntsaved;
  auto Itlength = sizeof(It) / sizeof(It[0]);
  if (Itlength % args->saveEvery > 0) {
    Ntsaved = ceil(Itlength / args->saveEvery);
  }

  if (args->verbose) {
    Logger::instance().info("CUASSolver.cpp: solve(): runtime = {}, time step = {}, Ntsaved = {} for saveEvery = {}.",
                            totaltime_secs, dt_secs, Ntsaved, args->saveEvery);
  }

  // TODO!! solution init (part of saving to netcdf, see original-python main: 272-274)
  // melt, creep and Q are supposed to be part of the solution class.
  PETScGrid melt(model->Ncols, model->Nrows);
  PETScGrid creep(model->Ncols, model->Nrows);
  PETScGrid Q(model->Ncols, model->Nrows);

  melt.setZero();
  creep.setZero();
  Q.setZero();

  PetscScalar currTime = 0.0;
  clock_t t;
  if (rank == 0) {
    t = clock();
  }

  // start
  // creating grids outside of loop to save time
  PETScGrid Se(Sp->getTotalNumOfCols(), Sp->getTotalNumOfRows());
  int size = model->Ncols * model->Nrows;
  PETScMat A(size, size);
  PetscScalar cavity_opening = 0;
  PETScGrid Teff(T->getTotalNumOfCols(), T->getTotalNumOfRows());
  PETScGrid TeffPowTexp(T->getTotalNumOfCols(), T->getTotalNumOfRows());

  for (int timeStep = 1; timeStep < Nt + 1; ++timeStep) {
    currTime += dt_secs;

    // time dependent forcing
    auto &currentQ = model->Q->getCurrentQ(currTime);

    // if (args->seaLevelForcing) {
    // TODO
    //}

    // copying T to Teff is done inside if and else
    // enable_unconfied
    if (!args->disableUnconfined) {
      enableUnconfined(Teff, TeffPowTexp, *Sp, *T_n, *K, *model->topg, *u_n, args->Texp, args->unconfSmooth,
                       args->layerThickness);
    } else {
      Sp->setZero();
      calculateTeffPowTexp(Teff, TeffPowTexp, *T, args->Texp);
    }

    calculateSeValues(Se, *Sp, *S);

    systemmatrix(A, b, model->Nrows, model->Ncols, Se, TeffPowTexp, model->dx, dt_secs, theta, *u_n, currentQ,
                 *dirichletValues, *dirichletMask);
    // solve the equation A*sol = b
    PETScSolver::solve(A, b, sol);

    u->setGlobalVecColMajor(sol);

    if (args->dochannels) {
      doChannels(melt, creep, *u_n, *gradMask, *T, *T_n, *model->pIce, *model->topg, *K, args->flowConstant, args->Texp,
                 args->roughnessFactor, args->noSmoothMelt, args->layerThickness, model->dx);
    } else {
      cavity_opening = 0;
      noChannels(melt, creep);
    }

    // switch pointers
    u_n.swap(u);
    T_n.swap(T);

    // we need solution.saveTimestep() for this to work
    /*if (timeStep % args->saveEvery == 0) {
      saveSolution(timeStep, *args, rank, *u, *u_n, *model, melt, cavity_opening);
    }*/
  }
  // end
  if (rank == 0) {
    t = clock() - t;
    Logger::instance().info("CUASSolver.cpp: solve(): computation took: {} seconds.", ((float)t) / CLOCKS_PER_SEC);
  }
}

}  // namespace CUAS
