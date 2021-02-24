#include "CUASSolver.h"

#include "CUASKernels.h"
#include "fill_matrix_coo.h"
#include "savetimestep.h"
#include "specialgradient.h"

#include "PETScMat.h"
#include "PETScSolver.h"
#include "PETScVec.h"

#include <cmath>
#include <memory>

// TODO should not be necessary in a solver use spdlog instead
#include <iostream>

namespace CUAS {

void solve(std::unique_ptr<PETScGrid> &u, std::unique_ptr<PETScGrid> &u_n, CUASModel &model, int const Nt,
           CUASArgs const &args, PetscScalar const totaltime_secs, PetscScalar const dt_secs) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //
  // SOLVER PREPARATION
  //
  PETScVec sol(model.Nrows * model.Ncols);

  PetscScalar It[Nt + 1];
  for (int i = 0; i < Nt + 1; ++i) {
    It[i] = i;
  }

  // auto u = std::make_unique<PETScGrid>(model.Ncols, model.Nrows);    // unknown u at new time level
  // auto u_n = std::make_unique<PETScGrid>(model.Ncols, model.Nrows);  // u at the previous time level

  // why is this not an arg?
  const int theta = 1;  // 1 means fully implicit, 0 means fully explicit, 0.5 is Crank-Nicholson

  PETScVec b(model.Nrows * model.Ncols);

  if (args.initialHead == "zero") {
    u_n->setZero();
  } else if (args.initialHead == "Nzero") {
    pressure2head(*u_n, *model.pIce, *model.topg, 0.0);
  } else if (args.initialHead == "topg") {
    u_n->copy(*model.topg);
  } else {
    // throw error
    std::cerr << "initialHead needs to be zero, Nzero or topg!" << std::endl;
    return;
  }

  //  if(args.restart){
  //    if(args.verbose){
  //      std::cout << "Read restart from file " << args.restart << std::endl;
  //    }
  // TODO
  //  }

  // save timedependent values
  PetscScalar Ntsaved;
  auto Itlength = sizeof(It) / sizeof(It[0]);
  if (Itlength % args.saveEvery > 0) {
    Ntsaved = ceil(Itlength / args.saveEvery);
  }

  if (args.verbose) {
    std::cout << "runtime = " << totaltime_secs << ", time step = " << dt_secs << ", Ntsaved = " << Ntsaved
              << " for saveEvery = " << args.saveEvery << std::endl;
  }

  // TODO!! solution init (part of saving to netcdf, see original-python main: 272-274)
  // melt, creep and Q are supposed to be part of the solution class.
  PETScGrid melt(model.Ncols, model.Nrows);
  PETScGrid creep(model.Ncols, model.Nrows);
  PETScGrid Q(model.Ncols, model.Nrows);

  melt.setZero();
  creep.setZero();
  Q.setZero();

  int time_current = 0;
  clock_t t;
  if (rank == 0) {
    t = clock();
  }

  // start
  // creating grids outside of loop to save time
  PETScGrid Se(model.Sp->getTotalNumOfCols(), model.Sp->getTotalNumOfRows());
  int size = model.Ncols * model.Nrows;
  PETScMat A(size, size);
  PetscScalar cavity_opening = 0;
  PETScGrid Teff(model.T->getTotalNumOfCols(), model.T->getTotalNumOfRows());
  PETScGrid TeffPowTexp(model.T->getTotalNumOfCols(), model.T->getTotalNumOfRows());

  for (int timeStep = 1; timeStep < Nt + 1; ++timeStep) {
    time_current += dt_secs;
    // TODO get_current_Q (part of time dependent forcing)
    PETScGrid &current_Q = *model.Q;  // get_current_Q(time_current);

    // if (args.seaLevelForcing) {
    // TODO
    //}

    PetscScalar bt = args.layerThickness;
    // copying T to Teff is done inside if and else
    // enable_unconfied
    if (!args.disableUnconfined) {
      enableUnconfined(Teff, TeffPowTexp, model, *u_n, args, bt);
    } else {
      model.Sp->setZero();
      calculateTeffPowTexp(Teff, TeffPowTexp, *model.T, args);
    }

    calculateSeValues(Se, *model.Sp, *model.S);

    fill_matrix_coo(A, b, model.Nrows, model.Ncols, Se, TeffPowTexp, model.dx, dt_secs, theta, *u_n, current_Q,
                    *model.dirichletValues, *model.dirichletMask);
    // solve the equation A*sol = b
    PETScSolver::solve(A, b, sol);

    u->setGlobalVecColMajor(sol);

    if (args.dochannels) {
      doChannels(melt, creep, *u_n, model, args, bt);
    } else {
      cavity_opening = 0;
      noChannels(melt, creep);
    }

    // switch pointers
    u_n.swap(u);
    model.T_n.swap(model.T);

    // we need solution.saveTimestep() for this to work
    if (timeStep % args.saveEvery == 0) {
      saveSolution(timeStep, args, rank, *u, *u_n, model, melt, cavity_opening);
    }
  }
  // end
  if (rank == 0) {
    t = clock() - t;
    std::cout << "computation took: " << ((float)t) / CLOCKS_PER_SEC << " seconds." << std::endl;
  }
}

}  // namespace CUAS
