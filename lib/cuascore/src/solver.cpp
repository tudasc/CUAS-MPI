#include "solver.h"
#include "CUASModel.h"
#include "PetscGrid.h"
#include "PetscMat.h"
#include "PetscSolver.h"
#include "PetscVec.h"
#include "fill_matrix_coo.h"
#include "helper.h"
#include "physicalConstants.h"
#include "specialgradient.h"

#include "petscdump.h"

#include <iostream>
#include <math.h>
#include <memory>

void solve(std::unique_ptr<PetscGrid> &u, std::unique_ptr<PetscGrid> &u_n, int const Nt, CUAS::CUASModel &model,
           CUAS::CUASArgs const &args, PetscScalar const totaltime_secs, PetscScalar const dt_secs) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //
  // SOLVER PREPARATION
  //
  PetscVec sol(model.Nrows * model.Ncols);

  PetscScalar It[Nt + 1];
  for (int i = 0; i < Nt + 1; ++i) {
    It[i] = i;
  }

  // auto u = std::make_unique<PetscGrid>(model.Ncols, model.Nrows);    // unknown u at new time level
  // auto u_n = std::make_unique<PetscGrid>(model.Ncols, model.Nrows);  // u at the previous time level

  // why is this not an arg?
  const int theta = 1;  // 1 means fully implicit, 0 means fully explicit, 0.5 is Crank-Nicholson

  PetscVec b(model.Nrows * model.Ncols);

  auto u_nGlobal = u_n->getAsGlobal2dArr();

  if (args.initialHead == "zero") {
    u_n->setZero();
  } else if (args.initialHead == "Nzero") {
    CUAS::pressure2head(*u_n, *model.p_ice, *model.topg, 0.0);
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
  int Itlength = sizeof(It) / sizeof(It[0]);
  if (Itlength % args.saveEvery > 0) {
    Ntsaved = ceil(Itlength / args.saveEvery);
  }

  if (args.verbose) {
    std::cout << "runtime = " << totaltime_secs << ", time step = " << dt_secs << ", Ntsaved = " << Ntsaved
              << " for saveEvery = " << args.saveEvery << std::endl;
  }

  // TODO!! solution init (part of saving to netcdf, see original-python main: 272-274)
  // melt, creep and Q are supposed to be part of the solution class.
  PetscGrid melt(model.Ncols, model.Nrows);
  PetscGrid creep(model.Ncols, model.Nrows);
  PetscGrid Q(model.Ncols, model.Nrows);

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
  PetscGrid Se(model.Sp->getTotalNumOfCols(), model.Sp->getTotalNumOfRows());
  int size = model.Ncols * model.Nrows;
  PetscMat A(size, size);
  PetscScalar cavity_opening = 0;
  PetscGrid Teff(model.T->getTotalNumOfCols(), model.T->getTotalNumOfRows());
  PetscGrid TeffPowTexp(model.T->getTotalNumOfCols(), model.T->getTotalNumOfRows());

  PetscSolver solver;
  for (int timeStep = 1; timeStep < Nt + 1; ++timeStep) {
    time_current += dt_secs;
    // TODO get_current_Q (part of time dependent forcing)
    PetscGrid &current_Q = *model.Q;  // get_current_Q(time_current);

    // if (args.seaLevelForcing) {
    // TODO
    //}

    // copying T to Teff is done inside if and else
    PetscScalar bt = args.layerThickness;
    // enable_unconfied
    if (!args.disableUnconfined) {
      auto u_nGlobal = u_n->getAsGlobal2dArr();
      auto topgGlobal = model.topg->getAsGlobal2dArr();
      auto TeffGlobal = Teff.getAsGlobal2dArr();
      auto TeffPowTexpGlobal = TeffPowTexp.getAsGlobal2dArr();
      auto KGlobal = model.K->getAsGlobal2dArr();
      auto T_nGlobal = model.T_n->getAsGlobal2dArr();
      auto SpGlobal = model.Sp->getAsGlobal2dArr();
      for (int j = 0; j < u_n->getLocalNumOfRows(); ++j) {
        for (int i = 0; i < u_n->getLocalNumOfCols(); ++i) {
          double psi = u_nGlobal[j][i] - topgGlobal[j][i];
          if (psi < 0) {
            psi = 0.01;
          }
          if (psi < bt) {
            TeffGlobal[j][i] = KGlobal[j][i] * psi;
          } else {
            TeffGlobal[j][i] = T_nGlobal[j][i];
          }
          TeffPowTexpGlobal[j][i] = pow(TeffGlobal[j][i], args.Texp);

          if (psi < (bt - args.unconfSmooth)) {
            SpGlobal[j][i] = 0.4;
          } else {
            SpGlobal[j][i] = 0;
          }
        }
      }
      u_n->restoreGlobal2dArr(u_nGlobal);
      model.topg->restoreGlobal2dArr(topgGlobal);
      Teff.setAsGlobal2dArr(TeffGlobal);
      TeffPowTexp.setAsGlobal2dArr(TeffPowTexpGlobal);
      model.T_n->restoreGlobal2dArr(T_nGlobal);
      model.Sp->setAsGlobal2dArr(SpGlobal);
      model.K->restoreGlobal2dArr(KGlobal);
    } else {
      model.Sp->setZero();

      // copy T to Teff
      auto TeffGlobal = Teff.getAsGlobal2dArr();
      auto TeffPowTexpGlobal = TeffPowTexp.getAsGlobal2dArr();
      auto TGlobal = model.T->getAsGlobal2dArr();
      for (int j = 0; j < Teff.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < Teff.getLocalNumOfCols(); ++i) {
          TeffGlobal[j][i] = TGlobal[j][i];
          TeffPowTexpGlobal[j][i] = pow(TeffGlobal[j][i], args.Texp);
        }
      }
      model.T->restoreGlobal2dArr(TGlobal);
      Teff.setAsGlobal2dArr(TeffGlobal);
      TeffPowTexp.setAsGlobal2dArr(TeffPowTexpGlobal);
    }

    auto SeGlobal = Se.getAsGlobal2dArr();
    auto SpGlobal = model.Sp->getAsGlobal2dArr();
    auto SGlobal = model.S->getAsGlobal2dArr();
    for (int j = 0; j < Se.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < Se.getLocalNumOfCols(); ++i) {
        SeGlobal[j][i] = SGlobal[j][i] + SpGlobal[j][i];  // effective storage
      }
    }
    Se.setAsGlobal2dArr(SeGlobal);
    model.Sp->restoreGlobal2dArr(SpGlobal);
    model.S->restoreGlobal2dArr(SGlobal);

    CUAS::fill_matrix_coo(A, b, model.Nrows, model.Ncols, Se, TeffPowTexp, model.dx, dt_secs, theta, *u_n, current_Q,
                          *model.dirichlet_values, *model.dirichlet_mask);

    // solve the equation A*sol = b
    solver.solve(A, b, sol);

    u->setGlobalVecColMajor(sol);

    if (args.dochannels) {
      auto maggrad2 = std::make_unique<PetscGrid>(u_n->getTotalNumOfCols(), u_n->getTotalNumOfRows());
      CUAS::gradient2(*u_n, *maggrad2, model.dx);

      auto maggrad2Global = maggrad2->getAsGlobal2dArr();
      auto grad_maskGlobal = model.grad_mask->getAsGlobal2dArr();
      for (int j = 0; j < u_n->getLocalNumOfRows(); ++j) {
        for (int i = 0; i < u_n->getLocalNumOfCols(); ++i) {
          if (maggrad2Global[j][i] == grad_maskGlobal[j][i]) {
            maggrad2Global[j][i] = 0;
          }
        }
      }
      maggrad2->setAsGlobal2dArr(maggrad2Global);
      model.grad_mask->restoreGlobal2dArr(grad_maskGlobal);

      auto creepGlobal = creep.getAsGlobal2dArr();
      auto TGlobal = model.T->getAsGlobal2dArr();
      auto p_iceGlobal = model.p_ice->getAsGlobal2dArr();
      auto u_nGlobal = u_n->getAsGlobal2dArr();
      auto topgGlobal = model.topg->getAsGlobal2dArr();
      PetscScalar multValue = 2 * args.flowConstant;
      PetscScalar onethird = 1.0 / 3.0;
      PetscScalar rhowater_gravity = RHO_WATER * GRAVITY;
      for (int j = 0; j < u_n->getLocalNumOfRows(); ++j) {
        for (int i = 0; i < u_n->getLocalNumOfCols(); ++i) {
          double p_w = (u_nGlobal[j][i] - topgGlobal[j][i] - bt) * rhowater_gravity;
          double N = p_iceGlobal[j][i] - p_w;
          creepGlobal[j][i] = multValue * pow(N * onethird, 3) * TGlobal[j][i];
        }
      }
      creep.setAsGlobal2dArr(creepGlobal);
      model.T->restoreGlobal2dArr(TGlobal);
      model.p_ice->restoreGlobal2dArr(p_iceGlobal);
      u_n->restoreGlobal2dArr(u_nGlobal);
      model.topg->restoreGlobal2dArr(topgGlobal);

      if (args.Texp != 1) {
        auto T_nExp = std::make_unique<PetscGrid>(u_n->getTotalNumOfCols(), u_n->getTotalNumOfRows());
        auto T_nExpGlobal = T_nExp->getAsGlobal2dArr();
        auto T_nGlobal = model.T_n->getAsGlobal2dArr();
        for (int j = 0; j < u_n->getLocalNumOfRows(); ++j) {
          for (int i = 0; i < u_n->getLocalNumOfCols(); ++i) {
            T_nExpGlobal[j][i] = pow(T_nGlobal[j][i], args.Texp);
          }
        }
        T_nExp->setAsGlobal2dArr(T_nExpGlobal);
        model.T_n->setAsGlobal2dArr(T_nGlobal);
        CUAS::compute_melt(melt, args.roughnessFactor, GRAVITY, RHO_WATER, *T_nExp, *model.K, *maggrad2, RHO_ICE,
                           LATENT_HEAT, bt);
      } else {
        CUAS::compute_melt(melt, args.roughnessFactor, GRAVITY, RHO_WATER, *model.T_n, *model.K, *maggrad2, RHO_ICE,
                           LATENT_HEAT, bt);
      }
      // TODO: if not args.noSmoothMelt:
    } else {
      cavity_opening = 0;
      auto meltGlobal = melt.getAsGlobal2dArr();
      auto creepGlobal = creep.getAsGlobal2dArr();
      for (int j = 0; j < u_n->getLocalNumOfRows(); ++j) {
        for (int i = 0; i < u_n->getLocalNumOfCols(); ++i) {
          meltGlobal[j][i] = 0;
          creepGlobal[j][i] = 0;
        }
      }
      melt.setAsGlobal2dArr(meltGlobal);
      creep.setAsGlobal2dArr(creepGlobal);
    }

    // switch pointers
    u_n.swap(u);
    model.T_n.swap(model.T);

    // we need solution.saveTimestep() for this to work
    // if (timeStep % args.saveEvery == 0) {
    //   if (args.verbose && rank == 0) {
    //     // TODO: implement time_display_units from helper to get time_scaling
    //     // std::cout << "time(" << i << "/" << Nt << ") = " << time_current/time_scaling << " (" << time_units << ")"
    //     <<
    //     // std::endl;
    //   }
    //   // eps_head: get highest differing value from u and u_n
    //   // eps_T: get highest differing value from T and T_n
    //   PetscScalar eps_head;
    //   PetscScalar eps_T;
    //   auto uGlobal = u->getAsGlobal2dArr();
    //   auto u_nGlobal = u_n->getAsGlobal2dArr();
    //   auto TGlobal = model.T->getAsGlobal2dArr();
    //   auto T_nGlobal = model.T_n->getAsGlobal2dArr();
    //   for (int j = 0; j < u->getLocalNumOfRows(); ++j) {
    //     for (int i = 0; i < u->getLocalNumOfCols(); ++i) {
    //       PetscScalar diffHead = abs(uGlobal[j][i] - u_nGlobal[j][i]);
    //       PetscScalar diffT = abs(TGlobal[j][i] - T_nGlobal[j][i]);
    //       if (diffHead > eps_head) {
    //         eps_head = diffHead;
    //       }
    //       if (diffT > eps_T) {
    //         eps_T = diffT;
    //       }
    //     }
    //   }
    //   u->restoreGlobal2dArr(uGlobal);
    //   u_n->restoreGlobal2dArr(u_nGlobal);
    //   model.T->restoreGlobal2dArr(TGlobal);
    //   model.T_n->restoreGlobal2dArr(T_nGlobal);

    //   auto meltGlobal = melt.getAsGlobal2dArr();
    //   for (int j = 0; j < melt.getLocalNumOfRows(); ++j) {
    //     for (int i = 0; i < melt.getLocalNumOfCols(); ++i) {
    //       meltGlobal[j][i] = meltGlobal[j][i] + cavity_opening;
    //     }
    //   }
    //   melt.setAsGlobal2dArr(meltGlobal);
    //   solution.saveTimestep();
    // }
  }
  // end
  if (rank == 0) {
    t = clock() - t;
    std::cout << "computation took: " << ((float)t) / CLOCKS_PER_SEC << " seconds." << std::endl;
  }
}
