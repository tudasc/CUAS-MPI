#ifndef CUAS_KERNELS_H
#define CUAS_KERNELS_H

#include "CUASArgs.h"
#include "CUASModel.h"
#include "helper.h"
#include "physicalConstants.h"
#include "specialgradient.h"

#include "PETScGrid.h"

namespace CUAS {

inline void binaryDialation(PetscGrid &output, PetscGrid const &input) {
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

  return;
}

inline void enableUnconfined(PetscGrid &Teff, PetscGrid &TeffPowTexp, CUASModel &model, PetscGrid const &u_n,
                             CUASArgs const &args, double const bt) {
  auto &u_nGlobal = u_n.getReadHandle();
  auto &T_nGlobal = model.T_n->getReadHandle();
  auto &topgGlobal = model.topg->getReadHandle();
  auto &KGlobal = model.K->getReadHandle();
  auto TeffGlobal = Teff.getWriteHandle();
  auto TeffPowTexpGlobal = TeffPowTexp.getWriteHandle();
  auto SpGlobalWrite = model.Sp->getWriteHandle();

  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      double psi = u_nGlobal(j, i) - topgGlobal(j, i);
      if (psi < 0) {
        psi = 0.01;
      }
      if (psi < bt) {
        TeffGlobal(j, i) = KGlobal(j, i) * psi;
      } else {
        // if (timeStep == 1) {
        //   std::cout << "t_nglobal : " << T_nGlobal(j, i) << std::endl;
        // }
        TeffGlobal(j, i) = T_nGlobal(j, i);
      }
      TeffPowTexpGlobal(j, i) = pow(TeffGlobal(j, i), args.Texp);

      if (psi < (bt - args.unconfSmooth)) {
        SpGlobalWrite(j, i) = 0.4;
      } else {
        SpGlobalWrite(j, i) = 0;
      }
    }
  }
}

inline void calculateTeffPowTexp(PetscGrid &Teff, PetscGrid &TeffPowTexp, PetscGrid const &T, CUASArgs const &args) {
  // copy T to Teff
  auto &TGlobal = T.getReadHandle();
  auto TeffGlobal = Teff.getWriteHandle();
  auto TeffPowTexpGlobal = TeffPowTexp.getWriteHandle();
  for (int j = 0; j < Teff.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Teff.getLocalNumOfCols(); ++i) {
      TeffGlobal(j, i) = TGlobal(j, i);
      TeffPowTexpGlobal(j, i) = pow(TeffGlobal(j, i), args.Texp);
    }
  }
}

inline void calculateSeValues(PetscGrid &Se, PetscGrid const &Sp, PetscGrid const &S) {  // extra scope for handles
  auto SeGlobal = Se.getWriteHandle();
  auto &SpGlobalRead = Sp.getReadHandle();
  auto &SGlobal = S.getReadHandle();
  for (int j = 0; j < Se.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Se.getLocalNumOfCols(); ++i) {
      SeGlobal(j, i) = SGlobal(j, i) + SpGlobalRead(j, i);  // effective storage
    }
  }
}

inline void doChannels(PetscGrid &melt, PetscGrid &creep, PetscGrid const &u_n, CUASModel const &model,
                       CUASArgs const &args, PetscScalar const bt) {
  PetscGrid maggrad2(u_n.getTotalNumOfCols(), u_n.getTotalNumOfRows());
  gradient2(maggrad2, u_n, model.dx);

  auto maggrad2Global = maggrad2.getWriteHandle();
  auto &grad_maskGlobal = model.grad_mask->getReadHandle();
  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      if (maggrad2Global(j, i) == grad_maskGlobal(j, i)) {
        maggrad2Global(j, i) = 0;
      }
    }
  }

  auto creepGlobal = creep.getWriteHandle();
  auto &TGlobal = model.T->getReadHandle();
  auto &p_iceGlobal = model.p_ice->getReadHandle();
  auto &u_nGlobal = u_n.getReadHandle();
  auto &topgGlobal = model.topg->getReadHandle();
  PetscScalar multValue = 2 * args.flowConstant;
  PetscScalar onethird = 1.0 / 3.0;
  PetscScalar rhowater_gravity = RHO_WATER * GRAVITY;
  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      double p_w = (u_nGlobal(j, i) - topgGlobal(j, i) - bt) * rhowater_gravity;
      double N = p_iceGlobal(j, i) - p_w;
      creepGlobal(j, i) = multValue * pow(N * onethird, 3) * TGlobal(j, i);
    }
  }

  if (args.Texp != 1) {
    PetscGrid T_nExp(u_n.getTotalNumOfCols(), u_n.getTotalNumOfRows());
    auto T_nExpGlobal = T_nExp.getWriteHandle();
    auto &T_nGlobal = model.T_n->getReadHandle();
    for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
        T_nExpGlobal(j, i) = pow(T_nGlobal(j, i), args.Texp);
      }
    }
    computeMelt(melt, args.roughnessFactor, GRAVITY, RHO_WATER, T_nExp, *model.K, maggrad2, RHO_ICE, LATENT_HEAT, bt);
  } else {
    computeMelt(melt, args.roughnessFactor, GRAVITY, RHO_WATER, *model.T_n, *model.K, maggrad2, RHO_ICE, LATENT_HEAT,
                bt);
  }
  // TODO: if not args.noSmoothMelt:
}

inline void noChannels(PetscGrid &melt, PetscGrid &creep) {
  melt.setZero();
  creep.setZero();
  // this is probably faster than setting zero
  /*auto meltGlobal = melt.getWriteHandle();
  auto creepGlobal = creep.getWriteHandle();
  for (int j = 0; j < melt.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < melt.getLocalNumOfCols(); ++i) {
      meltGlobal(j, i) = 0;
      creepGlobal(j, i) = 0;
    }
  }*/
}

}  // namespace CUAS

#endif
