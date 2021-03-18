#ifndef CUAS_KERNELS_H
#define CUAS_KERNELS_H

#include "CUASConstants.h"
#include "specialgradient.h"

#include "PETScGrid.h"

namespace CUAS {

/*
 * Converts hydraulic head to water pressure
 * bedElevation and seaLevel are relative to the mean sea level = 0.0
 */
inline void head2pressure(PETScGrid &result, PETScGrid const &head, PETScGrid const &bedElevation,
                          PetscScalar const seaLevel = 0.0) {
  if (!result.isCompatible(head) || !result.isCompatible(bedElevation)) {
    exit(1);
  }

  auto result2d = result.getWriteHandle();
  auto &head2d = head.getReadHandle();
  auto &bedElevation2d = bedElevation.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bedElevation2d(j, i) - seaLevel;
      result2d(j, i) = RHO_WATER * GRAVITY * (head2d(j, i) - effective_bed_elevation);
    }
  }
}

/*
 * Convert water pressure to hydraulic head
 * bed_elevation and sea_level are relative to the mean sea level = 0.0
 * pressure is ice pressure
 */
inline void pressure2head(PETScGrid &result, PETScGrid const &pressure, PETScGrid const &bedElevation,
                          PetscScalar const seaLevel = 0.0) {
  if (!result.isCompatible(pressure) || !result.isCompatible(bedElevation)) {
    exit(1);
  }

  auto result2d = result.getWriteHandle();
  auto &pressure2d = pressure.getReadHandle();
  auto &bedElevation2d = bedElevation.getReadHandle();

  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bedElevation2d(j, i) - seaLevel;
      result2d(j, i) = pressure2d(j, i) / (RHO_WATER * GRAVITY) + effective_bed_elevation;
    }
  }
}

/*
 * Compute ice overburden pressure
 */
inline void overburdenPressure(PETScGrid &result, PETScGrid const &thk) {
  if (!result.isCompatible(thk)) {
    exit(1);
  }

  auto result2d = result.getWriteHandle();
  auto &thk2d = thk.getReadHandle();
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      result2d(j, i) = thk2d(j, i) * RHO_ICE * GRAVITY;
    }
  }
}

/*
 * From Werder2013 / summers2018:
 * beta = (b r − b)/l r for b < b r , beta = 0 for b ≥ b r
 */
inline void cavityOpenB(PETScGrid &result, PetscScalar const beta, PetscScalar const v_b, PETScGrid const &K) {
  if (!result.isCompatible(K)) {
    exit(1);
  }

  auto resultGlobal = result.getWriteHandle();
  auto &KGlobal = K.getReadHandle();
  PetscScalar betaVb = beta * v_b;
  for (int j = 0; j < result.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result.getLocalNumOfCols(); ++i) {
      resultGlobal(j, i) = betaVb * KGlobal(j, i);
    }
  }
}

/*
 * like deFleurian2016
 */
inline void computeMelt(PETScGrid &melt, PetscScalar const r, PetscScalar const g, PetscScalar const rho_w,
                        PETScGrid const &T, PETScGrid const &K, PETScGrid const &gradh2, PetscScalar const rho_i,
                        PetscScalar const L, PetscScalar const bt) {
  if (!melt.isCompatible(T) || !melt.isCompatible(T) || !melt.isCompatible(K) || !melt.isCompatible(gradh2)) {
    exit(1);
  }

  auto resultGlobal = melt.getWriteHandle();
  auto &KGlobal = K.getReadHandle();
  auto &TGlobal = T.getReadHandle();
  auto &gradh2Global = gradh2.getReadHandle();
  const PetscScalar r_g_rhow = r * g * rho_w;
  const PetscScalar rhoi_L_inv = 1.0 / (rho_i * L);
  for (int j = 0; j < melt.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < melt.getLocalNumOfCols(); ++i) {
      resultGlobal(j, i) = r_g_rhow * TGlobal(j, i) * KGlobal(j, i) * gradh2Global(j, i) * rhoi_L_inv;
    }
  }
}

inline void binaryDialation(PETScGrid &output, PETScGrid const &input) {
  if (!output.isCompatible(input)) {
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

inline void enableUnconfined(PETScGrid &Teff, PETScGrid &TeffPowTexp, PETScGrid &Sp, PETScGrid const &T_n,
                             PETScGrid const &K, PETScGrid const &topg, PETScGrid const &u_n, PetscScalar const Texp,
                             PetscScalar const unconfSmooth, PetscScalar const bt) {
  if (!Teff.isCompatible(TeffPowTexp) || !Teff.isCompatible(u_n) || !Teff.isCompatible(T_n) ||
      !Teff.isCompatible(topg) || !Teff.isCompatible(K) || !Teff.isCompatible(Sp)) {
    exit(1);
  }

  auto &u_nGlobal = u_n.getReadHandle();
  auto &T_nGlobal = T_n.getReadHandle();
  auto &topgGlobal = topg.getReadHandle();
  auto &KGlobal = K.getReadHandle();
  auto TeffGlobal = Teff.getWriteHandle();
  auto TeffPowTexpGlobal = TeffPowTexp.getWriteHandle();
  auto SpGlobalWrite = Sp.getWriteHandle();

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
      TeffPowTexpGlobal(j, i) = pow(TeffGlobal(j, i), Texp);

      if (psi < (bt - unconfSmooth)) {
        SpGlobalWrite(j, i) = 0.4;
      } else {
        SpGlobalWrite(j, i) = 0;
      }
    }
  }
}

inline void calculateTeffPowTexp(PETScGrid &Teff, PETScGrid &TeffPowTexp, PETScGrid const &T, PetscScalar const Texp) {
  if (!Teff.isCompatible(TeffPowTexp) || !Teff.isCompatible(T)) {
    exit(1);
  }

  // copy T to Teff
  auto &TGlobal = T.getReadHandle();
  auto TeffGlobal = Teff.getWriteHandle();
  auto TeffPowTexpGlobal = TeffPowTexp.getWriteHandle();
  for (int j = 0; j < Teff.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Teff.getLocalNumOfCols(); ++i) {
      TeffGlobal(j, i) = TGlobal(j, i);
      TeffPowTexpGlobal(j, i) = pow(TeffGlobal(j, i), Texp);
    }
  }
}

inline void calculateSeValues(PETScGrid &Se, PETScGrid const &Sp, PETScGrid const &S) {  // extra scope for handles
  if (!Se.isCompatible(Sp) || !Se.isCompatible(S)) {
    exit(1);
  }

  auto SeGlobal = Se.getWriteHandle();
  auto &SpGlobalRead = Sp.getReadHandle();
  auto &SGlobal = S.getReadHandle();
  for (int j = 0; j < Se.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Se.getLocalNumOfCols(); ++i) {
      SeGlobal(j, i) = SGlobal(j, i) + SpGlobalRead(j, i);  // effective storage
    }
  }
}

inline void doChannels(PETScGrid &melt, PETScGrid &creep, PETScGrid const &u_n, PETScGrid const &gradMask,
                       PETScGrid const &T, PETScGrid const &T_n, PETScGrid const &pIce, PETScGrid const &topg,
                       PETScGrid const &K, PetscScalar const flowConstant, PetscScalar const Texp,
                       PetscScalar const roughnessFactor, PetscScalar const bt, PetscScalar const dx) {
  if (!melt.isCompatible(creep) || !melt.isCompatible(u_n)) {
    exit(1);
  }

  PETScGrid maggrad2(u_n.getTotalNumOfCols(), u_n.getTotalNumOfRows());
  gradient2(maggrad2, u_n, dx);

  auto maggrad2Global = maggrad2.getWriteHandle();
  auto &grad_maskGlobal = gradMask.getReadHandle();
  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      if (maggrad2Global(j, i) == grad_maskGlobal(j, i)) {
        maggrad2Global(j, i) = 0;
      }
    }
  }

  auto creepGlobal = creep.getWriteHandle();
  auto &TGlobal = T.getReadHandle();
  auto &p_iceGlobal = pIce.getReadHandle();
  auto &u_nGlobal = u_n.getReadHandle();
  auto &topgGlobal = topg.getReadHandle();
  PetscScalar multValue = 2 * flowConstant;
  PetscScalar onethird = 1.0 / 3.0;
  PetscScalar rhowater_gravity = RHO_WATER * GRAVITY;
  for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
      double p_w = (u_nGlobal(j, i) - topgGlobal(j, i) - bt) * rhowater_gravity;
      double N = p_iceGlobal(j, i) - p_w;
      creepGlobal(j, i) = multValue * pow(N * onethird, 3) * TGlobal(j, i);
    }
  }

  if (Texp != 1) {
    PETScGrid T_nExp(u_n.getTotalNumOfCols(), u_n.getTotalNumOfRows());
    auto T_nExpGlobal = T_nExp.getWriteHandle();
    auto &T_nGlobal = T_n.getReadHandle();
    for (int j = 0; j < u_n.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < u_n.getLocalNumOfCols(); ++i) {
        T_nExpGlobal(j, i) = pow(T_nGlobal(j, i), Texp);
      }
    }
    computeMelt(melt, roughnessFactor, GRAVITY, RHO_WATER, T_nExp, K, maggrad2, RHO_ICE, LATENT_HEAT, bt);
  } else {
    computeMelt(melt, roughnessFactor, GRAVITY, RHO_WATER, T_n, K, maggrad2, RHO_ICE, LATENT_HEAT, bt);
  }
  // TODO: if not args.noSmoothMelt:
}

inline void noChannels(PETScGrid &melt, PETScGrid &creep) {
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
