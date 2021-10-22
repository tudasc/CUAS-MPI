#include "CUASKernels.h"

#include "PETScGrid.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 9
#define GRID_SIZE_X 18
#define GRID_SIZE_Y 8

TEST(CUASKernelsTest, head2pressure) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid pressure(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid bedElevation(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid head(GRID_SIZE_Y, GRID_SIZE_X);

  PetscScalar seaLevel = 33.5;

  // setup
  {
    auto pressure2d = pressure.getWriteHandle();
    auto bedElevation2d = bedElevation.getWriteHandle();
    auto head2d = head.getWriteHandle();
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        pressure2d(j, i) = seaLevel * mpiRank - 2.3;
        bedElevation2d(j, i) = seaLevel * seaLevel - 31.3;
        head2d(j, i) = mpiRank * j + (i * 35);
      }
    }
  }

  // run head2pressure
  CUAS::head2pressure(pressure, head, bedElevation, seaLevel);

  // check
  {
    auto &pressure2d = pressure.getReadHandle();
    auto &bedElevation2d = bedElevation.getReadHandle();
    auto &head2d = head.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        double effective_bed_elevation = bedElevation2d(j, i) - seaLevel;
        // check result
        ASSERT_EQ(pressure2d(j, i), RHO_WATER * GRAVITY * (head2d(j, i) - effective_bed_elevation));

        // check for sideeffects
        ASSERT_EQ(bedElevation2d(j, i), seaLevel * seaLevel - 31.3);
        ASSERT_EQ(head2d(j, i), mpiRank * j + (i * 35));
      }
    }
  }
}

TEST(CUASKernelsTest, pressure2head) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid pressure(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid bedElevation(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid head(GRID_SIZE_Y, GRID_SIZE_X);

  PetscScalar seaLevel = 33.5;

  // setup
  {
    auto pressure2d = pressure.getWriteHandle();
    auto bedElevation2d = bedElevation.getWriteHandle();
    auto head2d = head.getWriteHandle();
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        pressure2d(j, i) = seaLevel * mpiRank - 2.3;
        bedElevation2d(j, i) = seaLevel * seaLevel - 31.3;
        head2d(j, i) = mpiRank * j + (i * 35);
      }
    }
  }

  // run pressure2head
  CUAS::pressure2head(head, pressure, bedElevation, seaLevel);

  // check
  {
    auto &pressure2d = pressure.getReadHandle();
    auto &bedElevation2d = bedElevation.getReadHandle();
    auto &head2d = head.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        double effective_bed_elevation = bedElevation2d(j, i) - seaLevel;
        // check result
        ASSERT_EQ(head2d(j, i), pressure2d(j, i) / (RHO_WATER * GRAVITY) + effective_bed_elevation);

        // check for sideeffects
        ASSERT_EQ(bedElevation2d(j, i), seaLevel * seaLevel - 31.3);
        ASSERT_EQ(pressure2d(j, i), seaLevel * mpiRank - 2.3);
      }
    }
  }
}

TEST(CUASKernelsTest, overburdenPressure) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid pressure(30, 12);
  PETScGrid thk(30, 12);

  PetscScalar seaLevel = 33.5;

  // setup
  {
    auto pressure2d = pressure.getWriteHandle();
    auto thk2d = thk.getWriteHandle();
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        pressure2d(j, i) = seaLevel * mpiRank - 2.3;
        thk2d(j, i) = mpiRank * j + (i * 35);
      }
    }
  }

  // run overburdenPressure
  CUAS::overburdenPressure(pressure, thk);

  // check
  {
    auto &pressure2d = pressure.getReadHandle();
    auto &thk2d = thk.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        // check result
        ASSERT_EQ(pressure2d(j, i), thk2d(j, i) * RHO_ICE * GRAVITY);

        // check for sideeffects
        ASSERT_EQ(thk2d(j, i), mpiRank * j + (i * 35));
      }
    }
  }
}

TEST(CUASKernelsTest, cavityOpenB) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid K(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid result(GRID_SIZE_Y, GRID_SIZE_X);

  PetscScalar beta = 33.5;
  PetscScalar v_b = 17.5;

  // setup
  {
    auto K2d = K.getWriteHandle();
    for (int j = 0; j < K.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < K.getLocalNumOfCols(); ++i) {
        K2d(j, i) = mpiRank * j + (i * 35);
      }
    }
  }

  // run pressure2head
  CUAS::cavityOpenB(result, beta, v_b, K);

  // check
  {
    auto &K2d = K.getReadHandle();
    auto &result2d = result.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < K.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < K.getLocalNumOfCols(); ++i) {
        // check result
        ASSERT_EQ(result2d(j, i), beta * v_b * K2d(j, i));

        // check for sideeffects
        ASSERT_EQ(K2d(j, i), mpiRank * j + (i * 35));
      }
    }
  }
}

TEST(CUASKernelsTest, computeMelt) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid K(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradh2(GRID_SIZE_X, GRID_SIZE_Y);

  PetscScalar r = 1.0;
  PetscScalar bt = 0.1;

  // init T values
  auto T2d = T.getWriteHandle();
  PetscScalar *Tarr[GRID_SIZE_Y];
  PetscScalar TarrValues[GRID_SIZE_X] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    Tarr[i] = TarrValues;
  }

  for (int i = 0; i < T.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < T.getLocalNumOfCols(); ++j) {
      T2d(i, j) = Tarr[T.getCornerY() + i][T.getCornerX() + j];
    }
  }
  T2d.setValues();

  // init K values
  auto K2d = K.getWriteHandle();

  PetscScalar *Karr[GRID_SIZE_Y];
  PetscScalar KarrValues[GRID_SIZE_X] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    Karr[i] = KarrValues;
  }

  for (int i = 0; i < K.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < K.getLocalNumOfCols(); ++j) {
      K2d(i, j) = Karr[T.getCornerY() + i][T.getCornerX() + j];
    }
  }
  K2d.setValues();

  // init gradh2
  auto grad2d = gradh2.getWriteHandle();

  PetscScalar *gradArr[GRID_SIZE_Y];
  PetscScalar gradArrBeginningEnd[GRID_SIZE_X] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  PetscScalar gradArrMiddle[GRID_SIZE_X] = {0,        0.008281, 0.008281, 0.008281, 0.008281, 0.008281,
                                            0.008281, 0.008281, 0.008281, 0.008281, 0.008281, 0.008281,
                                            0.008281, 0.008281, 0.008281, 0.008281, 0.008281, 0.008281};

  gradArr[0] = gradArrBeginningEnd;

  for (int i = 1; i <= 6; ++i) {
    gradArr[i] = gradArrMiddle;
  }
  gradArr[7] = gradArrBeginningEnd;

  for (int i = 0; i < gradh2.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradh2.getLocalNumOfCols(); ++j) {
      grad2d(i, j) = gradArr[gradh2.getCornerY() + i][gradh2.getCornerX() + j];
    }
  }
  grad2d.setValues();

  CUAS::computeMelt(result, r, GRAVITY, RHO_WATER, T, K, gradh2, RHO_ICE, LATENT_HEAT, bt);

  auto res2d = result.getReadHandle();

  PetscScalar *resArr[GRID_SIZE_Y];
  PetscScalar resArrBeginningEnd[GRID_SIZE_X] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  PetscScalar resArrMiddle[GRID_SIZE_X] = {0,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449,
                                           0.0000005345568862275449};

  resArr[0] = resArrBeginningEnd;

  for (int i = 1; i <= 6; ++i) {
    resArr[i] = resArrMiddle;
  }
  resArr[7] = resArrBeginningEnd;

  for (int i = 0; i < result.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < result.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(res2d(i, j), resArr[result.getCornerY() + i][result.getCornerX() + j]);
    }
  }
}

TEST(CUASKernelsTest, binaryDialation) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid noFlowMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradMask(GRID_SIZE_X, GRID_SIZE_Y);

  bool *nfMask[GRID_SIZE_Y + 2];
  bool nfBeginningEnd[GRID_SIZE_X + 2];
  for (int i = 0; i < GRID_SIZE_X + 2; ++i) {
    nfBeginningEnd[i] = true;
  }
  bool nfMiddle[GRID_SIZE_X + 2];
  nfMiddle[0] = true;
  for (int i = 1; i < GRID_SIZE_X + 2; ++i) {
    nfMiddle[i] = false;
  }

  nfMask[0] = nfBeginningEnd;

  for (int i = 1; i < (GRID_SIZE_Y + 1); ++i) {
    nfMask[i] = nfMiddle;
  }

  nfMask[GRID_SIZE_Y + 1] = nfBeginningEnd;

  auto nf2d = noFlowMask.getWriteHandleGhost();
  for (int i = 0; i < noFlowMask.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < noFlowMask.getLocalGhostNumOfCols(); ++j) {
      nf2d(i, j) = nfMask[noFlowMask.getCornerY() + i][noFlowMask.getCornerX() + j];
    }
  }

  nf2d.setValues();

  CUAS::binaryDialation(gradMask, noFlowMask);

  bool *gMask[GRID_SIZE_Y];
  bool gBeginningEnd[GRID_SIZE_X];
  for (int i = 0; i < GRID_SIZE_X; ++i) {
    gBeginningEnd[i] = true;
  }
  bool gMiddle[GRID_SIZE_X];
  gMiddle[0] = true;
  for (int i = 1; i < GRID_SIZE_X; ++i) {
    gMiddle[i] = false;
  }

  gMask[0] = gBeginningEnd;

  for (int i = 1; i < GRID_SIZE_Y - 1; ++i) {
    gMask[i] = gMiddle;
  }

  gMask[GRID_SIZE_Y - 1] = gBeginningEnd;

  auto grad2d = gradMask.getReadHandle();
  for (int i = 0; i < gradMask.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradMask.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(grad2d(i, j), gMask[gradMask.getCornerY() + i][gradMask.getCornerX() + j]);
    }
  }
}

TEST(CUASKernelsTest, enableUnconfined) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid Teff(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid TeffPowTexp(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T_n(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid K(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid Sp(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid topg(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid u_n(GRID_SIZE_X, GRID_SIZE_Y);
  PetscScalar bt = 0.1;
  PetscScalar unconfSmooth = 0.0;
  PetscScalar Texp = 1;

  // init Teff
  auto Teff2d = Teff.getWriteHandle();
  PetscScalar *TeffArr[GRID_SIZE_Y];
  PetscScalar TeffArrValues[GRID_SIZE_X] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                            0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    TeffArr[i] = TeffArrValues;
  }
  for (int i = 0; i < Teff.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < Teff.getLocalNumOfCols(); ++j) {
      Teff2d(i, j) = TeffArr[Teff.getCornerY() + i][Teff.getCornerX() + j];
    }
  }
  Teff2d.setValues();

  // init T_n
  auto T_n2d = T_n.getWriteHandle();

  PetscScalar *T_nArr[GRID_SIZE_Y];
  PetscScalar T_nArrValues[GRID_SIZE_X] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                           0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    T_nArr[i] = T_nArrValues;
  }

  for (int i = 0; i < T_n.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < T_n.getLocalNumOfCols(); ++j) {
      T_n2d(i, j) = T_nArr[T_n.getCornerY() + i][T_n.getCornerX() + j];
    }
  }
  T_n2d.setValues();

  // topg is zero

  // init K
  auto K2d = K.getWriteHandle();

  PetscScalar *Karr[GRID_SIZE_Y];
  PetscScalar KarrValues[GRID_SIZE_X] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    Karr[i] = KarrValues;
  }

  for (int i = 0; i < K.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < K.getLocalNumOfCols(); ++j) {
      K2d(i, j) = Karr[K.getCornerY() + i][K.getCornerX() + j];
    }
  }
  K2d.setValues();

  // Sp is zero

  // init u_n
  auto u_n2d = u_n.getWriteHandleGhost();
  PetscScalar *u_nArr[GRID_SIZE_Y + 2];
  PetscScalar u_nArrValues[GRID_SIZE_X + 2] = {1820, 1729, 1638, 1547, 1456, 1365, 1274, 1183, 1092, 1001,
                                               910,  819,  728,  637,  546,  455,  364,  273,  182,  91};
  for (int i = 0; i < (GRID_SIZE_Y + 2); ++i) {
    u_nArr[i] = u_nArrValues;
  }

  for (int i = 0; i < u_n.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < u_n.getLocalGhostNumOfCols(); ++j) {
      u_n2d(i, j) = u_nArr[u_n.getCornerY() + i][u_n.getCornerX() + j];
    }
  }
  u_n2d.setValues();

  CUAS::enableUnconfined(Teff, TeffPowTexp, Sp, T_n, K, topg, u_n, Texp, unconfSmooth, bt);

  // check results
  auto resTeff2d = Teff.getReadHandle();
  PetscScalar *resTeffArr[GRID_SIZE_Y];
  PetscScalar resTeffArrValues[GRID_SIZE_X] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                               0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    resTeffArr[i] = resTeffArrValues;
  }

  for (int i = 0; i < Teff.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < Teff.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(resTeff2d(i, j), resTeffArr[Teff.getCornerY() + i][Teff.getCornerX() + j]);
    }
  }
}

TEST(CUASKernelsTest, calculateTeffPowTexp) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid Teff(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid TeffPowTexp(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid TeffPowTexpRes(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  // In NoData Texp is 1 but for testing puposes 3 has been chosen
  PetscScalar Texp = 3;

  // init T
  auto T2d = T.getWriteHandle();
  PetscScalar *Tarr[GRID_SIZE_Y];
  PetscScalar TarrValues[GRID_SIZE_X] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    Tarr[i] = TarrValues;
  }

  for (int i = 0; i < T.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < T.getLocalNumOfCols(); ++j) {
      T2d(i, j) = Tarr[T.getCornerY() + i][T.getCornerX() + j];
    }
  }
  T2d.setValues();

  CUAS::calculateTeffPowTexp(Teff, TeffPowTexp, T, Texp);

  // init TeffPowTexpRes
  PetscScalar *TeffPowTexpResArr[GRID_SIZE_Y];
  PetscScalar TeffPowTexpResArrValues[GRID_SIZE_X] = {0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008,
                                                      0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    TeffPowTexpResArr[i] = TeffPowTexpResArrValues;
  }

  // compare results
  auto TeffPowOfFunction = TeffPowTexp.getReadHandle();
  for (int i = 0; i < TeffPowTexp.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < TeffPowTexp.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(TeffPowOfFunction(i, j),
                       TeffPowTexpResArr[TeffPowTexp.getCornerY() + i][TeffPowTexp.getCornerX() + j]);
    }
  }
}

TEST(CUASKernelsTest, calculateSeValues) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid S(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid Se(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid Sp(GRID_SIZE_X, GRID_SIZE_Y);

  // init S
  auto S2d = S.getWriteHandle();
  PetscScalar *Sarr[GRID_SIZE_Y];
  PetscScalar SarrValues[GRID_SIZE_X];
  for (int i = 0; i < GRID_SIZE_X; ++i) {
    SarrValues[i] = 0.0000982977696;
  }

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    Sarr[i] = SarrValues;
  }

  for (int i = 0; i < S.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < S.getLocalNumOfCols(); ++j) {
      S2d(i, j) = Sarr[S.getCornerY() + i][S.getCornerX() + j];
    }
  }
  S2d.setValues();

  // init Sp
  auto Sp2d = Sp.getWriteHandle();
  PetscScalar *SpArr[GRID_SIZE_Y];
  PetscScalar SpArrValues[GRID_SIZE_X];
  for (int i = 0; i < GRID_SIZE_X; ++i) {
    SpArrValues[i] = i;
  }

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    SpArr[i] = SpArrValues;
  }

  for (int i = 0; i < Sp.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < Sp.getLocalNumOfCols(); ++j) {
      Sp2d(i, j) = SpArr[Sp.getCornerY() + i][Sp.getCornerX() + j];
    }
  }
  Sp2d.setValues();

  CUAS::calculateSeValues(Se, Sp, S);

  // compare results
  PetscScalar *SeArr[GRID_SIZE_Y];
  PetscScalar SeValues[GRID_SIZE_X];

  for (int i = 0; i < GRID_SIZE_X; ++i) {
    SeValues[i] = (0.0000982977696) + i;
  }

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    SeArr[i] = SeValues;
  }

  auto Se2d = Se.getReadHandle();
  for (int i = 0; i < Se.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < Se.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(Se2d(i, j), SeArr[Se.getCornerY() + i][Se.getCornerX() + j]);
    }
  }
}

TEST(CUASKernelsTest, doChannels) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid u_n(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid melt(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid creep(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T_n(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid pIce(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid topg(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid K(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid noFlowMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid cavityOpening(GRID_SIZE_X, GRID_SIZE_Y);
  // init args
  PetscScalar flowConstant = 5e-25;
  PetscScalar roughnessFactor = 1.0;
  PetscScalar Texp = 1;
  PetscScalar noSmoothMelt = false;
  PetscScalar cavityBeta = 0.0005;
  PetscScalar basalVelocityIce = 1e-06;
  PetscScalar bt = 0.1;
  PetscScalar dx = 1000.0;
  PetscScalar dtSecs = 43200;
  PetscScalar tMin = 1e-07;
  PetscScalar tMax = 20.0;

  // init u_n
  auto u_n2d = u_n.getWriteHandleGhost();
  PetscScalar *u_nArr[GRID_SIZE_Y + 2];
  PetscScalar u_nArrValues[GRID_SIZE_X + 2] = {1820, 1729, 1638, 1547, 1456, 1365, 1274, 1183, 1092, 1001,
                                               910,  819,  728,  637,  546,  455,  364,  273,  182,  91};
  for (int i = 0; i < GRID_SIZE_Y + 2; ++i) {
    u_nArr[i] = u_nArrValues;
  }

  for (int i = 0; i < u_n.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < u_n.getLocalGhostNumOfCols(); ++j) {
      u_n2d(i, j) = u_nArr[u_n.getCornerY() + i][u_n.getCornerX() + j];
    }
  }
  u_n2d.setValues();

  // init T values
  auto T2d = T.getWriteHandle();
  PetscScalar *Tarr[GRID_SIZE_Y];
  PetscScalar TarrValues[GRID_SIZE_X] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    Tarr[i] = TarrValues;
  }

  for (int i = 0; i < T.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < T.getLocalNumOfCols(); ++j) {
      T2d(i, j) = Tarr[T.getCornerY() + i][T.getCornerX() + j];
    }
  }
  T2d.setValues();

  // init T_n
  auto T_n2d = T_n.getWriteHandle();

  PetscScalar *T_nArr[GRID_SIZE_Y];
  PetscScalar T_nArrValues[GRID_SIZE_X] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                           0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    T_nArr[i] = T_nArrValues;
  }

  for (int i = 0; i < T_n.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < T_n.getLocalNumOfCols(); ++j) {
      T_n2d(i, j) = T_nArr[T_n.getCornerY() + i][T_n.getCornerX() + j];
    }
  }
  T_n2d.setValues();

  // init K values
  auto K2d = K.getWriteHandle();

  PetscScalar *Karr[GRID_SIZE_Y];
  PetscScalar KarrValues[GRID_SIZE_X] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    Karr[i] = KarrValues;
  }

  for (int i = 0; i < K.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < K.getLocalNumOfCols(); ++j) {
      K2d(i, j) = Karr[K.getCornerY() + i][K.getCornerX() + j];
    }
  }
  K2d.setValues();

  // init noFlowMask
  auto noFlowMask2d = noFlowMask.getWriteHandleGhost();

  bool *nFMask[GRID_SIZE_Y + 2];
  bool nFBeginningEnd[GRID_SIZE_X + 2];
  for (int i = 0; i < GRID_SIZE_X + 2; ++i) {
    nFBeginningEnd[i] = true;
  }
  bool nFMiddle[GRID_SIZE_X + 2];
  nFMiddle[0] = true;
  for (int i = 1; i < GRID_SIZE_X + 2; ++i) {
    nFMiddle[i] = false;
  }

  nFMask[0] = nFBeginningEnd;

  for (int i = 1; i < GRID_SIZE_Y + 1; ++i) {
    nFMask[i] = nFMiddle;
  }

  nFMask[GRID_SIZE_Y + 1] = nFBeginningEnd;

  for (int i = 0; i < noFlowMask.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < noFlowMask.getLocalGhostNumOfCols(); ++j) {
      noFlowMask2d(i, j) = nFMask[noFlowMask.getCornerY() + i][noFlowMask.getCornerX() + j];
    }
  }
  noFlowMask2d.setValues();

  // set cavityOpening to zero
  cavityOpening.setZero();

  // init gradmask
  auto gradMask2d = gradMask.getWriteHandleGhost();

  PetscScalar *gradMaskArr[GRID_SIZE_Y + 2];
  PetscScalar gradMaskBeginningEnd[GRID_SIZE_X + 2];
  for (int i = 0; i < GRID_SIZE_X + 2; ++i) {
    gradMaskBeginningEnd[i] = 1;
  }

  PetscScalar gradMaskMiddle[GRID_SIZE_X + 2];
  gradMaskMiddle[0] = 1;
  gradMaskMiddle[1] = 1;
  for (int i = 2; i < GRID_SIZE_X + 2; ++i) {
    gradMaskMiddle[i] = 0;
  }

  gradMaskArr[0] = gradMaskBeginningEnd;
  gradMaskArr[1] = gradMaskBeginningEnd;
  for (int i = 2; i < GRID_SIZE_Y; ++i) {
    gradMaskArr[i] = gradMaskMiddle;
  }
  gradMaskArr[GRID_SIZE_Y] = gradMaskBeginningEnd;
  gradMaskArr[GRID_SIZE_Y + 1] = gradMaskBeginningEnd;

  for (int i = 0; i < gradMask.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < gradMask.getLocalGhostNumOfCols(); ++j) {
      gradMask2d(i, j) = gradMaskArr[gradMask.getCornerY() + i][gradMask.getCornerX() + j];
    }
  }
  gradMask2d.setValues();

  // init p_ice
  auto pIce2d = pIce.getWriteHandleGhost();
  PetscScalar *pIceArr[GRID_SIZE_Y + 2];
  PetscScalar pIceArrValues[GRID_SIZE_X + 2] = {17854200, 16961490, 16068780, 15176070, 14283360, 13390650, 12497940,
                                                11605230, 10712520, 9819810,  8927100,  8034390,  7141680,  6248970,
                                                5356260,  4463550,  3570840,  2678130,  1785420,  892710};
  for (int i = 0; i < GRID_SIZE_Y + 2; ++i) {
    pIceArr[i] = pIceArrValues;
  }

  for (int i = 0; i < pIce.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < pIce.getLocalGhostNumOfCols(); ++j) {
      pIce2d(i, j) = pIceArr[pIce.getCornerY() + i][pIce.getCornerX() + j];
    }
  }
  pIce2d.setValues();

  // topg is zero
  topg.setZero();

  CUAS::doChannels(melt, creep, u_n, gradMask, T, T_n, pIce, topg, K, noFlowMask, cavityOpening, flowConstant, Texp,
                   roughnessFactor, noSmoothMelt, cavityBeta, basalVelocityIce, tMin, tMax, bt, dx, dtSecs);

  // compare results
  PetscScalar *meltArr[GRID_SIZE_Y];
  PetscScalar meltFirstRow[GRID_SIZE_X] = {0,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312,
                                           0.00000006681961077844312};

  PetscScalar meltSecondRow[GRID_SIZE_X] = {
      0.00000006681961077844312, 0.0000004009176646706587, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004677372754491018, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004677372754491018, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004677372754491018, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004009176646706587};

  PetscScalar meltMiddleRow[GRID_SIZE_X] = {
      0.00000006681961077844312,      0.0000004677372754491018,       0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.00000046773727544910184};

  meltArr[0] = meltFirstRow;
  meltArr[1] = meltSecondRow;

  for (int i = 2; i < GRID_SIZE_Y - 2; ++i) {
    meltArr[i] = meltMiddleRow;
  }

  meltArr[GRID_SIZE_Y - 2] = meltSecondRow;
  meltArr[GRID_SIZE_Y - 1] = meltFirstRow;

  auto melt2d = melt.getReadHandle();
  for (int i = 0; i < melt.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < melt.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(melt2d(i, j), meltArr[melt.getCornerY() + i][melt.getCornerX() + j]);
    }
  }

  PetscScalar *creepArr[GRID_SIZE_Y];
  PetscScalar creepValues[GRID_SIZE_X];
  for (int i = 0; i < GRID_SIZE_X; ++i) {
    creepValues[i] = 0.0000000000000000069931566;
  }

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    creepArr[i] = creepValues;
  }

  auto creep2d = creep.getReadHandle();
  for (int i = 0; i < creep.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < creep.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(creep2d(i, j), creepArr[creep.getCornerY() + i][creep.getCornerX() + j]);
    }
  }

  PetscScalar *cavityArr[GRID_SIZE_Y];
  PetscScalar cavityMiddleRow[GRID_SIZE_X];
  for (int i = 0; i < GRID_SIZE_X; ++i) {
    cavityMiddleRow[i] = 0.000000005;
  }

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    cavityArr[i] = cavityMiddleRow;
  }

  auto cavity2d = cavityOpening.getReadHandle();
  for (int i = 0; i < cavityOpening.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < cavityOpening.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(cavity2d(i, j), cavityArr[cavityOpening.getCornerY() + i][cavityOpening.getCornerX() + j]);
    }
  }
  // check updated T
  PetscScalar *updatedTArr[GRID_SIZE_Y];
  PetscScalar updatedTFirstRow[GRID_SIZE_X] = {
      0.2002159999996979,  0.20310260718532666, 0.20310260718532666, 0.20310260718532666, 0.20310260718532666,
      0.20310260718532666, 0.20310260718532666, 0.20310260718532666, 0.20310260718532666, 0.20310260718532666,
      0.20310260718532666, 0.20310260718532666, 0.20310260718532666, 0.20310260718532666, 0.20310260718532666,
      0.20310260718532666, 0.20310260718532666, 0.20310260718532666};
  PetscScalar updatedTSecondRow[GRID_SIZE_X] = {
      0.20310260718532666, 0.21753564311347037, 0.2204222502990991, 0.2204222502990991, 0.2204222502990991,
      0.2204222502990991,  0.2204222502990991,  0.2204222502990991, 0.2204222502990991, 0.2204222502990991,
      0.2204222502990991,  0.2204222502990991,  0.2204222502990991, 0.2204222502990991, 0.2204222502990991,
      0.2204222502990991,  0.2204222502990991,  0.21753564311347037};
  PetscScalar updatedTThirdRow[GRID_SIZE_X] = {
      0.20310260718532666, 0.2204222502990991,  0.22330885748472784, 0.22330885748472784, 0.22330885748472784,
      0.22330885748472784, 0.22330885748472784, 0.22330885748472784, 0.22330885748472784, 0.22330885748472784,
      0.22330885748472784, 0.22330885748472784, 0.22330885748472784, 0.22330885748472784, 0.22330885748472784,
      0.22330885748472784, 0.22330885748472784, 0.22042225029909912};

  updatedTArr[0] = updatedTFirstRow;
  updatedTArr[1] = updatedTSecondRow;

  for (int i = 2; i < GRID_SIZE_Y - 2; ++i) {
    updatedTArr[i] = updatedTThirdRow;
  }

  updatedTArr[GRID_SIZE_Y - 2] = updatedTSecondRow;
  updatedTArr[GRID_SIZE_Y - 1] = updatedTFirstRow;

  auto updatedTGlobal = T.getReadHandle();
  for (int i = 0; i < T.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < T.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(updatedTGlobal(i, j), updatedTArr[T.getCornerY() + i][T.getCornerX() + j]);
    }
  }
}

// noChannels just sets both passed grids to zero
// TEST(CUASKernelsTest, noChannels) { ASSERT_EQ(mpiSize, MPI_SIZE); }

TEST(CUASKernelsTest, convolve) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid melt(GRID_SIZE_X, GRID_SIZE_Y);
  {
    auto melt2d = melt.getWriteHandle();

    PetscScalar *meltArr[GRID_SIZE_Y];
    PetscScalar meltArrBeginningEnd[GRID_SIZE_X] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    PetscScalar meltArrMiddle[GRID_SIZE_X] = {0,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449,
                                              0.0000005345568862275449};

    meltArr[0] = meltArrBeginningEnd;

    for (int i = 1; i <= (GRID_SIZE_Y - 2); ++i) {
      meltArr[i] = meltArrMiddle;
    }
    meltArr[GRID_SIZE_Y - 1] = meltArrBeginningEnd;

    for (int i = 0; i < melt.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < melt.getLocalNumOfCols(); ++j) {
        melt2d(i, j) = meltArr[melt.getCornerY() + i][melt.getCornerX() + j];
      }
    }
  }

  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);

  CUAS::convolveStar11411(melt, result);

  // compare results
  PetscScalar *resultArr[GRID_SIZE_Y];
  PetscScalar resultFirstRow[GRID_SIZE_X] = {
      0,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
      0.00000006681961077844312,
  };

  PetscScalar resultSecondRow[GRID_SIZE_X] = {
      0.00000006681961077844312, 0.0000004009176646706587, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004677372754491018, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004677372754491018, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004677372754491018, 0.0000004677372754491018, 0.0000004677372754491018,
      0.0000004677372754491018,  0.0000004009176646706587,
  };

  PetscScalar resultMiddleRow[GRID_SIZE_X] = {
      0.00000006681961077844312,      0.0000004677372754491018,       0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.0000005345568862275449275449,
      0.0000005345568862275449275449, 0.0000005345568862275449275449, 0.00000046773727544910184,
  };

  resultArr[0] = resultFirstRow;
  resultArr[1] = resultSecondRow;

  for (int i = 2; i <= GRID_SIZE_Y - 3; ++i) {
    resultArr[i] = resultMiddleRow;
  }

  resultArr[GRID_SIZE_Y - 2] = resultSecondRow;
  resultArr[GRID_SIZE_Y - 1] = resultFirstRow;

  {
    auto res2d = result.getReadHandle();

    for (int i = 0; i < result.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < result.getLocalNumOfCols(); ++j) {
        EXPECT_DOUBLE_EQ(res2d(i, j), resultArr[result.getCornerY() + i][result.getCornerX() + j]);
      }
    }
  }
}

TEST(CUASKernelsTest, clamp) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid input(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid result(GRID_SIZE_Y, GRID_SIZE_X);

  PetscScalar minimum = 1.0;
  PetscScalar maximum = 100.0;
  // setup
  {
    auto inputGlobal = input.getWriteHandle();
    for (int j = 0; j < input.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < input.getLocalNumOfCols(); ++i) {
        if (i % 3 == 0) {
          inputGlobal(j, i) = 110.0;
        } else if (i % 2 == 0) {
          inputGlobal(j, i) = -3.14;
        } else {
          inputGlobal(j, i) = 10;
        }
      }
    }
  }

  // clamp
  CUAS::clamp(input, minimum, maximum);

  // check
  {
    auto inputGlobal = input.getReadHandle();
    for (int j = 0; j < input.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < input.getLocalNumOfCols(); ++i) {
        ASSERT_EQ(inputGlobal(j, i) >= minimum, true);
        ASSERT_EQ(inputGlobal(j, i) <= maximum, true);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
