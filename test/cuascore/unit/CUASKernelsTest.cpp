/**
 * File: CUASKernelsTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASKernels.h"

#include "PETScGrid.h"
#include <cmath>
#include <cstdlib>

// #define TESTS_DUMP_NETCDF
#ifdef TESTS_DUMP_NETCDF
#include "NetCDFFile.h"
#endif

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 9
#define GRID_SIZE_X 18
#define GRID_SIZE_Y 8

TEST(CUASKernelsTest, headToPressure) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid pressure(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid bedElevation(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid head(GRID_SIZE_Y, GRID_SIZE_X);

  // We don't need to ensure 0 <= level <= waterLayerThickness here
  constexpr PetscScalar level = 33.5;

  // setup
  {
    auto pressure2d = pressure.getWriteHandle();
    auto bedElevation2d = bedElevation.getWriteHandle();
    auto head2d = head.getWriteHandle();
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        pressure2d(j, i) = level * mpiRank - 2.3;
        bedElevation2d(j, i) = level * level - 31.3;
        head2d(j, i) = mpiRank * j + (i * 35);
      }
    }
  }

  // run headToPressure
  CUAS::headToPressure(pressure, head, bedElevation, level);

  // check
  {
    auto &pressure2d = pressure.getReadHandle();
    auto &bedElevation2d = bedElevation.getReadHandle();
    auto &head2d = head.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        // check result
        ASSERT_EQ(pressure2d(j, i), RHO_WATER * GRAVITY * (head2d(j, i) - bedElevation2d(j, i) - level));
        // check for side effects
        ASSERT_EQ(bedElevation2d(j, i), level * level - 31.3);
        ASSERT_EQ(head2d(j, i), mpiRank * j + (i * 35));
      }
    }
  }
}

TEST(CUASKernelsTest, pressureToHead) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid pressure(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid bedElevation(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid head(GRID_SIZE_Y, GRID_SIZE_X);

  // We don't need to ensure 0 <= level <= waterLayerThickness here
  constexpr PetscScalar level = 33.5;

  // setup
  {
    auto pressure2d = pressure.getWriteHandle();
    auto bedElevation2d = bedElevation.getWriteHandle();
    auto head2d = head.getWriteHandle();
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        pressure2d(j, i) = level * mpiRank - 2.3;
        bedElevation2d(j, i) = level * level - 31.3;
        head2d(j, i) = mpiRank * j + (i * 35);
      }
    }
  }

  // run pressureToHead
  CUAS::pressureToHead(head, pressure, bedElevation, level);

  // check
  {
    auto &pressure2d = pressure.getReadHandle();
    auto &bedElevation2d = bedElevation.getReadHandle();
    auto &head2d = head.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        // check result
        ASSERT_EQ(head2d(j, i), pressure2d(j, i) / (RHO_WATER * GRAVITY) + bedElevation2d(j, i) + level);
        // check for side effects
        ASSERT_EQ(bedElevation2d(j, i), level * level - 31.3);
        ASSERT_EQ(pressure2d(j, i), level * mpiRank - 2.3);
      }
    }
  }
}

TEST(CUASKernelsTest, pressureToHead_ToPressure) {
  /*
   * Test if pressureToHead() and headToPressure() are the inverse operation
   * of each other. Test with default and given layer thickness as last argument to both
   * methods.
   */

  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid pressure(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid head(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid bedElevation(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid result(GRID_SIZE_Y, GRID_SIZE_X);

  // setup
  {
    auto pressure2d = pressure.getWriteHandle();
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        pressure2d(j, i) = 1.1e5 * mpiRank;
      }
    }
  }

  // run pressureToHead and headToPressure with defaults
  CUAS::pressureToHead(head, pressure, bedElevation);
  CUAS::headToPressure(result, head, bedElevation);

  // check initial pressure equals final pressure
  {
    auto &p_org = pressure.getReadHandle();
    auto &p_new = result.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        // check result
        ASSERT_NEAR(p_new(j, i), p_org(j, i), 1e-8);
      }
    }
  }

  // run pressureToHead and headToPressure with given level within aquifer
  constexpr PetscScalar z_w = 42.0;
  CUAS::pressureToHead(head, pressure, bedElevation, z_w);
  CUAS::headToPressure(result, head, bedElevation, z_w);

  // check initial pressure equals final pressure
  {
    auto &p_org = pressure.getReadHandle();
    auto &p_new = result.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
        // check result
        ASSERT_NEAR(p_new(j, i), p_org(j, i), 1e-8);
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

TEST(CUASKernelsTest, computeCavityOpening) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid basalVelocity(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid result(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid bndMask(GRID_SIZE_Y, GRID_SIZE_X);

  PetscScalar beta = 33.5;
  PetscScalar K = 17.5;

  // setup
  {
    auto vb = basalVelocity.getWriteHandle();
    for (int j = 0; j < basalVelocity.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < basalVelocity.getLocalNumOfCols(); ++i) {
        vb(j, i) = mpiRank * j + (i * 35);
      }
    }
  }

  CUAS::computeCavityOpening(result, beta, K, basalVelocity, bndMask);

  // check
  {
    auto &vb = basalVelocity.getReadHandle();
    auto &result2d = result.getReadHandle();
    // TODO: check ghost cells
    for (int j = 0; j < basalVelocity.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < basalVelocity.getLocalNumOfCols(); ++i) {
        // check result
        ASSERT_EQ(result2d(j, i), beta * K * vb(j, i));

        // check for sideeffects
        ASSERT_EQ(vb(j, i), mpiRank * j + (i * 35));
      }
    }
  }
}

TEST(CUASKernelsTest, computeMeltOpening) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid effTransNorth(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid effTransEast(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid effTransSouth(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid effTransWest(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid bndMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid topg(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid hydraulicHead(GRID_SIZE_X, GRID_SIZE_Y);

  constexpr PetscScalar r = 1.0;
  constexpr PetscScalar bt = 0.1;
  constexpr PetscScalar K = 10;
  constexpr PetscScalar dx = 1000.0;

  bndMask.setRealBoundary(DIRICHLET_FLAG);
  topg.setConst(0.0);

  // init T values
  T.setConst(0.2);
  effTransNorth.setConst(0.2);
  effTransEast.setConst(0.2);
  effTransSouth.setConst(0.2);
  effTransWest.setConst(0.2);

  // init hydraulicHead
  auto head = hydraulicHead.getWriteHandleGhost();
  PetscScalar *h_Arr[GRID_SIZE_Y + 2];
  PetscScalar h_ArrValues[GRID_SIZE_X + 2] = {1820, 1729, 1638, 1547, 1456, 1365, 1274, 1183, 1092, 1001,
                                              910,  819,  728,  637,  546,  455,  364,  273,  182,  91};
  for (int i = 0; i < GRID_SIZE_Y + 2; ++i) {
    h_Arr[i] = h_ArrValues;
  }

  for (int i = 0; i < hydraulicHead.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < hydraulicHead.getLocalGhostNumOfCols(); ++j) {
      head(i, j) = h_Arr[hydraulicHead.getCornerY() + i][hydraulicHead.getCornerX() + j];
    }
  }
  head.setValues();

  // CUAS::computeMeltOpeningDeprecated(result, r, K, T, gradh2);
  CUAS::computeMeltOpening(result, r, K, dx, effTransEast, effTransWest, effTransNorth, effTransSouth, hydraulicHead,
                           topg, bndMask);

  auto &res2d = result.getReadHandle();
  auto &mask = bndMask.getReadHandle();

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
      if (mask(i, j) == COMPUTE_FLAG) {
        EXPECT_DOUBLE_EQ(res2d(i, j), resArr[result.getCornerY() + i][result.getCornerX() + j]);
      }
    }
  }
}

TEST(CUASKernelsTest, binaryDilation) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid genericMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradMask(GRID_SIZE_X, GRID_SIZE_Y);

  // + 2 to account for the ghosts
  bool *nfMask[GRID_SIZE_Y + 2];
  bool nfBeginningEnd[GRID_SIZE_X + 2];
  for (int i = 0; i < GRID_SIZE_X + 2; ++i) {
    nfBeginningEnd[i] = true;
  }
  bool nfMiddle[GRID_SIZE_X + 2];
  nfMiddle[0] = true;  // this is outside the visible CUAS range in NetCDF
  for (int i = 1; i < GRID_SIZE_X + 2; ++i) {
    nfMiddle[i] = false;
  }

  // combine the rows
  nfMask[0] = nfBeginningEnd;
  for (int i = 1; i < (GRID_SIZE_Y + 1); ++i) {
    nfMask[i] = nfMiddle;
  }
  nfMask[GRID_SIZE_Y + 1] = nfBeginningEnd;

  auto mask = genericMask.getWriteHandleGhost();
  for (int i = 0; i < genericMask.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < genericMask.getLocalGhostNumOfCols(); ++j) {
      mask(i, j) = nfMask[genericMask.getCornerY() + i][genericMask.getCornerX() + j];
    }
  }
  mask.setValues();  // needed if in "{ ... }"?

  CUAS::binaryDilation(gradMask, genericMask);

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

  // combine the rows
  gMask[0] = gBeginningEnd;
  for (int i = 1; i < GRID_SIZE_Y - 1; ++i) {
    gMask[i] = gMiddle;
  }
  gMask[GRID_SIZE_Y - 1] = gBeginningEnd;

  auto &grad2d = gradMask.getReadHandle();
  for (int i = 0; i < gradMask.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradMask.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(grad2d(i, j), gMask[gradMask.getCornerY() + i][gradMask.getCornerX() + j]);
    }
  }
}

TEST(CUASKernelsTest, getFlux) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  // T = Teff = T_e = T_w = T_n = T_s = const.
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid head(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid fluxM(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid fluxX(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid fluxY(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid fluxM_analytical(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid fluxX_analytical(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid fluxY_analytical(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid bndMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradHeadSquared(GRID_SIZE_X, GRID_SIZE_Y);

  constexpr PetscScalar Tconst = 0.2;
  constexpr PetscScalar dx = 1.0;
  constexpr PetscScalar fillValue = 0.0;  // or NC_FILL_DOUBLE

  // initialize fields
  T.setConst(Tconst);
  fluxM.setConst(fillValue);
  fluxX.setConst(fillValue);
  fluxY.setConst(fillValue);
  fluxM_analytical.setConst(fillValue);
  fluxX_analytical.setConst(fillValue);
  fluxY_analytical.setConst(fillValue);

  bndMask.setRealBoundary(DIRICHLET_FLAG);

  // init head
  {
    auto cornerX = head.getCornerX();
    auto cornerY = head.getCornerY();
    auto h = head.getWriteHandle();
    auto fl = fluxM_analytical.getWriteHandle();
    auto flx = fluxX_analytical.getWriteHandle();
    auto fly = fluxY_analytical.getWriteHandle();
    auto gradh2 = gradHeadSquared.getWriteHandle();
    auto &mask = bndMask.getReadHandle();

    for (int row = 0; row < head.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < head.getLocalNumOfCols(); ++col) {
        const PetscScalar x = (cornerX + col) * dx;
        const PetscScalar y = (cornerY + row) * dx;
        h(row, col) = x * x + y * y;
        if (mask(row, col) == (PetscScalar)COMPUTE_FLAG) {
          // analytical solution:  flux = T^exp * |grad h|
          gradh2(row, col) = 4.0 * (x * x + y * y);
          fl(row, col) = Tconst * PetscSqrtScalar(gradh2(row, col));
          flx(row, col) = -Tconst * 2.0 * x;
          fly(row, col) = -Tconst * 2.0 * y;
        }
      }
    }
  }

  CUAS::getFlux(fluxM, fluxX, fluxY, bndMask, head, T, T, T, T, dx);

#ifdef TESTS_DUMP_NETCDF
  // Gets information about the currently running test.
  // Do NOT delete the returned object - it's managed by the UnitTest class.
  const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
  auto filename = std::string(test_info->test_suite_name())
                      .append(std::string("_-_"))
                      .append(std::string(test_info->name()))
                      .append(std::string(".nc"));
  CUAS::NetCDFFile file(filename, GRID_SIZE_X, GRID_SIZE_Y);
  file.defineGrid("T", LIMITED);
  file.defineGrid("head", LIMITED);
  file.defineGrid("fluxMagnitude", LIMITED);
  file.defineGrid("fluxXDir", LIMITED);
  file.defineGrid("fluxYDir", LIMITED);
  file.defineGrid("fluxMagnitude_analytical", LIMITED);
  file.defineGrid("fluxXDir_analytical", LIMITED);
  file.defineGrid("fluxYDir_analytical", LIMITED);
  file.write("T", T, 0);
  file.write("head", head, 0);
  file.write("fluxMagnitude", fluxM, 0);
  file.write("fluxXDir", fluxX, 0);
  file.write("fluxYDir", fluxY, 0);

  file.write("fluxMagnitude_analytical", fluxM_analytical, 0);
  file.write("fluxXDir_analytical", fluxX_analytical, 0);
  file.write("fluxYDir_analytical", fluxY_analytical, 0);

#endif

  // compare results
  {
    auto &f = fluxM.getReadHandle();
    auto &f_analytical = fluxM_analytical.getReadHandle();
    auto &mask = bndMask.getReadHandle();
    for (int row = 0; row < fluxM.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < fluxM.getLocalNumOfCols(); ++col) {
        if (mask(row, col) == COMPUTE_FLAG) {
          EXPECT_DOUBLE_EQ(f(row, col), f_analytical(row, col));
        } else {
          EXPECT_DOUBLE_EQ(f(row, col), 0.0);
        }
      }
    }
  }
}

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
    auto &res2d = result.getReadHandle();

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
    auto &inputGlobal = input.getReadHandle();
    for (int j = 0; j < input.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < input.getLocalNumOfCols(); ++i) {
        ASSERT_EQ(inputGlobal(j, i) >= minimum, true);
        ASSERT_EQ(inputGlobal(j, i) <= maximum, true);
      }
    }
  }
}

TEST(CUASKernelsTest, getEffectiveAquiferProperties) {
  /*
   * Note, we only test inside the boundary, because boundary conditions are applied elsewhere.
   */

  ASSERT_EQ(mpiSize, MPI_SIZE);

  // input
  PETScGrid head(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid transmissivity(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid topg(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid bndMask(GRID_SIZE_X, GRID_SIZE_Y);
  // output
  PETScGrid Teff(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid Seff(GRID_SIZE_X, GRID_SIZE_Y);

  // config confined-case
  constexpr PetscScalar specificStorage = 1e-3;  // unit: m^-1, selected for testing only; usually much lower
  constexpr PetscScalar specificYield = 0.3;     // unit: none, selected for testing only; usually 0.4
  constexpr PetscScalar Tinit = 10.0;            // unit: m^2/s
  constexpr PetscScalar hinit = 100.0;           // unit: m

  head.setConst(hinit);            // unit: m, > layer thickness->confined
  transmissivity.setConst(Tinit);  // unit: m^2/s
  topg.setConst(0.0);              // unit: m
  bndMask.setConst(COMPUTE_FLAG);
  bndMask.setRealBoundary(DIRICHLET_FLAG);

  // confined case: head > layer thickness
  {
    constexpr PetscScalar layerThickness = 0.1 * hinit;                  // unit: m
    constexpr PetscScalar unconfinedSmooth = 0.5 * layerThickness;       // not used in the confined case
    constexpr PetscScalar Seff_test = specificStorage * layerThickness;  // Ehlig & Halepaska 1976, eq. 7
    constexpr PetscScalar Teff_test = Tinit;                             // Ehlig & Halepaska 1976, eq. 5

    CUAS::getEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                        specificStorage, specificYield, unconfinedSmooth);

    // now test Teff
    auto &result1 = Teff.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), Teff_test) << "at i=" << i << ", j=" << j;
      }
    }
    // now test Seff
    auto &result2 = Seff.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), Seff_test) << "at i=" << i << ", j=" << j;
      }
    }
  }

  // no netcdf dump for the trivial confined case

  // unconfined case: head <= layer thickness and head < layer thickness - unconfinedSmooth
  {
    constexpr PetscScalar layerThickness = 2.0 * hinit;             // unit: m
    constexpr PetscScalar unconfinedSmooth = 0.5 * layerThickness;  // less than the layer thickness
    constexpr PetscScalar psi = hinit;
    constexpr PetscScalar Seff_test = specificStorage * layerThickness + specificYield;
    constexpr PetscScalar Teff_test = Tinit / layerThickness * psi;

    CUAS::getEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                        specificStorage, specificYield, unconfinedSmooth);

    // now test Teff
    auto &result1 = Teff.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), Teff_test) << "at i=" << i << ", j=" << j;
      }
    }
    // now test Seff
    auto &result2 = Seff.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), Seff_test) << "at i=" << i << ", j=" << j;
      }
    }
  }

  // unconfined case: head <= layer thickness and head > layer thickness - unconfinedSmooth
  {
    /*
     * Note, S'(h) = S'/d * (b - h) (Ehlig & Halepaska 1976, eq. 7), where S' = Sy in CUAS
     * We use b = 2 * h0, d = 1.5 * h0 and Sy = 0.3, and therefore we get:
     *       S'(h0) = Sy/(1.5*h0)*(2*h0-h0) = Sy/1.5 = 0.2
     *       Seff = Ss*b + S'(h0) = Ss * 2 * h0 + Sy/1.5 = 0.4
     */
    constexpr PetscScalar layerThickness = 2.0 * hinit;    // unit: m
    constexpr PetscScalar unconfinedSmooth = 1.5 * hinit;  // the head is now within the smoothing region
    constexpr PetscScalar psi = hinit;
    constexpr PetscScalar Seff_test = specificStorage * layerThickness + specificYield / 1.5;
    constexpr PetscScalar Teff_test = Tinit / layerThickness * psi;

    CUAS::getEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                        specificStorage, specificYield, unconfinedSmooth);

    // now test Teff
    auto &result1 = Teff.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), Teff_test) << "at i=" << i << ", j=" << j;
      }
    }
    // now test Seff
    auto &result2 = Seff.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), Seff_test) << "at i=" << i << ", j=" << j;
      }
    }
  }

  // head below bedrock: This should not happen in a simulation, but sometimes it happens and thus,
  // we need to consider this.
  {
    constexpr PetscScalar layerThickness = 2.0 * hinit;    // unit: m
    constexpr PetscScalar unconfinedSmooth = 1.5 * hinit;  // the head is now within the smoothing region
    constexpr PetscScalar psi = -10.0;                     // head below bedrock (topg=0)
    head.setConst(psi);                                    // the head is now below zero and below the smoothing region

    constexpr PetscScalar Seff_test =
        specificStorage * layerThickness + specificYield;  // max value that can be reached
    constexpr PetscScalar Teff_test = NOFLOW_VALUE;

    // fixme:
    CUAS::getEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                        specificStorage, specificYield, unconfinedSmooth);

    // now test Teff
    auto &result1 = Teff.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), Teff_test) << "at i=" << i << ", j=" << j;
      }
    }
    // now test Seff
    auto &result2 = Seff.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), Seff_test) << "at i=" << i << ", j=" << j;
      }
    }
  }

#ifdef TESTS_DUMP_NETCDF
  // Gets information about the currently running test.
  // Do NOT delete the returned object - it's managed by the UnitTest class.
  const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
  auto filename = std::string(test_info->test_suite_name())
                      .append(std::string("_-_"))
                      .append(std::string(test_info->name()))
                      .append(std::string(".nc"));

  CUAS::NetCDFFile file(filename, GRID_SIZE_X, GRID_SIZE_Y);
  file.defineGrid("bndMask", LIMITED);
  file.defineGrid("transmissivity", LIMITED);
  file.defineGrid("head", LIMITED);
  file.defineGrid("effective_transmissivity", LIMITED);
  file.defineGrid("effective_storativity", LIMITED);

  file.write("bndMask", bndMask, 0);
  file.write("transmissivity", transmissivity, 0);
  file.write("head", head, 0);
  file.write("effective_transmissivity", Teff, 0);
  file.write("effective_storativity", Seff, 0);

#endif
}

TEST(CUASKernelsTest, updateEffectiveAquiferProperties) {
  /*
   * Note, we only test inside the boundary, because boundary conditions are applied elsewhere.
   * This test re-uses code from the getEffectiveAquiferProperties test.
   */

  ASSERT_EQ(mpiSize, MPI_SIZE);

  // input
  PETScGrid head(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid transmissivity(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid topg(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid bndMask(GRID_SIZE_X, GRID_SIZE_Y);
  // defaults
  PETScGrid TeffDflt(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid SeffDflt(GRID_SIZE_X, GRID_SIZE_Y);
  // output
  PETScGrid Teff(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid Seff(GRID_SIZE_X, GRID_SIZE_Y);

  // config confined-case
  constexpr PetscScalar specificStorage = 1e-3;          // unit: m^-1, selected for testing only; usually much lower
  constexpr PetscScalar specificYield = 0.3;             // unit: none, selected for testing only; usually 0.4
  constexpr PetscScalar Tinit = 10.0;                    // unit: m^2/s
  constexpr PetscScalar hinit = 100.0;                   // unit: m
  constexpr PetscScalar layerThickness = 2.0 * hinit;    // unit: m
  constexpr PetscScalar unconfinedSmooth = 1.5 * hinit;  // the head is now within the smoothing region

  head.setConst(hinit);            // unit: m, > layer thickness->confined
  transmissivity.setConst(Tinit);  // unit: m^2/s
  topg.setConst(0.0);              // unit: m
  bndMask.setConst(COMPUTE_FLAG);
  bndMask.setRealBoundary(DIRICHLET_FLAG);

  // get the default values for this setup
  CUAS::getEffectiveAquiferProperties(SeffDflt, TeffDflt, transmissivity, head, topg, bndMask, layerThickness,
                                      specificStorage, specificYield, unconfinedSmooth);

  // case1: disableUnconfined true, and doAnyChannel true.
  {
    Teff.setConst(999.9);
    Seff.setConst(666.6);
    CUAS::updateEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                           specificStorage, specificYield, unconfinedSmooth, true, true);

    // now test Teff
    auto &result1 = Teff.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), Tinit) << "at i=" << i << ", j=" << j;
      }
    }
    // Seff should be untouched
    auto &result2 = Seff.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), 666.6) << "at i=" << i << ", j=" << j;
      }
    }
  }

  // case2: disableUnconfined true, and doAnyChannel false.
  {
    Teff.setConst(999.9);
    Seff.setConst(666.6);
    CUAS::updateEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                           specificStorage, specificYield, unconfinedSmooth, true, false);

    // Teff should be untouched
    auto &result1 = Teff.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), 999.9) << "at i=" << i << ", j=" << j;
      }
    }
    // Seff should also be untouched
    auto &result2 = Seff.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), 666.6) << "at i=" << i << ", j=" << j;
      }
    }
  }

  // case3: disableUnconfined false, and doAnyChannel true.
  {
    Teff.setConst(999.9);
    Seff.setConst(666.6);
    CUAS::updateEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                           specificStorage, specificYield, unconfinedSmooth, false, true);

    // Teff should have changed
    auto &result1 = Teff.getReadHandle();
    auto &default1 = TeffDflt.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), default1(i, j)) << "at i=" << i << ", j=" << j;
      }
    }
    // Seff should have changed
    auto &result2 = Seff.getReadHandle();
    auto &default2 = SeffDflt.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), default2(i, j)) << "at i=" << i << ", j=" << j;
      }
    }
  }

  // case4: disableUnconfined false, and doAnyChannel false. This should give the same results as in case 3.
  {
    Teff.setConst(999.9);
    Seff.setConst(666.6);
    CUAS::updateEffectiveAquiferProperties(Seff, Teff, transmissivity, head, topg, bndMask, layerThickness,
                                           specificStorage, specificYield, unconfinedSmooth, false, false);

    // Teff should have changed
    auto &result1 = Teff.getReadHandle();
    auto &default1 = TeffDflt.getReadHandle();
    for (int i = 1; i < Teff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Teff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result1(i, j), default1(i, j)) << "at i=" << i << ", j=" << j;
      }
    }
    // Seff should have changed
    auto &result2 = Seff.getReadHandle();
    auto &default2 = SeffDflt.getReadHandle();
    for (int i = 1; i < Seff.getLocalNumOfRows() - 1; ++i) {
      for (int j = 1; j < Seff.getLocalNumOfCols() - 1; ++j) {
        EXPECT_DOUBLE_EQ(result2(i, j), default2(i, j)) << "at i=" << i << ", j=" << j;
      }
    }
  }

#ifdef TESTS_DUMP_NETCDF
  // Gets information about the currently running test.
  // Do NOT delete the returned object - it's managed by the UnitTest class.
  const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
  auto filename = std::string(test_info->test_suite_name())
                      .append(std::string("_-_"))
                      .append(std::string(test_info->name()))
                      .append(std::string(".nc"));

  CUAS::NetCDFFile file(filename, GRID_SIZE_X, GRID_SIZE_Y);
  file.defineGrid("bndMask", LIMITED);
  file.defineGrid("transmissivity", LIMITED);
  file.defineGrid("head", LIMITED);
  file.defineGrid("effective_transmissivity", LIMITED);
  file.defineGrid("effective_storativity", LIMITED);

  file.write("bndMask", bndMask, 0);
  file.write("transmissivity", transmissivity, 0);
  file.write("head", head, 0);
  file.write("effective_transmissivity", Teff, 0);
  file.write("effective_storativity", Seff, 0);

#endif
}

TEST(CUASKernelsTest, blockInflow) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  // constexpr auto SHMIP_NX = 66;  // 122 would be the full domain
  constexpr auto SHMIP_NX = 10;  // 122 would be the full domain
  constexpr auto SHMIP_NY = 23;

  std::array<std::array<PetscScalar, SHMIP_NX>, SHMIP_NY> transmissivity_blockInflow = {
      {{1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 1e-14, 1, 1, 1, 1, 1, 1, 1, 1},
       {100, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {100, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {100, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 1e-14, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1, 1, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14}}};

  std::array<std::array<PetscScalar, SHMIP_NX>, SHMIP_NY> transmissivity_blockInflowInv = {
      {{1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 100, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 100, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 100, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {100, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {100, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {100, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 100, 1, 1, 1, 1, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 100, 1, 1, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 100, 1, 1, 1, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 100, 1},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14},
       {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14}}};

  // input
  PETScGrid hydraulicHead(SHMIP_NX, SHMIP_NY);
  PETScGrid hydraulicTransmissivity(SHMIP_NX, SHMIP_NY);
  PETScGrid bndMask(SHMIP_NX, SHMIP_NY);
  PETScGrid bedElevation(SHMIP_NX, SHMIP_NY);
  PETScGrid iceThickness(SHMIP_NX, SHMIP_NY);
  PETScGrid result_blockInflow(SHMIP_NX, SHMIP_NY);
  PETScGrid result_blockInflowInv(SHMIP_NX, SHMIP_NY);

  constexpr PetscScalar Twater = 100.0;

  // SHMIP valley glacier parameters
  // domain length
  constexpr auto xend = 6.0e3;
  constexpr auto yend = 550.;
  // surface parameters
  constexpr auto beta = 0.25;

  constexpr auto s2 = 100. / xend;
  constexpr auto min_thick = 1.0;
  // bed parameters
  constexpr auto g1 = .5e-6;
  constexpr auto alpha = 3.;
  constexpr auto bench_para = 300. / xend;
  constexpr auto para = bench_para;  // could also be bed_para

  {
    auto mask = bndMask.getWriteHandle();
    auto head = hydraulicHead.getWriteHandle();
    auto trans = hydraulicTransmissivity.getWriteHandle();
    auto topg = bedElevation.getWriteHandle();
    auto thk = iceThickness.getWriteHandle();
    auto cornerX = bndMask.getCornerX();
    auto cornerY = bndMask.getCornerY();

    for (int row = 0; row < bndMask.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < bndMask.getLocalNumOfCols(); ++col) {
        auto x = ((cornerX + col) * 50.0) - 50.0;  // add another line of points to the east
        auto y = ((cornerY + row) * 50.0) - yend;

        // surface
        auto surf = 100. * pow(x + 200., beta) + s2 * x - pow(2.0e10, beta) + min_thick;
        auto s_xend = 100. * pow(xend + 200., beta) + s2 * xend - 100. * pow(200., beta) + min_thick;

        // helper functions
        auto f_func = para * x + pow(x, 2) * (s_xend - para * 6.0e3) / pow(6.0e3, 2);
        auto f_Bench = bench_para * x + pow(x, 2) * (s_xend - bench_para * 6.0e3) / pow(6.0e3, 2);
        auto g_func = 0.5e-6 * pow(std::abs(y), 3);
        auto h_func = (5 - 4.5 * x / 6.0e3) * (surf - f_func) / (surf - f_Bench);
        // bed elevation
        topg(row, col) = f_func + g_func * h_func;
        // ice thickness
        thk(row, col) = std::max(surf - topg(row, col), 0.0);

        if (thk(row, col) > 0.0) {
          mask(row, col) = (PetscScalar)COMPUTE_FLAG;
          trans(row, col) = 1.0;  // arbitrary value
        } else {
          mask(row, col) = (PetscScalar)NOFLOW_FLAG;
          trans(row, col) = NOFLOW_VALUE;
        }
        head(row, col) = 0.9 * thk(row, col) + topg(row, col);
        // head(row, col) = topg(row, col) + 42.0;  // debug
      }
    }
  }

  // loop though again and set all points next to active cuas mask to outflow
  // we can't use mask with ghosted read and write, so use thk instead.
  {
    auto mask = bndMask.getWriteHandle();
    auto trans = hydraulicTransmissivity.getWriteHandle();
    auto &thk = iceThickness.getReadHandle();

    auto cornerX = bndMask.getCornerX();
    auto cornerY = bndMask.getCornerY();
    auto cornerXGhost = bndMask.getCornerXGhost();
    auto cornerYGhost = bndMask.getCornerYGhost();

    for (int row = 0; row < bndMask.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < bndMask.getLocalNumOfCols(); ++col) {
        if (mask(row, col) == (PetscScalar)NOFLOW_FLAG) {
          auto j = row + (cornerY - cornerYGhost);  // ghostLocalRow
          auto i = col + (cornerX - cornerXGhost);  // ghostLocalCol

          // neighborhood of 4 connected pixels (1-connectivity)
          if (thk(j, i + 1, GHOSTED) > 0.0 || thk(j, i - 1, GHOSTED) > 0.0 || thk(j + 1, i, GHOSTED) > 0.0 ||
              thk(j - 1, i, GHOSTED) > 0.0) {
            mask(row, col) = (PetscScalar)DIRICHLET_LAKE_FLAG;
            trans(row, col) = Twater;
          }
        }
      }
    }
  }

  result_blockInflow.copy(hydraulicTransmissivity);
  CUAS::blockInflow(result_blockInflow, bndMask, hydraulicHead, Twater);

  result_blockInflowInv.copy(hydraulicTransmissivity);
  CUAS::blockInflowInv(result_blockInflowInv, bndMask, hydraulicHead, Twater);

#ifdef TESTS_DUMP_NETCDF
  // Gets information about the currently running test.
  // Do NOT delete the returned object - it's managed by the UnitTest class.
  const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
  auto filename = std::string(test_info->test_suite_name())
                      .append(std::string("_-_"))
                      .append(std::string(test_info->name()))
                      .append(std::string(".nc"));

  CUAS::NetCDFFile file(filename, SHMIP_NX, SHMIP_NY);
  file.defineGrid("bndMask", LIMITED);
  file.defineGrid("transmissivity", LIMITED);
  file.defineGrid("head", LIMITED);
  file.defineGrid("topg", LIMITED);
  file.defineGrid("thk", LIMITED);
  file.defineGrid("transmissivity_blockInflow", LIMITED);
  file.defineGrid("transmissivity_blockInflowInv", LIMITED);
  file.write("bndMask", bndMask, 0);
  file.write("transmissivity", hydraulicTransmissivity, 0);
  file.write("head", hydraulicHead, 0);
  file.write("topg", bedElevation, 0);
  file.write("thk", iceThickness, 0);
  file.write("transmissivity_blockInflow", result_blockInflow, 0);
  file.write("transmissivity_blockInflowInv", result_blockInflowInv, 0);
#endif

  // check pixel by pixel
  {
    auto &result = result_blockInflow.getReadHandle();
    int cornerX = result_blockInflow.getCornerX();
    int cornerY = result_blockInflow.getCornerY();
    for (int row = 0; row < result_blockInflow.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < result_blockInflow.getLocalNumOfCols(); ++col) {
        ASSERT_EQ(result(row, col), transmissivity_blockInflow[cornerY + row][cornerX + col]);
      }
    }
  }

  // check pixel by pixel
  {
    auto &result = result_blockInflowInv.getReadHandle();
    int cornerX = result_blockInflowInv.getCornerX();
    int cornerY = result_blockInflowInv.getCornerY();
    for (int row = 0; row < result_blockInflowInv.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < result_blockInflowInv.getLocalNumOfCols(); ++col) {
        ASSERT_EQ(result(row, col), transmissivity_blockInflowInv[cornerY + row][cornerX + col]);
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
