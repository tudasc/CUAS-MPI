#include "CUASKernels.h"

#include "PETScGrid.h"

//#define TESTS_DUMP_NETCDF
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

  // Todo: update test to ensure 0 <= level <= waterLayerThickness
  PetscScalar level = 33.5;

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

        // check for sideeffects
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

  // Todo: update test to reflect changes in cuas-python
  PetscScalar seaLevel = 0.0;  // was 33.5

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

  // run pressureToHead
  CUAS::pressureToHead(head, pressure, bedElevation);

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
        ASSERT_EQ(head2d(j, i), pressure2d(j, i) / (RHO_WATER * GRAVITY) + bedElevation2d(j, i));

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

TEST(CUASKernelsTest, computeCavityOpening) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid basalVelocity(GRID_SIZE_Y, GRID_SIZE_X);
  PETScGrid result(GRID_SIZE_Y, GRID_SIZE_X);

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

  // run pressureToHead
  CUAS::computeCavityOpening(result, beta, K, basalVelocity);

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
  PETScGrid gradh2(GRID_SIZE_X, GRID_SIZE_Y);

  PetscScalar r = 1.0;
  PetscScalar bt = 0.1;
  PetscScalar K = 10;
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

  CUAS::computeMeltOpening(result, r, K, T, gradh2);

  auto &res2d = result.getReadHandle();

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

TEST(CUASKernelsTest, getFluxMagnitude) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid head(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid flux(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid flux_analytical(GRID_SIZE_X, GRID_SIZE_Y);

  PETScGrid gradHeadSquared(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradMask(GRID_SIZE_X, GRID_SIZE_Y);

  constexpr PetscScalar Tconst = 0.2;
  constexpr PetscScalar dx = 1.0;

  // initialize fields
  T.setConst(Tconst);
  flux.setConst(-9999.0);

  gradMask.setConst(0.0);
  gradMask.setRealBoundary(1.0);

  // init head
  {
    auto cornerX = head.getCornerX();
    auto cornerY = head.getCornerY();
    auto h = head.getWriteHandle();
    auto f = flux_analytical.getWriteHandle();
    auto gradh2 = gradHeadSquared.getWriteHandle();
    auto &gmask = gradMask.getReadHandle();
    for (int j = 0; j < head.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < head.getLocalNumOfCols(); ++i) {
        if (gmask(j, i) == 1.0) {  // gradMask has inverted meaning
          f(j, i) = 0.0;
        } else {
          const PetscScalar x = (cornerX + i) * dx;
          const PetscScalar y = (cornerY + j) * dx;
          h(j, i) = x * x + y * y;
          // analytical solution:  flux = T^exp * |grad h|
          gradh2(j, i) = 4.0 * (x * x + y * y);
          f(j, i) = Tconst * PetscSqrtScalar(gradh2(j, i));
        }
      }
    }
  }

  CUAS::getFluxMagnitude(flux, gradHeadSquared, T);

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
  file.defineGrid("gradMask", LIMITED);
  file.defineGrid("flux", LIMITED);
  file.defineGrid("flux_analytical", LIMITED);
  file.write("T", T, 0);
  file.write("head", head, 0);
  file.write("gradMask", gradMask, 0);
  file.write("flux", flux, 0);
  file.write("flux_analytical", flux_analytical, 0);
#endif

  // compare results

  auto &f = flux.getReadHandle();
  auto &f_analytical = flux.getReadHandle();
  auto &gmask = gradMask.getReadHandle();
  for (int j = 0; j < flux.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < flux.getLocalNumOfCols(); ++i) {
      if (gmask(j, i) == 1.0) {  // gradMask has inverted meaning
        EXPECT_DOUBLE_EQ(f(j, i), 0.0);
      } else {
        EXPECT_DOUBLE_EQ(f(j, i), f_analytical(j, i));
      }
    }
  }
}

TEST(CUASKernelsTest, doChannels) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid hydraulicHead(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid melt(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid creep(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid cavity(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T_n(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid topg(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid basalVelocityIce(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid rateFactorIce(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid bndMask(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid effectivePressure(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid icePressure(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid gradientHeadSquared(GRID_SIZE_X, GRID_SIZE_Y);

  // init args
  constexpr PetscScalar roughnessFactor = 1.0;
  constexpr PetscScalar noSmoothMelt = false;
  constexpr PetscScalar cavityBeta = 0.0005;
  constexpr PetscScalar layerThickness = 0.1;
  constexpr PetscScalar dx = 1000.0;
  constexpr PetscScalar dtSecs = 43200;
  constexpr PetscScalar Tmin = 1e-07;
  constexpr PetscScalar Tmax = 20.0;
  constexpr PetscScalar K = 10.0;

  basalVelocityIce.setConst(1e-06);
  rateFactorIce.setConst(5e-25);

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

  // init T values
  T.setConst(0.2);
  T_n.setConst(0.2);

  // init noFlowMask
  auto mask = bndMask.getWriteHandleGhost();

  PetscScalar *nFMask[GRID_SIZE_Y + 2];
  PetscScalar nFBeginningEnd[GRID_SIZE_X + 2];
  for (int i = 0; i < GRID_SIZE_X + 2; ++i) {
    nFBeginningEnd[i] = (PetscScalar)NOFLOW_FLAG;  // noFlowMask = true
  }
  PetscScalar nFMiddle[GRID_SIZE_X + 2];
  nFMiddle[0] = (PetscScalar)NOFLOW_FLAG;                   // noFlowMask = true
  nFMiddle[GRID_SIZE_X + 1] = (PetscScalar)DIRICHLET_FLAG;  // noFlowMask = false
  for (int i = 1; i < GRID_SIZE_X + 1; ++i) {
    nFMiddle[i] = (PetscScalar)COMPUTE_FLAG;  // noFlowMask = false
  }

  nFMask[0] = nFBeginningEnd;
  for (int i = 1; i < GRID_SIZE_Y + 1; ++i) {
    nFMask[i] = nFMiddle;
  }
  nFMask[GRID_SIZE_Y + 1] = nFBeginningEnd;

  for (int i = 0; i < bndMask.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < bndMask.getLocalGhostNumOfCols(); ++j) {
      mask(i, j) = nFMask[bndMask.getCornerY() + i][bndMask.getCornerX() + j];
    }
  }
  mask.setValues();

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
  auto pIce2d = icePressure.getWriteHandleGhost();
  PetscScalar *pIceArr[GRID_SIZE_Y + 2];
  PetscScalar pIceArrValues[GRID_SIZE_X + 2] = {17854200, 16961490, 16068780, 15176070, 14283360, 13390650, 12497940,
                                                11605230, 10712520, 9819810,  8927100,  8034390,  7141680,  6248970,
                                                5356260,  4463550,  3570840,  2678130,  1785420,  892710};
  for (int i = 0; i < GRID_SIZE_Y + 2; ++i) {
    pIceArr[i] = pIceArrValues;
  }

  for (int i = 0; i < icePressure.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < icePressure.getLocalGhostNumOfCols(); ++j) {
      pIce2d(i, j) = pIceArr[icePressure.getCornerY() + i][icePressure.getCornerX() + j];
    }
  }
  pIce2d.setValues();

  topg.setZero();  // topg is zero in this setup
  creep.setZero();
  melt.setZero();
  cavity.setZero();

  CUAS::headToEffectivePressure(effectivePressure, hydraulicHead, topg, icePressure, layerThickness);
  CUAS::computeCreepOpening(creep, rateFactorIce, effectivePressure, T_n);  // was creep closure in previous versions

  CUAS::getGradHeadSQR(gradientHeadSquared, hydraulicHead, dx, gradMask);
  CUAS::computeMeltOpening(melt, roughnessFactor, K, T_n, gradientHeadSquared);
  if (!noSmoothMelt) {
    PETScGrid tmp(melt.getTotalNumOfCols(), melt.getTotalNumOfRows());
    CUAS::convolveStar11411(melt, tmp);
    melt.copy(tmp);
  }

  CUAS::computeCavityOpening(cavity, cavityBeta, K, basalVelocityIce);

  // update
  CUAS::doChannels(T, T_n, creep, melt, cavity, bndMask, Tmin, Tmax, dtSecs);

  // Gets information about the currently running test.
  // Do NOT delete the returned object - it's managed by the UnitTest class.
  const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
  auto filename = std::string(test_info->test_suite_name())
                      .append(std::string("_-_"))
                      .append(std::string(test_info->name()))
                      .append(std::string(".nc"));
#ifdef TESTS_DUMP_NETCDF
  CUAS::NetCDFFile file(filename, GRID_SIZE_X, GRID_SIZE_Y);
  file.defineGrid("creep", LIMITED);
  file.defineGrid("melt", UNLIMITED);
  file.defineGrid("cavity", UNLIMITED);
  file.defineGrid("bndMask", LIMITED);
  file.defineGrid("T", LIMITED);
  file.defineGrid("T_n", LIMITED);
  file.write("creep", creep, 0);
  file.write("melt", melt, 0);
  file.write("cavity", cavity, 0);
  file.write("bndMask", bndMask, 0);
  file.write("T", T, 0);
  file.write("T_n", T_n, 0);
#endif

  //
  // melt opening
  //
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

  // tkleiner (22.02.2022): fails as we compute melt opening using the gradMask
  auto &melt2d = melt.getReadHandle();
  for (int i = 0; i < melt.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < melt.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(melt2d(i, j), meltArr[melt.getCornerY() + i][melt.getCornerX() + j])
          << "at i=" << i << ", j=" << j;
    }
  }

  //
  // creep opening
  //
  PetscScalar *creepArr[GRID_SIZE_Y];
  PetscScalar creepValues[GRID_SIZE_X];
  for (int i = 0; i < GRID_SIZE_X; ++i) {
    creepValues[i] = -0.0000000000000000069931566;  // we now have creep opening instead of creep closure
  }

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    creepArr[i] = creepValues;
  }

  // tkleiner (22.02.2022): ok, as we compute creep opening everywhere without using the bndMask
  auto &creep2d = creep.getReadHandle();
  for (int i = 0; i < creep.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < creep.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(creep2d(i, j), creepArr[creep.getCornerY() + i][creep.getCornerX() + j])
          << "at i=" << i << ", j=" << j;
    }
  }

  PetscScalar *cavityArr[GRID_SIZE_Y];
  PetscScalar cavityMiddleRow[GRID_SIZE_X];
  for (int i = 0; i < GRID_SIZE_X; ++i) {
    cavityMiddleRow[i] = 0.000000005;  //
  }

  for (int i = 0; i < GRID_SIZE_Y; ++i) {
    cavityArr[i] = cavityMiddleRow;
  }

  // tkleiner (22.02.2022): ok, as we compute cavity opening everywhere without using the bndMask
  auto &cavity2d = cavity.getReadHandle();
  for (int i = 0; i < cavity.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < cavity.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(cavity2d(i, j), cavityArr[cavity.getCornerY() + i][cavity.getCornerX() + j])
          << "at i=" << i << ", j=" << j;
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

  auto &updatedTGlobal = T.getReadHandle();
  for (int i = 0; i < T.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < T.getLocalNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(updatedTGlobal(i, j), updatedTArr[T.getCornerY() + i][T.getCornerX() + j])
          << "at i=" << i << ", j=" << j;
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

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
