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

// TODO: verify
TEST(CUASKernelsTest, computeMelt) {}

// TODO: verify
TEST(CUASKernelsTest, binaryDialation) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto nf_mask = std::make_unique<PETScGrid>(GRID_SIZE_X, GRID_SIZE_Y);
  auto grad_mask = std::make_unique<PETScGrid>(GRID_SIZE_X, GRID_SIZE_Y);

  // setup
  {
    grad_mask->setZero();
    auto nf2d = nf_mask->getWriteHandle();
    for (int i = 0; i < nf_mask->getLocalNumOfRows(); ++i) {
      for (int j = 0; j < nf_mask->getLocalNumOfCols(); ++j) {
        if (i == j) {
          nf2d(i, j) = true;
        }
      }
    }
  }

  CUAS::binaryDialation(*grad_mask, *nf_mask);

  // check
  {
    auto &nf_arr = nf_mask->getReadHandle();
    auto &grad_arr = grad_mask->getReadHandle();
    for (int i = 0; i < nf_mask->getLocalNumOfRows(); ++i) {
      for (int j = 0; j < nf_mask->getLocalNumOfCols(); ++j) {
        if (i == j) {
          ASSERT_EQ(nf_arr(i, j), true);
          ASSERT_EQ(grad_arr(i, j), true);
          if (i < nf_mask->getLocalNumOfRows() - 1)
            ASSERT_EQ(grad_arr(i + 1, j), true);
          if (i > 0)
            ASSERT_EQ(grad_arr(i - 1, j), true);
          if (j < nf_mask->getLocalNumOfCols() - 1)
            ASSERT_EQ(grad_arr(i, j + 1), true);
          if (j > 0)
            ASSERT_EQ(grad_arr(i, j - 1), true);
        }
      }
    }
  }
}

// TODO: What does this test check? Do we need it?
TEST(CUASKernelsTest, ghostCellsTest) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto nf_mask = std::make_unique<PETScGrid>(GRID_SIZE_X, GRID_SIZE_Y);
  auto grad_mask = std::make_unique<PETScGrid>(GRID_SIZE_X, GRID_SIZE_Y);

  auto nf2d = nf_mask->getWriteHandle();

  auto rows = nf_mask->getLocalNumOfRows();
  auto cols = nf_mask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (i == j) {
        nf2d(i, j) = true;
      }
    }
  }

  nf2d.setValues();
  grad_mask->setZero();

  CUAS::binaryDialation(*grad_mask, *nf_mask);

  auto &grad_arr = grad_mask->getReadHandle();

  auto rowsGhost = grad_mask->getLocalGhostNumOfRows();
  auto colsGhost = grad_mask->getLocalGhostNumOfCols();

  auto xGhost = grad_mask->getCornerXGhost();
  auto yGhost = grad_mask->getCornerYGhost();

  auto total_rows = grad_mask->getTotalNumOfRows();
  auto total_cols = grad_mask->getTotalNumOfCols();

  for (int i = 0; i < rowsGhost; ++i) {
    for (int j = 0; j < colsGhost; ++j) {
      if (yGhost == -1 && i == 0 || xGhost == -1 && j == 0 ||
          yGhost + rowsGhost - 1 == total_rows && i == rowsGhost - 1 ||
          xGhost + colsGhost - 1 == total_cols && j == colsGhost - 1) {
        ASSERT_EQ(grad_arr(i, j, GHOSTED), false);
      }
    }
  }
}

// TODO: verify
TEST(CUASKernelsTest, enableUnconfined) {}

// TODO: verify
TEST(CUASKernelsTest, calculateTeffPowTexp) {}

// TODO: verify
TEST(CUASKernelsTest, calculateSeValues) {}

// TODO: verify
TEST(CUASKernelsTest, doChannels) {}

// TODO: verify
TEST(CUASKernelsTest, noChannels) {}

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

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
