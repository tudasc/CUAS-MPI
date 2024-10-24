/**
 * File: initialHeadTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "initialHeadTest.h"
#include "CUASArgs.h"
#include "CUASKernels.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

TEST(initialHeadTest, Nopc) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "Nopc";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_EQ(head(row, col), topg(row, col));
    }
  }
}

TEST(initialHeadTest, Nopc_belowSL) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "Nopc";
  args.layerThickness = 0.0;
  args.verbose = true;

  model.topg->setConst(-10);  // deviate from default here

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_EQ(head(row, col), 0.0);  // positive preserving head
    }
  }
}

TEST(initialHeadTest, low) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "low";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_NEAR(head(row, col), 0.1 * THK0 + topg(row, col), 1e-6);
    }
  }
}

TEST(initialHeadTest, mid) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "mid";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_NEAR(head(row, col), 0.5 * THK0 + topg(row, col), 1e-6);
    }
  }
}

TEST(initialHeadTest, high) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "high";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_NEAR(head(row, col), 0.9 * THK0 + topg(row, col), 1e-6);
    }
  }
}

TEST(initialHeadTest, Nzero) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "Nzero";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_NEAR(head(row, col), THK0 + topg(row, col), 1e-6);
    }
  }
}

TEST(initialHeadTest, topg) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "topg";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_EQ(head(row, col), topg(row, col));
    }
  }
}

TEST(initialHeadTest, scalar) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "100.0";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_EQ(head(row, col), 100.0);
    }
  }
}

TEST(initialHeadTest, scalarScientific) {
  PETScGrid result(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid workSpace(GRID_SIZE_X, GRID_SIZE_Y);
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "1e3";
  args.layerThickness = 0.0;
  args.verbose = true;

  CUAS::setInitialHeadFromArgs(result, *model.bndMask, *model.topg, *model.pIce, args, workSpace);

  auto &head = result.getReadHandle();
  for (int row = 0; row < result.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < result.getLocalNumOfCols(); ++col) {
      ASSERT_EQ(head(row, col), 1e3);
    }
  }
}

int main(int argc, char *argv[]) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  result = RUN_ALL_TESTS();
  PetscFinalize();

  return result;
}
