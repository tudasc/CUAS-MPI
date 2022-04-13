#include "initialHeadTest.h"

#include "CUASArgs.h"
#include "CUASSolver.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

TEST(initialHeadTest, zero) {
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "zero";
  args.verbose = true;

  CUAS::CUASSolver solver(&model, &args);
  solver.setup();

  auto &head = solver.currHead->getReadHandle();
  for (int i = 0; i < solver.currHead->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < solver.currHead->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(head(i, j), 0.0);
    }
  }
}

TEST(initialHeadTest, Nzero) {
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "Nzero";
  args.verbose = true;

  CUAS::CUASSolver solver(&model, &args);
  solver.setup();

  auto &head = solver.currHead->getReadHandle();
  auto &pice = model.pIce->getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int i = 0; i < solver.currHead->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < solver.currHead->getLocalNumOfCols(); ++j) {
      ASSERT_NEAR(head(i, j), pice(i, j) / (RHO_WATER * GRAVITY) + topg(i, j), 1e-6);
    }
  }
}

TEST(initialHeadTest, topg) {
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "topg";
  args.verbose = true;

  CUAS::CUASSolver solver(&model, &args);
  solver.setup();

  auto &head = solver.currHead->getReadHandle();
  auto &topg = model.topg->getReadHandle();
  for (int i = 0; i < solver.currHead->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < solver.currHead->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(head(i, j), topg(i, j));
    }
  }
}

TEST(initialHeadTest, scalar) {
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "100.0";
  args.verbose = true;

  CUAS::CUASSolver solver(&model, &args);
  solver.setup();

  auto &head = solver.currHead->getReadHandle();
  for (int i = 0; i < solver.currHead->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < solver.currHead->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(head(i, j), 100.0);
    }
  }
}

TEST(initialHeadTest, scalarScientific) {
  auto pmodel = fillData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  args.initialHead = "1e3";
  args.verbose = true;

  CUAS::CUASSolver solver(&model, &args);
  solver.setup();

  auto &head = solver.currHead->getReadHandle();
  for (int i = 0; i < solver.currHead->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < solver.currHead->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(head(i, j), 1e3);
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
