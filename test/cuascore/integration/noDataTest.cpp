#include "noDataTest.h"

#include "../unit/fillNoData.h"

#include "CUASArgs.h"
#include "CUASSolver.h"
#include "Forcing/ConstantForcing.h"

#include "gtest/gtest.h"

#include <memory>

#define MPI_SIZE 9

int mpiRank;
int mpiSize;

int m_argc;
char **m_argv;

TEST(noDataTest, compareModelToPython) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  CUAS::CUASArgs args;
  CUAS::parseArgs(m_argc, m_argv, args);

  auto pmodel = fillNoData();
  auto &model = *pmodel;
  model.init();

  ASSERT_EQ(model.dx, 1000.0);

  // fill topg with only 0.
  std::array<PetscScalar, NX> zero;
  zero.fill(0.0);
  std::array<std::array<PetscScalar, NX>, NY> topg;
  topg.fill(zero);

  // fill with ones
  std::array<PetscScalar, NX> ones;
  ones.fill(1.0);
  std::array<std::array<PetscScalar, NX>, NY> bmelt;
  bmelt.fill(ones);

  PETScGrid topgPy(NX, NY);
  PETScGrid thkPy(NX, NY);
  PETScGrid p_icePy(NX, NY);
  PETScGrid bndPy(NX, NY);
  PETScGrid QPy(NX, NY);
  {
    auto thkPy2d = thkPy.getWriteHandle();
    auto topgPy2d = topgPy.getWriteHandle();
    auto p_icePy2d = p_icePy.getWriteHandle();
    auto bndPy2d = bndPy.getWriteHandle();
    auto QPy2d = QPy.getWriteHandle();

    int cornerX = thkPy.getCornerX();
    int cornerY = thkPy.getCornerY();
    // fill up py grids to compare with mpi
    for (int i = 0; i < thkPy.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < thkPy.getLocalNumOfCols(); ++j) {
        thkPy2d(i, j) = thk[cornerY + i][cornerX + j];
        topgPy2d(i, j) = topg[cornerY + i][cornerX + j];
        p_icePy2d(i, j) = p_ice[cornerY + i][cornerX + j];
        bndPy2d(i, j) = bnd_mask[cornerY + i][cornerX + j];
        QPy2d(i, j) = bmelt[cornerY + i][cornerX + j];
      }
    }
  }

  auto &thkPy2d = thkPy.getReadHandle();
  auto &topgPy2d = topgPy.getReadHandle();
  auto &p_icePy2d = p_icePy.getReadHandle();
  auto &bndPy2d = bndPy.getReadHandle();
  auto &QPy2d = QPy.getReadHandle();

  auto &topgGlob = model.topg->getReadHandle();
  for (int i = 0; i < model.topg->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < model.topg->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(topgGlob(i, j), topgPy2d(i, j));
    }
  }

  auto &thkGlob = model.thk->getReadHandle();
  for (int i = 0; i < model.thk->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < model.thk->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(thkGlob(i, j), thkPy2d(i, j));
    }
  }

  auto &p_iceGlob = model.pIce->getReadHandle();
  for (int i = 0; i < model.pIce->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < model.pIce->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(p_iceGlob(i, j), p_icePy2d(i, j));
    }
  }

  auto &bndGlob = model.bndMask->getReadHandle();
  for (int i = 0; i < model.bndMask->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < model.bndMask->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(bndGlob(i, j), bndPy2d(i, j));
    }
  }

  // take zero here for constant forcing
  auto &Q = model.Q->getCurrentQ(0);
  auto &QGlob = Q.getReadHandle();
  for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(QGlob(i, j), QPy2d(i, j) / SPY * args.supplyMultiplier);
    }
  }
}

TEST(noDataTest, solverComparisonDirect) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto pmodel = fillNoData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  CUAS::parseArgs(m_argc, m_argv, args);

  CUAS::CUASSolver solver(&model, &args);
  solver.setup();

  // auto n = 7300;
  auto n = 73;  // as T is const., we would only need a few time steps to pass
  CUAS::timeSecs dt_secs = 43200;
  CUAS::timeSecs totaltime_secs = dt_secs * n;

  auto timeSteps = CUAS::getTimeStepArray(0, totaltime_secs, dt_secs);

  solver.solve(timeSteps);

  PETScGrid uPy(NX, NY);
  PETScGrid u_nPy(NX, NY);
  {
    auto uPyH = uPy.getWriteHandle();
    auto u_nPyH = u_nPy.getWriteHandle();

    int cornerX = uPy.getCornerX();
    int cornerY = uPy.getCornerY();
    for (int i = 0; i < uPy.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < uPy.getLocalNumOfCols(); ++j) {
        uPyH(i, j) = u[cornerY + i][cornerX + j];
        u_nPyH(i, j) = u_n[cornerY + i][cornerX + j];
      }
    }
  }

  auto &uPyRH = uPy.getReadHandle();
  auto &u_nPyRH = u_nPy.getReadHandle();

  auto &uGRH = solver.nextHead->getReadHandle();
  auto &u_nGRH = solver.currHead->getReadHandle();
  for (int i = 0; i < solver.nextHead->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < solver.nextHead->getLocalNumOfCols(); ++j) {
      ASSERT_NEAR(uGRH(i, j), uPyRH(i, j), 0.6) << "at i=" << i << ", j=" << j;
      ASSERT_NEAR(u_nGRH(i, j), u_nPyRH(i, j), 0.6) << "at i=" << i << ", j=" << j;  // obsolete
    }
  }
}

TEST(noDataTest, solverComparison) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto pmodel = fillNoData();
  auto &model = *pmodel;
  model.init();

  CUAS::CUASArgs args;
  CUAS::parseArgs(m_argc, m_argv, args);

  CUAS::CUASSolver solver(&model, &args);
  solver.setup();

  auto n = 7300;
  CUAS::timeSecs dt_secs = 43200;
  CUAS::timeSecs totaltime_secs = dt_secs * n;

  auto timeSteps = CUAS::getTimeStepArray(0, totaltime_secs, dt_secs);

  solver.solve(timeSteps);

  PETScGrid uPy(NX, NY);
  PETScGrid u_nPy(NX, NY);
  {
    auto uPyH = uPy.getWriteHandle();
    auto u_nPyH = u_nPy.getWriteHandle();

    int cornerX = uPy.getCornerX();
    int cornerY = uPy.getCornerY();
    for (int i = 0; i < uPy.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < uPy.getLocalNumOfCols(); ++j) {
        uPyH(i, j) = u[cornerY + i][cornerX + j];
        u_nPyH(i, j) = u_n[cornerY + i][cornerX + j];
      }
    }
  }

  auto &uPyRH = uPy.getReadHandle();
  auto &u_nPyRH = u_nPy.getReadHandle();

  auto &uGRH = solver.nextHead->getReadHandle();
  auto &u_nGRH = solver.currHead->getReadHandle();
  for (int i = 0; i < solver.nextHead->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < solver.nextHead->getLocalNumOfCols(); ++j) {
      ASSERT_NEAR(uGRH(i, j), uPyRH(i, j), 0.6);
      ASSERT_NEAR(u_nGRH(i, j), u_nPyRH(i, j), 0.6);
    }
  }

  // check c/solution
  int low, high;
  VecGetOwnershipRange(solver.sol->getRaw(), &low, &high);
  int size = high - low;
  for (int i = 0; i < size; ++i) {
    double v;
    int p[] = {low + i};
    VecGetValues(solver.sol->getRaw(), 1, p, &v);
    ASSERT_NEAR((float)v, (float)c[low + i], 0.6);
  }
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  m_argc = argc;
  m_argv = argv;
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
