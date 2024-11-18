/**
 * File: noDataTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "noDataTest.h"

#include "../unit/fillNoData.h"

#include "CUASArgs.h"
#include "CUASSolver.h"
#include "Forcing/SteadyForcing.h"

// #define TESTS_DUMP_NETCDF
#ifdef TESTS_DUMP_NETCDF
#include "NetCDFFile.h"
#endif

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

#ifdef TESTS_DUMP_NETCDF
  {
    // Gets information about the currently running test.
    // Do NOT delete the returned object - it's managed by the UnitTest class.
    const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
    auto filename = std::string(test_info->test_suite_name())
                        .append(std::string("_-_"))
                        .append(std::string(test_info->name()))
                        .append(std::string(".nc"));
    CUAS::NetCDFFile file(filename, model.Ncols, model.Nrows);
    file.defineGrid("thk", LIMITED);
    file.defineGrid("topg", LIMITED);
    file.defineGrid("bnd_mask", LIMITED);
    file.defineGrid("bmelt", LIMITED);
    file.write("thk", *model.thk, 0);
    file.write("topg", *model.topg, 0);
    file.write("bnd_mask", *model.bndMask, 0);
    file.write("bmelt", model.getCurrentWaterSource(0), 0);  // m/s
  }
#endif

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

  auto &thkPy2d = thkPy.getReadHandle();
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
  auto &Q = model.getCurrentWaterSource(0);
  auto &QGlob = Q.getReadHandle();
  for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(QGlob(i, j), QPy2d(i, j) / SPY * args.supplyMultiplier);
    }
  }
}

TEST(noDataTest, solverComparison) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto pmodel = fillNoData();
  auto &model = *pmodel;
  model.init();

  int argc = 6;
  char arg0[] = "test";
  char arg1[] = "--verbose";
  char arg2[] = "--Tinit=0.2";                       // will be ignored if no channel evolution
  char arg3[] = "--specificStorage=0.000982977696";  // old cuas defaults: Ssmulti=1.0, Ss=0.000982977696
  char arg4[] = "--conductivity=2.0";                // ensure we have T = K*b = Tinit = 0.2"
  char arg5[] = "--layerThickness=0.1";              // old cuas default -> important to get S = Ss*b right
  char *argv[] = {arg0, arg1, arg2, arg3, arg4, arg5};
  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

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

#ifdef TESTS_DUMP_NETCDF
  {
    // Gets information about the currently running test.
    // Do NOT delete the returned object - it's managed by the UnitTest class.
    const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
    auto filename = std::string(test_info->test_suite_name())
                        .append(std::string("_-_"))
                        .append(std::string(test_info->name()))
                        .append(std::string(".nc"));
    CUAS::NetCDFFile file(filename, model.Ncols, model.Nrows);
    file.defineGrid("u", LIMITED);
    file.defineGrid("u_python", LIMITED);
    file.defineGrid("bndMask", LIMITED);
    file.defineGrid("T", LIMITED);

    file.write("u", *solver.nextHead, 0);
    file.write("u_python", uPy, 0);
    file.write("bndMask", *model.bndMask, 0);
    file.write("T", *solver.nextTransmissivity, 0);
  }
#endif

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

TEST(noDataTest, solverComparisonDirect) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto pmodel = fillNoData();
  auto &model = *pmodel;
  model.init();

  int argc = 7;
  char arg0[] = "test";
  char arg1[] = "--verbose";
  char arg2[] = "--Tinit=0.2";                       // will be ignored if no channel evolution
  char arg3[] = "--specificStorage=0.000982977696";  // old cuas defaults: Ssmulti=1.0, Ss=0.000982977696
  char arg4[] = "--conductivity=2.0";                // ensure we have T = K*b = Tinit = 0.2"
  char arg5[] = "--layerThickness=0.1";              // old cuas default -> important to get S = Ss*b right
  char arg6[] = "--directSolver";  // we need to use MUMPS to pass this test in 73 steps instead of 7300
  char *argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6};
  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

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

#ifdef TESTS_DUMP_NETCDF
  {
    // Gets information about the currently running test.
    // Do NOT delete the returned object - it's managed by the UnitTest class.
    const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
    auto filename = std::string(test_info->test_suite_name())
                        .append(std::string("_-_"))
                        .append(std::string(test_info->name()))
                        .append(std::string(".nc"));
    CUAS::NetCDFFile file(filename, model.Ncols, model.Nrows);
    file.defineGrid("u", LIMITED);
    file.defineGrid("u_python", LIMITED);
    file.defineGrid("bndMask", LIMITED);
    file.defineGrid("T", LIMITED);
    file.write("u", *solver.nextHead, 0);
    file.write("u_python", uPy, 0);
    file.write("bndMask", *model.bndMask, 0);
    file.write("T", *solver.nextTransmissivity, 0);
  }
#endif

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
