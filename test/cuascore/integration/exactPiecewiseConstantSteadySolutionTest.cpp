/**
 * File: exactPiecewiseConstantSteadySolutionTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "exactPiecewiseConstantSteadySolutionTest.h"
#include "CUASArgs.h"
#include "CUASSolver.h"
#include "PETScGrid.h"
#include "timeparse.h"  // CUAS::Time, ...
#include <algorithm>    // std::min, std::max
#include <cstdio>       // sprintf, snprintf

#ifdef TESTS_DUMP_NETCDF
#include "NetCDFFile.h"
#endif

#include "gtest/gtest.h"
#include <memory>

#define BUFLEN 80

int mpiRank;
int mpiSize;

int m_argc;
char **m_argv;

TEST(exactPiecewiseConstantSteadySolutionTest, transmissivityFuncXTest) {
#ifndef NDEBUG
  const std::vector<PetscScalar> xVec = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  for (auto &x : xVec) {
    CUAS_INFO_RANK0("x = {}, T(x) = {}", x, transmissivityFuncX(x))
  }
#endif

  // x outside to the left
  {
    auto xTest = layerBoundaries.front() - 0.666;
    auto T = transmissivityFuncX(xTest);
    ASSERT_EQ(NOFLOW_VALUE, T) << "at x=" << xTest;
  }

  // x outside to the right
  {
    auto xTest = layerBoundaries.back() + 0.666;
    auto T = transmissivityFuncX(xTest);
    ASSERT_EQ(NOFLOW_VALUE, T) << "at x=" << xTest;
  }
}

// const std::vector<PetscScalar> layerBoundaries = {0.0, 0.2, 0.5, 1};     // len M+1
// const std::vector<PetscScalar> layerTransmissivities = {0.2, 0.4, 4.0};  // len M
//   x = [0.0  0.1  0.2  0.3  0.4  0.5  0.6   0.7  0.8   0.9   1.   ]
//   h = [0.5  1.7  2.9  3.5  4.1  4.7  4.76  4.82 4.88  4.94  5.   ]
// num = [0.   0.5  1.   1.25 1.5  1.75 1.775 1.8  1.825 1.85  1.875]
// den = 1.875
TEST(exactPiecewiseConstantSteadySolutionTest, exactSolutionFuncXTest) {
#ifndef NDEBUG
  const std::vector<PetscScalar> xVec = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  for (auto &x : xVec) {
    CUAS_INFO_RANK0("x = {}, h_exact(x) = {}", x, exactSolutionFuncX(x))
  }
#endif

  // some known solution
  {
    // left BC
    auto xTest = layerBoundaries.front();
    ASSERT_EQ(head_x0, exactSolutionFuncX(xTest)) << "at x=" << xTest;
  }
  {
    // right BC
    auto xTest = layerBoundaries.back();
    ASSERT_EQ(head_xL, exactSolutionFuncX(xTest)) << "at x=" << xTest;
  }
}

TEST(exactPiecewiseConstantSteadySolutionTest, exactSolutionTest) {
  // command line arg
  // 0: testname,  1: nx, 2: nt, 3: dt, 4: filename
  int nx, ny, nt;
  CUAS::timeSecs dt;
  std::string fileName;

  ASSERT_LE(m_argc, 5);  // no more than 4 arguments allowed
  parseTestArgs(m_argc, m_argv, nx, nt, dt, fileName);

  ny = 11;  // y in [-5 * dx, 5 * dx]

  ASSERT_EQ(nx % 2, 1);
  ASSERT_EQ(ny % 2, 1);

  auto pmodel = fillModelData(nx, ny);
  auto &model = *pmodel;
  model.init();  // set dx, dy, and do some checks

  ASSERT_EQ(layerBoundaries.front(), model.xAxis.front());
  ASSERT_EQ(layerBoundaries.back(), model.xAxis.back());
  ASSERT_EQ(-5.0 * model.dx, model.yAxis.front());
  ASSERT_EQ(5.0 * model.dx, model.yAxis.back());

  // set-up command-line arguments for CUAS
  int argc = 10;
  char arg0[] = "test";
  char arg1[] = "--verbose";
  char arg2[] = "--directSolver";       //
  char arg3[] = "--outputSize=xlarge";  //
  char arg4[] = "--saveEvery=1";        //
  char arg5[BUFLEN];
  snprintf(arg5, BUFLEN, "--totaltime=%ld seconds", dt * nt);
  char arg6[BUFLEN];
  snprintf(arg6, BUFLEN, "--dt=%ld seconds", dt);
  char arg7[BUFLEN];
  snprintf(arg7, BUFLEN, "--specificStorage=%g", 1.0);
  char arg8[BUFLEN];
  snprintf(arg8, BUFLEN, "--layerThickness=%f", 1.0);
  char arg9[BUFLEN];
  snprintf(arg9, BUFLEN, "--disableUnconfined");  // we have 0.5m head at the left boundary -> this would be unconfined

  char *argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9};
  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

  // setup time from args
  CUAS::Time time;
  time.timeSteps = CUAS::getTimeStepArray(0, CUAS::parseTime(args.totaltime), CUAS::parseTime(args.dt));

  // setup solution handler for file io
  std::unique_ptr<CUAS::SolutionHandler> solutionHandler;
  if (!fileName.empty()) {
    solutionHandler = std::make_unique<CUAS::SolutionHandler>(fileName, model.Ncols, model.Nrows, args.outputSize);
    // fixme: The solution handler writes to fileName no matter what is given as args.output,
    //        but CUAS command-line output reports args.output
    args.output = fileName;  // for "cosmetic" reasons only
  }

  // CUAS::CUASSolver solver(&model, &args, solutionHandler.get());
  TestSolver solver(&model, &args, solutionHandler.get());
  solver.setup(pmodel->xAxis, pmodel->yAxis);
  solver.solve(time.timeSteps);

  // compute analytical solution
  PETScGrid exactHead(nx, ny);
  solver.getExactSolution(pmodel->xAxis, pmodel->yAxis, exactHead);

// Dump to netcdf for further inspection
#ifdef TESTS_DUMP_NETCDF
  {
    PETScGrid diff(nx, ny);
    {
      auto nRows = exactHead.getLocalNumOfRows();  // y-dir
      auto nCols = exactHead.getLocalNumOfCols();  // x-dir
      auto &hexact = exactHead.getReadHandle();
      auto &hcuas = solver.currHead->getReadHandle();
      auto d = diff.getWriteHandle();
      for (int j = 0; j < nRows; ++j) {    // y-dir
        for (int i = 0; i < nCols; ++i) {  // x-dir
          d(j, i) = hexact(j, i) - hcuas(j, i);
        }
      }
    }

    // Do NOT delete the returned object - it's managed by the UnitTest class.
    const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
    auto filename = std::string(test_info->test_suite_name())
                        .append(std::string("_-_"))
                        .append(std::string(test_info->name()))
                        .append(std::string(".nc"));
    CUAS::NetCDFFile file(filename, model.Ncols, model.Nrows);
    file.defineGrid("h_cuas", LIMITED);
    file.defineGrid("h_exact", LIMITED);
    file.defineGrid("h_diff", LIMITED);
    file.defineGrid("bndMask", LIMITED);
    file.write("h_cuas", *solver.currHead, 0);
    file.write("h_exact", exactHead, 0);
    file.write("bndMask", *model.bndMask, 0);
    file.write("h_diff", diff, 0);
  }
#endif

  // Report error norms for temporal, spatial or combined order of accuracy analysis.
  // Quite often in the literature authors report the error as L1 or L2 norm of the
  // difference, but this is not the correct error norm. Ferziger&Peric 2002 are very clear about this.
  // See e.g. "Richardson extrapolation"
  // Mean Absolute Error (MAE), Root Mean Squared Error (RMSE), MaxE
  const PetscScalar N = nx * ny;
  // fixme: We get large errors here, because at the northern and southern margin we have no-flow conditions
  //        and the head stays at the initial conditions!
  auto [L1, L2, Linf] = exactHead.getErrorNorms(*solver.currHead);
  CUAS_INFO_RANK0("Result: res={}, dt={}, MAE={}, RMSE={}, MaxE={}\n", pmodel->dx, dt, L1 / N, L2 / N, Linf)

  // compare on the grid
  {
    auto nRows = exactHead.getLocalNumOfRows();
    auto nCols = exactHead.getLocalNumOfCols();
    auto &hexact = exactHead.getReadHandle();
    auto &hcuas = solver.currHead->getReadHandle();
    auto &mask = model.bndMask->getReadHandle();
    for (int row = 0; row < nRows; ++row) {    // y-dir
      for (int col = 0; col < nCols; ++col) {  // x-dir
        // at the now-flow boundaries the head stays at the initial head
        if (mask(row, col) == COMPUTE_FLAG || mask(row, col) == DIRICHLET_FLAG) {
          ASSERT_NEAR(hexact(row, col), hcuas(row, col), 0.05) << "at row=" << row << ", col=" << col;
        }
      }
    }
  }

  // The steady-state flux must be constant, although we have a layered setup
  {
    auto nRows = solver.fluxMagnitude->getLocalNumOfRows();
    auto nCols = solver.fluxMagnitude->getLocalNumOfCols();
    auto &flux = solver.fluxMagnitude->getReadHandle();
    auto &mask = model.bndMask->getReadHandle();
    for (int row = 0; row < nRows; ++row) {    // y-dir
      for (int col = 0; col < nCols; ++col) {  // x-dir
        if (mask(row, col) == COMPUTE_FLAG) {
          // this works only for the given layerBoundaries and layerTransmissivities
          ASSERT_NEAR(flux(row, col), 2.4307900, 1e-6) << "at row=" << row << ", col=" << col;
        }
      }
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
