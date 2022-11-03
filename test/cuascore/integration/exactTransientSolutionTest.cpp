#include "exactTransientSolutionTest.h"
#include "CUASArgs.h"
#include "CUASSolver.h"
#include "timeparse.h"  // CUAS::Time, ...
#include <cstdio>       // sprintf, snprintf

#ifdef TESTS_DUMP_NETCDF
#include "NetCDFFile.h"
#endif

#include "gtest/gtest.h"

#include <memory>

#define MPI_SIZE 4
#define BUFLEN 80

int mpiRank;
int mpiSize;

int m_argc;
char **m_argv;

TEST(exactTransientSolutionTest, transientForcingTest) {
  // ASSERT_EQ(mpiSize, MPI_SIZE);

  int nx = 21;
  int ny = 11;
  PetscScalar res = 50000.0;

  // nx and ny must be odd for the test
  ASSERT_EQ(nx % 2, 1);
  ASSERT_EQ(ny % 2, 1);

  // setup time steps
  CUAS::Time time;
  CUAS::timeSecs dt = 86400;
  CUAS::timeSecs endTime = dt * 31;
  time.timeSteps = CUAS::getTimeStepArray(0, endTime, dt);

  // setup forcing
  auto forcing = std::make_unique<TransientForcing>(nx, ny, res);

  // sanity test the forcing: max(Q(x,y,t)) == Q_max(t), convert m/s to m/a
  const PetscScalar Lx = (nx - 1) * res;
  const PetscScalar Ly = (ny - 1) * res;
  for (auto &t_secs : time.timeSteps) {
    auto Q_max = max_Q_exact((PetscScalar)t_secs, Lx, Ly);
    auto maxQ = forcing->getCurrentQ(t_secs).getMax();
    ASSERT_DOUBLE_EQ(Q_max * SPY, maxQ * SPY) << "at t = " << t_secs << " seconds";
  }
}

TEST(exactTransientSolutionTest, exactSolutionTest) {
  // ASSERT_EQ(mpiSize, MPI_SIZE);

  // command line args
  // 0: testname,  1: res, 2: nx, 3: ny, 4: nt, 5: dt, 6: filename
  PetscScalar res;
  int nx, ny, nt;
  CUAS::timeSecs dt;
  std::string fileName;

  ASSERT_LE(m_argc, 7);  // no more than 6 arguments allowed
  parseTestArgs(m_argc, m_argv, res, nx, ny, nt, dt, fileName);
  CUAS_INFO_RANK0("{}:\n\tres={}, nx={}, ny={}, nt={}, dt={}, filename=<{}>", __PRETTY_FUNCTION__, res, nx, ny, nt, dt,
                  fileName);

  // for testing hmax and Qmax at Lx/2, Ly/2 we need a grid point at that location,
  // thus nx and ny must be odd
  ASSERT_EQ(nx % 2, 1);
  ASSERT_EQ(ny % 2, 1);

  auto pmodel = fillModelData(nx, ny, res);
  auto &model = *pmodel;
  model.init();  // set dx, dy, and does some checks

  // test x- and yAxis after model init
  const PetscScalar Lx = (nx - 1) * res;
  const PetscScalar Ly = (ny - 1) * res;
  ASSERT_EQ(Lx, model.xAxis.back() - model.xAxis.front());
  ASSERT_EQ(Ly, model.yAxis.back() - model.yAxis.front());

  // confined: the layer thickness must be below the initial head
  constexpr PetscScalar layerThickness = 0.5 * EXACT_TRANSIENT_SOLUTION_H0;

  // set-up command-line arguments for CUAS
  int argc = 12;
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
  snprintf(arg7, BUFLEN, "--initialHead=%.1f", EXACT_TRANSIENT_SOLUTION_H0);
  char arg8[BUFLEN];
  snprintf(arg8, BUFLEN, "--Tinit=%g", EXACT_TRANSIENT_SOLUTION_TRANSMISSIVITY);
  char arg9[BUFLEN];
  snprintf(arg9, BUFLEN, "--specificStorage=%g", EXACT_TRANSIENT_SOLUTION_STORATIVITY / layerThickness);
  char arg10[BUFLEN];
  snprintf(arg10, BUFLEN, "--conductivity=%g", EXACT_TRANSIENT_SOLUTION_TRANSMISSIVITY / layerThickness);
  char arg11[BUFLEN];
  snprintf(arg11, BUFLEN, "--layerThickness=%f", layerThickness);

  char *argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11};
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

  CUAS::CUASSolver solver(&model, &args, solutionHandler.get());
  solver.setup();
  solver.solve(time.timeSteps);

  // compute analytical solution
  PETScGrid exactHead(nx, ny);
  {
    // Note, for this 30day experiment, the last time steps are
    // ...
    // [2022-10-20 23:09:56.733] [console] [info] time(30/31) = 2592000 s, dt = 86400 s  <-- 30 days
    // [2022-10-20 23:09:56.738] [console] [info] time(31/31) = 2678400 s, dt = -9999 s  <-- update diagnostics
    //
    // That's because CUAS uses the n time.timeSteps to compute n-1 dt's,
    // and thus we only have n-1 solutions.
    // |---------|---------|--------- ... |---------|-------|
    // |<  dt1  >|<  dt2  >|<  dt3  > ... |< dtn-1 >|
    // t0        t1        t2         ... tn-1      tn
    auto t = (PetscScalar)time.timeSteps[time.timeSteps.size() - 2];  // second last element
    auto nRows = exactHead.getLocalNumOfRows();                       // y-dir
    auto nCols = exactHead.getLocalNumOfCols();                       // x-dir
    auto const cornerX = exactHead.getCornerX();
    auto const cornerY = exactHead.getCornerY();
    auto h = exactHead.getWriteHandle();
    for (int i = 0; i < nRows; ++i) {    // y-dir
      for (int j = 0; j < nCols; ++j) {  // x-dir
        auto x = (cornerX + j) * res;    // todo: not sure about x and y here
        auto y = (cornerY + i) * res;
        h(i, j) = head_exact(x, y, Lx, Ly, t);
      }
    }
  }

// Dump to netcdf for further inspection
#ifdef TESTS_DUMP_NETCDF
  PETScGrid diff(nx, ny);
  {
    // Do NOT delete the returned object - it's managed by the UnitTest class.
    const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();
    auto filename = std::string(test_info->test_suite_name())
                        .append(std::string("_-_"))
                        .append(std::string(test_info->name()))
                        .append(std::string(".nc"));
    CUAS::NetCDFFile file(filename, model.Ncols, model.Nrows);
    file.defineGrid("h_cuas", LIMITED);
    file.defineGrid("h_exact", LIMITED);
    file.defineGrid("bndMask", LIMITED);
    file.write("h_cuas", *solver.currHead, 0);
    file.write("h_exact", exactHead, 0);
    file.write("bndMask", *model.bndMask, 0);
  }
#endif

  // compare on the grid
  {
    auto nRows = exactHead.getLocalNumOfRows();
    auto nCols = exactHead.getLocalNumOfCols();
    auto &hexact = exactHead.getReadHandle();
    auto &hcuas = solver.currHead->getReadHandle();
    for (int i = 0; i < nRows; ++i) {    // y-dir
      for (int j = 0; j < nCols; ++j) {  // x-dir
        ASSERT_NEAR(hexact(i, j), hcuas(i, j), 1e-10) << "at i=" << i << ", j=" << j;
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
