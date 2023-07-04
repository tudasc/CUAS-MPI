/**
 * File: exactCircularSteadySolutionTest.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

/*
 * This setup is based on the example 5.3, page 101 in Galuzzo, 2011.
 *
 * Solve: div (T grad h) = -1, for S = Ss*b = 1
 *
 * Galluzzo, Benjamin Jason. "A finite-difference based approach to solving the subsurface fluid flow equation in
 * heterogeneous media." PhD (Doctor of Philosophy) thesis, University of Iowa, 2011.
 * https://doi.org/10.17077/etd.kbfbts1p
 */

#ifndef CUAS_EXACT_CIRCULAR_STEADY_SOLUTION_TEST_H
#define CUAS_EXACT_CIRCULAR_STEADY_SOLUTION_TEST_H

#include "CUASConstants.h"
#include "CUASModel.h"
#include "CUASSolver.h"
#include "Forcing/ConstantForcing.h"
#include "PETScGrid.h"
#include "cxxopts.hpp"
#include "timeparse.h"

//
class TestSolver : public CUAS::CUASSolver {
  using CUAS::CUASSolver::CUASSolver;  // c++11 makes all constructors of the CUAS::CUASSolver visible

 public:
  /** Implements initial conditions for head and transmissivity */
  void setup(const std::vector<PetscScalar> &xAxis, const std::vector<PetscScalar> &yAxis) {
    CUAS::CUASSolver::setup();

    // set initial transmissivity and head
    auto nRows = currHead->getLocalNumOfRows();
    auto nCols = currHead->getLocalNumOfCols();
    auto cornerX = currHead->getCornerX();
    auto cornerY = currHead->getCornerY();
    auto h0 = currHead->getWriteHandle();
    auto T0 = currTransmissivity->getWriteHandle();
    for (int j = 0; j < nRows; ++j) {
      for (int i = 0; i < nCols; ++i) {
        const PetscScalar xi = xAxis[cornerX + i];
        const PetscScalar yj = yAxis[cornerY + j];
        const auto r2 = xi * xi + yj * yj;
        h0(j, i) = (25.0 - r2) / 8.0;
        T0(j, i) = (r2 <= 1.0) ? 0.25 : 2.0;
      }
    }
  }

  /** Exact solution as in Galuzzo, 2011 (page 101) */
  void getExactSolution(const std::vector<PetscScalar> &xAxis, const std::vector<PetscScalar> &yAxis,
                        PETScGrid &result) {
    // todo: check nx and ny matches size of results and xAxis, yAxis

    auto nRows = result.getLocalNumOfRows();
    auto nCols = result.getLocalNumOfCols();
    auto cornerX = result.getCornerX();
    auto cornerY = result.getCornerY();
    auto head = result.getWriteHandle();
    for (int j = 0; j < nRows; ++j) {
      for (int i = 0; i < nCols; ++i) {
        const PetscScalar xi = xAxis[cornerX + i];
        const PetscScalar yj = yAxis[cornerY + j];
        const auto r2 = xi * xi + yj * yj;
        head(j, i) = (r2 <= 1.0) ? 4.0 - r2 : (25.0 - r2) / 8.0;
      }
    }
  }
};

/** Setup dummy domain on [-2.5, 2.5] x [-2.5, 2.5] and basal melt of 1 m/s water equivalent */
std::unique_ptr<CUAS::CUASModel> fillModelData(int nx, int ny) {
  if (nx != ny) {
    CUAS_ERROR("{}: Error nx != ny. Exiting.", __PRETTY_FUNCTION__);
    exit(1);
  }

  if (nx % 2 != 1) {
    CUAS_ERROR("{}: Error nx and ny must be odd. Exiting.", __PRETTY_FUNCTION__);
    exit(1);
  }

  const PetscScalar res = 5.0 / (PetscScalar)(nx - 1);
  auto n = std::floor(std::log2(res));
  if (std::pow(2, n) != res) {
    CUAS_ERROR("{}: Error dx, dy must be a power of two (dx = 2^(-n), n = 1, 2, ...). Exiting.", __PRETTY_FUNCTION__);
    exit(1);
  }

  auto pmodel = std::make_unique<CUAS::CUASModel>(nx, ny);
  auto &model = *pmodel;

  for (int i = 0; i < model.xAxis.size(); ++i) {
    auto v = -2.5 + i * res;
    model.xAxis[i] = v;
    model.yAxis[i] = v;
  }

  model.topg->setZero();
  model.thk->setConst(100.0);  // could be any positive value
  model.bndMask->setConst(COMPUTE_FLAG);
  model.bndMask->setGhostBoundary(DIRICHLET_FLAG);
  model.bndMask->setRealBoundary(DIRICHLET_FLAG);

  PETScGrid src(nx, ny);
  src.setConst(1.0);
  model.Q = std::make_unique<CUAS::ConstantForcing>(src);

  return pmodel;
}

void parseTestArgs(int argc, char **argv, int &nx, int &nt, CUAS::timeSecs &dt, std::string &fileName) {
  cxxopts::Options options("test");
  // clang-format off
  options
    .positional_help("[NX] [NT] [DT] [exactSteadySolutionTest.h]")
    .show_positional_help()
    .add_options()
      ("script", "The script file to execute", cxxopts::value<std::string>())
      ("server", "The server to execute on", cxxopts::value<std::string>())
      ("filenames", "The filename(s) to process", cxxopts::value<std::vector<std::string>>())
      ("nx", "number of nodes in x- and y-direction",
        cxxopts::value<int>()->default_value("21"))
      ("nt", "number of time steps",
        cxxopts::value<int>()->default_value("31"))
      ("dt", "time step length in seconds",
        cxxopts::value<CUAS::timeSecs>()->default_value("86400"))
      ("filename", "The server to execute on",
        cxxopts::value<std::string>());
  // clang-format on
  options.parse_positional({"nx", "nt", "dt", "filename"});

  cxxopts::ParseResult result = options.parse(argc, argv);
  nx = result["nx"].as<int>();
  nt = result["nt"].as<int>();
  dt = result["dt"].as<CUAS::timeSecs>();
  if (result.count("filename")) {
    fileName = result["filename"].as<std::string>();
  }  // no default
}

#endif  // CUAS_EXACT_CIRCULAR_STEADY_SOLUTION_TEST_H