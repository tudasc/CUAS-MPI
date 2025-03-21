/**
 * File: exactPiecewiseConstantSteadySolutionTest.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

/*
 * This setup is based on the example
 *
 * H. P. Langtangen and S. Linge, Finite Difference Computing with PDEs:
 * A Modern Software Approach. Springer International Publishing, 2017.
 * doi: 10.1007/978-3-319-55456-3.
 *
 *  We solve the steady solution for the 1D problem
 *       d/dx (T(x) dh/dx) = 0
 *  for a medium build of M layers with different material properties and
 *  with the layer boundaries b_0, \dots, b_M, where b_0 == 0 and b_M == Lx.
 *
 */

#ifndef CUAS_EXACT_PIECEWISE_CONSTANT_STEADY_SOLUTION_TEST_H
#define CUAS_EXACT_PIECEWISE_CONSTANT_STEADY_SOLUTION_TEST_H

#include "CUASConstants.h"
#include "CUASModel.h"
#include "CUASSolver.h"
#include "Forcing/SteadyForcing.h"
#include "PETScGrid.h"
#include "cxxopts.hpp"
#include "timeparse.h"
#include <cmath>  // std::floor, std::pow, std::log2

const std::vector<PetscScalar> layerBoundaries = {0.0, 0.2, 0.5, 1};     // len M+1
const std::vector<PetscScalar> layerTransmissivities = {0.2, 0.4, 4.0};  // len M
constexpr PetscScalar head_x0 = 0.5;                                     // left Dirichlet boundary condition
constexpr PetscScalar head_xL = 5.0;                                     // right Dirichlet boundary condition

/** piecewise constant function for transmissivity */
inline PetscScalar transmissivityFuncX(PetscScalar x) {
  auto numLayers = layerTransmissivities.size();
  assert(numLayers == layerBoundaries.size() - 1);

  // early exit for right margin
  if (x == layerBoundaries[numLayers]) {
    return layerTransmissivities[numLayers - 1];
  }
  // very inefficient way to find the actual layer and layer transmisivity
  for (auto i = 0; i < numLayers; ++i) {
    if ((layerBoundaries[i] <= x) && (x < layerBoundaries[i + 1])) {
      return layerTransmissivities[i];
    }
  }
  // should not be reached
  return NOFLOW_VALUE;
}

/**
 *
 */
inline PetscScalar exactSolutionFuncX(PetscScalar x) {
  auto numLayers = layerTransmissivities.size();
  assert(numLayers == layerBoundaries.size() - 1);

  // get layer number m

  // it points to the first element greater or equal to x
  auto it = std::lower_bound(layerBoundaries.begin(), layerBoundaries.end(), x);
  auto m = std::distance(layerBoundaries.begin(), it) - 1;
  m = (m < 0) ? 0 : m;

  // integrate through the first m-1 layers and then add the
  // contribution from the remaining part x - b_m into the m-th layer
  PetscScalar numerator = 0.0;
  for (int i = 0; i < m; ++i) {
    numerator += (layerBoundaries[i + 1] - layerBoundaries[i]) / transmissivityFuncX(layerBoundaries[i]);
  }
  // remaining part
  numerator += (x - layerBoundaries[m]) / transmissivityFuncX(layerBoundaries[m]);

  // the denominator in the exact formula is a constant for all x
  PetscScalar denominator = 0.0;
  for (int i = 0; i < numLayers; ++i) {
    denominator += (layerBoundaries[i + 1] - layerBoundaries[i]) / transmissivityFuncX(layerBoundaries[i]);
  }

  if (denominator > 0.0) {
    return head_x0 + (head_xL - head_x0) * numerator / denominator;
  } else {
    // should not happen!
    return head_x0;
  }
}

//
class TestSolver : public CUAS::CUASSolver {
  using CUAS::CUASSolver::CUASSolver;  // c++11 makes all constructors of the CUAS::CUASSolver visible

 public:
  /** Implements initial conditions for head and transmissivity */
  void setup(const std::vector<PetscScalar> &xAxis, const std::vector<PetscScalar> &yAxis) {
    CUAS::CUASSolver::setup();

    // set bnd_mask, initial transmissivity and head
    auto nRows = currHead->getLocalNumOfRows();
    auto nCols = currHead->getLocalNumOfCols();
    auto cornerX = currHead->getCornerX();
    auto cornerY = currHead->getCornerY();
    auto T0 = currTransmissivity->getWriteHandle();
    for (int row = 0; row < nRows; ++row) {
      for (int col = 0; col < nCols; ++col) {
        const PetscScalar xi = xAxis[cornerX + col];
        T0(row, col) = transmissivityFuncX(xi);
      }
    }
    currHead->setConst(0.5 * (head_x0 + head_xL));
    currHead->setRealBoundary(head_x0, PETScGrid::Direction::West);
    currHead->setRealBoundary(head_xL, PETScGrid::Direction::East);

    // no-flow transmissivity (the test should not depend on the solver to correct this)
    // the order is important for the corner points
    currTransmissivity->setRealBoundary(NOFLOW_VALUE, PETScGrid::Direction::North);
    currTransmissivity->setRealBoundary(NOFLOW_VALUE, PETScGrid::Direction::South);
    currTransmissivity->setRealBoundary(layerTransmissivities.back(), PETScGrid::Direction::East);
    currTransmissivity->setRealBoundary(layerTransmissivities.front(), PETScGrid::Direction::West);
  }

  /** TODO */
  void getExactSolution(const std::vector<PetscScalar> &xAxis, const std::vector<PetscScalar> &yAxis,
                        PETScGrid &result) {
    // todo: check nx and ny matches size of results and xAxis, yAxis

    auto nRows = result.getLocalNumOfRows();
    auto nCols = result.getLocalNumOfCols();
    auto cornerX = result.getCornerX();
    auto head = result.getWriteHandle();
    for (int row = 0; row < nRows; ++row) {
      for (int col = 0; col < nCols; ++col) {
        const PetscScalar xi = xAxis[cornerX + col];
        head(row, col) = exactSolutionFuncX(xi);
      }
    }
  }
};

/** Setup dummy domain on [0, 1] x [-Ly/2, Ly/2] and basal melt of 0 m/s water equivalent */
std::unique_ptr<CUAS::CUASModel> fillModelData(int nx, int ny) {
  if (nx % 2 != 1) {
    CUAS_ERROR("{}: Error nx and ny must be odd. Exiting.", __PRETTY_FUNCTION__);
    exit(1);
  }
  if (ny % 2 != 1) {
    CUAS_ERROR("{}: Error ny and ny must be odd. Exiting.", __PRETTY_FUNCTION__);
    exit(1);
  }

  auto pmodel = std::make_unique<CUAS::CUASModel>(nx, ny);
  auto &model = *pmodel;

  const PetscScalar Lx = layerBoundaries.back() - layerBoundaries.front();
  const PetscScalar res = Lx / (PetscScalar)(nx - 1);

  auto n = std::floor(std::log2(res));
  if (std::pow(2, n) != res) {
    CUAS_WARN_RANK0("{}: dx = {} found. A power of two (dx = 2^(-n), n = 1, 2, ...) is preferred.", __PRETTY_FUNCTION__,
                    res);
  }
  for (int i = 0; i < model.xAxis.size(); ++i) {
    model.xAxis[i] = layerBoundaries.front() + i * res;
  }

  const PetscScalar Ly = (ny - 1) * res;
  for (int i = 0; i < model.yAxis.size(); ++i) {
    model.yAxis[i] = -Ly / 2.0 + i * res;
  }

  model.topg->setZero();
  model.thk->setConst(100.0);  // could be any positive value

  model.bndMask->setConst(COMPUTE_FLAG);
  model.bndMask->setGhostBoundary(NOFLOW_FLAG);  // fixme: Why? (tkleiner)
  // the order is important for the corner points
  model.bndMask->setRealBoundary(NOFLOW_FLAG, PETScGrid::Direction::North);
  model.bndMask->setRealBoundary(NOFLOW_FLAG, PETScGrid::Direction::South);
  model.bndMask->setRealBoundary(DIRICHLET_FLAG, PETScGrid::Direction::East);
  model.bndMask->setRealBoundary(DIRICHLET_FLAG, PETScGrid::Direction::West);

  // no source term for this experiment
  PETScGrid src(nx, ny);
  src.setConst(0.0);
  model.setWaterSource(std::make_unique<CUAS::SteadyForcing>(src));

  return pmodel;
}

void parseTestArgs(int argc, char **argv, int &nx, int &nt, CUAS::timeSecs &dt, std::string &fileName) {
  cxxopts::Options options("test");
  // clang-format off
  options
    .positional_help("[NX] [NT] [DT] [exactPiecewiseConstantSteadySolutionTest.h]")
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

#endif  // CUAS_EXACT_PIECEWISE_CONSTANT_STEADY_SOLUTION_TEST_H