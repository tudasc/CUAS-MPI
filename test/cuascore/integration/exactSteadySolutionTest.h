#ifndef CUAS_EXACT_STEADY_SOLUTION_TEST_H
#define CUAS_EXACT_STEADY_SOLUTION_TEST_H

#include "CUASConstants.h"
#include "CUASModel.h"
#include "Forcing/Forcing.h"
#include "PETScGrid.h"
#include "cxxopts.hpp"
#include "timeparse.h"

#ifndef EXACT_STEADY_SOLUTION_STORATIVITY
#define EXACT_STEADY_SOLUTION_STORATIVITY 1.0e-6
// reach "steady-state" within 31 days, but with visible transient, 5.0e-4 is too high.
// #define EXACT_STEADY_SOLUTION_STORATIVITY 2.0e-4
#endif

#ifndef EXACT_STEADY_SOLUTION_TRANSMISSIVITY
#define EXACT_STEADY_SOLUTION_TRANSMISSIVITY 1.0
#endif

#ifndef EXACT_STEADY_SOLUTION_H0
#define EXACT_STEADY_SOLUTION_H0 10.0
#endif

#ifndef EXACT_STEADY_SOLUTION_HMAX
#define EXACT_STEADY_SOLUTION_HMAX 100.0
#endif

#ifndef EXACT_STEADY_SOLUTION_AFAC
#define EXACT_STEADY_SOLUTION_AFAC (16.0 * (EXACT_STEADY_SOLUTION_HMAX - EXACT_STEADY_SOLUTION_H0))
#endif

/** computes the exact solution for the head
 *
 * @return
 */
inline PetscScalar head_exact(PetscScalar const x, PetscScalar const y, PetscScalar const Lx, PetscScalar const Ly) {
  constexpr PetscScalar h0 = EXACT_STEADY_SOLUTION_H0;
  constexpr PetscScalar a = EXACT_STEADY_SOLUTION_AFAC;
  const auto c = Lx * Lx * Ly * Ly;
  return a * x * y * (Lx - x) * (Ly - y) / c + h0;
}

/** computes d^2 h / dx^2 based on the exact solution
 *
 * @return
 */
inline PetscScalar d2hdx2_exact([[maybe_unused]] PetscScalar const x, PetscScalar const y, PetscScalar const Lx,
                                PetscScalar const Ly) {
  constexpr PetscScalar a = EXACT_STEADY_SOLUTION_AFAC;
  const auto c = Lx * Lx * Ly * Ly;
  return -2.0 * a * y * (Ly - y) / c;
}

/** computes d^2 h / dy^2 based on the exact solution
 *
 * @return
 */
inline PetscScalar d2hdy2_exact(PetscScalar const x, [[maybe_unused]] PetscScalar const y, PetscScalar const Lx,
                                PetscScalar const Ly) {
  constexpr PetscScalar a = EXACT_STEADY_SOLUTION_AFAC;
  const auto c = Lx * Lx * Ly * Ly;
  return -2.0 * a * x * (Lx - x) / c;
}

/** Computes the time-dependent maximum head of the exact solution
 *
 * @return maxHead(t) in m/s
 */
inline PetscScalar max_head_exact() {
  constexpr PetscScalar h0 = EXACT_STEADY_SOLUTION_H0;
  constexpr PetscScalar a = EXACT_STEADY_SOLUTION_AFAC;
  return a / 16.0 + h0;
}

/** Computes the maximum melt of the exact solution
 *
 * @return maxQ in m/s
 */
inline PetscScalar max_Q_exact(PetscScalar const Lx, PetscScalar const Ly) {
  constexpr PetscScalar T = EXACT_STEADY_SOLUTION_TRANSMISSIVITY;
  constexpr PetscScalar a = EXACT_STEADY_SOLUTION_AFAC;
  const auto c = Lx * Lx * Ly * Ly;
  constexpr PetscScalar T_hmean = 2.0 * T * T / (2.0 * T + TINY);
  return 0.5 * a * T_hmean * ((Lx * Lx) + (Ly * Ly)) / c;
}

class SteadyForcing : public CUAS::Forcing {
 public:
  explicit SteadyForcing(int const nx, int const ny, PetscScalar const res, PetscScalar const multiplier = 1.0,
                         PetscScalar const offset = 0.0)
      : resolution(res), Lx(res * (nx - 1)), Ly(res * (ny - 1)) {
    currQ = std::make_unique<PETScGrid>(nx, ny);

    if (multiplier != 1.0) {
      SteadyForcing::applyMultiplier(multiplier);
    }
    if (offset != 0.0) {
      SteadyForcing::applyOffset(offset);
    }
  }
  SteadyForcing(SteadyForcing &) = delete;
  SteadyForcing(SteadyForcing &&) = delete;

  PETScGrid const &getCurrentQ(CUAS::timeSecs currTime) override {
    {
      // some abbreviations
      constexpr PetscScalar T = EXACT_STEADY_SOLUTION_TRANSMISSIVITY;
      // CUAS uses the harmonic mean to build the solution matrix.
      // See lib/cuascore/src/systemmatrix.cpp and thus,
      // T_hmean = 2.0 * T1 * T2 / (T1 + T2 + TINY) =  2*T^2 / (2*T + TINY) < T_exact
      constexpr PetscScalar T_hmean = 2.0 * T * T / (2.0 * T + TINY);

      // compute source term Q on the grid
      auto nRows = currQ->getLocalNumOfRows();  // y-dir
      auto nCols = currQ->getLocalNumOfCols();  // x-dir
      auto const cornerX = currQ->getCornerX();
      auto const cornerY = currQ->getCornerY();
      auto currQWrite = currQ->getWriteHandle();
      for (int i = 0; i < nRows; ++i) {    // y-dir
        for (int j = 0; j < nCols; ++j) {  // x-dir
          auto x = (cornerX + j) * resolution;
          auto y = (cornerY + i) * resolution;
          auto d2hdx2 = d2hdx2_exact(x, y, Lx, Ly);
          auto d2hdy2 = d2hdy2_exact(x, y, Lx, Ly);
          currQWrite(i, j) = -T_hmean * (d2hdx2 + d2hdy2);  // T_hmean instead of T
        }
      }
    }
    return *currQ;
  }

 private:
  std::unique_ptr<PETScGrid> currQ;
  PetscScalar resolution;
  PetscScalar Lx, Ly;

  void applyMultiplier(PetscScalar multiplier) override {
    if (multiplier != 1.0) {
      currQ->applyMultiplier(multiplier);
    }
  }

  void applyOffset(PetscScalar offset) override {
    if (offset != 0.0) {
      currQ->applyOffset(offset);
    }
  }
};

std::unique_ptr<CUAS::CUASModel> fillModelData(int nx, int ny, PetscScalar res) {
  auto pmodel = std::make_unique<CUAS::CUASModel>(nx, ny);
  auto &model = *pmodel;

  // initialize x and y like in python: x = np.arange(20) * 1000.0
  for (int i = 0; i < model.xAxis.size(); ++i) {
    model.xAxis[i] = i * res;
  }
  for (int i = 0; i < model.yAxis.size(); ++i) {
    model.yAxis[i] = i * res;
  }

  model.topg->setZero();
  model.thk->setConst(100.0);  // could be any positive value
  model.bndMask->setConst(COMPUTE_FLAG);
  model.bndMask->setGhostBoundary(DIRICHLET_FLAG);
  model.bndMask->setRealBoundary(DIRICHLET_FLAG);
  model.Q = std::make_unique<SteadyForcing>(nx, ny, res);

  return pmodel;
}

void parseTestArgs(int argc, char **argv, PetscScalar &res, int &nx, int &ny, int &nt, CUAS::timeSecs &dt,
                   std::string &fileName) {
  cxxopts::Options options("test");
  // clang-format off
  options
    .positional_help("[RES] [NX] [NY] [NT] [DT] [FILENAME]")
    .show_positional_help()
    .add_options()
      ("script", "The script file to execute", cxxopts::value<std::string>())
      ("server", "The server to execute on", cxxopts::value<std::string>())
      ("filenames", "The filename(s) to process", cxxopts::value<std::vector<std::string>>())
      ("res", "horizontal resolution in meter",
        cxxopts::value<PetscScalar>()->default_value("50000.0"))
      ("nx", "number of nodes in x-direction",
        cxxopts::value<int>()->default_value("21"))
      ("ny", "number of nodes in y-direction",
        cxxopts::value<int>()->default_value("11"))
      ("nt", "number of time steps",
        cxxopts::value<int>()->default_value("31"))
      ("dt", "time step length in seconds",
        cxxopts::value<CUAS::timeSecs>()->default_value("86400"))
      ("filename", "The server to execute on",
        cxxopts::value<std::string>());
  // clang-format on
  options.parse_positional({"res", "nx", "ny", "nt", "dt", "filename"});

  cxxopts::ParseResult result = options.parse(argc, argv);
  res = result["res"].as<PetscScalar>();
  nx = result["nx"].as<int>();
  ny = result["ny"].as<int>();
  nt = result["nt"].as<int>();
  dt = result["dt"].as<CUAS::timeSecs>();
  if (result.count("filename")) {
    fileName = result["filename"].as<std::string>();
  }  // no default
}

#endif