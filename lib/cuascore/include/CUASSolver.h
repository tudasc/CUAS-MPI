#ifndef CUAS_SOLVER_H
#define CUAS_SOLVER_H

#include "CUASArgs.h"
#include "CUASModel.h"
#include "SolutionHandler.h"
#include "timeparse.h"

#include "PETScGrid.h"

namespace CUAS {

class CUASSolver {
 public:
  explicit CUASSolver(CUASModel *model, CUASArgs const *const args, CUAS::SolutionHandler *solutionHandler = nullptr)
      : model(model), args(args), solutionHandler(solutionHandler), numOfCols(model->Ncols), numOfRows(model->Nrows) {
    nextHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    currHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    nextTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    currTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    gradMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    dirichletValues = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    // rate factor from flow law (ice rheology)
    rateFactorIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    // basal velocity of ice
    basalVelocityIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    sol = std::make_unique<PETScVector>(numOfCols * numOfRows);
  }
  CUASSolver(CUASSolver &) = delete;
  CUASSolver(CUASSolver &&) = delete;

  void solve(std::vector<CUAS::timeSecs> &timeSteps);

  void setup();

 public:
  std::unique_ptr<PETScGrid> nextHead;  // unknown head at new time level
  std::unique_ptr<PETScGrid> currHead;  // head at the current time level
  std::unique_ptr<PETScGrid> nextTransmissivity;
  std::unique_ptr<PETScGrid> currTransmissivity;
  std::unique_ptr<PETScVector> sol;

 private:
  int const numOfCols;
  int const numOfRows;

  std::unique_ptr<PETScGrid> gradMask;
  std::unique_ptr<PETScGrid> dirichletValues;

  std::unique_ptr<PETScGrid> rateFactorIce;
  std::unique_ptr<PETScGrid> basalVelocityIce;

  CUASModel *const model;
  CUASArgs const *const args;
  CUAS::SolutionHandler *const solutionHandler;
};

}  // namespace CUAS

#endif
