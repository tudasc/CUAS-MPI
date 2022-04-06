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
      : model(model), args(args), solutionHandler(solutionHandler) {
    int numOfCols = model->Ncols;
    int numOfRows = model->Nrows;

    nextHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    currHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    S = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    Sp = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    K = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    nextTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    currTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    gradMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    dirichletValues = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    sol = std::make_unique<PETScVector>(numOfCols * numOfRows);
  }
  CUASSolver(CUASSolver &) = delete;
  CUASSolver(CUASSolver &&) = delete;

  void solve(std::vector<CUAS::timeSecs> &timeSteps);

  void setup();

 public:
  std::unique_ptr<PETScGrid> nextHead;  // unknown head at new time level
  std::unique_ptr<PETScGrid> currHead;  // head at the current time level
  std::unique_ptr<PETScVector> sol;

 private:
  std::unique_ptr<PETScGrid> S;
  std::unique_ptr<PETScGrid> Sp;
  std::unique_ptr<PETScGrid> K;
  std::unique_ptr<PETScGrid> nextTransmissivity;
  std::unique_ptr<PETScGrid> currTransmissivity;
  std::unique_ptr<PETScGrid> gradMask;
  std::unique_ptr<PETScGrid> dirichletValues;

  CUASModel *const model;
  CUASArgs const *const args;
  CUAS::SolutionHandler *const solutionHandler;
};

}  // namespace CUAS

#endif
