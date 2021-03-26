#ifndef CUAS_SOLVER_H
#define CUAS_SOLVER_H

#include "CUASArgs.h"
#include "CUASModel.h"

#include "PETScGrid.h"

namespace CUAS {

class CUASSolver {
 public:
  explicit CUASSolver(CUASModel *model, CUASArgs const *const args) : model(model), args(args) {
    int numOfCols = model->Ncols;
    int numOfRows = model->Nrows;

    S = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    Sp = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    K = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    T = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    T_n = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    noFlowMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    gradMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    seaLevelForcingMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    dirichletMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    dirichletValues = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  }
  CUASSolver(CUASSolver &) = delete;
  CUASSolver(CUASSolver &&) = delete;

  void solve(std::unique_ptr<PETScGrid> &u, std::unique_ptr<PETScGrid> &u_n, int const Nt,
             PetscScalar const totaltime_secs, PetscScalar const dt_secs);

  void setup();

 private:
  std::unique_ptr<PETScGrid> S;
  std::unique_ptr<PETScGrid> Sp;
  std::unique_ptr<PETScGrid> K;
  std::unique_ptr<PETScGrid> T;
  std::unique_ptr<PETScGrid> T_n;
  std::unique_ptr<PETScGrid> noFlowMask;
  std::unique_ptr<PETScGrid> gradMask;
  std::unique_ptr<PETScGrid> seaLevelForcingMask;
  std::unique_ptr<PETScGrid> dirichletMask;
  std::unique_ptr<PETScGrid> dirichletValues;

  CUASModel *const model;
  CUASArgs const *const args;
};

}  // namespace CUAS

#endif
