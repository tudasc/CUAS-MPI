#ifndef CUAS_MODEL_H
#define CUAS_MODEL_H

#include "PETScGrid.h"

#include <memory>
#include <vector>

namespace CUAS {

class CUASModel {
 public:
  CUASModel(int numOfCols, int numOfRows);
  //~CUASModel();

  // x = cols, y = rows
  int const Ncols, Nrows;
  PetscScalar dx, dy;
  std::vector<PetscScalar> cols, rows;
  std::unique_ptr<PETScGrid> usurf;
  std::unique_ptr<PETScGrid> topg;
  std::unique_ptr<PETScGrid> thk;
  std::unique_ptr<PETScGrid> bndMask;
  std::unique_ptr<PETScGrid> Q;
  std::unique_ptr<PETScGrid> pIce;

  // grids for setup
  std::unique_ptr<PETScGrid> S = nullptr;
  std::unique_ptr<PETScGrid> Sp = nullptr;
  std::unique_ptr<PETScGrid> K = nullptr;
  std::unique_ptr<PETScGrid> T = nullptr;
  std::unique_ptr<PETScGrid> T_n = nullptr;
  std::unique_ptr<PETScGrid> noFlowMask = nullptr;
  std::unique_ptr<PETScGrid> gradMask = nullptr;
  std::unique_ptr<PETScGrid> seaLevelForcingMask = nullptr;
  std::unique_ptr<PETScGrid> dirichletMask = nullptr;
  std::unique_ptr<PETScGrid> dirichletValues = nullptr;

  // type not sure; if NULL: assume constant forcing
  void *time_forcing = nullptr;

  void init();
};

}  // namespace CUAS

#endif
