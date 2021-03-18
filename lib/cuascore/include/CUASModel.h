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

  // type not sure; if NULL: assume constant forcing
  void *time_forcing = nullptr;

  void init();
};

}  // namespace CUAS

#endif
