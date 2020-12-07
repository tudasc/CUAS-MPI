#ifndef CUAS_MODEL_H
#define CUAS_MODEL_H

#include "PetscGrid.h"

namespace CUAS {

class CUASModel {
 public:
  CUASModel(int numOfCols, int numOfRows);
  ~CUASModel();

  // x = cols, y = rows
  int const Ncols, Nrows;
  PetscScalar *cols, *rows;
  PetscScalar dx;
  PetscGrid *usurf;
  PetscGrid *topg;
  PetscGrid *thk;
  PetscGrid *bnd_mask;
  PetscGrid *Q;
  PetscGrid *p_ice;
  // type not sure; if NULL: assume constant forcing
  void *time_forcing = nullptr;

  void init();
};

}  // namespace CUAS

#endif
