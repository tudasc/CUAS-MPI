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
  std::vector<PetscScalar> cols, rows;
  PetscScalar dx;
  std::unique_ptr<PetscGrid> usurf;
  std::unique_ptr<PetscGrid> topg;
  std::unique_ptr<PetscGrid> thk;
  std::unique_ptr<PetscGrid> bnd_mask;
  std::unique_ptr<PetscGrid> Q;
  std::unique_ptr<PetscGrid> p_ice;

  // grids for setup
  std::unique_ptr<PetscGrid> S;
  std::unique_ptr<PetscGrid> Sp;
  std::unique_ptr<PetscGrid> K;
  std::unique_ptr<PetscGrid> T;
  std::unique_ptr<PetscGrid> T_n;
  std::unique_ptr<PetscGrid> no_flow_mask;
  std::unique_ptr<PetscGrid> grad_mask;
  std::unique_ptr<PetscGrid> sea_level_forcing_mask;
  std::unique_ptr<PetscGrid> dirichlet_mask;
  std::unique_ptr<PetscGrid> dirichlet_values;

  // type not sure; if NULL: assume constant forcing
  void *time_forcing = nullptr;

  void init();
};

}  // namespace CUAS

#endif
