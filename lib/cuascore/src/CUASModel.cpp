#include "CUASModel.h"

#include "physicalConstants.h"

namespace CUAS {

CUASModel::CUASModel(int numOfCols, int numOfRows) : Ncols(numOfCols), Nrows(numOfRows) {
  usurf = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  topg = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  thk = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  bnd_mask = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  Q = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  p_ice = std::make_unique<PetscGrid>(numOfCols, numOfRows);

  // grids for setup
  S = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  Sp = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  K = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  T = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  T_n = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  no_flow_mask = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  grad_mask = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  sea_level_forcing_mask = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  dirichlet_mask = std::make_unique<PetscGrid>(numOfCols, numOfRows);
  dirichlet_values = std::make_unique<PetscGrid>(numOfCols, numOfRows);
}

// call init before using Model but after thk has been set
void CUASModel::init() {
  dx = cols[1] - cols[0];

  // p_ice = thk * RHO_ICE * GRAVITY (python)
  auto &thk2d = thk->getReadHandle();
  auto p_iceGlobal = p_ice->getWriteHandle();
  for (int j = 0; j < thk->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < thk->getLocalNumOfCols(); ++i) {
      p_iceGlobal(j, i) = thk2d(j, i) * RHO_ICE * GRAVITY;
    }
  }
}

}  // namespace CUAS
