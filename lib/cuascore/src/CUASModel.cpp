#include "CUASModel.h"

#include "physicalConstants.h"

namespace CUAS {

CUASModel::CUASModel(int numOfCols, int numOfRows) : Ncols(numOfCols), Nrows(numOfRows) {
  usurf = new PetscGrid(numOfCols, numOfRows);
  topg = new PetscGrid(numOfCols, numOfRows);
  thk = new PetscGrid(numOfCols, numOfRows);
  bnd_mask = new PetscGrid(numOfCols, numOfRows);
  Q = new PetscGrid(numOfCols, numOfRows);
  p_ice = new PetscGrid(numOfCols, numOfRows);

  // grids for setup
  S = new PetscGrid(numOfCols, numOfRows);
  Sp = new PetscGrid(numOfCols, numOfRows);
  K = new PetscGrid(numOfCols, numOfRows);
  T = new PetscGrid(numOfCols, numOfRows);
  T_n = new PetscGrid(numOfCols, numOfRows);
  no_flow_mask = new PetscGrid(numOfCols, numOfRows);
  grad_mask = new PetscGrid(numOfCols, numOfRows);
  sea_level_forcing_mask = new PetscGrid(numOfCols, numOfRows);
  dirichlet_mask = new PetscGrid(numOfCols, numOfRows);
  dirichlet_values = new PetscGrid(numOfCols, numOfRows);
}

// call init before using Model but after thk has been set
void CUASModel::init() {
  dx = cols[1] - cols[0];

  // p_ice = thk * RHO_ICE * GRAVITY (python)
  PetscScalar **thk2d = thk->getAsGlobal2dArr();
  PetscScalar **p_iceGlobal = p_ice->getAsGlobal2dArr();
  for (int j = 0; j < thk->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < thk->getLocalNumOfCols(); ++i) {
      p_iceGlobal[j][i] = thk2d[j][i] * RHO_ICE * GRAVITY;
    }
  }
  p_ice->setAsGlobal2dArr(p_iceGlobal);
  thk->restoreGlobal2dArr(thk2d);
}

CUASModel::~CUASModel() {
  delete usurf;
  delete topg;
  delete thk;
  delete bnd_mask;
  delete Q;
  delete p_ice;

  // delete grids for setup
  delete S;
  delete Sp;
  delete K;
  delete T;
  delete T_n;
  delete no_flow_mask;
  delete grad_mask;
  delete sea_level_forcing_mask;
  delete dirichlet_mask;
  delete dirichlet_values;
}

}  // namespace CUAS
