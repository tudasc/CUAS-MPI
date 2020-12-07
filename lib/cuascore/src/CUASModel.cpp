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
}

}  // namespace CUAS
