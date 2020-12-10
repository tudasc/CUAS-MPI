#ifndef CUAS_FILLMODEL_H
#define CUAS_FILLMODEL_H

#include "CUASModel.h"
#include "physicalConstants.h"

void fillNoData(CUAS::CUASModel &model) {
  PetscScalar x[20];
  PetscScalar y[10];
  // initialize x and y like in python: x = np.arange(20) * 1000.0
  for (int i = 0; i < 20; i++) {
    int val = (i + 1) * 1000;
    x[i] = val;
    if (i < 10) {
      y[i] = val;
    }
  }

  model.cols = x;
  model.rows = y;

  PetscGrid *usurf = model.usurf;
  PetscGrid *thk = model.thk;
  PetscScalar **usurfLocal2d = usurf->getAsGlobal2dArr();
  PetscScalar **thkGlobal2d = thk->getAsGlobal2dArr();
  int numOfCols = usurf->getLocalNumOfCols();
  int numOfRows = usurf->getLocalNumOfRows();
  int cornerX = usurf->getCornerX();
  int cornerY = usurf->getCornerY();

  PetscScalar val = 2000.0;
  for (int j = 0; j < numOfRows; ++j) {
    for (int i = 0; i < numOfCols; ++i) {
      usurfLocal2d[j][i] = val - (cornerX + i) * 100;
      thkGlobal2d[j][i] = val - (cornerX + i) * 100;
    }
  }
  usurf->setAsGlobal2dArr(usurfLocal2d);
  thk->setAsGlobal2dArr(thkGlobal2d);

  model.topg->setZero();

  // bnd-mask: just last row -> DIRICHLET_Flag
  // first row + first&last col -> NOFLOW_Flag
  // inside: 0
  auto bnd_mask = model.bnd_mask;
  PetscScalar **bnd_maskLocal2d = bnd_mask->getAsGlobal2dArr();
  numOfCols = bnd_mask->getLocalNumOfCols();
  numOfRows = bnd_mask->getLocalNumOfRows();
  int totalNumOfCols = bnd_mask->getTotalNumOfCols();
  int totalNumOfRows = bnd_mask->getTotalNumOfRows();
  cornerX = bnd_mask->getCornerX();
  cornerY = bnd_mask->getCornerY();

  for (int j = 0; j < numOfRows; ++j) {
    for (int i = 0; i < numOfCols; ++i) {
      if (cornerY + j == 0) {
        bnd_maskLocal2d[j][i] = NOFLOW_FLAG;
      } else if (cornerY + j == totalNumOfRows - 1) {
        bnd_maskLocal2d[j][i] = NOFLOW_FLAG;
      } else if (cornerX + i == 0) {
        bnd_maskLocal2d[j][i] = NOFLOW_FLAG;
      } else if (cornerX + i == totalNumOfCols - 1) {
        bnd_maskLocal2d[j][i] = DIRICHLET_FLAG;
      } else {
        bnd_maskLocal2d[j][i] = 0;
      }
    }
  }
  bnd_mask->setAsGlobal2dArr(bnd_maskLocal2d);

  PetscGrid *bmelt = model.Q;
  bmelt->setConst(1);

  model.time_forcing = NULL;
}

#endif
