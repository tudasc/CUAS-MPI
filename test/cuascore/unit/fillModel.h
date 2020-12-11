#ifndef CUAS_FILLMODEL_H
#define CUAS_FILLMODEL_H

#include "CUASModel.h"
#include "physicalConstants.h"

void fillNoData(CUAS::CUASModel &model) {
  model.cols.resize(20);
  model.rows.resize(10);
  // initialize x and y like in python: x = np.arange(20) * 1000.0
  for (int i = 0; i < 20; i++) {
    int val = (i + 1) * 1000;
    model.cols[i] = val;
    if (i < 10) {
      model.rows[i] = val;
    }
  }

  auto usurfLocal2d = model.usurf->getAsGlobal2dArr();
  auto thkGlobal2d = model.thk->getAsGlobal2dArr();
  int numOfCols = model.usurf->getLocalNumOfCols();
  int numOfRows = model.usurf->getLocalNumOfRows();
  int cornerX = model.usurf->getCornerX();
  int cornerY = model.usurf->getCornerY();

  PetscScalar val = 2000.0;
  for (int j = 0; j < numOfRows; ++j) {
    for (int i = 0; i < numOfCols; ++i) {
      usurfLocal2d[j][i] = val - (cornerX + i) * 100;
      thkGlobal2d[j][i] = val - (cornerX + i) * 100;
    }
  }
  model.usurf->setAsGlobal2dArr(usurfLocal2d);
  model.thk->setAsGlobal2dArr(thkGlobal2d);

  model.topg->setZero();

  // bnd-mask: just last row -> DIRICHLET_Flag
  // first row + first&last col -> NOFLOW_Flag
  // inside: 0
  auto bnd_maskLocal2d = model.bnd_mask->getAsGlobal2dArr();
  numOfCols = model.bnd_mask->getLocalNumOfCols();
  numOfRows = model.bnd_mask->getLocalNumOfRows();
  int totalNumOfCols = model.bnd_mask->getTotalNumOfCols();
  int totalNumOfRows = model.bnd_mask->getTotalNumOfRows();
  cornerX = model.bnd_mask->getCornerX();
  cornerY = model.bnd_mask->getCornerY();

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
  model.bnd_mask->setAsGlobal2dArr(bnd_maskLocal2d);

  auto &bmelt = *model.Q;
  bmelt.setConst(1);

  model.time_forcing = NULL;
}

#endif
