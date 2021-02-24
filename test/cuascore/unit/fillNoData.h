#ifndef CUAS_FILLMODEL_H
#define CUAS_FILLMODEL_H

#include "CUASConstants.h"
#include "CUASModel.h"

#define NODATA_COLS 20
#define NODATA_ROWS 10

std::unique_ptr<CUAS::CUASModel> fillNoData() {
  auto pmodel = std::make_unique<CUAS::CUASModel>(NODATA_COLS, NODATA_ROWS);
  auto &model = *pmodel;

  // initialize x and y like in python: x = np.arange(20) * 1000.0
  for (int i = 0; i < model.cols.size(); ++i) {
    model.cols[i] = i * 1000;
  }
  for (int i = 0; i < model.rows.size(); ++i) {
    model.rows[i] = i * 1000;
  }

  {
    auto usurf2d = model.usurf->getWriteHandleGhost();
    auto thk2d = model.thk->getWriteHandleGhost();
    auto const cornerX = model.usurf->getCornerX();
    for (int j = 0; j < model.usurf->getLocalGhostNumOfRows(); ++j) {
      for (int i = 0; i < model.usurf->getLocalGhostNumOfCols(); ++i) {
        usurf2d(j, i) = 2000.0 - (cornerX + i) * 100;
        thk2d(j, i) = 2000.0 - (cornerX + i) * 100;
      }
    }
  }

  model.topg->setZero();

  // bnd-mask: just last col -> DIRICHLET_Flag
  // first col + first&last row -> NOFLOW_Flag
  // inside: 0
  // first set total boundary to NOFLOW_FLAG
  model.bndMask->setGlobalBoundariesConst(NOFLOW_FLAG);
  // second set last col, do not override corners
  if (model.bndMask->getCornerXGhost() + model.bndMask->getLocalGhostNumOfCols() ==
      model.bndMask->getTotalGhostNumOfCols() - 1) {
    auto bndMask2d = model.bndMask->getWriteHandleGhost();
    for (int i = 0; i < model.bndMask->getLocalGhostNumOfRows(); ++i) {
      bndMask2d(i, model.bndMask->getLocalGhostNumOfCols() - 1) = DIRICHLET_FLAG;
    }
    // reset corner
    if (model.bndMask->getCornerY() == 0) {
      bndMask2d(0, 7) = NOFLOW_FLAG;
    }
    if (model.bndMask->getCornerYGhost() + model.bndMask->getLocalGhostNumOfRows() ==
        model.bndMask->getTotalGhostNumOfRows() - 1) {
      bndMask2d(5, 7) = NOFLOW_FLAG;
    }
  }

  // TODO bmelt vs Q
  auto &bmelt = *model.Q;
  bmelt.setConst(1);

  return pmodel;
}

#endif
