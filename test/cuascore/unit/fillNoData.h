/**
 * File: fillNoData.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_FILLMODEL_H
#define CUAS_FILLMODEL_H

#include "CUASConstants.h"
#include "CUASModel.h"
#include "Forcing/SteadyForcing.h"

#define NODATA_COLS 20
#define NODATA_ROWS 10

std::unique_ptr<CUAS::CUASModel> fillNoData() {
  auto pmodel = std::make_unique<CUAS::CUASModel>(NODATA_COLS, NODATA_ROWS);
  auto &model = *pmodel;

  // initialize x and y like in python: x = np.arange(20) * 1000.0
  for (int i = 0; i < model.xAxis.size(); ++i) {
    model.xAxis[i] = i * 1000;
  }
  for (int i = 0; i < model.yAxis.size(); ++i) {
    model.yAxis[i] = i * 1000;
  }

  {
    auto thk2d = model.thk->getWriteHandle();
    auto const cornerX = model.thk->getCornerX();
    for (int j = 0; j < model.thk->getLocalNumOfRows(); ++j) {
      for (int i = 0; i < model.thk->getLocalNumOfCols(); ++i) {
        thk2d(j, i) = 2000.0 - (cornerX + i) * 100;
      }
    }
  }

  model.topg->setZero();

  // bnd-mask: just last col -> DIRICHLET_Flag
  // first col + first&last row -> NOFLOW_Flag
  // inside: 0
  // first set total boundary to NOFLOW_FLAG
  model.bndMask->setGhostBoundary(NOFLOW_FLAG);
  model.bndMask->setRealBoundary(NOFLOW_FLAG);
  // second set last col, do not override corners
  if (model.bndMask->getCornerXGhost() + model.bndMask->getLocalGhostNumOfCols() ==
      model.bndMask->getTotalGhostNumOfCols() - 1) {
    auto bndMask2d = model.bndMask->getWriteHandleGhost();
    for (int i = 0; i < model.bndMask->getLocalGhostNumOfRows(); ++i) {
      if (model.bndMask->getCornerYGhost() + i <= 0 ||
          model.bndMask->getCornerYGhost() + i >= model.bndMask->getTotalGhostNumOfRows() - 3) {
        continue;
      }
      bndMask2d(i, model.bndMask->getLocalGhostNumOfCols() - 1) = DIRICHLET_FLAG;
      bndMask2d(i, model.bndMask->getLocalGhostNumOfCols() - 2) = DIRICHLET_FLAG;
    }
  }

  // TODO bmelt vs Q
  PETScGrid bmelt(NODATA_COLS, NODATA_ROWS);
  bmelt.setConst(1);
  model.setWaterSource(std::make_unique<CUAS::SteadyForcing>(bmelt, 1.0 / SPY));

  return pmodel;
}

#endif
