#ifndef INITIALHEADTEST_H
#define INITIALHEADTEST_H

#include "CUASConstants.h"
#include "CUASModel.h"

#ifndef GRID_SIZE_X
#define GRID_SIZE_X 9
#endif

#ifndef GRID_SIZE_Y
#define GRID_SIZE_Y 5
#endif

#ifndef THK
#define THK 1000.0
#endif

#ifndef TOPG
#define TOPG 10.0
#endif

std::unique_ptr<CUAS::CUASModel> fillData() {
  auto pmodel = std::make_unique<CUAS::CUASModel>(GRID_SIZE_X, GRID_SIZE_Y);
  auto &model = *pmodel;

  // initialize x and y like in python: x = np.arange(20) * 1000.0
  for (int i = 0; i < model.xAxis.size(); ++i) {
    model.xAxis[i] = i * 1000;
  }
  for (int i = 0; i < model.yAxis.size(); ++i) {
    model.yAxis[i] = i * 1000;
  }

  model.topg->setConst(TOPG);
  model.thk->setConst(THK);
  model.bndMask->setConst(DIRICHLET_FLAG);  // all dirichlet --> nothing to solve for

  return pmodel;
}

#endif
