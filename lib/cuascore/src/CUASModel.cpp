#include "CUASModel.h"

#include "CUASConstants.h"

#include "Logger.h"

namespace CUAS {

CUASModel::CUASModel(int numOfCols, int numOfRows) : Ncols(numOfCols), Nrows(numOfRows) {
  xAxis.resize(Ncols);
  yAxis.resize(Nrows);

  topg = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  thk = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  bndMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  Q = nullptr;
  pIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);
}

// call init before using Model but after thk has been set
void CUASModel::init() {
  dx = xAxis[1] - xAxis[0];
  dy = yAxis[1] - yAxis[0];
  if (dx != dy) {
    CUAS_ERROR("CUASModel.cpp: init(): dx and dy are not equal. Exiting.");
    exit(1);
  }

  // the check if dy == 0 is just for readability. The previous implies that if dx == 0 -> dy == 0
  if (dx == 0 || dy == 0) {
    CUAS_ERROR("CUASModel.cpp: init(): dx and dy are zero. Exiting.");
    exit(1);
  }

  for (int i = 0; i < xAxis.size() - 1; ++i) {
    if (std::fabs(xAxis[i + 1] - xAxis[i] - dx) > std::numeric_limits<double>::epsilon()) {
      CUAS_ERROR("CUASModel.cpp: init(): The values of xAxis are not evenly spaced across all time steps. Exiting.");
      exit(1);
    }
  }

  for (int i = 0; i < yAxis.size() - 1; ++i) {
    if (std::fabs(yAxis[i + 1] - yAxis[i] - dy) > std::numeric_limits<double>::epsilon()) {
      CUAS_ERROR("CUASModel.cpp: init(): The values of yAxis are not evenly spaced across all time steps. Exiting.");
      exit(1);
    }
  }

  // pIce = thk * RHO_ICE * GRAVITY (python)
  auto &thk2d = thk->getReadHandle();
  auto pIce2d = pIce->getWriteHandleGhost();
  for (int j = 0; j < thk->getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < thk->getLocalGhostNumOfCols(); ++i) {
      pIce2d(j, i) = thk2d(j, i, GHOSTED) * RHO_ICE * GRAVITY;
    }
  }
}

}  // namespace CUAS
