#include "CUASModel.h"

#include "CUASConstants.h"

namespace CUAS {

CUASModel::CUASModel(int numOfCols, int numOfRows) : Ncols(numOfCols), Nrows(numOfRows) {
  cols.resize(Ncols);
  rows.resize(Nrows);

  usurf = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  topg = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  thk = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  bndMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  Q = nullptr;
  pIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);
}

// call init before using Model but after thk has been set
void CUASModel::init() {
  dx = cols[1] - cols[0];
  dy = rows[1] - rows[0];
  if (dx != dy)
    exit(1);

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
