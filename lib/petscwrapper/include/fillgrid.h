#ifndef CUAS_FILLGRID_H
#define CUAS_FILLGRID_H

#include "PETScGrid.h"

inline void fillGlobalIndicesBlocked(PETScGrid &grid) {
  auto indexHandle = grid.getWriteHandle();
  auto startIndex = (grid.getCornerY() * grid.getTotalNumOfCols()) + (grid.getCornerX() * grid.getLocalNumOfRows());
  for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
      indexHandle(j, i) = startIndex;
      startIndex++;
    }
  }
}

#endif
