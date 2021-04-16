#include "specialgradient.h"

#include "Logger.h"

#include <cmath>

namespace CUAS {

void gradient2(PETScGrid &gradient, PETScGrid const &input, PetscScalar const dx) {
  if (!gradient.isCompatible(input) || dx == 0) {
    Logger::instance().error("specialgradient.cpp: gradient2(): PETScGrids are not compatible or dx is 0. Exiting.");
    exit(1);
  }
  auto dx_inv = 1.0 / dx;
  auto &localh = input.getReadHandle();  // local
  auto grad = gradient.getWriteHandle();
  for (int i = 0; i < gradient.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient.getLocalNumOfCols(); ++j) {
      int indexI = i + 1;
      int indexJ = j + 1;
      PetscScalar squaredDiffX =
          pow((localh(indexI, indexJ, GHOSTED) - localh(indexI - 1, indexJ, GHOSTED)) * dx_inv, 2) +
          pow((localh(indexI + 1, indexJ, GHOSTED) - localh(indexI, indexJ, GHOSTED)) * dx_inv, 2);
      PetscScalar squaredDiffY =
          pow((localh(indexI, indexJ, GHOSTED) - localh(indexI, indexJ - 1, GHOSTED)) * dx_inv, 2) +
          pow((localh(indexI, indexJ + 1, GHOSTED) - localh(indexI, indexJ, GHOSTED)) * dx_inv, 2);
      grad(i, j) = 0.5 * (squaredDiffX + squaredDiffY);
    }
  }
}

void gradient2_central(PETScGrid &gradient, PETScGrid const &input, PetscScalar const dx) {
  if (!gradient.isCompatible(input) || dx == 0) {
    Logger::instance().error(
        "specialgradient.cpp: gradient2_central(): PETScGrids are not compatible or dx is 0. Exiting.");
    exit(1);
  }
  auto &localh = input.getReadHandle();  // local
  auto grad = gradient.getWriteHandle();
  for (int i = 0; i < gradient.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient.getLocalNumOfCols(); ++j) {
      int indexI = i + 1;
      int indexJ = j + 1;
      PetscScalar grad_x = (localh(indexI + 1, indexJ, GHOSTED) - localh(indexI - 1, indexJ, GHOSTED)) / (2 * dx);
      PetscScalar grad_y = (localh(indexI, indexJ + 1, GHOSTED) - localh(indexI, indexJ - 1, GHOSTED)) / (2 * dx);
      grad(i, j) = pow(grad_x, 2) + pow(grad_y, 2);
    }
  }
}

}  // namespace CUAS
