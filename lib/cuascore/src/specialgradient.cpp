#include "specialgradient.h"

#include <cmath>

namespace CUAS {

void gradient2(PetscGrid &gradient, PetscGrid const &h, PetscScalar const dx) {
  gradient.setZero();
  auto &localh = h.getReadHandle();  // local
  auto grad = gradient.getWriteHandle();
  for (int i = 0; i < gradient.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient.getLocalNumOfCols(); ++j) {
      int indexI = i + 1;
      int indexJ = j + 1;
      PetscScalar grad_x = 0.5 * (pow((localh(indexI, indexJ, GHOSTED) - localh(indexI - 1, indexJ, GHOSTED)) / dx, 2) +
                                  pow((localh(indexI + 1, indexJ, GHOSTED) - localh(indexI, indexJ, GHOSTED)) / dx, 2));
      PetscScalar grad_y = 0.5 * (pow((localh(indexI, indexJ, GHOSTED) - localh(indexI, indexJ - 1, GHOSTED)) / dx, 2) +
                                  pow((localh(indexI, indexJ + 1, GHOSTED) - localh(indexI, indexJ, GHOSTED)) / dx, 2));
      grad(i, j) = grad_x + grad_y;
    }
  }
}

void gradient2_central(PetscGrid &gradient, PetscGrid const &h, PetscScalar const dx) {
  gradient.setZero();
  auto &localh = h.getReadHandle();  // local
  auto grad = gradient.getWriteHandle();
  for (int i = 0; i < gradient.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient.getLocalNumOfCols(); ++j) {
      int indexI = i + 1;
      int indexJ = j + 1;
      PetscScalar grad_x = (localh(indexI + 1, indexJ, true) - localh(indexI - 1, indexJ, true)) / (2 * dx);
      PetscScalar grad_y = (localh(indexI, indexJ + 1, true) - localh(indexI, indexJ - 1, true)) / (2 * dx);
      grad(i, j) = pow(grad_x, 2) + pow(grad_y, 2);
    }
  }
}

}  // namespace CUAS
