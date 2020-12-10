#include "specialgradient.h"

#include <cmath>

namespace CUAS {

void gradient2(PetscGrid &h, PetscGrid &gradient, PetscScalar const dx) {
  gradient.setZero();
  PetscScalar **localh = h.getAsLocal2dArr();
  PetscScalar **grad = gradient.getAsGlobal2dArr();
  for (int i = 0; i < gradient.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient.getLocalNumOfCols(); ++j) {
      int indexI = i + 1;
      int indexJ = j + 1;
      double grad_x = 0.5 * (pow((localh[indexI][indexJ] - localh[indexI - 1][indexJ]) / dx, 2) +
                             pow((localh[indexI + 1][indexJ] - localh[indexI][indexJ]) / dx, 2));
      double grad_y = 0.5 * (pow((localh[indexI][indexJ] - localh[indexI][indexJ - 1]) / dx, 2) +
                             pow((localh[indexI][indexJ + 1] - localh[indexI][indexJ]) / dx, 2));
      grad[i][j] = grad_x + grad_y;
    }
  }
  gradient.setAsGlobal2dArr(grad);
}

void gradient2_central(PetscGrid &h, PetscGrid &gradient, PetscScalar const dx) {
  gradient.setZero();
  PetscScalar **localh = h.getAsLocal2dArr();
  PetscScalar **grad = gradient.getAsGlobal2dArr();
  for (int i = 0; i < gradient.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient.getLocalNumOfCols(); ++j) {
      int indexI = i + 1;
      int indexJ = j + 1;
      double grad_x = (localh[indexI + 1][indexJ] - localh[indexI - 1][indexJ]) / (2 * dx);
      double grad_y = (localh[indexI][indexJ + 1] - localh[indexI][indexJ - 1]) / (2 * dx);
      grad[i][j] = pow(grad_x, 2) + pow(grad_y, 2);
    }
  }
  gradient.setAsGlobal2dArr(grad);
}

}  // namespace CUAS
