#include "specialgradient.h"

#include <cmath>

namespace CUAS {

void gradient2(PetscGrid *const h, PetscGrid *const gradient, PetscScalar const dx) {
  PetscScalar grad_arr_x, grad_arr_y;
  int cols = h->getTotalNumOfCols();
  int rows = h->getTotalNumOfRows();
  PetscScalar **grad_arr = gradient->getAsGlobal2dArr();
  PetscScalar **h_arr = h->getAsGlobal2dArr();
  // write content of grid h to the gradient grid
  for (int i = 0; i < gradient->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient->getLocalNumOfCols(); ++j) {
      grad_arr[i][j] = h_arr[i][j];
    }
  }
  PetscScalar **localh = h->getAsLocal2dArr();
  // chose startindex for the 2d-arr of the gradient. The logic is needed because the starting-indices in the final
  // gradient differ on the edges
  int index_cols;
  int index_rows;
  (h->getCornerY() == 0) ? index_rows = 1 : index_rows = 0;
  (h->getCornerX() == 0) ? index_cols = 1 : index_cols = 0;
  int old_index_cols = index_cols;
  int iter_rows = h->getLocalGhostNumOfRows() - 1;
  int iter_cols = h->getLocalGhostNumOfCols() - 1;
  for (int i = 1; i < iter_rows; ++i) {
    for (int j = 1; j < iter_cols; ++j) {
      grad_arr_x =
          (0.5 * (pow(((localh[i][j] - localh[i - 1][j]) / dx), 2) + pow(((localh[i + 1][j] - localh[i][j]) / dx), 2)));
      grad_arr_y =
          (0.5 * (pow(((localh[i][j] - localh[i][j - 1]) / dx), 2) + pow(((localh[i][j + 1] - localh[i][j]) / dx), 2)));
      // write the calculated value to the gradient
      grad_arr[index_rows][index_cols] = grad_arr_x + grad_arr_y;
      ++index_cols;
    }
    index_cols = old_index_cols;
    ++index_rows;
  }
  // restore
  h->restoreLocal2dArr(localh);
  gradient->setAsGlobal2dArr(grad_arr);
}

void gradient2_central(PetscGrid *const h, PetscGrid *const gradient, PetscScalar const dx) {
  // apart from the calculation of grad_arr_x/y this function does the same as gradient2(...);
  PetscScalar grad_arr_x, grad_arr_y;
  int cols = h->getTotalNumOfCols();
  int rows = h->getTotalNumOfRows();
  PetscScalar **grad_arr = gradient->getAsGlobal2dArr();
  PetscScalar **h_arr = h->getAsGlobal2dArr();
  for (int i = 0; i < gradient->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < gradient->getLocalNumOfCols(); ++j) {
      grad_arr[i][j] = h_arr[i][j];
    }
  }
  PetscScalar **localh = h->getAsLocal2dArr();
  int index_cols;
  int index_rows;
  (h->getCornerY() == 0) ? index_rows = 1 : index_rows = 0;
  (h->getCornerX() == 0) ? index_cols = 1 : index_cols = 0;
  int old_index_cols = index_cols;
  int iter_rows = h->getLocalGhostNumOfRows() - 1;
  int iter_cols = h->getLocalGhostNumOfCols() - 1;
  for (int i = 1; i < iter_rows; ++i) {
    for (int j = 1; j < iter_cols; ++j) {
      grad_arr_x = (localh[i + 1][j] - localh[i - 1][j]) / (2 * dx);
      grad_arr_y = (localh[i][j + 1] - localh[i][j - 1]) / (2 * dx);
      grad_arr_x = pow(grad_arr_x, 2);
      grad_arr_y = pow(grad_arr_y, 2);
      grad_arr[index_rows][index_cols] = grad_arr_x + grad_arr_y;
      ++index_cols;
    }
    index_cols = old_index_cols;
    ++index_rows;
  }
  // restore
  h->restoreLocal2dArr(localh);
  gradient->setAsGlobal2dArr(grad_arr);
}

}  // namespace CUAS
