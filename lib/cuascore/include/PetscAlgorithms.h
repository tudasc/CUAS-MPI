#ifndef PETSCALGORITHMS_H
#define PETSCALGORITHMS_H

#include "PetscGrid.h"

inline void binaryDialation(PetscGrid &no_flow_mask, PetscGrid &grad_mask) {
  auto nf_arr = no_flow_mask.getAsLocal2dArr();
  auto grad_arr = grad_mask.getAsGlobal2dArr();

  int index_cols = 1;
  int index_rows = 1;

  auto iter_rows = grad_mask.getLocalNumOfRows();
  auto iter_cols = grad_mask.getLocalNumOfCols();

  for (int i = 0; i < iter_rows; ++i) {
    for (int j = 0; j < iter_cols; ++j) {
      if (nf_arr[index_rows][index_cols] == true)
        grad_arr[i][j] = true;
      if (nf_arr[index_rows + 1][index_cols] == true)
        grad_arr[i][j] = true;
      if (nf_arr[index_rows - 1][index_cols] == true)
        grad_arr[i][j] = true;
      if (nf_arr[index_rows][index_cols + 1] == true)
        grad_arr[i][j] = true;
      if (nf_arr[index_rows][index_cols - 1] == true)
        grad_arr[i][j] = true;
      ++index_cols;
    }
    index_cols = 1;
    ++index_rows;
  }

  grad_mask.setAsGlobal2dArr(grad_arr);
  no_flow_mask.restoreLocal2dArr(nf_arr);

  return;
}

#endif
