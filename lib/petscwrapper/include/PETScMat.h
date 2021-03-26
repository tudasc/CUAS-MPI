#ifndef CUAS_PETSCMAT_H
#define CUAS_PETSCMAT_H

#include "petsc.h"

/*
 * Matrices in PETSc are stored row major
 * access values[j * cols + i]
 * keep in mind that we ues Compressed Sparse Row (CSR) structure
 * therefore the following is only an abstraction
 *
 * +------------------>(cols, i, n)
 * | 0  1  2  3  4  5  6
 * | 7  8  9 10 11 12 13
 * |
 * v
 * (rows, j, m)
 *
 * +------------------>(cols, i, n)
 * | process 0
 * | process 1
 * | process 2
 * v
 * (rows, j, m)
 */

class PETScMat {
 public:
  // constructor
  explicit PETScMat(int numOfRows, int numOfCols) : rows(numOfRows), cols(numOfCols) {
    MatCreate(PETSC_COMM_WORLD, &mat);
    MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, numOfRows, numOfCols);
    MatSetFromOptions(mat);
    MatSetUp(mat);
  }
  PETScMat(PETScMat &) = delete;
  PETScMat(PETScMat &&) = delete;
  ~PETScMat() { MatDestroy(&mat); }

  // getter
  Mat getRaw() { return mat; }
  int getCols() const { return cols; }
  int getRows() const { return rows; }
  // setter
  void setValue(int row, int col, PetscScalar val) { MatSetValue(mat, row, col, val, INSERT_VALUES); }
  // utils
  void assemble() {
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  }

 private:
  Mat mat;
  const int cols;
  const int rows;

  friend class PETScSolver;
};

#endif
