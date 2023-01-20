/**
 * File: PETScMatrix.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PETSCMAT_H
#define CUAS_PETSCMAT_H

#include "petsc.h"

/*
 * Matrices in PETSc are stored row major
 * access values[j * xAxis + i]
 * keep in mind that we ues Compressed Sparse Row (CSR) structure
 * therefore the following is only an abstraction
 *
 * +------------------>(cols, xAxis, i, n)
 * | 0  1  2  3  4  5  6
 * | 7  8  9 10 11 12 13
 * |
 * v
 * (rows, yAxis, j, m)
 *
 * +------------------>(cols, xAxis, i, n)
 * | process 0
 * | process 1
 * | process 2
 * v
 * (rows, yAxis, j, m)
 */

class PETScMatrix {
 public:
  // constructor
  explicit PETScMatrix(int numOfRows, int numOfCols) : nRows(numOfRows), nCols(numOfCols) {
    MatCreate(PETSC_COMM_WORLD, &mat);
    MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, numOfRows, numOfCols);
    MatSetFromOptions(mat);
    MatSetUp(mat);
  }
  explicit PETScMatrix(Mat mat) : mat(mat) { MatGetSize(mat, &nRows, &nCols); }
  PETScMatrix(PETScMatrix &) = delete;
  PETScMatrix(PETScMatrix &&) = delete;
  ~PETScMatrix() { MatDestroy(&mat); }

  // getter
  Mat getRaw() { return mat; }
  int getNumberOfCols() const { return nCols; }
  int getNumberOfRows() const { return nRows; }
  // setter
  void setValue(int row, int col, PetscScalar val) { MatSetValue(mat, row, col, val, INSERT_VALUES); }
  void setZero() { MatZeroEntries(mat); }
  // utils
  void assemble() {
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  }

 private:
  Mat mat;
  int nCols;
  int nRows;

  friend class PETScSolver;
};

#endif
