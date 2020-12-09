#ifndef CUAS_PETSCMAT_H
#define CUAS_PETSCMAT_H

#include "petsc.h"

class PetscMat {
 private:
  Mat mat;
  const int cols;
  const int rows;

 public:
  // constructor
  PetscMat(int numOfRows, int numOfCols);
  // getter
  Mat getPetscRaw() { return mat; }
  int getCols() const { return cols; }
  int getRows() const { return rows; }
  // setter
  void setValue(int row, int col, PetscScalar val);
  // utils
  void assemble();

  ~PetscMat();

  friend class PetscSolver;
};

#endif
