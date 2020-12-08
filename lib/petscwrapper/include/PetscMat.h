#ifndef CUAS_PETSCMAT_H
#define CUAS_PETSCMAT_H

#include "petsc.h"

class PetscMat {
 private:
  Mat M;
  const int cols;
  const int rows;

 public:
  // constructor
  PetscMat(int numOfCols, int numOfRows);
  // getter
  Mat getPetscRaw() { return M; }
  int getCols() const { return cols; }
  int getRows() const { return rows; }
  // setter
  void setValue(int row, int col, PetscScalar val);
  // utils
  void assemble();

  ~PetscMat();
};

#endif
