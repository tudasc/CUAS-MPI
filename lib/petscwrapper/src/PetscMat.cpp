#include "PetscMat.h"

#include <iostream>

PetscMat::PetscMat(int numOfRows, int numOfCols) : rows(numOfRows), cols(numOfCols) {
  MatCreate(PETSC_COMM_WORLD, &M);
  MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, numOfRows, numOfCols);
  MatSetFromOptions(M);
  MatSetUp(M);
}

void PetscMat::setValue(int row, int col, PetscScalar val) { MatSetValue(M, row, col, val, INSERT_VALUES); }

void PetscMat::assemble() {
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
}

PetscMat::~PetscMat() { MatDestroy(&M); }
