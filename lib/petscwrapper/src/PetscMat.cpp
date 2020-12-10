#include "PetscMat.h"

PetscMat::PetscMat(int numOfRows, int numOfCols) : rows(numOfRows), cols(numOfCols) {
  MatCreate(PETSC_COMM_WORLD, &mat);
  MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, numOfRows, numOfCols);
  MatSetFromOptions(mat);
  MatSetUp(mat);
}

void PetscMat::setValue(int row, int col, PetscScalar val) { MatSetValue(mat, row, col, val, INSERT_VALUES); }

void PetscMat::assemble() {
  MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
}

PetscMat::~PetscMat() { MatDestroy(&mat); }
