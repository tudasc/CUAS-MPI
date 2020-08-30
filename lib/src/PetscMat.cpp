#include "PetscMat.h"

#include "petsc.h"

#include <iostream>

PetscMat::PetscMat(int numOfCols, int numOfRows) : cols(numOfCols), rows(numOfRows) {
  MatCreate(PETSC_COMM_WORLD, &M);
  MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, numOfCols, numOfRows);
  MatSetFromOptions(M);
  MatSetUp(M);
}

void PetscMat::setValue(int row, int col, PetscScalar val) { MatSetValue(M, row, col, val, INSERT_VALUES); }

void PetscMat::assemble() {
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
}

void PetscMat::viewGlobal() const { MatView(M, PETSC_VIEWER_STDOUT_WORLD); }

void PetscMat::viewLocal() const { MatView(M, PETSC_VIEWER_STDOUT_SELF); }

PetscMat::~PetscMat() { MatDestroy(&M); }
