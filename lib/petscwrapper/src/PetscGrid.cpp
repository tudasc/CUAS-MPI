#include "PetscGrid.h"

#include <iostream>

PetscGrid::PetscGrid(int numOfCols, int numOfRows) : totalNumOfCols(numOfCols), totalNumOfRows(numOfRows) {
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, numOfCols, numOfRows,
               PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &dm);

  DMSetFromOptions(dm);
  DMSetUp(dm);

  DMCreateGlobalVector(dm, &global);
  DMCreateLocalVector(dm, &local);

  DMDAGetCorners(dm, &cornerX, &cornerY, NULL, &localNumOfCols, &localNumOfRows, NULL);
  DMDAGetGhostCorners(dm, &cornerXGhost, &cornerYGhost, NULL, &localGhostNumOfCols, &localGhostNumOfRows, NULL);
}

void PetscGrid::setGlobalVecAndUpdate(Vec global) {
  global = global;
  DMGlobalToLocalBegin(dm, global, INSERT_VALUES, local);
  DMGlobalToLocalEnd(dm, global, INSERT_VALUES, local);
}

PetscScalar **PetscGrid::getAsLocal2dArr() {
  PetscScalar **values;
  VecGetArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &values);
  return values;
}

void PetscGrid::restoreLocal2dArr(PetscScalar **values) {
  VecRestoreArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &values);
}

PetscScalar **PetscGrid::getAsGlobal2dArr() {
  PetscScalar **values;
  VecGetArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
  return values;
}

void PetscGrid::restoreGlobal2dArr(PetscScalar **values) {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
}

void PetscGrid::setAsGlobal2dArr(PetscScalar **globalValues) {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &globalValues);
  DMGlobalToLocal(dm, global, INSERT_VALUES, local);
}

void PetscGrid::setConst(PetscScalar value) {
  VecSet(local, value);
  DMLocalToGlobal(dm, local, INSERT_VALUES, global);
}

int PetscGrid::countNonZero() const {
  PetscScalar **gridArr2d;
  VecGetArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &gridArr2d);
  int count = 0;
  for (int j = 0; j < localNumOfRows; ++j) {
    for (int i = 0; i < localNumOfCols; ++i) {
      if (gridArr2d[j][i] != 0) {
        count++;
      }
    }
  }
  MPIU_Allreduce(&count, &count, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &gridArr2d);
  return count;
}

PetscGrid::~PetscGrid() {
  DMRestoreLocalVector(dm, &local);
  DMRestoreGlobalVector(dm, &global);
  DMDestroy(&dm);
}
