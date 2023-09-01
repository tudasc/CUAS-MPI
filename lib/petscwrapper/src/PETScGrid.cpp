/**
 * File: PETScGrid.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "PETScGrid.h"

PETScGrid::PETScGrid(int numOfCols, int numOfRows, PetscScalar boundaryValue)
    : totalNumOfCols(numOfCols),
      totalNumOfRows(numOfRows),
      totalGhostNumOfCols(numOfCols + 2),
      totalGhostNumOfRows(numOfRows + 2),
      readHandle(this) {
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, totalNumOfCols,
               totalNumOfRows, PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr, &dm);

  DMSetFromOptions(dm);
  DMSetUp(dm);

  DMCreateGlobalVector(dm, &global);
  DMCreateLocalVector(dm, &local);

  DMDAGetCorners(dm, &cornerX, &cornerY, nullptr, &localNumOfCols, &localNumOfRows, nullptr);
  DMDAGetGhostCorners(dm, &cornerXGhost, &cornerYGhost, nullptr, &localGhostNumOfCols, &localGhostNumOfRows, nullptr);

  // create values-Pointer to access using the handle
  VecGetArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
  VecGetArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &valuesGhosted);

  setGhostBoundary(boundaryValue);
}

PETScGrid::~PETScGrid() {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
  VecRestoreArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &valuesGhosted);
  VecDestroy(&local);
  VecDestroy(&global);
  DMDestroy(&dm);
}

PetscScalar PETScGrid::getMaxAbsDiff(PETScGrid const &sub) const {
  if (!isCompatible(sub)) {
    CUAS_ERROR("PETScGrid.cpp: copy: input is not compatible. Exiting.")
    exit(1);
  }
  PetscScalar result = 0.0;
  Vec diff;
  DMGetGlobalVector(dm, &diff);
  VecCopy(global, diff);
  VecAXPY(diff, -1, sub.global);
  VecNorm(diff, NORM_INFINITY, &result);
  DMRestoreGlobalVector(dm, &diff);
  return result;
}

std::array<PetscScalar, 3> PETScGrid::getErrorNorms(PETScGrid const &sub) const {
  if (!isCompatible(sub)) {
    CUAS_ERROR("{}: input is not compatible. Exiting.", __PRETTY_FUNCTION__);
    exit(1);
  }
  PetscScalar L1, L2, Linf;
  Vec diff;
  DMCreateGlobalVector(dm, &diff);
  VecCopy(global, diff);
  VecAXPY(diff, -1, sub.global);
  VecNorm(diff, NORM_1, &L1);
  VecNorm(diff, NORM_2, &L2);
  VecNorm(diff, NORM_INFINITY, &Linf);
  return {L1, L2, Linf};
}

PetscScalar PETScGrid::getMax() const {
  PetscScalar result;
  VecMax(global, PETSC_NULL, &result);
  return result;
}

void PETScGrid::setConst(PetscScalar value) {
  VecSet(local, value);
  DMLocalToGlobal(dm, local, INSERT_VALUES, global);
}

void PETScGrid::copy(PETScGrid const &input) {
  if (!isCompatible(input)) {
    CUAS_ERROR("PETScGrid.cpp: copy: input is not compatible. Exiting.")
    exit(1);
  }

  VecCopy(input.local, local);
  DMLocalToGlobal(dm, local, INSERT_VALUES, global);
}

void PETScGrid::copyGlobal(PETScGrid const &input) {
  if (!isCompatible(input)) {
    CUAS_ERROR("PETScGrid.cpp: copy: input is not compatible. Exiting.")
    exit(1);
  }

  VecCopy(input.global, global);
  DMGlobalToLocal(dm, global, INSERT_VALUES, local);
}

bool PETScGrid::isOnGhostBoundary(int row, int col) const {
  return cornerYGhost == -1 && row == 0 || cornerXGhost == -1 && col == 0 ||
         cornerYGhost + localGhostNumOfRows == totalGhostNumOfRows - 1 && row == localGhostNumOfRows - 1 ||
         cornerXGhost + localGhostNumOfCols == totalGhostNumOfCols - 1 && col == localGhostNumOfCols - 1;
}

bool PETScGrid::isOnRealBoundary(int row, int col) const {
  return cornerYGhost == -1 && row == 0 || cornerXGhost == -1 && col == 0 ||
         cornerYGhost + localGhostNumOfRows == totalGhostNumOfRows - 1 && row == localNumOfRows - 1 ||
         cornerXGhost + localGhostNumOfCols == totalGhostNumOfCols - 1 && col == localNumOfCols - 1;
}

void PETScGrid::setGhostBoundary(PetscScalar value) {
  auto handle = getWriteHandleGhost();
  for (int row = 0; row < localGhostNumOfRows; ++row) {
    for (int col = 0; col < localGhostNumOfCols; ++col) {
      if (isOnGhostBoundary(row, col)) {
        handle(row, col) = value;
      }
    }
  }
}

void PETScGrid::setRealBoundary(PetscScalar value) {
  auto handle = getWriteHandle();
  for (int row = 0; row < localNumOfRows; ++row) {
    for (int col = 0; col < localNumOfCols; ++col) {
      if (isOnRealBoundary(row, col)) {
        handle(row, col) = value;
      }
    }
  }
}

void PETScGrid::findAndReplaceGhostBoundary(PetscScalar oldValue, PetscScalar newValue) {
  auto handle = getWriteHandleGhost();
  for (int row = 0; row < localGhostNumOfRows; ++row) {
    for (int col = 0; col < localGhostNumOfCols; ++col) {
      if (isOnGhostBoundary(row, col)) {
        if (handle(row, col) == oldValue) {
          handle(row, col) = newValue;
        }
      }
    }
  }
}

void PETScGrid::findAndReplaceRealBoundary(PetscScalar oldValue, PetscScalar newValue) {
  auto handle = getWriteHandle();
  for (int row = 0; row < localNumOfRows; ++row) {
    for (int col = 0; col < localNumOfCols; ++col) {
      if (isOnRealBoundary(row, col)) {
        if (handle(row, col) == oldValue) {
          handle(row, col) = newValue;
        }
      }
    }
  }
}

/*void PETScGrid::setGlobalVecColMajor(PETScVector &globalVec, bool ghosted) {
  // the following lines are copying the globalVec Vector to each processor.
  // could probably be optimized to get the specific needed values.

  if (ghosted) {
    if (globalVec.getSize() != totalGhostNumOfCols * totalGhostNumOfRows) {
      CUAS_ERROR("PETScGrid.cpp: setGlobalVecColMajor(): sizes do not fit! Exiting.")
      exit(1);
    }

    VecScatter vs;
    Vec accessVector;
    VecScatterCreateToAll(globalVec.getRaw(), &vs, &accessVector);
    VecScatterBegin(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar *vecArr;
    VecGetArrayRead(accessVector, &vecArr);

    auto global2d = getWriteHandleGhost();
    for (int i = 0; i < localGhostNumOfRows; ++i) {
      int indexI = getCornerY() + i;
      for (int j = 0; j < localGhostNumOfCols; ++j) {
        int indexJ = getCornerX() + j;
        global2d(i, j) = vecArr[totalGhostNumOfRows * indexJ + indexI];
      }
    }

    VecRestoreArrayRead(accessVector, &vecArr);
    VecDestroy(&accessVector);
    VecScatterDestroy(&vs);

  } else {
    if (globalVec.getSize() != totalNumOfCols * totalNumOfRows) {
      CUAS_ERROR("PETScGrid.cpp: setGlobalVecColMajor(): sizes do not fit! Exiting.")
      exit(1);
    }

    VecScatter vs;
    Vec accessVector;
    VecScatterCreateToAll(globalVec.getRaw(), &vs, &accessVector);
    VecScatterBegin(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar *vecArr;
    VecGetArrayRead(accessVector, &vecArr);

    auto global2d = getWriteHandle();
    for (int i = 0; i < localNumOfRows; ++i) {
      int indexI = getCornerY() + i;
      for (int j = 0; j < localNumOfCols; ++j) {
        int indexJ = getCornerX() + j;
        global2d(i, j) = vecArr[totalNumOfRows * indexJ + indexI];
      }
    }

    VecRestoreArrayRead(accessVector, &vecArr);
    VecDestroy(&accessVector);
    VecScatterDestroy(&vs);
  }
}*/

/*void PETScGrid::setGlobalVecRowMajor(PETScVector &globalVec, bool ghosted) {
  // the following lines are copying the globalVec Vector to each processor.
  // could probably be optimized to get the specific needed values.

  if (ghosted) {
    if (globalVec.getSize() != totalGhostNumOfCols * totalGhostNumOfRows) {
      CUAS_ERROR("PETScGrid.cpp: setGlobalVecRowMajor(): sizes do not fit! Exiting.")
      exit(1);
    }

    VecScatter vs;
    Vec accessVector;
    VecScatterCreateToAll(globalVec.getRaw(), &vs, &accessVector);
    VecScatterBegin(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar *vecArr;
    VecGetArrayRead(accessVector, &vecArr);

    auto global2d = getWriteHandleGhost();
    for (int i = 0; i < localGhostNumOfRows; ++i) {
      int indexI = getCornerY() + i;
      for (int j = 0; j < localGhostNumOfCols; ++j) {
        int indexJ = getCornerX() + j;
        global2d(i, j) = vecArr[totalGhostNumOfCols * indexI + indexJ];
      }
    }

    VecRestoreArrayRead(accessVector, &vecArr);
    VecDestroy(&accessVector);
    VecScatterDestroy(&vs);

  } else {
    if (globalVec.getSize() != totalNumOfCols * totalNumOfRows) {
      CUAS_ERROR("PETScGrid.cpp: setGlobalVecRowMajor(): sizes do not fit! Exiting.")
      exit(1);
    }

    VecScatter vs;
    Vec accessVector;
    VecScatterCreateToAll(globalVec.getRaw(), &vs, &accessVector);
    VecScatterBegin(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(vs, globalVec.getRaw(), accessVector, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar *vecArr;
    VecGetArrayRead(accessVector, &vecArr);

    auto global2d = getWriteHandle();
    for (int i = 0; i < localNumOfRows; ++i) {
      int indexI = getCornerY() + i;
      for (int j = 0; j < localNumOfCols; ++j) {
        int indexJ = getCornerX() + j;
        global2d(i, j) = vecArr[totalNumOfCols * indexI + indexJ];
      }
    }

    VecRestoreArrayRead(accessVector, &vecArr);
    VecDestroy(&accessVector);
    VecScatterDestroy(&vs);
  }
}*/

/*PetscScalar **PETScGrid::getAsLocal2dArr() {
  PetscScalar **values;
  VecGetArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &values);
  return values;
}*/

/*void PETScGrid::restoreLocal2dArr(PetscScalar **values) {
  VecRestoreArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &values);
}*/

/*PetscScalar **PETScGrid::getAsGlobal2dArr() {
  PetscScalar **values;
  VecGetArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
  return values;
}*/

/*void PETScGrid::restoreGlobal2dArr(PetscScalar **values) {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
}*/

/*void PETScGrid::setAsGlobal2dArr(PetscScalar **globalValues) {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &globalValues);
  DMGlobalToLocal(dm, global, INSERT_VALUES, local);
}*/

// int PETScGrid::countNonZero() const {
//   PetscScalar **gridArr2d;
//   VecGetArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &gridArr2d);
//   int count = 0;
//   for (int j = 0; j < localNumOfRows; ++j) {
//     for (int i = 0; i < localNumOfCols; ++i) {
//       if (gridArr2d[j][i] != 0) {
//         count++;
//       }
//     }
//   }
//   MPIU_Allreduce(&count, &count, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
//   VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &gridArr2d);
//   return count;
// }
