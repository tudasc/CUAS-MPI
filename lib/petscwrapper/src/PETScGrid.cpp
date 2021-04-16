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

  setGlobalBoundariesConst(boundaryValue);
}

void PETScGrid::setGlobalBoundariesConst(PetscScalar value) {
  // TODO use WriteHandleGhost
  // we write to ghost-cells here
  PetscScalar **local2d;
  VecGetArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &local2d);

  for (int i = 0; i < localGhostNumOfRows; ++i) {
    for (int j = 0; j < localGhostNumOfCols; ++j) {
      if (cornerYGhost == -1 && i == 0 || cornerXGhost == -1 && j == 0 ||
          cornerYGhost + localGhostNumOfRows == totalGhostNumOfRows - 1 && i == localGhostNumOfRows - 1 ||
          cornerXGhost + localGhostNumOfCols == totalGhostNumOfCols - 1 && j == localGhostNumOfCols - 1) {
        local2d[i][j] = value;
      }
    }
  }
  VecRestoreArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &local2d);
  DMLocalToGlobal(dm, local, INSERT_VALUES, global);
}

void PETScGrid::setInnerBoundariesConst(PetscScalar value) {
  // TODO use WriteHandle
  PetscScalar **global2d;
  VecGetArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &global2d);

  for (int i = 0; i < localNumOfRows; ++i) {
    for (int j = 0; j < localNumOfCols; ++j) {
      if (cornerYGhost == -1 && i == 0 || cornerXGhost == -1 && j == 0 ||
          cornerYGhost + localGhostNumOfRows == totalGhostNumOfRows - 1 && i == localNumOfRows - 1 ||
          cornerXGhost + localGhostNumOfCols == totalGhostNumOfCols - 1 && j == localNumOfCols - 1) {
        global2d[i][j] = value;
      }
    }
  }
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &global2d);
  DMGlobalToLocal(dm, global, INSERT_VALUES, local);
}

void PETScGrid::setGlobalVecColMajor(PETScVec &globalVec, bool ghosted) {
  // the following lines are copying the globalVec Vector to each processor.
  // could probably be optimized to get the specific needed values.

  if (ghosted) {
    if (globalVec.getSize() != totalGhostNumOfCols * totalGhostNumOfRows) {
      Logger::instance().error("PETScGrid.cpp: setGlobalVecColMajor(): sizes do not fit! Exiting.");
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
      Logger::instance().error("PETScGrid.cpp: setGlobalVecColMajor(): sizes do not fit! Exiting.");
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
}

void PETScGrid::setGlobalVecRowMajor(PETScVec &globalVec, bool ghosted) {
  // the following lines are copying the globalVec Vector to each processor.
  // could probably be optimized to get the specific needed values.

  if (ghosted) {
    if (globalVec.getSize() != totalGhostNumOfCols * totalGhostNumOfRows) {
      Logger::instance().error("PETScGrid.cpp: setGlobalVecRowMajor(): sizes do not fit! Exiting.");
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
      Logger::instance().error("PETScGrid.cpp: setGlobalVecRowMajor(): sizes do not fit! Exiting.");
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
}

PetscScalar **PETScGrid::getAsLocal2dArr() {
  PetscScalar **values;
  VecGetArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &values);
  return values;
}

void PETScGrid::restoreLocal2dArr(PetscScalar **values) {
  VecRestoreArray2d(local, localGhostNumOfRows, localGhostNumOfCols, 0, 0, &values);
}

PetscScalar **PETScGrid::getAsGlobal2dArr() {
  PetscScalar **values;
  VecGetArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
  return values;
}

void PETScGrid::restoreGlobal2dArr(PetscScalar **values) {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
}

void PETScGrid::setAsGlobal2dArr(PetscScalar **globalValues) {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &globalValues);
  DMGlobalToLocal(dm, global, INSERT_VALUES, local);
}

void PETScGrid::setConst(PetscScalar value) {
  VecSet(local, value);
  DMLocalToGlobal(dm, local, INSERT_VALUES, global);
}

void PETScGrid::copy(PETScGrid const &input) {
  if (!isCompatible(input)) {
    Logger::instance().error("PETScGrid.cpp: copy: input is not compatible. Exiting.");
    exit(1);
  }

  VecCopy(input.local, local);
  DMLocalToGlobal(dm, local, INSERT_VALUES, global);
}

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

PETScGrid::~PETScGrid() {
  VecRestoreArray2d(global, localNumOfRows, localNumOfCols, 0, 0, &values);
  DMRestoreLocalVector(dm, &local);
  DMRestoreGlobalVector(dm, &global);
  DMDestroy(&dm);
}
