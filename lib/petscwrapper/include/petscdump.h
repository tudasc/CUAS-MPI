#ifndef CUAS_PETSCDUMP_H
#define CUAS_PETSCDUMP_H

#include "PetscGrid.h"
#include "PetscMat.h"
#include "PetscVec.h"

#include "petsc.h"

#include <iostream>

inline void dump(PetscGrid &grid, bool showGhostCells = true) {
  if (showGhostCells) {
    PetscScalar **gridArr2d = grid.getAsLocal2dArr();

    int size, rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for (int proc = 0; proc < size; ++proc) {
      MPI_Barrier(PETSC_COMM_WORLD);
      if (rank == proc) {
        std::cout << "process num: " << rank << std::endl;
        for (int j = 0; j < grid.getLocalGhostNumOfRows(); ++j) {
          for (int i = 0; i < grid.getLocalGhostNumOfCols(); ++i) {
            std::cout << gridArr2d[j][i];
          }
          std::cout << std::endl;
        }
        std::cout << "------------" << std::endl;
      }
    }
    grid.restoreLocal2dArr(gridArr2d);
  } else {
    PetscScalar **gridArr2d = grid.getAsGlobal2dArr();
    int size, rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for (int proc = 0; proc < size; ++proc) {
      MPI_Barrier(PETSC_COMM_WORLD);
      if (rank == proc) {
        std::cout << "process num: " << rank << std::endl;
        for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
          for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
            std::cout << gridArr2d[j][i];
          }
          std::cout << std::endl;
        }
        std::cout << "------------" << std::endl;
      }
    }
    grid.restoreGlobal2dArr(gridArr2d);
  }
}

inline void dump(PetscMat &mat, bool showGlobal = true) {
  // mat.assemble();
  if (showGlobal) {
    MatView(mat.getPetscRaw(), PETSC_VIEWER_STDOUT_WORLD);
  } else {
    MatView(mat.getPetscRaw(), PETSC_VIEWER_STDOUT_SELF);
  }
}

inline void dump(PetscVec &vec) { VecView(vec.getPetscRaw(), PETSC_VIEWER_STDOUT_WORLD); }

#endif
