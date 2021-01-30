#ifndef CUAS_PETSCDUMP_H
#define CUAS_PETSCDUMP_H

#include "PetscGrid.h"
#include "PetscMat.h"
#include "PetscVec.h"

#include "petsc.h"

#include <iostream>

inline void dump(PetscGrid const &grid, bool showGhostCells = true) {
  if (showGhostCells) {
    auto &gridArr2d = grid.getReadHandle();  // local

    int size, rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for (int proc = 0; proc < size; ++proc) {
      MPI_Barrier(PETSC_COMM_WORLD);
      if (rank == proc) {
        std::cout << "process num: " << rank << std::endl;
        for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
          for (int j = 0; j < grid.getLocalGhostNumOfCols(); ++j) {
            std::cout << gridArr2d(i, j, GHOSTED) << " ";
          }
          std::cout << std::endl;
        }
        std::cout << "------------" << std::endl;
      }
    }
  } else {
    auto &gridArr2d = grid.getReadHandle();
    int size, rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for (int proc = 0; proc < size; ++proc) {
      MPI_Barrier(PETSC_COMM_WORLD);
      if (rank == proc) {
        std::cout << "process num: " << rank << std::endl;
        for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
          for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
            std::cout << gridArr2d(i, j);
          }
          std::cout << std::endl;
        }
        std::cout << "------------" << std::endl;
      }
    }
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
