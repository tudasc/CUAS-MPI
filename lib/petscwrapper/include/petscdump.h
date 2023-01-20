/**
 * File: petscdump.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PETSCDUMP_H
#define CUAS_PETSCDUMP_H

#include "PETScGrid.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

#include "petsc.h"

#include <chrono>
#include <iostream>
#include <thread>

inline void dump(PETScGrid const &grid, bool showGhostCells = true) {
  if (showGhostCells) {
    auto &gridArr2d = grid.getReadHandle();  // local

    int size, rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for (int proc = 0; proc < size; ++proc) {
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
      MPI_Barrier(PETSC_COMM_WORLD);
      // sleep to finalize output
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  } else {
    auto &gridArr2d = grid.getReadHandle();
    int size, rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for (int proc = 0; proc < size; ++proc) {
      if (rank == proc) {
        std::cout << "process num: " << rank << std::endl;
        for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
          for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
            std::cout << gridArr2d(i, j) << " ";
          }
          std::cout << std::endl;
        }
        std::cout << "------------" << std::endl;
      }
      MPI_Barrier(PETSC_COMM_WORLD);
      // sleep to finalize output
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  }
}

inline void dump(PETScMatrix &mat, bool showGlobal = true) {
  // mat.assemble();
  if (showGlobal) {
    MatView(mat.getRaw(), PETSC_VIEWER_STDOUT_WORLD);
  } else {
    MatView(mat.getRaw(), PETSC_VIEWER_STDOUT_SELF);
  }
}

inline void dump(PETScVector &vec) { VecView(vec.getRaw(), PETSC_VIEWER_STDOUT_WORLD); }

#endif
