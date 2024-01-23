/**
 * File: timing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_TIMING_H
#define CUAS_TIMING_H

#include "mpi.h"

clock_t beginTime;

void beginSolverTiming() {
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // start
  if (rank == 0) {
    beginTime = clock();
  }
}

void endSolverTiming() {
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (rank == 0) {
    auto elapsedTime = clock() - beginTime;
    CUAS_INFO_RANK0("CUASSolver.cpp: solve(): computation took: {} seconds.", ((float)elapsedTime) / CLOCKS_PER_SEC)
  }
}

#endif  // CUAS_TIMING_H
