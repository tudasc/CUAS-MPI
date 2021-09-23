#ifndef CUAS_PETSCSOLVER_H
#define CUAS_PETSCSOLVER_H

#include "PETScMatrix.h"
#include "PETScVector.h"

#include "petsc.h"

class PETScSolver {
 public:
  inline static void solve(PETScMatrix const &A, PETScVector const &b, PETScVector &solution) {
    KSP ksp;
    PC pc;

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A.mat, A.mat);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, b.vec, solution.vec);

    KSPDestroy(&ksp);
  }
};

#endif
