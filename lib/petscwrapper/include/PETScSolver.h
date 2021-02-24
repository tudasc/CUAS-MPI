#ifndef CUAS_PETSCSOLVER_H
#define CUAS_PETSCSOLVER_H

#include "PETScMat.h"
#include "PETScVec.h"

#include "petsc.h"

class PETScSolver {
 public:
  inline static void solve(PETScMat const &A, PETScVec const &b, PETScVec &solution) {
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
