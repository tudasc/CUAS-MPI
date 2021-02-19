#ifndef CUAS_PETSCSOLVER_H
#define CUAS_PETSCSOLVER_H

#include "PETScMat.h"
#include "PETScVec.h"

#include "petsc.h"

class PetscSolver {
 public:
  inline void solve(PetscMat const &A, PetscVec const &b, PetscVec &solution) {
    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    Mat mat_A = A.mat;
    KSPSetOperators(ksp, mat_A, mat_A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetFromOptions(ksp);
    Vec solVec = solution.vec;
    KSPSolve(ksp, b.vec, solVec);
  }
};

#endif