#ifndef CUAS_PETSCSOLVER_H
#define CUAS_PETSCSOLVER_H

#include "Logger.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

#include "petsc.h"

class PETScSolver {
 public:
  inline static void solve(PETScMatrix const &A, PETScVector const &b, PETScVector &solution, bool verbose = false) {
    KSP ksp;
    PC pc;

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A.mat, A.mat);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetFromOptions(ksp);

    KSPSolve(ksp, b.vec, solution.vec);
    {
      KSPConvergedReason reason;
      KSPGetConvergedReason(ksp, &reason);
      if (reason < 0) {
        // diverged
        CUAS_ERROR_RANK0("PETScSolver: solve(): failed to converge (KSP reason {}).", KSPConvergedReasons[reason]);
        exit(1);
      } else {
        if (verbose) {
          // converged
          PetscInt its;     // number of iterations to solution
          PetscReal rnorm;  // residual norm
          KSPGetIterationNumber(ksp, &its);
          KSPGetResidualNorm(ksp, &rnorm);
          CUAS_INFO_RANK0("PETScSolver: solve(): converged (KSP reason {}) in {} iterations (residual norm {}).",
                          KSPConvergedReasons[reason], (int)its, (double)rnorm);
        }
      }
    }

    KSPDestroy(&ksp);
  }
};

#endif
