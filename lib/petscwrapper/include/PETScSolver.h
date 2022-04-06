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

  /** Solve with MUMPS direct solver
   *
   * Some options are taken from the ISSM toolkits file for MUMPS
   *    -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps
   *    -mat_mumps_icntl_14 120 -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2
   * For more options see https://petsc.org/main/docs/manualpages/Mat/MATSOLVERMUMPS.html
   */
  inline static void solveDirectMUMPS(PETScMatrix const &A, PETScVector const &b, PETScVector &solution,
                                      bool verbose = false) {
    KSP ksp;
    PC pc;

    KSPCreate(PETSC_COMM_WORLD, &ksp);   // creates the default KSP context
    KSPSetOperators(ksp, A.mat, A.mat);  // PETSc >= 3.5.0

    // KSPPREONLY applies only the preconditioner as MUMPS being a direct solver
    KSPSetType(ksp, KSPPREONLY);

    // PREONLY solver need zero initial guess
    KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

    // getting the associated preconditioner context
    KSPGetPC(ksp, &pc);

    // set preconditioner (PC):  PCJACOBI, PCCHOLESKY (if symmetric), PCLU (if non-symmetric)
    PCSetType(pc, PCLU);  // fixme: PCCHOLESKY should work, but fails for 2D Poisson solver. Why?

    // tell the solver to perform the LU factorization and get the
    PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);  // PETSc >= 3.9.0

    PetscOptionsSetValue(PETSC_NULL, "-mat_mumps_icntl_14", "120");  // needed?
    PetscOptionsSetValue(PETSC_NULL, "-mat_mumps_icntl_28", "2");    // parallel analysis
    PetscOptionsSetValue(PETSC_NULL, "-mat_mumps_icntl_29", "2");    // parallel ordering: 2 = parmetis

    // set the command line options provided by the user to override the defaults
    KSPSetFromOptions(ksp);

    // solve the linear system, number of iteration will be 1 and rnorm = 0.0
    KSPSolve(ksp, b.vec, solution.vec);
    {
      KSPConvergedReason reason;
      KSPGetConvergedReason(ksp, &reason);
      if (reason < 0) {
        // diverged
        CUAS_ERROR("PETScSolver: solve(): failed to converge (KSP reason {}).", KSPConvergedReasons[reason]);
        exit(1);
      }
    }

    KSPDestroy(&ksp);
  }
};

#endif
