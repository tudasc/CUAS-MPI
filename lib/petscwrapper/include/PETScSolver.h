/**
 * File: PETScSolver.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PETSCSOLVER_H
#define CUAS_PETSCSOLVER_H

#include "Logger.h"
#include "PETScGrid.h"
#include "PETScMatrix.h"

#include "petsc.h"

struct PETScSolverConvergenceInformation {
  int numberOfIterations;  // to solution
  double residualNorm;
  std::string reason;  // the reason why the `KSP` iteration was stopped
                       // as in PETSC src/ksp/ksp/interface/dlregisksp.c
};

class PETScSolver {
 public:
  inline static PETScSolverConvergenceInformation solve(PETScMatrix const &A, PETScGrid const &bGrid,
                                                        PETScGrid &solGrid) {
    KSP ksp;
    PC pc;
    PETScSolverConvergenceInformation result = {0, 0.0, ""};

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A.mat, A.mat);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetFromOptions(ksp);

    KSPSolve(ksp, bGrid.global, solGrid.global);
    {
      KSPConvergedReason reason;
      PetscInt numberOfIterations;  // to solution
      PetscReal residualNorm;
      KSPGetConvergedReason(ksp, &reason);
      if (reason < 0) {
        // diverged
        CUAS_ERROR_RANK0("PETScSolver: solve(): failed to converge (KSP reason {}).", KSPConvergedReasons[reason])
        exit(1);
      }

      KSPGetIterationNumber(ksp, &numberOfIterations);
      KSPGetResidualNorm(ksp, &residualNorm);

      // fill in return values
      result.numberOfIterations = static_cast<int>(numberOfIterations);
      result.residualNorm = static_cast<double>(residualNorm);
      result.reason = std::string(KSPConvergedReasons[reason]);
    }

    KSPDestroy(&ksp);
    return result;
  }

  /** Solve with MUMPS direct solver
   *
   * Some options are taken from the ISSM toolkits file for MUMPS
   *    -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps
   *    -mat_mumps_icntl_14 120 -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2
   * For more options see https://petsc.org/main/docs/manualpages/Mat/MATSOLVERMUMPS.html
   */
  inline static PETScSolverConvergenceInformation solveDirectMUMPS(PETScMatrix const &A, PETScGrid const &bGrid,
                                                                   PETScGrid &solGrid) {
    KSP ksp;
    PC pc;
    PETScSolverConvergenceInformation result = {1, 0.0, ""};

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

    // solve the linear system, number of iteration will be 1 and residualNorm = 0.0
    KSPSolve(ksp, bGrid.global, solGrid.global);
    {
      KSPConvergedReason reason;
      KSPGetConvergedReason(ksp, &reason);
      if (reason < 0) {
        // diverged
        CUAS_ERROR("PETScSolver: solve(): failed to converge (KSP reason {}).", KSPConvergedReasons[reason])
        exit(1);
      }
      result.reason = KSPConvergedReasons[reason];
    }

    KSPDestroy(&ksp);
    return result;
  }
};

#endif
