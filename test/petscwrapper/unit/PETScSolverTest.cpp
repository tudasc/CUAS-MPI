/**
 * File: PETScSolverTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "PETScSolver.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 4
#define VEC_SIZE 4
#define VERBOSE true

/*TEST(PetscSolverTest, solve) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto A = std::make_unique<PETScMatrix>(VEC_SIZE, VEC_SIZE);
  auto b = std::make_unique<PETScVector>(VEC_SIZE);
  auto s = std::make_unique<PETScVector>(VEC_SIZE);

  A->setValue(0, 0, 2);
  A->setValue(1, 1, 2);
  A->setValue(2, 2, 2);
  A->setValue(3, 3, 2);
  A->assemble();

  b->setValue(mpiRank, mpiRank);
  b->assemble();

  PETScSolver::solve(*A, *b, *s, VERBOSE);

  PetscScalar v;
  int p[] = {mpiRank};
  VecGetValues(s->getRaw(), 1, p, &v);
  ASSERT_EQ(v, mpiRank / 2.0);
}

#if (PETSC_HAVE_MUMPS == 1) && (PETSC_HAVE_PARMETIS == 1)
TEST(PetscSolverTest, solveMUMPS) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto A = std::make_unique<PETScMatrix>(VEC_SIZE, VEC_SIZE);
  auto b = std::make_unique<PETScVector>(VEC_SIZE);
  auto s = std::make_unique<PETScVector>(VEC_SIZE);

  A->setValue(0, 0, 2);
  A->setValue(1, 1, 2);
  A->setValue(2, 2, 2);
  A->setValue(3, 3, 2);
  A->assemble();

  b->setValue(mpiRank, mpiRank);
  b->assemble();

  PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_type", "preonly");
  PetscOptionsSetValue(PETSC_NULLPTR, "-pc_type", "lu");
  PetscOptionsSetValue(PETSC_NULLPTR, "-pc_factor_mat_solver_type", "mumps");
  PetscOptionsSetValue(PETSC_NULLPTR, "-mat_mumps_icntl_14", "120");
  PetscOptionsSetValue(PETSC_NULLPTR, "-mat_mumps_icntl_28", "2");
  PetscOptionsSetValue(PETSC_NULLPTR, "-mat_mumps_icntl_29", "2");
  PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_error_if_not_converged", nullptr);
  PETScSolver::solve(*A, *b, *s, VERBOSE);

  PetscScalar v;
  int p[] = {mpiRank};
  VecGetValues(s->getRaw(), 1, p, &v);
  ASSERT_EQ(v, mpiRank / 2.0);
}
#endif

#if (PETSC_HAVE_MUMPS == 1) && (PETSC_HAVE_PARMETIS == 1)
// This is used to compare with a regular PETScSolver::solve()
// with PETSc args set from options
TEST(PetscSolverTest, solveDirectMUMPS) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto A = std::make_unique<PETScMatrix>(VEC_SIZE, VEC_SIZE);
  auto b = std::make_unique<PETScVector>(VEC_SIZE);
  auto s = std::make_unique<PETScVector>(VEC_SIZE);

  A->setValue(0, 0, 2);
  A->setValue(1, 1, 2);
  A->setValue(2, 2, 2);
  A->setValue(3, 3, 2);
  A->assemble();

  b->setValue(mpiRank, mpiRank);
  b->assemble();

  PETScSolver::solveDirectMUMPS(*A, *b, *s, VERBOSE);

  PetscScalar v;
  int p[] = {mpiRank};
  VecGetValues(s->getRaw(), 1, p, &v);
  ASSERT_EQ(v, mpiRank / 2.0);
}
#endif*/

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
