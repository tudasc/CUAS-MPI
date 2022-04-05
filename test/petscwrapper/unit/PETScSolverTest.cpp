#include "PETScSolver.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 4
#define VEC_SIZE 4
#define VERBOSE true

TEST(PetscSolverTest, solve) {
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

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
