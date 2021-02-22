#include "PETScMat.h"
#include "PETScSolver.h"
#include "PETScVec.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 4
#define VEC_SIZE 4

TEST(PetscSolverTest, solve) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  PetscSolver solver;
  auto A = std::make_unique<PetscMat>(VEC_SIZE, VEC_SIZE);
  auto b = std::make_unique<PetscVec>(VEC_SIZE);
  auto s = std::make_unique<PetscVec>(VEC_SIZE);

  A->setValue(0, 0, 2);
  A->setValue(1, 1, 2);
  A->setValue(2, 2, 2);
  A->setValue(3, 3, 2);
  A->assemble();

  b->setValue(mpiRank, mpiRank);
  b->assemble();

  solver.solve(*A, *b, *s);

  PetscScalar v;
  int p[] = {mpiRank};
  VecGetValues(s->getPetscRaw(), 1, p, &v);
  ASSERT_EQ(v, mpiRank / 2.0);
}

int main(int argc, char *argv[]) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  result = RUN_ALL_TESTS();
  PetscFinalize();
  return 0;

  return result;
}