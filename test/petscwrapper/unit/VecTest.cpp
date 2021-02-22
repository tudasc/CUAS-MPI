#include "PETScVec.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 5
#define VEC_SIZE 20

TEST(PetscVecTest, size) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PetscVec>(size);
  ASSERT_EQ(vec->getSize(), size);
  int petscsize;
  VecGetSize(vec->getPetscRaw(), &petscsize);
  ASSERT_EQ(petscsize, size);
}

TEST(PetscVecTest, localsize) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  // ASSERT_EQ(size % mpiSize, 0);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PetscVec>(size);
  int petsclocalsize;
  VecGetLocalSize(vec->getPetscRaw(), &petsclocalsize);
  ASSERT_EQ(petsclocalsize, size / mpiSize);
}

TEST(PetscVecTest, initialvalues) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PetscVec>(size);

  PetscScalar v;
  int p[] = {mpiRank * (size / mpiSize) + 1};
  VecGetValues(vec->getPetscRaw(), 1, p, &v);
  ASSERT_EQ(v, 0);
}

TEST(PetscVecTest, setconst) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PetscVec>(size);

  vec->setConst(5.3);
  PetscScalar v;
  int p[] = {mpiRank * (size / mpiSize) + 1};
  VecGetValues(vec->getPetscRaw(), 1, p, &v);
  ASSERT_EQ(v, 5.3);

  vec->setZero();
  VecGetValues(vec->getPetscRaw(), 1, p, &v);
  ASSERT_EQ(v, 0);
}

TEST(PetscVecTest, setvalue) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PetscVec>(size);
  if (mpiRank == 1) {
    // vec->setValue(15, 2.3);
    vec->setValue(6, 2.7);
  }
  if (mpiRank == 4) {
    vec->setValue(16, 4.3);
    vec->setValue(17, 4.3);
    vec->setValue(15, 4.3);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  if (mpiRank == 3) {
    PetscScalar v;
    int p[] = {15};
    VecGetValues(vec->getPetscRaw(), 1, p, &v);
    ASSERT_EQ(v, 0);
  }
  if (mpiRank == 1) {
    PetscScalar v;
    int p[] = {6};
    VecGetValues(vec->getPetscRaw(), 1, p, &v);
    ASSERT_EQ(v, 2.7);
  }

  vec->assemble();

  if (mpiRank == 3) {
    PetscScalar v;
    int p[] = {15};
    VecGetValues(vec->getPetscRaw(), 1, p, &v);
    ASSERT_EQ(v, 4.3);
  }
  if (mpiRank == 1) {
    PetscScalar v;
    int p[] = {6};
    VecGetValues(vec->getPetscRaw(), 1, p, &v);
    ASSERT_EQ(v, 2.7);
  }
}

int main(int argc, char *argv[]) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  result = RUN_ALL_TESTS();
  PetscFinalize();

  return result;
}
