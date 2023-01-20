/**
 * File: PETScVectorTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "PETScVector.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 5
#define VEC_SIZE 20

TEST(PETScVectorTest, size) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PETScVector>(size);
  ASSERT_EQ(vec->getSize(), size);
  int petscsize;
  VecGetSize(vec->getRaw(), &petscsize);
  ASSERT_EQ(petscsize, size);
}

TEST(PETScVectorTest, localsize) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  // ASSERT_EQ(size % mpiSize, 0);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PETScVector>(size);
  int petsclocalsize;
  VecGetLocalSize(vec->getRaw(), &petsclocalsize);
  ASSERT_EQ(petsclocalsize, size / mpiSize);
}

TEST(PETScVectorTest, initialvalues) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PETScVector>(size);

  PetscScalar v;
  int p[] = {mpiRank * (size / mpiSize) + 1};
  VecGetValues(vec->getRaw(), 1, p, &v);
  ASSERT_EQ(v, 0);
}

TEST(PETScVectorTest, setconst) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PETScVector>(size);

  vec->setConst(5.3);
  PetscScalar v;
  int p[] = {mpiRank * (size / mpiSize) + 1};
  VecGetValues(vec->getRaw(), 1, p, &v);
  ASSERT_EQ(v, 5.3);

  vec->setZero();
  VecGetValues(vec->getRaw(), 1, p, &v);
  ASSERT_EQ(v, 0);
}

TEST(PETScVectorTest, setvalue) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int size = VEC_SIZE;
  auto vec = std::make_unique<PETScVector>(size);
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
    VecGetValues(vec->getRaw(), 1, p, &v);
    ASSERT_EQ(v, 0);
  }
  if (mpiRank == 1) {
    PetscScalar v;
    int p[] = {6};
    VecGetValues(vec->getRaw(), 1, p, &v);
    ASSERT_EQ(v, 2.7);
  }

  vec->assemble();

  if (mpiRank == 3) {
    PetscScalar v;
    int p[] = {15};
    VecGetValues(vec->getRaw(), 1, p, &v);
    ASSERT_EQ(v, 4.3);
  }
  if (mpiRank == 1) {
    PetscScalar v;
    int p[] = {6};
    VecGetValues(vec->getRaw(), 1, p, &v);
    ASSERT_EQ(v, 2.7);
  }
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
