#include "PETScMatrix.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 5
#define MAT_COLS 20
#define MAT_ROWS 25

TEST(PETScMatrixTest, size) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = MAT_COLS;
  const int rows = MAT_ROWS;
  auto mat = std::make_unique<PETScMatrix>(rows, cols);
  ASSERT_EQ(mat->getNumberOfCols(), cols);
  ASSERT_EQ(mat->getNumberOfRows(), rows);

  int petscCols;
  int petscRows;
  MatGetSize(mat->getRaw(), &petscRows, &petscCols);
  ASSERT_EQ(petscCols, cols);
  ASSERT_EQ(petscRows, rows);
}

// mpiSize needs to be a divisor of cols and rows
TEST(PETScMatrixTest, localsize) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  // ASSERT_EQ(size % mpiSize, 0);

  const int cols = MAT_COLS;
  const int rows = MAT_ROWS;
  auto mat = std::make_unique<PETScMatrix>(rows, cols);
  int petsclocalcols;
  int petsclocalrows;
  MatGetLocalSize(mat->getRaw(), &petsclocalrows, &petsclocalcols);
  ASSERT_EQ(petsclocalrows, rows / mpiSize);
  ASSERT_EQ(petsclocalcols, cols / mpiSize);
}

TEST(PETScMatrixTest, getColsRows) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = MAT_COLS;
  const int rows = MAT_ROWS;
  auto mat = std::make_unique<PETScMatrix>(rows, cols);

  int colsFromMat = mat->getNumberOfCols();
  ASSERT_EQ(cols, colsFromMat);

  int rowsFromMat = mat->getNumberOfRows();
  ASSERT_EQ(rows, rowsFromMat);
}

TEST(PETScMatrixTest, setvalue) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = MAT_COLS;
  const int rows = MAT_ROWS;
  auto mat = std::make_unique<PETScMatrix>(rows, cols);
  mat->assemble();

  if (mpiRank == 0) {
    mat->setValue(18, 19, 2.7);
  }
  if (mpiRank == 2) {
    mat->setValue(9, 2, 1.337);
  }
  if (mpiRank == 3) {
    mat->setValue(24, 19, 13.37);
  }
  if (mpiRank == 4) {
    mat->setValue(1, 2, 133.7);
  }

  mat->assemble();

  // to access the matrix, the value needs to be in the processes share of
  // the matrix
  // seems like petsc stores the matrix row-wise
  // 0: 0-4
  // 1: 5-9
  // 2: 10-14
  // 3: 15-19
  // 4: 20-24
  if (mpiRank == 3) {
    PetscScalar v;
    int p1[] = {18};
    int p2[] = {19};
    MatGetValues(mat->getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 2.7);
  }
  if (mpiRank == 1) {
    PetscScalar v;
    int p1[] = {9};
    int p2[] = {2};
    MatGetValues(mat->getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 1.337);
  }
  if (mpiRank == 4) {
    PetscScalar v;
    int p1[] = {24};
    int p2[] = {19};
    MatGetValues(mat->getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 13.37);
  }
  if (mpiRank == 0) {
    PetscScalar v;
    int p1[] = {1};
    int p2[] = {2};
    MatGetValues(mat->getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 133.7);
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
