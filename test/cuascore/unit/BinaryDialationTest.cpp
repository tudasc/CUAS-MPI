#include "PetscAlgorithms.h"
#include "PetscGrid.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 9
#define GRID_SIZE_X 20
#define GRID_SIZE_Y 15

TEST(PetscAlgorithmsTest, dialationtest) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto nf_mask = std::make_unique<PetscGrid>(GRID_SIZE_X, GRID_SIZE_Y);
  auto grad_mask = std::make_unique<PetscGrid>(GRID_SIZE_X, GRID_SIZE_Y);

  auto nf2d = nf_mask->getAsGlobal2dArr();

  auto rows = nf_mask->getLocalNumOfRows();
  auto cols = nf_mask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (i == j) {
        nf2d[i][j] = true;
      }
    }
  }

  nf_mask->setAsGlobal2dArr(nf2d);
  grad_mask->setZero();

  binaryDialation(*nf_mask, *grad_mask);

  auto nf_arr = nf_mask->getAsGlobal2dArr();
  auto grad_arr = grad_mask->getAsGlobal2dArr();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (i == j) {
        ASSERT_EQ(nf_arr[i][j], true);
        ASSERT_EQ(grad_arr[i][j], true);
        if (i < rows - 1)
          ASSERT_EQ(grad_arr[i + 1][j], true);
        if (i > 0)
          ASSERT_EQ(grad_arr[i - 1][j], true);
        if (j < cols - 1)
          ASSERT_EQ(grad_arr[i][j + 1], true);
        if (j > 0)
          ASSERT_EQ(grad_arr[i][j - 1], true);
      }
    }
  }

  nf_mask->restoreGlobal2dArr(nf_arr);
  grad_mask->restoreGlobal2dArr(grad_arr);
}

TEST(PetscAlgorithmsTest, ghostCellsTest) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto nf_mask = std::make_unique<PetscGrid>(GRID_SIZE_X, GRID_SIZE_Y);
  auto grad_mask = std::make_unique<PetscGrid>(GRID_SIZE_X, GRID_SIZE_Y);

  auto nf2d = nf_mask->getAsGlobal2dArr();

  auto rows = nf_mask->getLocalNumOfRows();
  auto cols = nf_mask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (i == j) {
        nf2d[i][j] = true;
      }
    }
  }

  nf_mask->setAsGlobal2dArr(nf2d);
  grad_mask->setZero();

  binaryDialation(*nf_mask, *grad_mask);

  auto grad_arr = grad_mask->getAsLocal2dArr();

  auto rowsGhost = grad_mask->getLocalGhostNumOfRows();
  auto colsGhost = grad_mask->getLocalGhostNumOfCols();

  auto xGhost = grad_mask->getCornerXGhost();
  auto yGhost = grad_mask->getCornerYGhost();

  auto total_rows = grad_mask->getTotalNumOfRows();
  auto total_cols = grad_mask->getTotalNumOfCols();

  for (int i = 0; i < rowsGhost; ++i) {
    for (int j = 0; j < colsGhost; ++j) {
      if (yGhost == -1 && i == 0 || xGhost == -1 && j == 0 ||
          yGhost + rowsGhost - 1 == total_rows && i == rowsGhost - 1 ||
          xGhost + colsGhost - 1 == total_cols && j == colsGhost - 1) {
        ASSERT_EQ(grad_arr[i][j], false);
      }
    }
  }

  grad_mask->restoreLocal2dArr(grad_arr);
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
