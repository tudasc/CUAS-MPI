#include "specialgradient.h"

#include "PetscGrid.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 9

TEST(GradientTest, result) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto grid = std::make_unique<PetscGrid>(20, 20);
  grid->setZero();

  auto myglobarr = grid->getWriteHandle();
  for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
      myglobarr(i, j) = j;
    }
  }
  myglobarr.setValues();

  auto gradient1 = std::make_unique<PetscGrid>(20, 20);
  CUAS::gradient2(*grid, *gradient1, 1.0);
  auto &grad_arr_1 = gradient1->getReadHandle();
  ASSERT_EQ(grad_arr_1(3, 3), 1);
  if (gradient1->getCornerX() == 0) {
    ASSERT_EQ(grad_arr_1(0, 0), 0.5);
  }
  if (gradient1->getCornerX() == 0 && gradient1->getCornerY() == 0) {
    ASSERT_EQ(grad_arr_1(0, 0), 0.5);
    ASSERT_EQ(grad_arr_1(0, 6), 36.5);
  }
}

TEST(GradientTest, ghostCells) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto grid = std::make_unique<PetscGrid>(20, 20);
  grid->setZero();

  auto myglobarr = grid->getWriteHandle();
  for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
      myglobarr(i, j) = j;
    }
  }

  myglobarr.setValues();

  auto gradient1 = std::make_unique<PetscGrid>(20, 20);
  CUAS::gradient2(*grid, *gradient1, 1.0);

  int localGhostNumOfRows = gradient1->getLocalGhostNumOfRows();
  int localGhostNumOfCols = gradient1->getLocalGhostNumOfCols();

  int cornerYGhost = gradient1->getCornerYGhost();
  int cornerXGhost = gradient1->getCornerXGhost();

  int totalNumOfRows = gradient1->getTotalNumOfRows();
  int totalNumOfCols = gradient1->getTotalNumOfCols();

  auto &local2d = gradient1->getReadHandle();  // local
  for (int i = 0; i < localGhostNumOfRows; ++i) {
    for (int j = 0; j < localGhostNumOfCols; ++j) {
      if (cornerYGhost == -1 && i == 0 || cornerXGhost == -1 && j == 0 ||
          cornerYGhost + localGhostNumOfRows - 1 == totalNumOfRows && i == localGhostNumOfRows - 1 ||
          cornerXGhost + localGhostNumOfCols - 1 == totalNumOfCols && j == localGhostNumOfCols - 1) {
        ASSERT_EQ(local2d(i, j, GHOSTED), 0);
      }
    }
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
