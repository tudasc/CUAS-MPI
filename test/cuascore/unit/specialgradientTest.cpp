#include "specialgradient.h"

#include "PETScGrid.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 6
#define GRID_SIZE 10

TEST(GradientTest, result) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto grid = std::make_unique<PETScGrid>(GRID_SIZE, GRID_SIZE);
  auto gradient = std::make_unique<PETScGrid>(GRID_SIZE, GRID_SIZE);
  grid->setZero();

  {
    auto myglobarr = grid->getWriteHandle();
    for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
      for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
        myglobarr(i, j) = (i + 3) * j;
      }
    }
  }

  CUAS::gradient2(*gradient, *grid, 2.0);

  {
    auto &gradHandle = gradient->getReadHandle();
    switch (mpiRank) {
      case 0: {
        ASSERT_EQ(gradHandle(0, 0, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 1, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 2, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 3, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 4, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 5, GHOSTED), 0);
        ASSERT_EQ(gradHandle(1, 0, GHOSTED), 0);
        ASSERT_EQ(gradHandle(1, 1, GHOSTED), 1.125);
        ASSERT_EQ(gradHandle(1, 2, GHOSTED), 3.5);
        ASSERT_EQ(gradHandle(1, 3, GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(1, 4, GHOSTED), 13.5);
        ASSERT_EQ(gradHandle(1, 5, GHOSTED), 39.125);
        ASSERT_EQ(gradHandle(2, 0, GHOSTED), 0);
        ASSERT_EQ(gradHandle(2, 1, GHOSTED), 2);
        ASSERT_EQ(gradHandle(2, 2, GHOSTED), 4.25);
        ASSERT_EQ(gradHandle(2, 3, GHOSTED), 5);
        ASSERT_EQ(gradHandle(2, 4, GHOSTED), 6.25);
        ASSERT_EQ(gradHandle(2, 5, GHOSTED), 38);
        ASSERT_EQ(gradHandle(3, 0, GHOSTED), 0);
        ASSERT_EQ(gradHandle(3, 1, GHOSTED), 3.125);
        ASSERT_EQ(gradHandle(3, 2, GHOSTED), 6.5);
        ASSERT_EQ(gradHandle(3, 3, GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(3, 4, GHOSTED), 8.5);
        ASSERT_EQ(gradHandle(3, 5, GHOSTED), 57.125);
        /*ASSERT_EQ(gradHandle(4, 0, GHOSTED), 0);
        ASSERT_EQ(gradHandle(4, 1, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 2, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 3, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 4, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 5, GHOSTED), );*/

        ASSERT_EQ(gradHandle(0, 0, NONE_GHOSTED), 1.125);
        ASSERT_EQ(gradHandle(0, 1, NONE_GHOSTED), 3.5);
        ASSERT_EQ(gradHandle(0, 2, NONE_GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(0, 3, NONE_GHOSTED), 13.5);
        ASSERT_EQ(gradHandle(1, 0, NONE_GHOSTED), 2);
        ASSERT_EQ(gradHandle(1, 1, NONE_GHOSTED), 4.25);
        ASSERT_EQ(gradHandle(1, 2, NONE_GHOSTED), 5);
        ASSERT_EQ(gradHandle(1, 3, NONE_GHOSTED), 6.25);
        ASSERT_EQ(gradHandle(2, 0, NONE_GHOSTED), 3.125);
        ASSERT_EQ(gradHandle(2, 1, NONE_GHOSTED), 6.5);
        ASSERT_EQ(gradHandle(2, 2, NONE_GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(2, 3, NONE_GHOSTED), 8.5);
        break;
      }
      case 1: {
        ASSERT_EQ(gradHandle(0, 0, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 1, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 2, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 3, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 4, GHOSTED), 0);
        ASSERT_EQ(gradHandle(0, 5, GHOSTED), 0);
        ASSERT_EQ(gradHandle(1, 0, GHOSTED), 39.125);
        ASSERT_EQ(gradHandle(1, 1, GHOSTED), 19.125);
        ASSERT_EQ(gradHandle(1, 2, GHOSTED), 3.5);
        ASSERT_EQ(gradHandle(1, 3, GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(1, 4, GHOSTED), 13.5);
        ASSERT_EQ(gradHandle(1, 5, GHOSTED), 39.125);
        ASSERT_EQ(gradHandle(2, 0, GHOSTED), 38);
        ASSERT_EQ(gradHandle(2, 1, GHOSTED), 34);
        ASSERT_EQ(gradHandle(2, 2, GHOSTED), 4.25);
        ASSERT_EQ(gradHandle(2, 3, GHOSTED), 5);
        ASSERT_EQ(gradHandle(2, 4, GHOSTED), 6.25);
        ASSERT_EQ(gradHandle(2, 5, GHOSTED), 38);
        ASSERT_EQ(gradHandle(3, 0, GHOSTED), 57.125);
        ASSERT_EQ(gradHandle(3, 1, GHOSTED), 53.125);
        ASSERT_EQ(gradHandle(3, 2, GHOSTED), 6.5);
        ASSERT_EQ(gradHandle(3, 3, GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(3, 4, GHOSTED), 8.5);
        ASSERT_EQ(gradHandle(3, 5, GHOSTED), 57.125);
        /*ASSERT_EQ(gradHandle(4, 0, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 1, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 2, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 3, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 4, GHOSTED), );
        ASSERT_EQ(gradHandle(4, 5, GHOSTED), 0);*/

        ASSERT_EQ(gradHandle(0, 0, NONE_GHOSTED), 19.125);
        ASSERT_EQ(gradHandle(0, 1, NONE_GHOSTED), 3.5);
        ASSERT_EQ(gradHandle(0, 2, NONE_GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(0, 3, NONE_GHOSTED), 13.5);
        ASSERT_EQ(gradHandle(1, 0, NONE_GHOSTED), 34);
        ASSERT_EQ(gradHandle(1, 1, NONE_GHOSTED), 4.25);
        ASSERT_EQ(gradHandle(1, 2, NONE_GHOSTED), 5);
        ASSERT_EQ(gradHandle(1, 3, NONE_GHOSTED), 6.25);
        ASSERT_EQ(gradHandle(2, 0, NONE_GHOSTED), 53.125);
        ASSERT_EQ(gradHandle(2, 1, NONE_GHOSTED), 6.5);
        ASSERT_EQ(gradHandle(2, 2, NONE_GHOSTED), 7.25);
        ASSERT_EQ(gradHandle(2, 3, NONE_GHOSTED), 8.5);
        break;
      }
        // TODO complete
    }
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
