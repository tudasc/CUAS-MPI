#include "fillgrid.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 4
#define GRID_COLS 7
#define GRID_ROWS 5

TEST(fillgridTest, fillGlobalIndicesBlocked) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid grid(GRID_COLS, GRID_ROWS);

  fillGlobalIndicesBlocked(grid);

  {
    auto &handle = grid.getReadHandle();
    switch (mpiRank) {
      case 0: {
        ASSERT_EQ(handle(0, 0), 0);
        ASSERT_EQ(handle(0, 1), 1);
        ASSERT_EQ(handle(0, 2), 2);
        ASSERT_EQ(handle(0, 3), 3);
        ASSERT_EQ(handle(1, 0), 4);
        ASSERT_EQ(handle(1, 1), 5);
        ASSERT_EQ(handle(1, 2), 6);
        ASSERT_EQ(handle(1, 3), 7);
        ASSERT_EQ(handle(2, 0), 8);
        ASSERT_EQ(handle(2, 1), 9);
        ASSERT_EQ(handle(2, 2), 10);
        ASSERT_EQ(handle(2, 3), 11);
        break;
      }
      case 1: {
        ASSERT_EQ(handle(0, 0), 12);
        ASSERT_EQ(handle(0, 1), 13);
        ASSERT_EQ(handle(0, 2), 14);
        ASSERT_EQ(handle(1, 0), 15);
        ASSERT_EQ(handle(1, 1), 16);
        ASSERT_EQ(handle(1, 2), 17);
        ASSERT_EQ(handle(2, 0), 18);
        ASSERT_EQ(handle(2, 1), 19);
        ASSERT_EQ(handle(2, 2), 20);
        break;
      }
      case 2: {
        ASSERT_EQ(handle(0, 0), 21);
        ASSERT_EQ(handle(0, 1), 22);
        ASSERT_EQ(handle(0, 2), 23);
        ASSERT_EQ(handle(0, 3), 24);
        ASSERT_EQ(handle(1, 0), 25);
        ASSERT_EQ(handle(1, 1), 26);
        ASSERT_EQ(handle(1, 2), 27);
        ASSERT_EQ(handle(1, 3), 28);
        break;
      }
      case 3: {
        ASSERT_EQ(handle(0, 0), 29);
        ASSERT_EQ(handle(0, 1), 30);
        ASSERT_EQ(handle(0, 2), 31);
        ASSERT_EQ(handle(1, 0), 32);
        ASSERT_EQ(handle(1, 1), 33);
        ASSERT_EQ(handle(1, 2), 34);
        break;
      }
      default: {
      }
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
