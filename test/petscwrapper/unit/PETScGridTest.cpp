#include "PETScGrid.h"
#include "PETScVec.h"

#include "gtest/gtest.h"

#include <memory>

int mpiRank;
int mpiSize;

#define MPI_SIZE 4
#define GRID_COLS 12
#define GRID_ROWS 16

TEST(PETScGridTest, size) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PETScGrid>(cols, rows);
  ASSERT_EQ(grid->getTotalNumOfCols(), cols);
  ASSERT_EQ(grid->getTotalNumOfRows(), rows);
  ASSERT_EQ(grid->getTotalGhostNumOfCols(), cols + 2);
  ASSERT_EQ(grid->getTotalGhostNumOfRows(), rows + 2);
}

// mpiSize needs to be a divisor of cols and rows
TEST(PETScGridTest, localsize) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  // ASSERT_EQ(size % mpiSize, 0);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PETScGrid>(cols, rows);

  // we have 4 blocks of 7*5 elements
  // this leads to 15*10 elements + 2 cols and 2 rows of ghost cells
  ASSERT_EQ(grid->getLocalNumOfRows(), rows / 2);
  ASSERT_EQ(grid->getLocalNumOfCols(), cols / 2);
  ASSERT_EQ(grid->getLocalGhostNumOfRows(), rows / 2 + 2);
  ASSERT_EQ(grid->getLocalGhostNumOfCols(), cols / 2 + 2);
}

TEST(PETScGridTest, setandgetGlobal) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PETScGrid>(cols, rows);
  auto glob = grid->getWriteHandle();
  for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
      glob(i, j) = mpiRank;
    }
  }
  glob.setValues();

  // TODO check ghost cells
  auto &glob2 = grid->getReadHandle();
  for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(glob2(i, j), mpiRank);
    }
  }
}

TEST(PETScGridTest, constandzero) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PETScGrid>(cols, rows);

  grid->setConst(5);

  auto &loc = grid->getReadHandle();  // local
  for (int i = 0; i < grid->getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalGhostNumOfCols(); ++j) {
      ASSERT_EQ(loc(i, j, GHOSTED), 5);
    }
  }

  grid->setZero();

  for (int i = 0; i < grid->getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalGhostNumOfCols(); ++j) {
      ASSERT_EQ(loc(i, j, GHOSTED), 0);
    }
  }
}

TEST(PETScGridTest, checkCorners) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PETScGrid>(cols, rows);

  int x = grid->getCornerX();
  int xGhost = grid->getCornerXGhost();

  int y = grid->getCornerY();
  int yGhost = grid->getCornerYGhost();

  if (mpiRank == 0) {
    ASSERT_EQ(xGhost, -1);
    ASSERT_EQ(yGhost, -1);
    ASSERT_EQ(x, 0);
    ASSERT_EQ(y, 0);
  }
  if (mpiRank == 1) {
    ASSERT_EQ(xGhost, 5);
    ASSERT_EQ(yGhost, -1);
    ASSERT_EQ(x, 6);
    ASSERT_EQ(y, 0);
  }
  if (mpiRank == 2) {
    ASSERT_EQ(xGhost, -1);
    ASSERT_EQ(yGhost, 7);
    ASSERT_EQ(x, 0);
    ASSERT_EQ(y, 8);
  }
  if (mpiRank == 3) {
    ASSERT_EQ(xGhost, 5);
    ASSERT_EQ(yGhost, 7);
    ASSERT_EQ(x, 6);
    ASSERT_EQ(y, 8);
  }
}

TEST(PETScGridTest, boundaryTest) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PETScGrid>(cols, rows, 13.37);
  int counter = 0;

  // grid->setAsGlobal2dArr(grid->getAsGlobal2dArr());

  auto &loc = grid->getReadHandle();  // local
  for (int i = 0; i < grid->getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalGhostNumOfCols(); ++j) {
      if (mpiRank == 0 && (i == 0 || j == 0)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        ++counter;
      }
      if (mpiRank == 1 && (i == 0 || j == grid->getLocalGhostNumOfCols() - 1)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        ++counter;
      }
      if (mpiRank == 2 && (i == grid->getLocalGhostNumOfRows() - 1 || j == 0)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        ++counter;
      }
      if (mpiRank == 3 && (i == grid->getLocalGhostNumOfRows() - 1 || j == grid->getLocalGhostNumOfCols() - 1)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        ++counter;
      }
    }
  }
  ASSERT_EQ(counter, 17);
}

TEST(PETScGridTest, copyTest) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid_1 = std::make_unique<PETScGrid>(cols, rows, 5);

  grid_1->setConst(5);

  auto grid_2 = std::make_unique<PETScGrid>(cols, rows, 10);

  grid_2->setConst(10);

  if (grid_1->copy(*grid_2)) {
    FAIL() << "Did not expect an error!";
  }

  // test if all values have been written to the other grid

  auto &loc = grid_1->getReadHandle();
  for (int i = 0; i < grid_1->getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < grid_1->getLocalGhostNumOfCols(); ++j) {
      ASSERT_EQ(loc(i, j, GHOSTED), 10);
    }
  }

  // test if an error is thrown when the grids have different sizes

  auto grid_3 = std::make_unique<PETScGrid>(cols + 1, rows, 15);

  grid_3->setConst(15);

  if (int error = grid_1->copy(*grid_3)) {
    EXPECT_EQ(error, 1);
  }
}

TEST(PETScGridTest, handleTest) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid grid(217, 13);

  // test handles with explicit setValues()
  auto w = grid.getWriteHandle();
  for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
      w(i, j) = i * grid.getLocalNumOfRows() + j;
    }
  }
  w.setValues();
  auto &r = grid.getReadHandle();
  for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(r(i, j), i * grid.getLocalNumOfRows() + j);
    }
  }
  // test handles with out of scope
  {
    auto w2 = grid.getWriteHandle();
    for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
        w2(i, j) = (i * grid.getLocalNumOfRows() + j) * 5;
      }
    }
  }

  for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(r(i, j), (i * grid.getLocalNumOfRows() + j) * 5);
    }
  }
}

TEST(PETScGridTest, setGlobalVecColMajor) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 6;
  const int rows = 4;

  PETScGrid grid(cols, rows);
  grid.setZero();

  {
    PETScVec inputNoneGhosted(cols * rows);

    if (mpiRank == 0) {
      for (int i = 0; i < inputNoneGhosted.getSize(); ++i) {
        inputNoneGhosted.setValue(i, i);
      }
    }

    inputNoneGhosted.assemble();

    grid.setGlobalVecColMajor(inputNoneGhosted, NONE_GHOSTED);

    auto &handle = grid.getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(handle(0, 0, GHOSTED), 0);
      ASSERT_EQ(handle(0, 1, GHOSTED), 0);
      ASSERT_EQ(handle(0, 2, GHOSTED), 0);
      ASSERT_EQ(handle(0, 3, GHOSTED), 0);
      ASSERT_EQ(handle(1, 0, GHOSTED), 0);
      ASSERT_EQ(handle(1, 1, GHOSTED), 0);
      ASSERT_EQ(handle(1, 2, GHOSTED), 4);
      ASSERT_EQ(handle(1, 3, GHOSTED), 8);
      ASSERT_EQ(handle(2, 0, GHOSTED), 0);
      ASSERT_EQ(handle(2, 1, GHOSTED), 1);
      ASSERT_EQ(handle(2, 2, GHOSTED), 5);
      ASSERT_EQ(handle(2, 3, GHOSTED), 9);
    }
  }

  {
    PETScVec inputGhosted((cols + 2) * (rows + 2));

    if (mpiRank == 0) {
      for (int i = 0; i < inputGhosted.getSize(); ++i) {
        inputGhosted.setValue(i, i);
      }
    }

    inputGhosted.assemble();

    grid.setGlobalVecColMajor(inputGhosted, GHOSTED);

    auto &handle = grid.getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(handle(0, 0, GHOSTED), 0);
      ASSERT_EQ(handle(0, 1, GHOSTED), 6);
      ASSERT_EQ(handle(0, 2, GHOSTED), 12);
      ASSERT_EQ(handle(0, 3, GHOSTED), 18);
      ASSERT_EQ(handle(1, 0, GHOSTED), 1);
      ASSERT_EQ(handle(1, 1, GHOSTED), 7);
      ASSERT_EQ(handle(1, 2, GHOSTED), 13);
      ASSERT_EQ(handle(1, 3, GHOSTED), 19);
      ASSERT_EQ(handle(2, 0, GHOSTED), 2);
      ASSERT_EQ(handle(2, 1, GHOSTED), 8);
      ASSERT_EQ(handle(2, 2, GHOSTED), 14);
      ASSERT_EQ(handle(2, 3, GHOSTED), 20);
    }
  }
}

TEST(PETScGridTest, setGlobalVecRowMajor) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 6;
  const int rows = 4;

  PETScGrid grid(cols, rows);
  grid.setZero();

  {
    PETScVec inputNoneGhosted(cols * rows);

    if (mpiRank == 0) {
      for (int i = 0; i < inputNoneGhosted.getSize(); ++i) {
        inputNoneGhosted.setValue(i, i);
      }
    }

    inputNoneGhosted.assemble();

    grid.setGlobalVecRowMajor(inputNoneGhosted, NONE_GHOSTED);

    auto &handle = grid.getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(handle(0, 0, GHOSTED), 0);
      ASSERT_EQ(handle(0, 1, GHOSTED), 0);
      ASSERT_EQ(handle(0, 2, GHOSTED), 0);
      ASSERT_EQ(handle(0, 3, GHOSTED), 0);
      ASSERT_EQ(handle(1, 0, GHOSTED), 0);
      ASSERT_EQ(handle(1, 1, GHOSTED), 0);
      ASSERT_EQ(handle(1, 2, GHOSTED), 1);
      ASSERT_EQ(handle(1, 3, GHOSTED), 2);
      ASSERT_EQ(handle(2, 0, GHOSTED), 0);
      ASSERT_EQ(handle(2, 1, GHOSTED), 6);
      ASSERT_EQ(handle(2, 2, GHOSTED), 7);
      ASSERT_EQ(handle(2, 3, GHOSTED), 8);
    }
  }

  {
    PETScVec inputGhosted((cols + 2) * (rows + 2));

    if (mpiRank == 0) {
      for (int i = 0; i < inputGhosted.getSize(); ++i) {
        inputGhosted.setValue(i, i);
      }
    }

    inputGhosted.assemble();

    grid.setGlobalVecRowMajor(inputGhosted, GHOSTED);

    auto &handle = grid.getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(handle(0, 0, GHOSTED), 0);
      ASSERT_EQ(handle(0, 1, GHOSTED), 1);
      ASSERT_EQ(handle(0, 2, GHOSTED), 2);
      ASSERT_EQ(handle(0, 3, GHOSTED), 3);
      ASSERT_EQ(handle(1, 0, GHOSTED), 8);
      ASSERT_EQ(handle(1, 1, GHOSTED), 9);
      ASSERT_EQ(handle(1, 2, GHOSTED), 10);
      ASSERT_EQ(handle(1, 3, GHOSTED), 11);
      ASSERT_EQ(handle(2, 0, GHOSTED), 16);
      ASSERT_EQ(handle(2, 1, GHOSTED), 17);
      ASSERT_EQ(handle(2, 2, GHOSTED), 18);
      ASSERT_EQ(handle(2, 3, GHOSTED), 19);
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
