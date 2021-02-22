#include "PETScGrid.h"
#include "PETScVec.h"

#include "gtest/gtest.h"

#include <memory>

#include "petscdump.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 4
#define GRID_COLS 12
#define GRID_ROWS 16

TEST(PetscGridTest, size) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PetscGrid>(cols, rows);
  ASSERT_EQ(grid->getTotalNumOfCols(), cols);
  ASSERT_EQ(grid->getTotalNumOfRows(), rows);
}

// mpiSize needs to be a divisor of cols and rows
TEST(PetscGridTest, localsize) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  // ASSERT_EQ(size % mpiSize, 0);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PetscGrid>(cols, rows);
  int localcols;
  int localrows;
  localcols = grid->getLocalNumOfCols();
  localrows = grid->getLocalNumOfRows();

  // we have 4 blocks of 8*6 elements
  ASSERT_EQ(localrows, rows / 2);
  ASSERT_EQ(localcols, cols / 2);
}

TEST(PetscGridTest, setandgetGlobal) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PetscGrid>(cols, rows);
  auto glob = grid->getWriteHandle();
  for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
      glob(i, j) = mpiRank;
    }
  }
  glob.setValues();

  auto &glob2 = grid->getReadHandle();
  for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
      ASSERT_EQ(glob2(i, j), mpiRank);
    }
  }
}

TEST(PetscGridTest, constandzero) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PetscGrid>(cols, rows);

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

TEST(PetscGridTest, checkCorners) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PetscGrid>(cols, rows);

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

TEST(PetscGridTest, boundaryTest) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid = std::make_unique<PetscGrid>(cols, rows, 13.37);
  int counter = 0;

  // grid->setAsGlobal2dArr(grid->getAsGlobal2dArr());

  auto &loc = grid->getReadHandle();  // local
  for (int i = 0; i < grid->getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalGhostNumOfCols(); ++j) {
      if (mpiRank == 0 && (i == 0 || j == 0)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        counter++;
      }
      if (mpiRank == 1 && (i == 0 || j == grid->getLocalGhostNumOfCols() - 1)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        counter++;
      }
      if (mpiRank == 2 && (i == grid->getLocalGhostNumOfRows() - 1 || j == 0)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        counter++;
      }
      if (mpiRank == 3 && (i == grid->getLocalGhostNumOfRows() - 1 || j == grid->getLocalGhostNumOfCols() - 1)) {
        ASSERT_EQ(loc(i, j, GHOSTED), 13.37);
        counter++;
      }
    }
  }
  ASSERT_EQ(counter, 17);
}

TEST(PetscMatTest, copyTest) {
  // ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  auto grid_1 = std::make_unique<PetscGrid>(cols, rows, 5);

  grid_1->setConst(5);

  auto grid_2 = std::make_unique<PetscGrid>(cols, rows, 10);

  grid_2->setConst(10);

  if (grid_1->copy(*grid_2.get())) {
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

  auto grid_3 = std::make_unique<PetscGrid>(cols + 1, rows, 15);

  grid_3->setConst(15);

  if (int error = grid_1->copy(*grid_3.get())) {
    EXPECT_EQ(error, 1);
  }
}

TEST(PetscGridTest, setGlobalVecTest) {
  // ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = GRID_COLS;
  const int rows = GRID_ROWS;
  PetscGrid grid(cols, rows);

  grid.setConst(5);

  auto &loc = grid.getReadHandle();  // local
  for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
    for (int j = 0; j < grid.getLocalGhostNumOfCols(); ++j) {
      ASSERT_EQ(loc(i, j, GHOSTED), 5);
    }
  }

  PetscVec global(cols * rows);
  global.setConst(3);
  grid.setGlobalVecColMajor(global);

  // auto loc2 = grid.getAsGlobal2dArr();
  for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(loc(i, j), 3);
    }
  }

  for (int i = 0; i < cols * rows; ++i) {
    global.setValue(i, i);
  }
  global.assemble();

  grid.setGlobalVecColMajor(global);
  // dump(grid);

  for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
    int indexI = grid.getCornerY() + i;
    for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
      int indexJ = grid.getCornerX() + j;
      ASSERT_EQ(loc(i, j), rows * indexJ + indexI);
    }
  }
}

TEST(PetscGridTest, handleTest) {
  PetscGrid grid(217, 13);

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
