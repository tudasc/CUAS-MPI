#include "PETScGrid.h"
#include "PETScVector.h"

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
  grid_1->copy(*grid_2);

  SUCCEED();

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

  ASSERT_EXIT(grid_1->copy(*grid_3), ::testing::ExitedWithCode(1),
              "PETScGrid.cpp: copy: input is not compatible. Exiting.");
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

/*TEST(PETScGridTest, setGlobalVecColMajor) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 6;
  const int rows = 4;

  PETScGrid grid(cols, rows);
  grid.setZero();

  {
    PETScVector inputNoneGhosted(cols * rows);

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
    PETScVector inputGhosted((cols + 2) * (rows + 2));

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
}*/

/*TEST(PETScGridTest, setGlobalVecRowMajor) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 6;
  const int rows = 4;

  PETScGrid grid(cols, rows);
  grid.setZero();

  {
    PETScVector inputNoneGhosted(cols * rows);

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
    PETScVector inputGhosted((cols + 2) * (rows + 2));

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
}*/

TEST(PETScGridTest, setGhostBoundary) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 7;
  const int rows = 5;

  PETScGrid grid(cols, rows);
  grid.setGhostBoundary(5);

  {
    auto &gridHandle = grid.getReadHandle();

    // check unchanged
    for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
        ASSERT_EQ(gridHandle(i, j), 0.0) << "(row(i): " << i << ", col(j): " << j << ")";
      }
    }

    // check boundary
    for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalGhostNumOfCols(); ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0 || grid.getCornerXGhost() == -1 && j == 0 ||
            grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1 &&
                i == grid.getLocalGhostNumOfRows() - 1 ||
            grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1 &&
                j == grid.getLocalGhostNumOfCols() - 1) {
          ASSERT_EQ(gridHandle(i, j, GHOSTED), 5.0) << "(row(i): " << i << ", col(j): " << j << ")";
        }
      }
    }
  }
}

TEST(PETScGridTest, findAndReplaceGhostBoundary) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 7;
  const int rows = 5;

  PETScGrid grid(cols, rows);
  grid.setGhostBoundary(5);

  // set some values 3
  {
    auto gridHandle = grid.getWriteHandleGhost();
    auto startJ = (grid.getCornerXGhost() == -1) ? 1 : 0;
    auto endJ = (grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1)
                    ? grid.getLocalGhostNumOfCols() - 1
                    : grid.getLocalGhostNumOfCols();
    for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
      for (int j = startJ; j < endJ; ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0) {
          gridHandle(i, j) = 3;
        }
      }
    }
  }

  grid.findAndReplaceGhostBoundary(5, 7);

  {
    auto &gridHandle = grid.getReadHandle();

    // check unchanged
    for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
        ASSERT_EQ(gridHandle(i, j), 0.0) << "(row(i): " << i << ", col(j): " << j << ")";
      }
    }
    auto startJ = (grid.getCornerXGhost() == -1) ? 1 : 0;
    auto endJ = (grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1)
                    ? grid.getLocalGhostNumOfCols() - 1
                    : grid.getLocalGhostNumOfCols();
    for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
      for (int j = startJ; j < endJ; ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0)
          ASSERT_EQ(gridHandle(i, j, GHOSTED), 3.0) << "(row(i): " << i << ", col(j): " << j << ")";
      }
    }

    // check boundary
    for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalGhostNumOfCols(); ++j) {
        if (grid.getCornerXGhost() == -1 && j == 0 ||
            grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1 &&
                i == grid.getLocalGhostNumOfRows() - 1 ||
            grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1 &&
                j == grid.getLocalGhostNumOfCols() - 1) {
          ASSERT_EQ(gridHandle(i, j, GHOSTED), 7.0) << "(row(i): " << i << ", col(j): " << j << ")";
        }
      }
    }
  }
}

TEST(PETScGridTest, setRealBoundary) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 7;
  const int rows = 5;

  PETScGrid grid(cols, rows);
  grid.setRealBoundary(5);

  {
    auto &gridHandle = grid.getReadHandle();

    // check unchanged
    auto startI = (grid.getCornerYGhost() == -1) ? 1 : 0;
    auto startJ = (grid.getCornerXGhost() == -1) ? 1 : 0;
    auto endI = (grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1)
                    ? grid.getLocalNumOfRows() - 1
                    : grid.getLocalNumOfRows();
    auto endJ = (grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1)
                    ? grid.getLocalNumOfCols() - 1
                    : grid.getLocalNumOfCols();
    for (int i = startI; i < endI; ++i) {
      for (int j = startJ; j < endJ; ++j) {
        ASSERT_EQ(gridHandle(i, j), 0.0) << "(row(i): " << i << ", col(j): " << j << ")";
      }
    }
    for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalGhostNumOfCols(); ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0 || grid.getCornerXGhost() == -1 && j == 0 ||
            grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1 &&
                i == grid.getLocalGhostNumOfRows() - 1 ||
            grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1 &&
                j == grid.getLocalGhostNumOfCols() - 1) {
          ASSERT_EQ(gridHandle(i, j, GHOSTED), 0.0) << "(row(i): " << i << ", col(j): " << j << ")";
        }
      }
    }

    // check changed inner boundary ring
    for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0 || grid.getCornerXGhost() == -1 && j == 0 ||
            grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1 &&
                i == grid.getLocalNumOfRows() - 1 ||
            grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1 &&
                j == grid.getLocalNumOfCols() - 1) {
          ASSERT_EQ(gridHandle(i, j), 5.0) << "(row(i): " << i << ", col(j): " << j << ")";
        }
      }
    }
  }
}

TEST(PETScGridTest, findAndReplaceRealBoundary) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 7;
  const int rows = 5;

  PETScGrid grid(cols, rows);
  grid.setRealBoundary(5);

  // set some values 3
  {
    auto gridHandle = grid.getWriteHandle();
    auto startJ = (grid.getCornerXGhost() == -1) ? 1 : 0;
    auto endJ = (grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1)
                    ? grid.getLocalNumOfCols() - 1
                    : grid.getLocalNumOfCols();
    for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
      for (int j = startJ; j < endJ; ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0) {
          gridHandle(i, j) = 3;
        }
      }
    }
  }

  grid.findAndReplaceRealBoundary(5, 7);

  {
    auto &gridHandle = grid.getReadHandle();

    // check unchanged
    auto startI = (grid.getCornerYGhost() == -1) ? 1 : 0;
    auto startJ = (grid.getCornerXGhost() == -1) ? 1 : 0;
    auto endI = (grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1)
                    ? grid.getLocalNumOfRows() - 1
                    : grid.getLocalNumOfRows();
    auto endJ = (grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1)
                    ? grid.getLocalNumOfCols() - 1
                    : grid.getLocalNumOfCols();
    for (int i = startI; i < endI; ++i) {
      for (int j = startJ; j < endJ; ++j) {
        ASSERT_EQ(gridHandle(i, j), 0.0) << "(row(i): " << i << ", col(j): " << j << ")";
      }
    }
    for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalGhostNumOfCols(); ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0 || grid.getCornerXGhost() == -1 && j == 0 ||
            grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1 &&
                i == grid.getLocalGhostNumOfRows() - 1 ||
            grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1 &&
                j == grid.getLocalGhostNumOfCols() - 1) {
          ASSERT_EQ(gridHandle(i, j, GHOSTED), 0.0) << "(row(i): " << i << ", col(j): " << j << ")";
        }
      }
    }
    for (int i = 0; i < grid.getLocalGhostNumOfRows(); ++i) {
      for (int j = startJ; j < endJ; ++j) {
        if (grid.getCornerYGhost() == -1 && i == 0)
          ASSERT_EQ(gridHandle(i, j), 3.0) << "(row(i): " << i << ", col(j): " << j << ")";
      }
    }

    // check changed inner boundary ring
    for (int i = 0; i < grid.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < grid.getLocalNumOfCols(); ++j) {
        if (grid.getCornerXGhost() == -1 && j == 0 ||
            grid.getCornerYGhost() + grid.getLocalGhostNumOfRows() == grid.getTotalGhostNumOfRows() - 1 &&
                i == grid.getLocalNumOfRows() - 1 ||
            grid.getCornerXGhost() + grid.getLocalGhostNumOfCols() == grid.getTotalGhostNumOfCols() - 1 &&
                j == grid.getLocalNumOfCols() - 1) {
          ASSERT_EQ(gridHandle(i, j), 7.0) << "(row(i): " << i << ", col(j): " << j << ")";
        }
      }
    }
  }
}

TEST(PETScGridTest, getMaxAbsDiff) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 5;
  const int rows = 3;

  PETScGrid grid1(cols, rows);
  PETScGrid grid2(cols, rows);
  {
    auto grid1Handle = grid1.getWriteHandle();
    auto grid2Handle = grid2.getWriteHandle();
    for (int j = 0; j < grid1.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < grid1.getLocalNumOfCols(); ++i) {
        grid1Handle(j, i) = (i + 1) * (j + 1) * (-1 * mpiRank - 1);
        grid2Handle(j, i) = (i - j) * 10;
      }
    }
  }

  auto maxAbsDiff = grid1.getMaxAbsDiff(grid2);
  ASSERT_EQ(maxAbsDiff, 29.0);

  maxAbsDiff = grid2.getMaxAbsDiff(grid1);
  ASSERT_EQ(maxAbsDiff, 29.0);

  maxAbsDiff = grid1.getMaxAbsDiff(grid1);
  ASSERT_EQ(maxAbsDiff, 0.0);

  maxAbsDiff = grid2.getMaxAbsDiff(grid2);
  ASSERT_EQ(maxAbsDiff, 0.0);

  // check if both grids have not changed
  {
    auto &grid1Handle = grid1.getReadHandle();
    auto &grid2Handle = grid2.getReadHandle();
    for (int j = 0; j < grid1.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < grid1.getLocalNumOfCols(); ++i) {
        ASSERT_EQ(grid1Handle(j, i), (i + 1) * (j + 1) * (-1 * mpiRank - 1));
        ASSERT_EQ(grid2Handle(j, i), (i - j) * 10);
      }
    }
  }
}

TEST(PETScGridTest, getErrorNorms) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  const int cols = 5;
  const int rows = 3;

  PETScGrid grid1(cols, rows);
  PETScGrid grid2(cols, rows);
  {
    auto grid1Handle = grid1.getWriteHandle();
    auto grid2Handle = grid2.getWriteHandle();
    for (int j = 0; j < grid1.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < grid1.getLocalNumOfCols(); ++i) {
        grid1Handle(j, i) = (i + 1) * (j + 1) * (-1 * mpiRank - 1);
        grid2Handle(j, i) = (i - j) * 10;
      }
    }
  }

  auto [L1, L2, Linf] = grid1.getErrorNorms(grid2);
  CUAS_INFO_RANK0("{}: L1={}, L2={}, Linf={}\n", __PRETTY_FUNCTION__, L1, L2, Linf);

  ASSERT_EQ(L1, 164.0);
  // L2 = sqrt(1^2+12^2+23^2+8^2+4^2+16^2+2^2+14^2+6^2+8^2+3^2+16^2+29^2+4^2+18^2)
  ASSERT_NEAR(L2, 52.4976189936, 1e-6);
  ASSERT_EQ(Linf, 29.0);  // same as in maxAbsDiff()

  // check if both grids have not changed
  {
    auto &grid1Handle = grid1.getReadHandle();
    auto &grid2Handle = grid2.getReadHandle();
    for (int j = 0; j < grid1.getLocalNumOfRows(); ++j) {
      for (int i = 0; i < grid1.getLocalNumOfCols(); ++i) {
        ASSERT_EQ(grid1Handle(j, i), (i + 1) * (j + 1) * (-1 * mpiRank - 1));
        ASSERT_EQ(grid2Handle(j, i), (i - j) * 10);
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
