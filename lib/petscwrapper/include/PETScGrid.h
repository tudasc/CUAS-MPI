/**
 * File: PETScGrid.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PETSCGRID_H
#define CUAS_PETSCGRID_H

#include "petsc.h"

#include "Logger.h"
#include "PETScVector.h"

#include <array>

#define GHOSTED true
#define NONE_GHOSTED false

/*
 * We interpret grids in row major
 * access values[i * xAxis + j]
 * access 2darray[i][j]
 *
 * grids are stored row-wise in the global vector
 * +-------------------->(m, x, first dimension, cols, xAxis, j)
 * |  0  1  4  5  8  9
 * |  2  3  6  7 10 11
 * | 12 13 16 17 20 21
 * | 14 15 18 19 22 23
 * |
 * v
 * (n, y, second dimension, rows, yAxis, i)
 *
 * petsc does split the vector along the first dimension (m) of the grid first
 * they are assigned to processes in the following order:
 * +---------------------------------->(m, x, first dimension, cols, xAxis, j)
 * | process 0   process 1   process 2
 * | process 4   process 2   process 5
 * |
 * v
 * (n, y, second dimension, rows, yAxis, i)
 */

class PETScGrid {
 public:
  /** Use ReadHandle to only read values of the grid.
   *
   * If ghosted: GhostValues are being read as well.
   */
  struct ReadHandle {
    PETScGrid const *const grid;
    explicit ReadHandle(PETScGrid const *grid) : grid(grid){};
    // disable copy of ReadHandle to ensure that getReadHandle is only used to get a reference and avoid copies
    ReadHandle(const ReadHandle &) = delete;
    ReadHandle &operator=(ReadHandle const &) = delete;
    ReadHandle(const ReadHandle &&) = delete;
    ReadHandle &operator=(ReadHandle const &&) = delete;
    ~ReadHandle() = default;

    PetscScalar operator()(int row, int col, bool ghosted = NONE_GHOSTED) const {
      PetscScalar result;
      if (ghosted) {
        if (row < 0 || col < 0 || row >= grid->getLocalGhostNumOfRows() || col >= grid->getLocalGhostNumOfCols()) {
          CUAS_ERROR("PETScGrid.h: ReadHandle: Access out of range. Exiting.")
          exit(1);
        }
        result = grid->valuesGhosted[row][col];
      } else {
        if (row < 0 || col < 0 || row >= grid->getLocalNumOfRows() || col >= grid->getLocalNumOfCols()) {
          CUAS_ERROR("PETScGrid.h: ReadHandle: Access out of range. Exiting.")
          exit(1);
        }
        result = grid->values[row][col];
      }
      return result;
    };

    PetscScalar const *const *getRaw() const { return grid->values; }
    PetscScalar const *const *getRawGhosted() const { return grid->valuesGhosted; }
  };

  /** Use WriteHandle to read and write values of the grid.
   *
   * To explicitly set values now use setValues().
   * Otherwise the values are being set automatically when WriteHandle is out of scope.
   * Using WriteHandle you can not write to ghost-cells.
   */
  struct WriteHandle {
    PETScGrid *const grid;
    explicit WriteHandle(PETScGrid *grid) : grid(grid){};

    PetscScalar &operator()(int row, int col) {
      if (row < 0 || col < 0 || row >= grid->getLocalNumOfRows() || col >= grid->getLocalNumOfCols()) {
        CUAS_ERROR("PETScGrid.h: WriteHandle: Access out of range. Exiting.")
        exit(1);
      }
      return grid->values[row][col];
    };

    // only use if you want to set values before destruction of WriteHandle
    // Please keep in mind that this might cause double communication
    void setValues() { DMGlobalToLocal(grid->dm, grid->global, INSERT_VALUES, grid->local); };

    ~WriteHandle() { setValues(); }

    PetscScalar *const *getRaw() { return grid->values; }
  };

  struct WriteHandleGhost {
    PETScGrid *const grid;
    explicit WriteHandleGhost(PETScGrid *grid) : grid(grid){};

    PetscScalar &operator()(int row, int col) {
      if (row < 0 || col < 0 || row >= grid->getLocalGhostNumOfRows() || col >= grid->getLocalGhostNumOfCols()) {
        CUAS_ERROR("PETScGrid.h: WriteHandleGhost: Access out of range. Exiting.")
        exit(1);
      }
      return grid->valuesGhosted[row][col];
    };

    // only use if you want to set values before destruction of WriteHandle
    // Please keep in mind that this might cause double communication
    void setValues() { DMLocalToGlobal(grid->dm, grid->local, INSERT_VALUES, grid->global); };

    ~WriteHandleGhost() { setValues(); }

    PetscScalar *const *getRaw() { return grid->valuesGhosted; }
  };

  // creates a Grid of numOfCols times numOfRows,
  // sets value of most outer cells to boundaryValue
  explicit PETScGrid(int numOfCols, int numOfRows, PetscScalar boundaryValue = 0);
  PETScGrid(PETScGrid &) = delete;
  PETScGrid(PETScGrid &&) = delete;
  ~PETScGrid();

  ReadHandle const &getReadHandle() const { return readHandle; }

  WriteHandle getWriteHandle() { return WriteHandle(this); }

  WriteHandleGhost getWriteHandleGhost() { return WriteHandleGhost(this); }

  // max(abs(this - sub)) this is L_inf error norm (PETSc: NORM_INFINITY)
  PetscScalar getMaxAbsDiff(PETScGrid const &sub) const;

  // returns NORM_1, NORM_2 and NORM_INFINITY
  std::array<PetscScalar, 3> getErrorNorms(PETScGrid const &sub) const;

  // maximum value of the grid
  PetscScalar getMax() const;

  // sets each entry of grid to value
  void setConst(PetscScalar value);

  // zeroes each entry of grid
  void setZero() { setConst(0); }

  void applyOffset(PetscScalar offset) {
    auto fWrite = getWriteHandle();
    auto rows = getLocalNumOfRows();
    auto cols = getLocalNumOfCols();
    for (int row = 0; row < rows; ++row) {
      for (int col = 0; col < cols; ++col) {
        fWrite(row, col) = fWrite(row, col) + offset;
      }
    }
  }

  void applyMultiplier(PetscScalar multiplier) {
    auto fWrite = getWriteHandle();
    auto rows = getLocalNumOfRows();
    auto cols = getLocalNumOfCols();
    for (int row = 0; row < rows; ++row) {
      for (int col = 0; col < cols; ++col) {
        fWrite(row, col) = fWrite(row, col) * multiplier;
      }
    }
  }

  // copys the content of one grid to another
  void copy(PETScGrid const &input);
  void copyGlobal(PETScGrid const &input);

  bool isCompatible(PETScGrid const &grid) const {
    return totalNumOfCols == grid.totalNumOfCols && totalNumOfRows == grid.totalNumOfRows &&
           localNumOfCols == grid.localNumOfCols && localNumOfRows == grid.localNumOfRows &&
           localGhostNumOfCols == grid.localGhostNumOfCols && localGhostNumOfRows == grid.localGhostNumOfRows &&
           totalGhostNumOfCols == grid.totalGhostNumOfCols && totalGhostNumOfRows == grid.totalGhostNumOfRows;
  }

  bool isOnGhostBoundary(int row, int col) const;
  bool isOnRealBoundary(int row, int col) const;

  // sets the outer boundaries (ghost-cells) of the grid
  void setGhostBoundary(PetscScalar value);
  void setRealBoundary(PetscScalar value);
  void findAndReplaceGhostBoundary(PetscScalar oldValue, PetscScalar newValue);
  void findAndReplaceRealBoundary(PetscScalar oldValue, PetscScalar newValue);

  // sets the grids values from globalVec using the column major layout
  // void setGlobalVecColMajor(PETScVector &globalVec, bool ghosted = NONE_GHOSTED);
  // sets the grids values from globalVec using the column major layout
  // void setGlobalVecRowMajor(PETScVector &globalVec, bool ghosted = NONE_GHOSTED);

  // int countNonZero() const;

  /* getter */

  // Vec getLocalVec() { return local; }
  // Vec getGlobalVec() { return global; }

  int getLocalNumOfCols() const { return localNumOfCols; }
  int getLocalNumOfRows() const { return localNumOfRows; }
  int getLocalGhostNumOfCols() const { return localGhostNumOfCols; }
  int getLocalGhostNumOfRows() const { return localGhostNumOfRows; }
  int getTotalNumOfRows() const { return totalNumOfRows; }
  int getTotalNumOfCols() const { return totalNumOfCols; }
  int getTotalGhostNumOfRows() const { return totalGhostNumOfRows; }
  int getTotalGhostNumOfCols() const { return totalGhostNumOfCols; }

  int getCornerX() const { return cornerX; }
  int getCornerXGhost() const { return cornerXGhost; }
  int getCornerY() const { return cornerY; }
  int getCornerYGhost() const { return cornerYGhost; }

  // DM getDM() { return dm; }
 public:
  DM dm;

  Vec local;
  Vec global;

 private:
  PetscScalar **values;
  PetscScalar **valuesGhosted;

  const ReadHandle readHandle;

  const int totalNumOfCols;
  const int totalNumOfRows;

  const int totalGhostNumOfCols;
  const int totalGhostNumOfRows;

  int localNumOfCols;
  int localNumOfRows;

  int localGhostNumOfCols;
  int localGhostNumOfRows;

  int cornerX, cornerY, cornerXGhost, cornerYGhost;

  // indices: from 0 to LocalGhostWidth
  // you MUST call restoreLocal2dArr afterwards!
  // returns a 2d Array of the processes share of data
  // including ghostcells of neighbour processes
  // PetscScalar **getAsLocal2dArr();

  // indices: from 0 to LocalNumOfCols
  // you MUST call setAsGlobal2dArr or restoreGlobal2dArr afterwards!
  // returns a 2d Array of the processes share of data
  // PetscScalar **getAsGlobal2dArr();

  // call getAsGlobal2dArr before; parameter: the changed return values.
  // sets the whole grid to given 2d values that have been obtained
  // using getAsGlobal2dArr()
  // void setAsGlobal2dArr(PetscScalar **values);

  // setting local arrays is not possible because of
  // possible race conditions when writing in ghost-cells.
  // restores the values from getAsLocal2dArr()
  // void restoreLocal2dArr(PetscScalar **values);

  // use only if you do not want to change values.
  // restores the values from getAsGlobal2dArr()
  // void restoreGlobal2dArr(PetscScalar **values);
};

#endif
