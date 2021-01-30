#ifndef CUAS_PETSCGRID_H
#define CUAS_PETSCGRID_H

#include "petsc.h"

#include "PetscVec.h"

#include <iostream>

#define GHOSTED true
#define NONE_GHOSTED false

/*
 * We interpret grids in row major
 * access values[i * cols + j]
 * access 2darray[i][j]
 *
 * grids are stored row-wise in the global vector
 * +-------------------->(m, x, first dimension, cols, j)
 * |  0  1  4  5  8  9
 * |  2  3  6  7 10 11
 * | 12 13 16 17 20 21
 * | 14 15 18 19 22 23
 * |
 * v
 * (n, y, second dimension, rows, i)
 *
 * petsc does split the vector along the first dimension (m) of the grid first
 * they are assigned to processes in the following order:
 * +---------------------------------->(m, x, first dimension, cols, j)
 * | process 0   process 1   process 2
 * | process 4   process 2   process 5
 * |
 * v
 * (n, y, second dimension, rows, i)
 */
class PetscGrid {
 public:
  /*
   * Use ReadHandle to only read values of the grid.
   * if ghosted: GhostValues are being read aswell.
   */
  struct ReadHandle {
    PetscGrid const *grid;
    ReadHandle(PetscGrid const *grid) : grid(grid){};
    PetscScalar operator()(int i, int j, bool ghosted = NONE_GHOSTED) const {
      if (ghosted) {
        return grid->valuesGhosted[i][j];
      } else {
        return grid->values[i][j];
      }
    };
  };

  /*
   * Use WriteHandle to read and write values of the grid.
   * To explicitly set values now use setValues().
   * Otherwise the values are being set automatically when WriteHandle is out of scope.
   * Using WriteHandle you can not write to ghost-cells.
   */
  struct WriteHandle {
    PetscGrid *grid;
    WriteHandle(PetscGrid *grid) : grid(grid){};

    PetscScalar &operator()(int i, int j) { return grid->values[i][j]; };

    // only use if you want to set values before WriteHandle
    void setValues() { DMGlobalToLocal(grid->dm, grid->global, INSERT_VALUES, grid->local); };

    ~WriteHandle() { DMGlobalToLocal(grid->dm, grid->global, INSERT_VALUES, grid->local); }
  };

  // creates a Grid of numOfCols times numOfRows,
  // sets value of most outer cells to boundaryValue
  PetscGrid(int numOfCols, int numOfRows, PetscScalar boundaryValue = 0);
  ~PetscGrid();

  ReadHandle const &getReadHandle() const { return readHandle; }

  WriteHandle getWriteHandle() { return WriteHandle(this); }

  // indices: from 0 to LocalNumOfCols
  // you MUST call setAsGlobal2dArr or restoreGlobal2dArr afterwards!
  // returns a 2d Array of the processes share of data
  PetscScalar **getAsGlobal2dArr();

  // call getAsGlobal2dArr before; parameter: the changed return values.
  // sets the whole grid to given 2d values that have been obtained
  // using getAsGlobal2dArr()
  void setAsGlobal2dArr(PetscScalar **values);

  // setting local arrays is not &&e because of
  // possible race conditions when writing in ghost-cells.
  // restores the values from getAsLocal2dArr()
  void restoreLocal2dArr(PetscScalar **values);

  // use only if you do not want to change values.
  // restores the values from getAsGlobal2dArr()
  void restoreGlobal2dArr(PetscScalar **values);

  // sets the grids values from globalVec using the column major layout
  void setGlobalVecColMajor(PetscVec &globalVec);

  // sets each entry of grid to value
  void setConst(PetscScalar value);

  // zeroes each entry of grid
  void setZero() { setConst(0); }

  // copys the content of one grid to another
  int copy(PetscGrid const &input);

  // int countNonZero() const;

  /* getter */

  // Vec getLocalVec() { return local; }
  // Vec getGlobalVec() { return global; }

  int getLocalNumOfCols() const { return localNumOfCols; }
  int getLocalGhostNumOfCols() const { return localGhostNumOfCols; }
  int getLocalNumOfRows() const { return localNumOfRows; }
  int getLocalGhostNumOfRows() const { return localGhostNumOfRows; }
  int getTotalNumOfRows() const { return totalNumOfRows; }
  int getTotalNumOfCols() const { return totalNumOfCols; }

  int getCornerX() const { return cornerX; }
  int getCornerXGhost() const { return cornerXGhost; }
  int getCornerY() const { return cornerY; }
  int getCornerYGhost() const { return cornerYGhost; }

  // DM getDM() { return dm; }
 private:
  DM dm;

  Vec local;
  Vec global;

  PetscScalar **values;
  PetscScalar **valuesGhosted;

  // is this the way to go? compiler complains without it
  ReadHandle readHandle;

  const int totalNumOfCols;
  const int totalNumOfRows;

  int localNumOfCols;
  int localNumOfRows;

  int localGhostNumOfCols;
  int localGhostNumOfRows;

  int cornerX, cornerY, cornerXGhost, cornerYGhost;

  // indices: from 0 to LocalGhostWidth
  // you MUST call restoreLocal2dArr afterwards!
  // returns a 2d Array of the processes share of data
  // including ghostcells of neighbour processes.
  // normally this method is not needed, use Read/WriteHandles instead
  PetscScalar **getAsLocal2dArr();

  // sets the outer boundaries (ghost-cells) of the grid
  void setGlobalBoundariesConst(PetscScalar value);
};

#endif
