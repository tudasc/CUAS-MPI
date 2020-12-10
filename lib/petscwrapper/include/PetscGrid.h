#ifndef CUAS_PETSCGRID_H
#define CUAS_PETSCGRID_H

#include "petsc.h"

class PetscGrid {
 public:
  // creates a Grid of numOfCols times numOfRows,
  // sets value of most outer cells to boundaryValue
  PetscGrid(int numOfCols, int numOfRows, PetscScalar boundaryValue = 0);
  ~PetscGrid();

  // indices: from 0 to LocalGhostWidth
  // you MUST call restoreLocal2dArr afterwards!
  // returns a 2d Array of the processes share of data
  // including ghostcells of neighbour processes
  PetscScalar **getAsLocal2dArr();

  // indices: from 0 to LocalNumOfCols
  // you MUST call setAsGlobal2dArr or restoreGlobal2dArr afterwards!
  // returns a 2d Array of the processes share of data
  PetscScalar **getAsGlobal2dArr();

  // call getAsGlobal2dArr before; parameter: the changed return values.
  // sets the whole grid to given 2d values that have been obtained
  // using getAsGlobal2dArr()
  void setAsGlobal2dArr(PetscScalar **values);

  // setting local arrays is not possible because of
  // possible race conditions when writing in ghost-cells.
  // restores the values from getAsLocal2dArr()
  void restoreLocal2dArr(PetscScalar **values);

  // use only if you do not want to change values.
  // restores the values from getAsGlobal2dArr()
  void restoreGlobal2dArr(PetscScalar **values);

  // void setGlobalVecAndUpdate(Vec glob);

  // sets each entry of grid to value
  void setConst(PetscScalar value);

  // zeroes each entry of grid
  void setZero() { setConst(0); }

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

  const int totalNumOfCols;
  const int totalNumOfRows;

  int localNumOfCols;
  int localNumOfRows;

  int localGhostNumOfCols;
  int localGhostNumOfRows;

  int cornerX, cornerY, cornerXGhost, cornerYGhost;

  // sets the outer boundaries (ghost-cells) of the grid
  void setGlobalBoundariesConst(PetscScalar value);
};

#endif
