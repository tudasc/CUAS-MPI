#ifndef CUAS_PETSCGRID_H
#define CUAS_PETSCGRID_H

#include "petsc.h"

class PetscGrid {
 public:
  PetscGrid(int numOfCols, int numOfRows);
  ~PetscGrid();

  // indices: from 0 to LocalGhostWidth
  // you MUST call restoreLocal2dArr afterwards!
  PetscScalar **getAsLocal2dArr();

  // indices: from 0 to LocalNumOfCols
  // you MUST call setAs2dArr or restoreGlobal2dArr afterwards!
  PetscScalar **getAsGlobal2dArr();

  // call getAsGlobal2dArr before; parameter: changed return values.
  void setAsGlobal2dArr(PetscScalar **values);

  void restoreLocal2dArr(PetscScalar **values);
  void restoreGlobal2dArr(PetscScalar **values);

  void setGlobalVecAndUpdate(Vec glob);

  void setConst(PetscScalar value);
  void setZero() { setConst(0); }

  int countNonZero() const;

  /* getter */

  Vec getLocalVec() { return local; }
  Vec getGlobalVec() { return global; }

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

  DM getDM() { return dm; }

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
};

#endif
