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

  Vec getLocalVec() { return local; }
  Vec getGlobalVec() { return global; }

  void setGlobalVecAndUpdate(Vec glob);

  void viewGridWithGhost() const;
  void viewGridNoGhost() const;

  int countNonZero() const;

  int getLocalNumOfCols() const { return localNumOfCols; }
  int getLocalGhostNumOfCols() const { return localNumOfCols; }
  int getLocalNumOfRows() const { return localNumOfRows; }
  int getLocalGhostNumOfRows() const { return localGhostNumOfRows; }
  int getTotalNumOfRows() const { return totalNumOfRows; }
  int getTotalNumOfCols() const { return totalNumOfCols; }

  DM getDM() { return dm; }

  int getCornerX() const { return cornerX; }
  int getCornerY() const { return cornerY; }

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
