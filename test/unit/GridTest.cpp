#include "PetscGrid.h"

#include <iostream>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscGrid *grid = new PetscGrid(10, 10);

  Vec myglobal = grid->getGlobalVec();
  VecSet(myglobal, 0);
  grid->setGlobalVecAndUpdate(myglobal);

  PetscScalar **mylocalarr, **myglobarr;
  mylocalarr = grid->getAsLocal2dArr();
  myglobarr = grid->getAsGlobal2dArr();
  for (int j = 0; j < grid->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < grid->getLocalNumOfCols(); ++i) {
      myglobarr[j][i] = mylocalarr[j][i] + rank;
    }
  }
  grid->setAsGlobal2dArr(myglobarr);
  grid->restoreLocal2dArr(mylocalarr);

  // ignore - just to test restore
  myglobarr = grid->getAsGlobal2dArr();
  grid->restoreGlobal2dArr(myglobarr);

  // only makes sense to view one at a time.
  grid->viewGridWithGhost();
  // grid->viewGridNoGhost();

  int nonZero = grid->countNonZero();
  if (rank == 0) {
    std::cout << nonZero << std::endl;
  }

  // test getters from header
  DM myGrid = grid->getDM();
  DMView(myGrid, PETSC_VIEWER_STDOUT_WORLD);
  if (rank == 0) {
    VecView(grid->getLocalVec(), PETSC_VIEWER_STDOUT_SELF);
    std::cout << grid->getLocalNumOfCols() << grid->getLocalGhostNumOfCols() << grid->getLocalNumOfRows()
              << grid->getLocalGhostNumOfRows() << grid->getTotalNumOfRows() << grid->getTotalNumOfCols()
              << grid->getCornerX() << grid->getCornerY() << std::endl;
  }

  delete grid;
  PetscFinalize();
  return 0;
}
