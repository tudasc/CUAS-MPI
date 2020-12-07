#include "PetscGrid.h"

#include <iostream>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscGrid *grid = new PetscGrid(10, 10);

  grid->setZero();

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

  // they should all be the rank
  grid->viewGridWithGhost();
  grid->viewGridNoGhost();

  // all values should be 5
  grid->setConst(5);
  grid->viewGridNoGhost();
  grid->viewGridWithGhost();

  // all values should be 0
  grid->setZero();
  grid->viewGridNoGhost();
  grid->viewGridWithGhost();

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
