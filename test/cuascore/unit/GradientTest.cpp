#include "specialgradient.h"

#include "petscdump.h"

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscGrid *grid = new PetscGrid(5, 5);
  Vec myglobal = grid->getGlobalVec();
  VecSet(myglobal, 0);
  grid->setGlobalVecAndUpdate(myglobal);
  PetscScalar **mylocalarr, **myglobarr;
  mylocalarr = grid->getAsLocal2dArr();
  myglobarr = grid->getAsGlobal2dArr();
  for (int i = 0; i < grid->getLocalNumOfRows(); ++i) {
    for (int j = 0; j < grid->getLocalNumOfCols(); ++j) {
      myglobarr[i][j] = rank + 3;
    }
  }
  grid->setAsGlobal2dArr(myglobarr);
  grid->restoreLocal2dArr(mylocalarr);
  PetscGrid *gradient_1 = new PetscGrid(5, 5);
  CUAS::gradient2(*grid, *gradient_1, 1.0);
  PetscGrid *gradient_2 = new PetscGrid(5, 5);
  CUAS::gradient2_central(*grid, *gradient_2, 1.0);
  dump(*gradient_1, false);
  dump(*gradient_2, false);
  delete gradient_1;
  delete gradient_2;
  delete grid;
  PetscFinalize();
  return 0;
}
