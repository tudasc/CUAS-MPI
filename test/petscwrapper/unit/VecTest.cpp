#include "PetscVec.h"

#include "petscdump.h"

#include <iostream>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscVec *vec = new PetscVec(5);
  vec->setZero();
  vec->setValue(0, 7);
  vec->assemble();
  dump(*vec);
  vec->setConst(10);
  dump(*vec);
  delete vec;
  PetscFinalize();
  return 0;
}
