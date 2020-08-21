#include "PetscVec.h"

#include "petsc.h"
#include <iostream>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscVec *vec = new PetscVec(5);
  vec->zero();
  vec->setValue(0, 7);
  vec->assemble();
  vec->view();
  delete vec;
  PetscFinalize();
  return 0;
}
