#include "PetscMat.h"

#include "petscdump.h"

#include <iostream>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscMat *mat = new PetscMat(5, 5);

  int numberCols = mat->getCols();
  int numberRows = mat->getRows();
  std::cout << "Rows, Cols" << numberRows << " " << numberCols << std::endl;
  for (int i = 0; i < 5; i++) {
    mat->setValue(i, i, 100);
  }
  mat->assemble();
  dump(*mat, true);
  mat->setValue(3, 2, 500);
  mat->assemble();
  dump(*mat, true);
  mat->assemble();
  //dump(*mat, false);
  delete mat;
}
