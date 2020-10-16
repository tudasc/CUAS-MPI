#include "helper.h"
#include "PetscGrid.h"

#include "petsc.h"

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscGrid *pressure = new PetscGrid(200, 200);
  PetscGrid *bed_elevation = new PetscGrid(200, 200);
  PetscGrid *head = new PetscGrid(200, 200);
  PetscScalar sea_level = 33.5;
  PetscScalar **pressure2d = pressure->getAsGlobal2dArr();
  PetscScalar **bed_elevation2d = bed_elevation->getAsGlobal2dArr();
  PetscScalar **head2d = head->getAsGlobal2dArr();
  for (int j = 0; j < pressure->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < pressure->getLocalNumOfCols(); ++i) {
      pressure2d[j][i] = sea_level * rank - 2.3;
      bed_elevation2d[j][i] = sea_level * sea_level - 31.3;
      head2d[j][i] = rank * j + (i * 35);
    }
  }
  pressure->setAsGlobal2dArr(pressure2d);
  bed_elevation->setAsGlobal2dArr(bed_elevation2d);
  head->setAsGlobal2dArr(head2d);
  CUAS::pressure2head(head, pressure, bed_elevation, sea_level);
  CUAS::pressure2head(pressure, head, bed_elevation, sea_level);
  CUAS::overburdenPressure(pressure, head);
  delete pressure;
  delete bed_elevation;
  delete head;
  PetscFinalize();
  return 0;
}
