#include "fillModel.h"

#include "petsc.h"

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);

  CUAS::CUASModel *model = new CUAS::CUASModel(20, 10);
  fillNoData(model);
  model->init();

  delete model;
  PetscFinalize();
  return 0;
}
