#include "setup.h"
#include "fillModel.h"
#include "parseCxxopts.h"

#include "PetscGrid.h"
#include "petscdump.h"

#include <iostream>
#include <memory>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);

  {
    auto model = std::make_unique<CUAS::CUASModel>(10, 10);
    fillNoData(*model);
    model->init();

    CUAS::CUASArgs args;
    CUAS::parseArgs(argc, argv, args);

    CUAS::setup(*model, args);

    std::cout << "SP" << std::endl;
    dump(*model->Sp, false);

    std::cout << "K" << std::endl;
    dump(*model->K, false);

    std::cout << "T" << std::endl;
    dump(*model->T, false);

    std::cout << "S" << std::endl;
    dump(*model->S, false);

    std::cout << "no_flow_mask" << std::endl;
    dump(*model->no_flow_mask, false);

    std::cout << "T_n" << std::endl;
    dump(*model->T_n, false);

    std::cout << "sea_level_forcing_mask" << std::endl;
    dump(*model->sea_level_forcing_mask, false);

    std::cout << "grad_mask" << std::endl;
    dump(*model->grad_mask, false);

    std::cout << "dirichlet_mask" << std::endl;
    dump(*model->dirichlet_mask, false);

    std::cout << "dirichlet_values" << std::endl;
    dump(*model->dirichlet_values, false);

    std::cout << "Q from model" << std::endl;
    dump(*model->Q, false);

    std::cout << "bnd_mask from model" << std::endl;
    dump(*model->bnd_mask, false);
  }

  PetscFinalize();

  return 0;
}
