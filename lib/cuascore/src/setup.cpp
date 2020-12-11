#include "setup.h"

#include "PetscAlgorithms.h"
#include "helper.h"
#include "physicalConstants.h"

#include "PetscGrid.h"

namespace CUAS {

void setup(CUASModel &model, CUASArgs const &args) {
  model.Sp->setZero();
  model.K->setConst(args.conductivity);

  auto bt = args.layerThickness;

  model.T->setConst(0.2);

  model.S->setConst(Ss * bt * args.Ssmulti);

  auto global_mask = model.no_flow_mask->getAsGlobal2dArr();
  auto K_arr = model.K->getAsGlobal2dArr();
  auto T_arr = model.T->getAsGlobal2dArr();
  auto mask = model.bnd_mask->getAsGlobal2dArr();

  auto rows = model.no_flow_mask->getLocalNumOfRows();
  auto cols = model.no_flow_mask->getLocalNumOfCols();

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (mask[i][j] == (PetscScalar)NOFLOW_FLAG) {
        global_mask[i][j] = true;
        K_arr[i][j] = NOFLOW_VALUE;
        T_arr[i][j] = NOFLOW_VALUE;
      } else {
        global_mask[i][j] = false;
      }
    }
  }

  model.no_flow_mask->setAsGlobal2dArr(global_mask);
  model.K->setAsGlobal2dArr(K_arr);

  auto T_n_arr = model.T_n->getAsGlobal2dArr();

  rows = model.T->getLocalNumOfRows();
  cols = model.T->getLocalNumOfCols();

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      T_n_arr[i][j] = T_arr[i][j];
    }
  }

  model.T_n->setAsGlobal2dArr(T_n_arr);
  model.T->setAsGlobal2dArr(T_arr);

  auto Q_arr = model.Q->getAsGlobal2dArr();

  rows = model.Q->getLocalNumOfRows();
  cols = model.Q->getLocalNumOfCols();

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      Q_arr[i][j] = (Q_arr[i][j]) / SPY * args.supplyMultiplier;
    }
  }

  model.Q->setAsGlobal2dArr(Q_arr);

  auto global_dir_mask = model.dirichlet_mask->getAsGlobal2dArr();

  rows = model.dirichlet_mask->getLocalNumOfRows();
  cols = model.dirichlet_mask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask[i][j] == (PetscScalar)DIRICHLET_FLAG) {
        global_dir_mask[i][j] = true;
      } else {
        global_dir_mask[i][j] = false;
      }
    }
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask[i][j] == (PetscScalar)DIRICHLET_LAKE_FLAG) {
        global_dir_mask[i][j] = true;
      }
    }
  }

  global_mask = model.no_flow_mask->getAsGlobal2dArr();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      global_dir_mask[i][j] = global_dir_mask[i][j] || global_mask[i][j];
    }
  }

  model.dirichlet_mask->setAsGlobal2dArr(global_dir_mask);
  model.no_flow_mask->restoreGlobal2dArr(global_mask);

  CUAS::pressure2head(*model.dirichlet_values, *model.p_ice, *model.topg, 0.0);

  auto global_sea_mask = model.sea_level_forcing_mask->getAsGlobal2dArr();

  rows = model.sea_level_forcing_mask->getLocalNumOfRows();
  cols = model.sea_level_forcing_mask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask[i][j] == (PetscScalar)DIRICHLET_FLAG) {
        global_sea_mask[i][j] = true;
      } else {
        global_sea_mask[i][j] = false;
      }
    }
  }

  model.sea_level_forcing_mask->setAsGlobal2dArr(global_sea_mask);
  model.bnd_mask->restoreGlobal2dArr(mask);

  binaryDialation(*model.no_flow_mask, *model.grad_mask);

  // TODO: Time dependent forcing (Interpolierung)
  //
  // TODO: Time dependent tidal forcing
}

}  // namespace CUAS