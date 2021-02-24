#include "CUASSetup.h"

#include "CUASConstants.h"
#include "CUASKernels.h"

#include "PETScGrid.h"

namespace CUAS {

void setup(CUASModel &model, CUASArgs const &args) {
  model.Sp->setZero();
  model.K->setConst(args.conductivity);

  auto bt = args.layerThickness;

  model.T->setConst(0.2);

  model.S->setConst(Ss * bt * args.Ssmulti);

  auto global_mask = model.noFlowMask->getWriteHandle();
  auto K_arr = model.K->getWriteHandle();
  auto T_arr = model.T->getWriteHandle();
  auto &mask = model.bndMask->getReadHandle();

  auto rows = model.noFlowMask->getLocalNumOfRows();
  auto cols = model.noFlowMask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)NOFLOW_FLAG) {
        global_mask(i, j) = true;
        K_arr(i, j) = NOFLOW_VALUE;
        T_arr(i, j) = NOFLOW_VALUE;
      } else {
        global_mask(i, j) = false;
      }
    }
  }

  // needed for function call later on
  global_mask.setValues();
  T_arr.setValues();

  model.T_n->copy(*model.T);

  auto Q_arr = model.Q->getWriteHandle();

  rows = model.Q->getLocalNumOfRows();
  cols = model.Q->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      Q_arr(i, j) = (Q_arr(i, j)) / SPY * args.supplyMultiplier;
    }
  }

  auto global_dir_mask = model.dirichletMask->getWriteHandle();

  rows = model.dirichletMask->getLocalNumOfRows();
  cols = model.dirichletMask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)DIRICHLET_FLAG) {
        global_dir_mask(i, j) = true;
      } else {
        global_dir_mask(i, j) = false;
      }
    }
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)DIRICHLET_LAKE_FLAG) {
        global_dir_mask(i, j) = true;
      }
    }
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      global_dir_mask(i, j) = global_dir_mask(i, j) || global_mask(i, j);
    }
  }

  pressure2head(*model.dirichletValues, *model.pIce, *model.topg, 0.0);

  auto global_sea_mask = model.seaLevelForcingMask->getWriteHandle();

  rows = model.seaLevelForcingMask->getLocalNumOfRows();
  cols = model.seaLevelForcingMask->getLocalNumOfCols();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (mask(i, j) == (PetscScalar)DIRICHLET_FLAG) {
        global_sea_mask(i, j) = true;
      } else {
        global_sea_mask(i, j) = false;
      }
    }
  }

  binaryDialation(*model.gradMask, *model.noFlowMask);

  // TODO: Time dependent forcing (Interpolierung)
  //
  // TODO: Time dependent tidal forcing
}

}  // namespace CUAS
