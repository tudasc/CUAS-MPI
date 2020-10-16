#include "CUASModel.h"
#include "helper.h"
#include "PetscGrid.h"
#include "fillModel.h"
#include "fill_matrix_coo.h"

#include "petsc.h"

#include <math.h>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  CUAS::CUASModel *model = new CUAS::CUASModel(20, 10);
  fillNoData(model);
  // not sure if init is needed
  model->init();
  PetscMat *matToBeFilled = new PetscMat(200, 200);
  PetscGrid *Se = new PetscGrid(20, 10);
  PetscScalar **SeGlobal = Se->getAsGlobal2dArr();
  for (int j = 0; j < Se->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Se->getLocalNumOfCols(); ++i) {
      SeGlobal[j][i] = 9.82977696 * pow(10, -5);  // normally: physconst.Ss * args.layerthickness * args.Ssmulti
    }
  }
  Se->setAsGlobal2dArr(SeGlobal);
  PetscGrid *noFlow_mask = new PetscGrid(20, 10);
  PetscScalar **noFlow_maskGlobal = noFlow_mask->getAsGlobal2dArr();
  PetscScalar **bnd_mask2d = model->bnd_mask->getAsGlobal2dArr();
  for (int j = 0; j < noFlow_mask->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < noFlow_mask->getLocalNumOfCols(); ++i) {
      if (bnd_mask2d[j][i] == NOFLOW_FLAG) {
        noFlow_maskGlobal[j][i] = true;
      } else {
        noFlow_maskGlobal[j][i] = false;
      }
    }
  }
  noFlow_mask->setAsGlobal2dArr(noFlow_maskGlobal);
  noFlow_maskGlobal = noFlow_mask->getAsGlobal2dArr();
  PetscGrid *Teff = new PetscGrid(20, 10);
  PetscScalar **TeffGlobal = Teff->getAsGlobal2dArr();
  for (int j = 0; j < Teff->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Teff->getLocalNumOfCols(); ++i) {
      if (noFlow_maskGlobal[j][i]) {
        TeffGlobal[j][i] = NOFLOW_VALUE;
      } else {
        TeffGlobal[j][i] = 0.2;
      }
    }
  }
  Teff->setAsGlobal2dArr(TeffGlobal);
  // u_n = select_initial_head... --> standard: nzero --> pressure2head
  PetscGrid *u_n = new PetscGrid(20, 10);
  CUAS::pressure2head(u_n, model->p_ice, model->topg, 0.0);
  PetscGrid *dirichlet_mask = new PetscGrid(20, 10);
  PetscScalar **dirichlet_maskGlobal = dirichlet_mask->getAsGlobal2dArr();
  bool dirich_mask1;
  for (int j = 0; j < dirichlet_mask->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < dirichlet_mask->getLocalNumOfCols(); ++i) {
      if (bnd_mask2d[j][i] == DIRICHLET_FLAG) {
        dirich_mask1 = true;
      } else {
        dirich_mask1 = false;
      }
      dirich_mask1 = dirich_mask1 || (bnd_mask2d[j][i] == DIRICHLET_LAKE_FLAG);
      dirich_mask1 = dirich_mask1 || noFlow_maskGlobal[j][i];
      dirichlet_maskGlobal[j][i] = dirich_mask1;
    }
  }
  dirichlet_mask->setAsGlobal2dArr(dirichlet_maskGlobal);
  model->bnd_mask->restoreGlobal2dArr(bnd_mask2d);
  noFlow_mask->restoreGlobal2dArr(noFlow_maskGlobal);
  // stolen from output from nodata test
  PetscGrid *WertfuerQ = new PetscGrid(20, 10);
  PetscScalar **WertfuerQGlobal = WertfuerQ->getAsGlobal2dArr();
  for (int j = 0; j < WertfuerQ->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < WertfuerQ->getLocalNumOfCols(); ++i) {
      WertfuerQGlobal[j][i] = 3.17 * pow(10, -8);  // normally: (args.conductivity * args.layerthickness)^args.exponent
    }
  }
  WertfuerQ->setAsGlobal2dArr(WertfuerQGlobal);
  PetscScalar dt = 43200;
  PetscScalar theta = 1;
  // different values
  PetscGrid *dirichlet_values = new PetscGrid(20, 10);
  CUAS::pressure2head(dirichlet_values, model->p_ice, model->topg, 0.0);
  PetscVec *b = new PetscVec(model->Ncols * model->Nrows);
  CUAS::fill_matrix_coo(matToBeFilled, b, model->Nrows, model->Nrows, Se, Teff, model->dx, dt, theta, u_n, WertfuerQ,
                        dirichlet_values, dirichlet_mask);
  matToBeFilled->viewGlobal();
  b->view();
  PetscFinalize();
  return 0;
}
