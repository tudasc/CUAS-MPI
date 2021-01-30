#include "CUASModel.h"
#include "PetscGrid.h"
#include "fillModel.h"
#include "fill_matrix_coo.h"
#include "helper.h"
#include "petscdump.h"

#include <math.h>

#include "gtest/gtest.h"

#include <iostream>
int mpiRank;
int mpiSize;

#define MPI_SIZE 4

TEST(fillMatrixTest, result) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  CUAS::CUASModel model(20, 10);
  fillNoData(model);
  // not sure if init is needed
  model.init();
  PetscMat matToBeFilled(200, 200);

  PetscGrid Se(20, 10);
  Se.setConst(9.82977696 * pow(10, -5));  // normally: physconst.Ss * args.layerthickness * args.Ssmulti
  PetscGrid noFlow_mask(20, 10);
  auto noFlow_maskGlobal = noFlow_mask.getAsGlobal2dArr();
  auto bnd_mask2d = model.bnd_mask->getAsGlobal2dArr();
  for (int j = 0; j < noFlow_mask.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < noFlow_mask.getLocalNumOfCols(); ++i) {
      if (bnd_mask2d[j][i] == NOFLOW_FLAG) {
        noFlow_maskGlobal[j][i] = true;
      } else {
        noFlow_maskGlobal[j][i] = false;
      }
    }
  }
  noFlow_mask.setAsGlobal2dArr(noFlow_maskGlobal);

  noFlow_maskGlobal = noFlow_mask.getAsGlobal2dArr();
  PetscGrid Teff(20, 10);
  auto TeffGlobal = Teff.getAsGlobal2dArr();
  for (int j = 0; j < Teff.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Teff.getLocalNumOfCols(); ++i) {
      if (noFlow_maskGlobal[j][i]) {
        TeffGlobal[j][i] = NOFLOW_VALUE;
      } else {
        TeffGlobal[j][i] = 0.2;
      }
    }
  }
  Teff.setAsGlobal2dArr(TeffGlobal);

  // u_n = select_initial_head... --> standard: nzero --> pressure2head
  PetscGrid u_n(20, 10);
  CUAS::pressure2head(u_n, *model.p_ice, *model.topg, 0.0);
  PetscGrid dirichlet_mask(20, 10);
  auto dirichlet_maskGlobal = dirichlet_mask.getAsGlobal2dArr();
  bool dirich_mask1;
  for (int j = 0; j < dirichlet_mask.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < dirichlet_mask.getLocalNumOfCols(); ++i) {
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
  dirichlet_mask.setAsGlobal2dArr(dirichlet_maskGlobal);
  model.bnd_mask->restoreGlobal2dArr(bnd_mask2d);
  noFlow_mask.restoreGlobal2dArr(noFlow_maskGlobal);
  // stolen from output of nodata test
  PetscGrid ValueQ(20, 10);
  ValueQ.setConst(3.17 * pow(10, -8));

  PetscScalar dt = 43200;
  PetscScalar theta = 1;
  // different values
  PetscGrid dirichlet_values(20, 10);
  CUAS::pressure2head(dirichlet_values, *model.p_ice, *model.topg, 0.0);
  PetscVec b(model.Ncols * model.Nrows);
  CUAS::fill_matrix_coo(matToBeFilled, b, model.Nrows, model.Ncols, Se, Teff, model.dx, dt, theta, u_n, ValueQ,
                        dirichlet_values, dirichlet_mask);

  // check single values (compared to original python)
  PetscScalar v;
  int p1[] = {11};
  int p2[] = {1};

  if (mpiRank == 0) {
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -8.789614337456268e-12);

    p1[0] = 11;
    p2[0] = 10;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -8.789614337456268e-12);

    p1[0] = 11;
    p2[0] = 11;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 176.79233913932018);

    p1[0] = 11;
    p2[0] = 21;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -87.89616956965129);

    p1[0] = 13;
    p2[0] = 3;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -8.789614337456268e-12);

    p1[0] = 13;
    p2[0] = 12;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -87.89616956965129);

    p1[0] = 13;
    p2[0] = 13;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 264.6885087089627);

    p1[0] = 0;
    p2[0] = 0;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 1.0);

    p1[0] = 9;
    p2[0] = 9;
    MatGetValues(matToBeFilled.getPetscRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 1.0);
  }
}

int main(int argc, char *argv[]) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  result = RUN_ALL_TESTS();
  PetscFinalize();

  return result;
}
