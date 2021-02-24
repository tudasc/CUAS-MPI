#include "CUASKernels.h"
#include "CUASModel.h"
#include "fillNoData.h"
#include "fill_matrix_coo.h"

#include "PETScGrid.h"

#include "gtest/gtest.h"

#include <math.h>

int mpiRank;
int mpiSize;

#define MPI_SIZE 4

// TODO: verify
/*TEST(fillMatrixTest, result) {
  ASSERT_EQ(mpiSize, MPI_SIZE);
  auto pmodel = fillNoData();
  auto &model = *pmodel;
  // not sure if init is needed
  model.init();
  PETScMat matToBeFilled(200, 200);

  PETScGrid Se(20, 10);
  Se.setConst(9.82977696 * pow(10, -5));  // normally: physconst.Ss * args.layerthickness * args.Ssmulti
  PETScGrid noFlow_mask(20, 10);
  auto noFlow_maskGlobal = noFlow_mask.getWriteHandle();
  auto &bnd_mask2d = model.bndMask->getReadHandle();
  for (int j = 0; j < noFlow_mask.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < noFlow_mask.getLocalNumOfCols(); ++i) {
      if (bnd_mask2d(j, i) == NOFLOW_FLAG) {
        noFlow_maskGlobal(j, i) = true;
      } else {
        noFlow_maskGlobal(j, i) = false;
      }
    }
  }
  noFlow_maskGlobal.setValues();

  auto &noFlow_maskRead = noFlow_mask.getReadHandle();
  PETScGrid Teff(20, 10);
  auto TeffGlobal = Teff.getWriteHandle();
  for (int j = 0; j < Teff.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < Teff.getLocalNumOfCols(); ++i) {
      if (noFlow_maskRead(j, i)) {
        TeffGlobal(j, i) = NOFLOW_VALUE;
      } else {
        TeffGlobal(j, i) = 0.2;
      }
    }
  }
  TeffGlobal.setValues();

  // u_n = select_initial_head... --> standard: nzero --> pressure2head
  PETScGrid u_n(20, 10);
  CUAS::pressure2head(u_n, *model.pIce, *model.topg, 0.0);
  PETScGrid dirichlet_mask(20, 10);
  auto dirichlet_maskGlobal = dirichlet_mask.getWriteHandle();
  bool dirich_mask1;
  for (int j = 0; j < dirichlet_mask.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < dirichlet_mask.getLocalNumOfCols(); ++i) {
      if (bnd_mask2d(j, i) == DIRICHLET_FLAG) {
        dirich_mask1 = true;
      } else {
        dirich_mask1 = false;
      }
      dirich_mask1 = dirich_mask1 || (bnd_mask2d(j, i) == DIRICHLET_LAKE_FLAG);
      dirich_mask1 = dirich_mask1 || noFlow_maskGlobal(j, i);
      dirichlet_maskGlobal(j, i) = dirich_mask1;
    }
  }
  dirichlet_maskGlobal.setValues();
  // stolen from output of nodata test
  PETScGrid ValueQ(20, 10);
  ValueQ.setConst(3.17 * pow(10, -8));

  PetscScalar dt = 43200;
  PetscScalar theta = 1;
  // different values
  PETScGrid dirichlet_values(20, 10);
  CUAS::pressure2head(dirichlet_values, *model.pIce, *model.topg, 0.0);
  PETScVec b(model.Ncols * model.Nrows);
  CUAS::fill_matrix_coo(matToBeFilled, b, model.Nrows, model.Ncols, Se, Teff, model.dx, dt, theta, u_n, ValueQ,
                        dirichlet_values, dirichlet_mask);

  // check single values (compared to original python)
  PetscScalar v;
  int p1[] = {11};
  int p2[] = {1};

  if (mpiRank == 0) {
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -8.789614337456268e-12);

    p1[0] = 11;
    p2[0] = 10;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -8.789614337456268e-12);

    p1[0] = 11;
    p2[0] = 11;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 176.79233913932018);

    p1[0] = 11;
    p2[0] = 21;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -87.89616956965129);

    p1[0] = 13;
    p2[0] = 3;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -8.789614337456268e-12);

    p1[0] = 13;
    p2[0] = 12;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, -87.89616956965129);

    p1[0] = 13;
    p2[0] = 13;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 264.6885087089627);

    p1[0] = 0;
    p2[0] = 0;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 1.0);

    p1[0] = 9;
    p2[0] = 9;
    MatGetValues(matToBeFilled.getRaw(), 1, p1, 1, p2, &v);
    ASSERT_EQ(v, 1.0);
  }
}*/

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
