/**
 * File: systemmatrixTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "systemmatrixTest.h"

#include "CUASConstants.h"
#include "Forcing/SteadyForcing.h"
#include "fillgrid.h"
#include "systemmatrix.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

TEST(fillMatrixTest, randomValues) {
  if (mpiSize != 1) {
    return;
  }

  ASSERT_EQ(mpiSize, 1);

  int dx = 10;
  int dtsecs = 40278;
  int theta = 3;
  PETScGrid SeGrid(NX, NY);
  PETScGrid TeffPowGrid(NX, NY);
  PETScGrid u_nGrid(NX, NY);
  PETScGrid QGrid(NX, NY);
  PETScGrid dirichValGrid(NX, NY);
  PETScGrid bndMaskGrid(NX, NY);
  {
    auto SeWrite = SeGrid.getWriteHandle();
    auto TeffPowWrite = TeffPowGrid.getWriteHandle();
    auto u_nWrite = u_nGrid.getWriteHandle();
    auto QWrite = QGrid.getWriteHandle();
    auto dirichValWrite = dirichValGrid.getWriteHandle();
    auto bndMaskWrite = bndMaskGrid.getWriteHandle();

    for (int i = 0; i < NY; ++i) {
      for (int j = 0; j < NX; ++j) {
        SeWrite(i, j) = Se[i][j];
        TeffPowWrite(i, j) = TeffPow[i][j];
        u_nWrite(i, j) = u_n[i][j];
        QWrite(i, j) = current_Q[i][j];
        dirichValWrite(i, j) = dirich_val[i][j];
        bndMaskWrite(i, j) = bnd_mask[i][j];
      }
    }
  }

  PETScGrid globalIndices(NX, NY);
  fillGlobalIndicesBlocked(globalIndices);

  std::unique_ptr<CUAS::Forcing> forcing = std::make_unique<CUAS::SteadyForcing>(QGrid);

  PETScGrid gridB(NX, NY);
  Mat petscA;
  DMCreateMatrix(gridB.dm, &petscA);
  PETScMatrix matA(petscA);

  CUAS::systemmatrix(matA, gridB, SeGrid, TeffPowGrid, dx, dtsecs, theta, u_nGrid, forcing->getCurrent(0),
                     dirichValGrid, bndMaskGrid, globalIndices);

  // compare nnz of Matrix A with nnz of python
  {
    std::vector<double> APetscNonZero;
    {
      PetscScalar v;
      for (int i = 0; i < matA.getNumberOfCols(); ++i) {
        for (int j = 0; j < matA.getNumberOfRows(); ++j) {
          MatGetValues(matA.getRaw(), 1, &i, 1, &j, &v);
          if (v != 0) {
            APetscNonZero.push_back(v);
          }
        }
      }
      ASSERT_EQ(ANonZeroPython.size(), APetscNonZero.size());
    }

    // compare values of Matrix A
    /*for (int i = 0; i < ANonZeroPython.size(); ++i) {
      ASSERT_NEAR(ANonZeroPython[i], APetscNonZero[i], 0.1);
    }*/
  }

  // compare results of Vector b with result of python
  /*{
    PetscScalar val;
    for (int i = 0; i < VecB.getSize(); ++i) {
      VecGetValues(VecB.getRaw(), 1, &i, &val);
      ASSERT_NEAR(val, b[i], 0.1);
    }
  }*/

  // check equality of sum(A*b)
  {
    PETScGrid productGrid(NX, NY);
    MatMult(petscA, gridB.global, productGrid.global);
    PetscScalar sum;
    VecSum(productGrid.global, &sum);
    ASSERT_DOUBLE_EQ(sum, 6708180545601.590820);
  }
}

TEST(fillMatrixTest, randomValues3x3) {
  constexpr int gridedge = 3;

  int dx = 10;
  int dtsecs = 40278;
  int theta = 3;
  PETScGrid SeGrid(gridedge, gridedge);
  PETScGrid TeffPowGrid(gridedge, gridedge);
  PETScGrid u_nGrid(gridedge, gridedge);
  PETScGrid QGrid(gridedge, gridedge);
  PETScGrid dirichValGrid(gridedge, gridedge);
  PETScGrid bndMaskGrid(gridedge, gridedge);

  {
    int cornerX = SeGrid.getCornerX();
    int cornerY = SeGrid.getCornerY();

    auto SeWrite = SeGrid.getWriteHandle();
    auto TeffPowWrite = TeffPowGrid.getWriteHandle();
    auto u_nWrite = u_nGrid.getWriteHandle();
    auto QWrite = QGrid.getWriteHandle();
    auto dirichValWrite = dirichValGrid.getWriteHandle();
    auto bndMaskWrite = bndMaskGrid.getWriteHandle();

    for (int row = 0; row < SeGrid.getLocalNumOfRows(); ++row) {
      for (int col = 0; col < SeGrid.getLocalNumOfCols(); ++col) {
        SeWrite(row, col) = Se[cornerY + row][cornerX + col];
        TeffPowWrite(row, col) = TeffPow[cornerY + row][cornerX + col];
        u_nWrite(row, col) = u_n[cornerY + row][cornerX + col];
        QWrite(row, col) = current_Q[cornerY + row][cornerX + col];
        dirichValWrite(row, col) = dirich_val[cornerY + row][cornerX + col];
        bndMaskWrite(row, col) = bnd_mask[cornerY + row][cornerX + col];
      }
    }
  }

  bndMaskGrid.setRealBoundary(NOFLOW_FLAG);

  PETScGrid globalIndices(gridedge, gridedge);
  fillGlobalIndicesBlocked(globalIndices);

  std::unique_ptr<CUAS::Forcing> forcing = std::make_unique<CUAS::SteadyForcing>(QGrid);

  PETScGrid gridB(gridedge, gridedge);
  Mat petscA;
  DMCreateMatrix(gridB.dm, &petscA);
  PETScMatrix matA(petscA);

  CUAS::systemmatrix(matA, gridB, SeGrid, TeffPowGrid, dx, dtsecs, theta, u_nGrid, forcing->getCurrent(0),
                     dirichValGrid, bndMaskGrid, globalIndices);

  // check equality of sum(A*b)
  {
    PETScGrid product(gridedge, gridedge);
    MatMult(petscA, gridB.global, product.global);
    PetscScalar sum;
    VecSum(product.global, &sum);
    ASSERT_DOUBLE_EQ(sum, -1339558223.324310);
  }
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
