#include "systemmatrixTest.h"

#include "Forcing/ConstantForcing.h"
#include "systemmatrix.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 1

TEST(fillMatrixTest, randomValues) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  int dx = 10;
  int dtsecs = 40278;
  int theta = 3;
  PETScMatrix MatA(NY * NX, NY * NX);
  PETScVector VecB(NY * NX);
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
  // using SPY here because the dumped Q was from currentQ after applying forcing
  std::unique_ptr<CUAS::Forcing> forcing = std::make_unique<CUAS::ConstantForcing>(QGrid);
  CUAS::systemmatrix(MatA, VecB, NY, NX, SeGrid, TeffPowGrid, dx, dtsecs, theta, u_nGrid, forcing->getCurrentQ(),
                     dirichValGrid, bndMaskGrid);

  std::vector<double> APetscNonZero;

  // compare nnz of Matrix A
  {
    PetscScalar v;
    for (int i = 0; i < MatA.getNumberOfCols(); ++i) {
      for (int j = 0; j < MatA.getNumberOfRows(); ++j) {
        MatGetValues(MatA.getRaw(), 1, &i, 1, &j, &v);
        if (v != 0) {
          APetscNonZero.push_back(v);
        }
      }
    }
    ASSERT_EQ(ANonZeroPython.size(), APetscNonZero.size());
  }

  // compare values of Matrix A
  for (int i = 0; i < ANonZeroPython.size(); ++i) {
    ASSERT_NEAR(ANonZeroPython[i], APetscNonZero[i], 0.1);
  }

  // compare results of Vector b
  {
    PetscScalar val;
    int p;
    for (int i = 0; i < VecB.getSize(); ++i) {
      VecGetValues(VecB.getRaw(), 1, &i, &val);
      ASSERT_NEAR(val, b[i], 0.1);
    }
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
