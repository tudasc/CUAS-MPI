#include "CUASFile.h"
#include "Logger.h"
#include "ModelReader.h"
#include "SolutionHandler.h"
#include "fillNoData.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define FILENAME "test.nc"

TEST(readTest, writeToNetcdf) {
  PETScGrid u(NODATA_COLS, NODATA_ROWS);
  auto u2d = u.getWriteHandle();
  for (int i = 0; i < u.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < u.getLocalNumOfCols(); ++j) {
      u2d(i, j) = mpiRank + 222;
    }
  }

  PETScGrid u_n(NODATA_COLS, NODATA_ROWS);
  auto u_n2d = u_n.getWriteHandle();
  for (int i = 0; i < u_n.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < u_n.getLocalNumOfCols(); ++j) {
      u_n2d(i, j) = mpiRank + 333;
    }
  }

  PETScGrid melt(NODATA_COLS, NODATA_ROWS);
  auto melt2d = melt.getWriteHandle();
  for (int i = 0; i < melt.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < melt.getLocalNumOfCols(); ++j) {
      melt2d(i, j) = mpiRank + 1 * 3.14;
    }
  }

  PetscScalar cavityOpening = 10.0;
  int Nt = 100;
  int saveEvery = 20;
  CUAS::CUASArgs args;
  auto model = fillNoData();
  CUAS::SolutionHandler handler(FILENAME, Nt, saveEvery, NODATA_COLS, NODATA_ROWS, mpiRank);

  for (int i = 1; i < Nt + 1; ++i) {
    if (i % saveEvery == 0) {
      handler.saveSolution(i, args, mpiRank, u, u_n, *model, melt, cavityOpening);
    }
  }

  // test results
  CUAS::CUASFile file(FILENAME, 'r');

  // test u
  PETScGrid uResult(NODATA_COLS, NODATA_ROWS);
  file.read("u20", uResult);
  auto uResult2d = uResult.getReadHandle();
  auto uActual = u.getReadHandle();

  for (int i = 0; i < uResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < uResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(uResult2d(i, j), uActual(i, j));
    }
  }

  // test u_n
  PETScGrid u_nResult(NODATA_COLS, NODATA_ROWS);
  file.read("u_n20", u_nResult);
  auto u_nResult2d = u_nResult.getReadHandle();
  auto u_nActual = u_n.getReadHandle();

  for (int i = 0; i < u_nResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < u_nResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(u_nResult2d(i, j), u_nActual(i, j));
    }
  }

  // test model
  // test usurf
  PETScGrid usurfResult(NODATA_COLS, NODATA_ROWS);
  file.read("usurf20", usurfResult);
  auto usurfResult2d = usurfResult.getReadHandle();
  auto usurfActual = model->usurf->getReadHandle();

  for (int i = 0; i < usurfResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < usurfResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(usurfResult2d(i, j), usurfActual(i, j));
    }
  }

  // test topg
  PETScGrid topgResult(NODATA_COLS, NODATA_ROWS);
  file.read("topg20", topgResult);
  auto topgResult2d = topgResult.getReadHandle();
  auto topgActual = model->topg->getReadHandle();

  for (int i = 0; i < topgResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < topgResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(topgResult2d(i, j), topgActual(i, j));
    }
  }

  // test thk
  PETScGrid thkResult(NODATA_COLS, NODATA_ROWS);
  file.read("thk20", thkResult);
  auto thkResult2d = thkResult.getReadHandle();
  auto thkActual = model->thk->getReadHandle();

  for (int i = 0; i < thkResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < thkResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(thkResult2d(i, j), thkActual(i, j));
    }
  }

  // test bndMask
  PETScGrid bndMaskResult(NODATA_COLS, NODATA_ROWS);
  file.read("bndMask20", bndMaskResult);
  auto bndMaskResult2d = bndMaskResult.getReadHandle();
  auto bndMaskActual = model->bndMask->getReadHandle();

  for (int i = 0; i < bndMaskResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < bndMaskResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(bndMaskResult2d(i, j), bndMaskActual(i, j));
    }
  }

  // test pIce
  PETScGrid pIceResult(NODATA_COLS, NODATA_ROWS);
  file.read("pIce20", pIceResult);
  auto pIceResult2d = pIceResult.getReadHandle();
  auto pIceActual = model->pIce->getReadHandle();

  for (int i = 0; i < pIceResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < pIceResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(pIceResult2d(i, j), pIceActual(i, j));
    }
  }

  // test melt
  PETScGrid meltResult(NODATA_COLS, NODATA_ROWS);
  file.read("melt20", meltResult);
  auto meltResult2d = meltResult.getReadHandle();
  auto meltActual = melt.getReadHandle();

  for (int i = 0; i < meltResult.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < meltResult.getLocalNumOfCols(); ++j) {
      ASSERT_EQ(meltResult2d(i, j), meltActual(i, j));
    }
  }

  // test cavityOpening
  PetscScalar cavityResult;
  file.read("cavity_opening20", cavityResult);
  ASSERT_EQ(cavityResult, cavityOpening);
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
