#include "fillNoData.h"

#include "CUASConstants.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 6

//#define NODATA_COLS 20
//#define NODATA_ROWS 10

TEST(fillNoDataTest, fillNoData) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  auto model = fillNoData();
  model->init();

  ASSERT_EQ(model->Ncols, NODATA_COLS);
  ASSERT_EQ(model->Nrows, NODATA_ROWS);

  ASSERT_EQ(model->dx, 1000.0);
  ASSERT_EQ(model->dy, 1000.0);

  ASSERT_EQ(model->cols[0], 0.0);
  ASSERT_EQ(model->cols[1], 1000.0);
  ASSERT_EQ(model->cols[13], 13000.0);

  ASSERT_EQ(model->rows[0], 0.0);
  ASSERT_EQ(model->rows[1], 1000.0);
  ASSERT_EQ(model->rows[7], 7000.0);

  ASSERT_TRUE(model->usurf->isCompatible(*model->topg));
  ASSERT_TRUE(model->usurf->isCompatible(*model->thk));
  ASSERT_TRUE(model->usurf->isCompatible(*model->bndMask));
  ASSERT_TRUE(model->usurf->isCompatible(*model->Q));
  ASSERT_TRUE(model->usurf->isCompatible(*model->pIce));

  ASSERT_EQ(model->usurf->getLocalNumOfRows(), 4);
  ASSERT_EQ(model->usurf->getLocalNumOfCols(), 6);
  ASSERT_EQ(model->usurf->getLocalGhostNumOfRows(), 6);
  ASSERT_EQ(model->usurf->getLocalGhostNumOfCols(), 8);
  ASSERT_EQ(model->usurf->getTotalNumOfRows(), 8);
  ASSERT_EQ(model->usurf->getTotalNumOfCols(), 18);
  ASSERT_EQ(model->usurf->getTotalGhostNumOfRows(), 10);
  ASSERT_EQ(model->usurf->getTotalGhostNumOfCols(), 20);

  {
    auto &usurfHandle = model->usurf->getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(usurfHandle(0, 0, GHOSTED), 2000);
      ASSERT_EQ(usurfHandle(0, 1, GHOSTED), 1900);
      ASSERT_EQ(usurfHandle(0, 2, GHOSTED), 1800);
      ASSERT_EQ(usurfHandle(1, 0, GHOSTED), 2000);
      ASSERT_EQ(usurfHandle(1, 1, GHOSTED), 1900);
      ASSERT_EQ(usurfHandle(1, 2, GHOSTED), 1800);
      ASSERT_EQ(usurfHandle(2, 0, GHOSTED), 2000);
      ASSERT_EQ(usurfHandle(2, 1, GHOSTED), 1900);
      ASSERT_EQ(usurfHandle(2, 2, GHOSTED), 1800);
    } else if (mpiRank == 1) {
      ASSERT_EQ(usurfHandle(0, 0, GHOSTED), 1400);
      ASSERT_EQ(usurfHandle(0, 1, GHOSTED), 1300);
      ASSERT_EQ(usurfHandle(0, 2, GHOSTED), 1200);
      ASSERT_EQ(usurfHandle(1, 0, GHOSTED), 1400);
      ASSERT_EQ(usurfHandle(1, 1, GHOSTED), 1300);
      ASSERT_EQ(usurfHandle(1, 2, GHOSTED), 1200);
      ASSERT_EQ(usurfHandle(2, 0, GHOSTED), 1400);
      ASSERT_EQ(usurfHandle(2, 1, GHOSTED), 1300);
      ASSERT_EQ(usurfHandle(2, 2, GHOSTED), 1200);
    } else if (mpiRank == 5) {
      ASSERT_EQ(usurfHandle(3, 5, GHOSTED), 300);
      ASSERT_EQ(usurfHandle(3, 6, GHOSTED), 200);
      ASSERT_EQ(usurfHandle(3, 7, GHOSTED), 100);
      ASSERT_EQ(usurfHandle(4, 5, GHOSTED), 300);
      ASSERT_EQ(usurfHandle(4, 6, GHOSTED), 200);
      ASSERT_EQ(usurfHandle(4, 7, GHOSTED), 100);
      ASSERT_EQ(usurfHandle(5, 5, GHOSTED), 300);
      ASSERT_EQ(usurfHandle(5, 6, GHOSTED), 200);
      ASSERT_EQ(usurfHandle(5, 7, GHOSTED), 100);
    }
  }

  {
    auto &topgHandle = model->topg->getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(topgHandle(0, 0, GHOSTED), 0);
      ASSERT_EQ(topgHandle(0, 1, GHOSTED), 0);
      ASSERT_EQ(topgHandle(0, 2, GHOSTED), 0);
      ASSERT_EQ(topgHandle(1, 0, GHOSTED), 0);
      ASSERT_EQ(topgHandle(1, 1, GHOSTED), 0);
      ASSERT_EQ(topgHandle(1, 2, GHOSTED), 0);
      ASSERT_EQ(topgHandle(2, 0, GHOSTED), 0);
      ASSERT_EQ(topgHandle(2, 1, GHOSTED), 0);
      ASSERT_EQ(topgHandle(2, 2, GHOSTED), 0);
    } else if (mpiRank == 1) {
      ASSERT_EQ(topgHandle(0, 0, GHOSTED), 0);
      ASSERT_EQ(topgHandle(0, 1, GHOSTED), 0);
      ASSERT_EQ(topgHandle(0, 2, GHOSTED), 0);
      ASSERT_EQ(topgHandle(1, 0, GHOSTED), 0);
      ASSERT_EQ(topgHandle(1, 1, GHOSTED), 0);
      ASSERT_EQ(topgHandle(1, 2, GHOSTED), 0);
      ASSERT_EQ(topgHandle(2, 0, GHOSTED), 0);
      ASSERT_EQ(topgHandle(2, 1, GHOSTED), 0);
      ASSERT_EQ(topgHandle(2, 2, GHOSTED), 0);
    } else if (mpiRank == 5) {
      ASSERT_EQ(topgHandle(3, 5, GHOSTED), 0);
      ASSERT_EQ(topgHandle(3, 6, GHOSTED), 0);
      ASSERT_EQ(topgHandle(3, 7, GHOSTED), 0);
      ASSERT_EQ(topgHandle(4, 5, GHOSTED), 0);
      ASSERT_EQ(topgHandle(4, 6, GHOSTED), 0);
      ASSERT_EQ(topgHandle(4, 7, GHOSTED), 0);
      ASSERT_EQ(topgHandle(5, 5, GHOSTED), 0);
      ASSERT_EQ(topgHandle(5, 6, GHOSTED), 0);
      ASSERT_EQ(topgHandle(5, 7, GHOSTED), 0);
    }
  }

  {
    auto &thkHandle = model->thk->getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(thkHandle(0, 0, GHOSTED), 2000);
      ASSERT_EQ(thkHandle(0, 1, GHOSTED), 1900);
      ASSERT_EQ(thkHandle(0, 2, GHOSTED), 1800);
      ASSERT_EQ(thkHandle(1, 0, GHOSTED), 2000);
      ASSERT_EQ(thkHandle(1, 1, GHOSTED), 1900);
      ASSERT_EQ(thkHandle(1, 2, GHOSTED), 1800);
      ASSERT_EQ(thkHandle(2, 0, GHOSTED), 2000);
      ASSERT_EQ(thkHandle(2, 1, GHOSTED), 1900);
      ASSERT_EQ(thkHandle(2, 2, GHOSTED), 1800);
    } else if (mpiRank == 1) {
      ASSERT_EQ(thkHandle(0, 0, GHOSTED), 1400);
      ASSERT_EQ(thkHandle(0, 1, GHOSTED), 1300);
      ASSERT_EQ(thkHandle(0, 2, GHOSTED), 1200);
      ASSERT_EQ(thkHandle(1, 0, GHOSTED), 1400);
      ASSERT_EQ(thkHandle(1, 1, GHOSTED), 1300);
      ASSERT_EQ(thkHandle(1, 2, GHOSTED), 1200);
      ASSERT_EQ(thkHandle(2, 0, GHOSTED), 1400);
      ASSERT_EQ(thkHandle(2, 1, GHOSTED), 1300);
      ASSERT_EQ(thkHandle(2, 2, GHOSTED), 1200);
    } else if (mpiRank == 5) {
      ASSERT_EQ(thkHandle(3, 5, GHOSTED), 300);
      ASSERT_EQ(thkHandle(3, 6, GHOSTED), 200);
      ASSERT_EQ(thkHandle(3, 7, GHOSTED), 100);
      ASSERT_EQ(thkHandle(4, 5, GHOSTED), 300);
      ASSERT_EQ(thkHandle(4, 6, GHOSTED), 200);
      ASSERT_EQ(thkHandle(4, 7, GHOSTED), 100);
      ASSERT_EQ(thkHandle(5, 5, GHOSTED), 300);
      ASSERT_EQ(thkHandle(5, 6, GHOSTED), 200);
      ASSERT_EQ(thkHandle(5, 7, GHOSTED), 100);
    }
  }

  {
    auto &bndMaskHandle = model->bndMask->getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(bndMaskHandle(0, 0, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(0, 1, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(0, 2, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(1, 0, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(1, 1, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(1, 2, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(2, 0, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(2, 1, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(2, 2, GHOSTED), 0);
    } else if (mpiRank == 2) {
      ASSERT_EQ(bndMaskHandle(0, 5, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(0, 6, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(0, 7, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(1, 5, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(1, 6, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(1, 7, GHOSTED), DIRICHLET_FLAG);
      ASSERT_EQ(bndMaskHandle(2, 5, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(2, 6, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(2, 7, GHOSTED), DIRICHLET_FLAG);
    } else if (mpiRank == 5) {
      ASSERT_EQ(bndMaskHandle(3, 5, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(3, 6, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(3, 7, GHOSTED), DIRICHLET_FLAG);
      ASSERT_EQ(bndMaskHandle(4, 5, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(4, 6, GHOSTED), 0);
      ASSERT_EQ(bndMaskHandle(4, 7, GHOSTED), DIRICHLET_FLAG);
      ASSERT_EQ(bndMaskHandle(5, 5, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(5, 6, GHOSTED), NOFLOW_FLAG);
      ASSERT_EQ(bndMaskHandle(5, 7, GHOSTED), NOFLOW_FLAG);
    }
  }

  // TODO difference between Q and bmelt
  {
    auto &QHandle = model->Q->getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(QHandle(0, 0, GHOSTED), 1);
      ASSERT_EQ(QHandle(0, 1, GHOSTED), 1);
      ASSERT_EQ(QHandle(0, 2, GHOSTED), 1);
      ASSERT_EQ(QHandle(1, 0, GHOSTED), 1);
      ASSERT_EQ(QHandle(1, 1, GHOSTED), 1);
      ASSERT_EQ(QHandle(1, 2, GHOSTED), 1);
      ASSERT_EQ(QHandle(2, 0, GHOSTED), 1);
      ASSERT_EQ(QHandle(2, 1, GHOSTED), 1);
      ASSERT_EQ(QHandle(2, 2, GHOSTED), 1);
    } else if (mpiRank == 1) {
      ASSERT_EQ(QHandle(0, 0, GHOSTED), 1);
      ASSERT_EQ(QHandle(0, 1, GHOSTED), 1);
      ASSERT_EQ(QHandle(0, 2, GHOSTED), 1);
      ASSERT_EQ(QHandle(1, 0, GHOSTED), 1);
      ASSERT_EQ(QHandle(1, 1, GHOSTED), 1);
      ASSERT_EQ(QHandle(1, 2, GHOSTED), 1);
      ASSERT_EQ(QHandle(2, 0, GHOSTED), 1);
      ASSERT_EQ(QHandle(2, 1, GHOSTED), 1);
      ASSERT_EQ(QHandle(2, 2, GHOSTED), 1);
    } else if (mpiRank == 5) {
      ASSERT_EQ(QHandle(3, 5, GHOSTED), 1);
      ASSERT_EQ(QHandle(3, 6, GHOSTED), 1);
      ASSERT_EQ(QHandle(3, 7, GHOSTED), 1);
      ASSERT_EQ(QHandle(4, 5, GHOSTED), 1);
      ASSERT_EQ(QHandle(4, 6, GHOSTED), 1);
      ASSERT_EQ(QHandle(4, 7, GHOSTED), 1);
      ASSERT_EQ(QHandle(5, 5, GHOSTED), 1);
      ASSERT_EQ(QHandle(5, 6, GHOSTED), 1);
      ASSERT_EQ(QHandle(5, 7, GHOSTED), 1);
    }
  }

  {
    auto &pIceHandle = model->pIce->getReadHandle();
    if (mpiRank == 0) {
      ASSERT_EQ(pIceHandle(0, 0, GHOSTED), 2000 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(0, 1, GHOSTED), 1900 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(0, 2, GHOSTED), 1800 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(1, 0, GHOSTED), 2000 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(1, 1, GHOSTED), 1900 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(1, 2, GHOSTED), 1800 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(2, 0, GHOSTED), 2000 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(2, 1, GHOSTED), 1900 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(2, 2, GHOSTED), 1800 * RHO_ICE * GRAVITY);
    } else if (mpiRank == 1) {
      ASSERT_EQ(pIceHandle(0, 0, GHOSTED), 1400 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(0, 1, GHOSTED), 1300 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(0, 2, GHOSTED), 1200 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(1, 0, GHOSTED), 1400 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(1, 1, GHOSTED), 1300 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(1, 2, GHOSTED), 1200 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(2, 0, GHOSTED), 1400 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(2, 1, GHOSTED), 1300 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(2, 2, GHOSTED), 1200 * RHO_ICE * GRAVITY);
    } else if (mpiRank == 5) {
      ASSERT_EQ(pIceHandle(3, 5, GHOSTED), 300 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(3, 6, GHOSTED), 200 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(3, 7, GHOSTED), 100 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(4, 5, GHOSTED), 300 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(4, 6, GHOSTED), 200 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(4, 7, GHOSTED), 100 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(5, 5, GHOSTED), 300 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(5, 6, GHOSTED), 200 * RHO_ICE * GRAVITY);
      ASSERT_EQ(pIceHandle(5, 7, GHOSTED), 100 * RHO_ICE * GRAVITY);
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
