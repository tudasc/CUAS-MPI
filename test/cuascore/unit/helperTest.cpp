#include "helper.h"
#include "physicalConstants.h"

#include "PetscGrid.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 9

TEST(HelperTest, pressure2head) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PetscGrid pressure(30, 12);
  PetscGrid bed_elevation(30, 12);
  PetscGrid head(30, 12);

  PetscScalar sea_level = 33.5;

  // setup
  auto pressure2d = pressure.getWriteHandle();
  auto bed_elevation2d = bed_elevation.getWriteHandle();
  auto head2d = head.getWriteHandle();
  for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
      pressure2d(j, i) = sea_level * mpiRank - 2.3;
      bed_elevation2d(j, i) = sea_level * sea_level - 31.3;
      head2d(j, i) = mpiRank * j + (i * 35);
    }
  }
  pressure2d.setValues();
  bed_elevation2d.setValues();
  head2d.setValues();

  // run pressure2head
  CUAS::pressure2head(head, pressure, bed_elevation, sea_level);

  auto &pressure2d2 = pressure.getReadHandle();
  auto &bed_elevation2d2 = bed_elevation.getReadHandle();
  auto &head2d2 = head.getReadHandle();
  for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d2(j, i) - sea_level;
      // check result
      ASSERT_EQ(head2d(j, i), pressure2d2(j, i) / (RHO_WATER * GRAVITY) + effective_bed_elevation);

      // check for sideeffects
      ASSERT_EQ(bed_elevation2d2(j, i), sea_level * sea_level - 31.3);
      ASSERT_EQ(pressure2d2(j, i), sea_level * mpiRank - 2.3);
    }
  }
  CUAS::head2pressure(pressure, head, bed_elevation, sea_level);
  CUAS::overburdenPressure(pressure, head);
}

TEST(HelperTest, head2pressure) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PetscGrid pressure(30, 12);
  PetscGrid bed_elevation(30, 12);
  PetscGrid head(30, 12);

  PetscScalar sea_level = 33.5;

  // setup
  auto pressure2d = pressure.getWriteHandle();
  auto bed_elevation2d = bed_elevation.getWriteHandle();
  auto head2d = head.getWriteHandle();
  for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
      pressure2d(j, i) = sea_level * mpiRank - 2.3;
      bed_elevation2d(j, i) = sea_level * sea_level - 31.3;
      head2d(j, i) = mpiRank * j + (i * 35);
    }
  }
  pressure2d.setValues();
  bed_elevation2d.setValues();
  head2d.setValues();

  // run head2pressure
  CUAS::head2pressure(pressure, head, bed_elevation, sea_level);

  auto &pressure2d2 = pressure.getReadHandle();
  auto &bed_elevation2d2 = bed_elevation.getReadHandle();
  auto &head2d2 = head.getReadHandle();
  for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d2(j, i) - sea_level;
      // check result
      ASSERT_EQ(pressure2d2(j, i), RHO_WATER * GRAVITY * (head2d(j, i) - effective_bed_elevation));

      // check for sideeffects
      ASSERT_EQ(bed_elevation2d2(j, i), sea_level * sea_level - 31.3);
      ASSERT_EQ(head2d2(j, i), mpiRank * j + (i * 35));
    }
  }
  CUAS::overburdenPressure(pressure, head);
}

TEST(HelperTest, overburdenPressure) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PetscGrid pressure(30, 12);
  PetscGrid head(30, 12);

  PetscScalar sea_level = 33.5;

  // setup
  auto pressure2d = pressure.getWriteHandle();
  auto head2d = head.getWriteHandle();
  for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
      pressure2d(j, i) = sea_level * mpiRank - 2.3;
      head2d(j, i) = mpiRank * j + (i * 35);
    }
  }
  pressure2d.setValues();
  pressure2d.setValues();

  // run overburdenPressure
  CUAS::overburdenPressure(pressure, head);

  auto &pressure2d2 = pressure.getReadHandle();
  auto &head2d2 = head.getReadHandle();
  for (int j = 0; j < pressure.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < pressure.getLocalNumOfCols(); ++i) {
      // check result
      ASSERT_EQ(pressure2d2(j, i), head2d2(j, i) * RHO_ICE * GRAVITY);

      // check for sideeffects
      ASSERT_EQ(head2d2(j, i), mpiRank * j + (i * 35));
    }
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
