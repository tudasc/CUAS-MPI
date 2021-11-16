#include "ModelReader.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define GRID_SIZE_X 5
#define GRID_SIZE_Y 3

TEST(ModelReaderTest, constructor) { CUAS::ModelReader reader("testmodel.nc"); }

TEST(ModelReaderTest, fillModelFromNetcdf) {
  CUAS::ModelReader reader("testmodel.nc");
  auto model = reader.fillModelFromNetcdf();
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
