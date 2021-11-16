#include "NetCDFFile.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define GRID_SIZE_X 9
#define GRID_SIZE_Y 5

TEST(NetCDFFileTest, define) {
  CUAS::NetCDFFile file("define.nc", GRID_SIZE_X, GRID_SIZE_Y);

  file.defineScalar("limitedScalar", LIMITED);
  file.defineScalar("unlimitedScalar", UNLIMITED);
  file.defineVectorX("vectorX");
  file.defineVectorY("vectorY");
  file.defineGrid("limitedGrid", LIMITED);
  file.defineGrid("unlimitedGrid", UNLIMITED);
}

TEST(NetCDFFileTest, addAttributes) {
  CUAS::NetCDFFile file("attributes.nc", GRID_SIZE_X, GRID_SIZE_Y);

  file.defineScalar("limitedScalar", LIMITED);
  file.addAttributeToVariable("limitedScalar", "testAttribute1", "testvalue1");
  file.addAttributeToVariable("limitedScalar", "testAttribute2", "testvalue2");
  file.defineScalar("unlimitedScalar", UNLIMITED);
  file.addAttributeToVariable("unlimitedScalar", "testAttribute1", "testvalue1");
  file.addAttributeToVariable("unlimitedScalar", "testAttribute2", "testvalue2");
  file.defineVectorX("vectorX");
  file.addAttributeToVariable("vectorX", "testAttribute1", "testvalue1");
  file.addAttributeToVariable("vectorX", "testAttribute2", "testvalue2");
  file.defineVectorY("vectorY");
  file.addAttributeToVariable("vectorY", "testAttribute1", "testvalue1");
  file.addAttributeToVariable("vectorY", "testAttribute2", "testvalue2");
  file.defineGrid("limitedGrid", LIMITED);
  file.addAttributeToVariable("limitedGrid", "testAttribute1", "testvalue1");
  file.addAttributeToVariable("limitedGrid", "testAttribute2", "testvalue2");
  file.defineGrid("unlimitedGrid", UNLIMITED);
  file.addAttributeToVariable("unlimitedGrid", "testAttribute1", "testvalue1");
  file.addAttributeToVariable("unlimitedGrid", "testAttribute2", "testvalue2");
}

TEST(NetCDFFileTest, addGlobalAttributes) {
  CUAS::NetCDFFile file("globalAttributes.nc", GRID_SIZE_X, GRID_SIZE_Y);

  file.addGlobalAttribute("globalText1", "testvalue1");
  file.addGlobalAttribute("globalText2", "testvalue2");
  file.addGlobalAttribute("globalBool", true);
  file.addGlobalAttribute("globalScalar", 1.23456);
  file.addGlobalAttribute("globalInt", 42);
}

TEST(NetCDFFileTest, write) {
  CUAS::NetCDFFile file("output.nc", GRID_SIZE_X, GRID_SIZE_Y);

  file.defineScalar("limitedScalar", LIMITED);
  file.defineScalar("unlimitedScalar", UNLIMITED);
  file.defineVectorX("vectorX");
  file.defineVectorY("vectorY");
  file.defineVectorY("vectorSTD");
  file.defineGrid("limitedGrid", LIMITED);
  file.defineGrid("unlimitedGrid", UNLIMITED);

  {
    // Test Scalars
    file.write("limitedScalar", 5);
    // overwrite value of limitedScalar
    file.write("limitedScalar", 3);
    file.write("unlimitedScalar", 5, 0);
    file.write("unlimitedScalar", 3, 1);
    file.write("unlimitedScalar", 1, 2);
  }

  {
    // Test Vectors
    PETScVector vecX(GRID_SIZE_X);
    int rangeBegin, rangeEnd;
    vecX.getOwnershipRange(rangeBegin, rangeEnd);
    for (int i = rangeBegin; i < rangeEnd; ++i) {
      vecX.setValue(i, i * 1.1);
    }
    vecX.assemble();
    file.write("vectorX", vecX);
  }

  {
    PETScVector vecY(GRID_SIZE_Y);
    int rangeBegin, rangeEnd;
    vecY.getOwnershipRange(rangeBegin, rangeEnd);
    for (int i = rangeBegin; i < rangeEnd; ++i) {
      vecY.setValue(i, i * 1.1);
    }
    vecY.assemble();
    file.write("vectorY", vecY);
  }

  {
    std::vector<PetscScalar> vecSTD(GRID_SIZE_Y);
    for (auto &entry : vecSTD) {
      entry = 10;
    }
    file.write("vectorSTD", vecSTD);
  }

  {
    // Test Grids
    PETScGrid grid(GRID_SIZE_X, GRID_SIZE_Y);
    {
      auto cornerX = grid.getCornerX();
      auto cornerY = grid.getCornerY();
      auto handle = grid.getWriteHandle();

      for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
          auto index = (cornerY + j) * grid.getTotalNumOfCols() + cornerX + i;
          handle(j, i) = index * 1.1;
        }
      }
    }

    file.write("limitedGrid", grid, 0);

    file.write("unlimitedGrid", grid, 0);

    {
      auto cornerX = grid.getCornerX();
      auto cornerY = grid.getCornerY();
      auto handle = grid.getWriteHandle();

      for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
          auto index = (cornerY + j) * grid.getTotalNumOfCols() + cornerX + i;
          handle(j, i) = index * 1.11;
        }
      }
    }
    file.write("unlimitedGrid", grid, 1);
    {
      auto cornerX = grid.getCornerX();
      auto cornerY = grid.getCornerY();
      auto handle = grid.getWriteHandle();

      for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
          auto index = (cornerY + j) * grid.getTotalNumOfCols() + cornerX + i;
          handle(j, i) = index * 1.111;
        }
      }
    }
    file.write("unlimitedGrid", grid, 2);
  }
}

TEST(NetCDFFileTest, readLimitedFile) {
  // TODO
  // EXPECT_TRUE(false);
}

TEST(NetCDFFileTest, readUnlimitedFile) {
  CUAS::NetCDFFile file("output.nc", 'r');

  {
    PetscScalar v;
    file.read("limitedScalar", v);
    EXPECT_DOUBLE_EQ(v, 3);
    file.read("unlimitedScalar", v);
    EXPECT_DOUBLE_EQ(v, 1);
  }

  {
    PETScVector vecX(GRID_SIZE_X);
    file.read("vectorX", vecX);
    int rangeBegin, rangeEnd;
    vecX.getOwnershipRange(rangeBegin, rangeEnd);
    int localLength = rangeEnd - rangeBegin;
    std::vector<PetscInt> index;
    index.resize(localLength);
    for (int i = 0; i < localLength; ++i) {
      index[i] = rangeBegin + i;
    }
    std::vector<PetscScalar> values;
    values.resize(localLength);
    VecGetValues(vecX.getRaw(), localLength, index.data(), values.data());
    for (int i = 0; i < localLength; ++i) {
      EXPECT_DOUBLE_EQ(values[i], 1.1 * (rangeBegin + i));
    }
  }

  {
    PETScVector vecY(GRID_SIZE_Y);
    file.read("vectorY", vecY);
    int rangeBegin, rangeEnd;
    vecY.getOwnershipRange(rangeBegin, rangeEnd);
    int localLength = rangeEnd - rangeBegin;
    std::vector<PetscInt> index;
    index.resize(localLength);
    for (int i = 0; i < localLength; ++i) {
      index[i] = rangeBegin + i;
    }
    std::vector<PetscScalar> values;
    values.resize(localLength);
    VecGetValues(vecY.getRaw(), localLength, index.data(), values.data());
    for (int i = 0; i < localLength; ++i) {
      EXPECT_DOUBLE_EQ(values[i], 1.1 * (rangeBegin + i));
    }
  }

  {
    PETScGrid grid(GRID_SIZE_X, GRID_SIZE_Y);
    file.read("limitedGrid", grid);
    {
      auto cornerX = grid.getCornerX();
      auto cornerY = grid.getCornerY();
      auto handle = grid.getReadHandle();
      for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
          auto index = (cornerY + j) * grid.getTotalNumOfCols() + cornerX + i;
          EXPECT_DOUBLE_EQ(handle(j, i), index * 1.1);
        }
      }
    }

    file.read("unlimitedGrid", grid);
    {
      auto cornerX = grid.getCornerX();
      auto cornerY = grid.getCornerY();
      auto handle = grid.getReadHandle();
      for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
          auto index = (cornerY + j) * grid.getTotalNumOfCols() + cornerX + i;
          EXPECT_DOUBLE_EQ(handle(j, i), index * 1.111);
        }
      }
    }
  }
}

/*TEST(NetCDFFileTest, readTimeStepArray) {
  // TODO
}*/

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
