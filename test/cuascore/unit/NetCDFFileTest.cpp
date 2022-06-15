#include "NetCDFFile.h"
#include "timeparse.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define GRID_SIZE_X 9
#define GRID_SIZE_Y 5
#define TIME_LEN 20

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

  file.addGlobalAttribute("globalText1", std::string("testvalue1"));
  file.addGlobalAttribute("globalText2", std::string("testvalue2"));
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
      auto &handle = grid.getReadHandle();
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
      auto &handle = grid.getReadHandle();
      for (int j = 0; j < grid.getLocalNumOfRows(); ++j) {
        for (int i = 0; i < grid.getLocalNumOfCols(); ++i) {
          auto index = (cornerY + j) * grid.getTotalNumOfCols() + cornerX + i;
          EXPECT_DOUBLE_EQ(handle(j, i), index * 1.111);
        }
      }
    }
  }
}

TEST(NetCDFFileTest, readTimeArrayOfDifferentTypes) {
  // Note, NC_LONG and NC_INT are actually the same data type (signed 4 byte integer)
  // See netcdf.h

  std::string fileName = "readTimeArray.nc";

  // CUAS data type for time, this is also "long" but could change in the future
  std::vector<CUAS::timeSecs> time(TIME_LEN);
  // default types
  std::vector<double> time_double(TIME_LEN);
  std::vector<float> time_float(TIME_LEN);
  std::vector<long> time_long(TIME_LEN);
  std::vector<int> time_int(TIME_LEN);
  // unsupported, but should work
  std::vector<short> time_short(TIME_LEN);
  std::vector<unsigned int> time_uint(TIME_LEN);
  std::vector<unsigned short> time_ushort(TIME_LEN);

  for (short i = 0; i < TIME_LEN; ++i) {
    time[i] = (CUAS::timeSecs)i;
    time_double[i] = (double)i;
    time_float[i] = (float)i;
    time_long[i] = (long)i;
    time_int[i] = (int)i;
    time_short[i] = i;
    time_uint[i] = (unsigned int)i;
    time_ushort[i] = (unsigned short)i;
  }

  // generate a valid input file with different types for time
  {
    int fileId, dimId, varId, tmp;
    nc_create_par(fileName.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId);
    nc_def_dim(fileId, "time", TIME_LEN, &dimId);
    // NetCDFFile::NetCDFFile(name, 'r') needs x- and y-dims
    nc_def_dim(fileId, "x", GRID_SIZE_X, &tmp);  // dummy dim without data
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &tmp);  // dummy dim without data

    // default types
    nc_def_var(fileId, "time_double", NC_DOUBLE, 1, &dimId, &tmp);
    nc_def_var(fileId, "time_float", NC_FLOAT, 1, &dimId, &tmp);
    nc_def_var(fileId, "time_long", NC_LONG, 1, &dimId, &tmp);
    nc_def_var(fileId, "time_int", NC_INT, 1, &dimId, &tmp);
    // unsupported, but should work
    nc_def_var(fileId, "time_short", NC_SHORT, 1, &dimId, &tmp);
    nc_def_var(fileId, "time_uint", NC_UINT, 1, &dimId, &tmp);
    nc_def_var(fileId, "time_ushort", NC_USHORT, 1, &dimId, &tmp);

    nc_enddef(fileId);

    // now write the data
    nc_inq_varid(fileId, "time_double", &varId);
    nc_put_var_double(fileId, varId, time_double.data());
    nc_inq_varid(fileId, "time_float", &varId);
    nc_put_var_float(fileId, varId, time_float.data());
    nc_inq_varid(fileId, "time_long", &varId);
    nc_put_var_long(fileId, varId, time_long.data());
    nc_inq_varid(fileId, "time_int", &varId);
    nc_put_var_int(fileId, varId, time_int.data());

    nc_inq_varid(fileId, "time_short", &varId);
    nc_put_var_short(fileId, varId, time_short.data());
    nc_inq_varid(fileId, "time_uint", &varId);
    nc_put_var_uint(fileId, varId, time_uint.data());
    nc_inq_varid(fileId, "time_ushort", &varId);
    nc_put_var_ushort(fileId, varId, time_ushort.data());
    // done
    nc_close(fileId);
  }

  // open test file for reading
  CUAS::NetCDFFile file(fileName, 'r');
  auto nt = file.getDimLength("time");
  ASSERT_EQ(nt, TIME_LEN);

  std::vector<std::string> varNames = {"time_double", "time_float", "time_long",  "time_int",
                                       "time_short",  "time_uint",  "time_ushort"};

  for (auto &name : varNames) {
    std::vector<CUAS::timeSecs> time_in(nt);
    // std::cout << name << std::endl;
    file.read(name, time_in);
    for (int i = 0; i < nt; ++i) {
      EXPECT_EQ(time[i], time_in[i]);
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
