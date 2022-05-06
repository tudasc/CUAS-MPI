#include "CUASArgs.h"
#include "ModelReader.h"
#include "NetCDFFile.h"

#include "gtest/gtest.h"
#include <iterator>
#include <vector>

int mpiRank;
int mpiSize;

#define GRID_SIZE_X 9
#define GRID_SIZE_Y 5
#define TIME 20
#define CORRUPTED_TIME 30
#define CORRUPTED_DIM 2

TEST(TimeForcingSucessTest, readSuccessThreeDims) {
  std::string timeForcingFile = "timeForcingReadSuccessThreeDims.nc";
  // create test file that is correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds;
    dimIds.reserve(3);
    nc_def_dim(fileId, "time", TIME, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[2]);
    int varIdTime;
    if (int retval = nc_def_var(fileId, "time", NC_DOUBLE, 1, &dimIds[0], &varIdTime)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 3, dimIds.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    std::vector<PetscScalar> time(TIME);
    for (int i = 0; i < TIME; ++i) {
      time[i] = i;
    }

    // create 3d Matrix
    PetscScalar bmeltTime[TIME][GRID_SIZE_Y][GRID_SIZE_X];
    for (int t = 0; t < TIME; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        for (int x = 0; x < GRID_SIZE_X; ++x) {
          bmeltTime[t][y][x] = t;
        }
      }
    }
    if (int retval = nc_put_var_double(fileId, varIdTime, time.data())) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  auto isTimeForcing = CUAS::ModelReader::isTimeDependentField(timeForcingFile, "bmelt");
  EXPECT_TRUE(isTimeForcing);
}

TEST(TimeForcingSucessTest, readSuccessTwoDims) {
  std::string timeForcingFile = "timeForcingReadSuccessTwoDims.nc";
  // create test file that is correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds;
    dimIds.reserve(2);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[0]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[1]);

    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 2, dimIds.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    // create 2d Matrix
    PetscScalar bmeltTime[GRID_SIZE_Y][GRID_SIZE_X];
    for (int y = 0; y < GRID_SIZE_Y; ++y) {
      for (int x = 0; x < GRID_SIZE_X; ++x) {
        bmeltTime[y][x] = y;
      }
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  auto isTimeForcing = CUAS::ModelReader::isTimeDependentField(timeForcingFile, "bmelt");
  EXPECT_FALSE(isTimeForcing);
}

TEST(TimeForcingDeathTest, fieldAndTimeDifferentDimensions) {
  std::string timeForcingFile = "timeForcingFieldAndTimeDifferentDimensions.nc";
  // create test file that is not correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds(3);
    nc_def_dim(fileId, "time", TIME, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[2]);

    // this id vector causes a problem because the time dimension of the variable bmelt is now different than the time
    // dimension of the variable time
    std::vector<int> dimIdsCorrupted(3);
    nc_def_dim(fileId, "time_corrupted", CORRUPTED_TIME, &dimIdsCorrupted[0]);
    // copy the y and x dimensions from dimIds into the corrupted version
    std::copy(dimIds.begin() + 1, dimIds.end(), dimIdsCorrupted.begin() + 1);

    int varIdTime;
    if (int retval = nc_def_var(fileId, "time", NC_DOUBLE, 1, &dimIds[0], &varIdTime)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineVector(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 3, dimIdsCorrupted.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    std::vector<PetscScalar> time(TIME);
    for (int i = 0; i < TIME; ++i) {
      time[i] = i;
    }

    // create 3d Matrix
    PetscScalar bmeltTime[CORRUPTED_TIME][GRID_SIZE_Y][GRID_SIZE_X];
    for (int t = 0; t < CORRUPTED_TIME; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        for (int x = 0; x < GRID_SIZE_X; ++x) {
          bmeltTime[t][y][x] = t;
        }
      }
    }
    if (int retval = nc_put_var_double(fileId, varIdTime, time.data())) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  CUAS::ModelReader reader(timeForcingFile);
  // check that isTimeDependentField fails. The error message is arbitrary
  EXPECT_EXIT(reader.isTimeDependentField(timeForcingFile, "bmelt"), testing::ExitedWithCode(1), ".*");
}

TEST(TimeForcingDeathTest, fileHasNoFieldName) {
  std::string timeForcingFile = "timeForcingFileHasNoFieldName.nc";
  // create test file that is not correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds(3);
    nc_def_dim(fileId, "time", TIME, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[2]);

    int varIdTime;
    if (int retval = nc_def_var(fileId, "time", NC_DOUBLE, 1, &dimIds[0], &varIdTime)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineVector(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    int varIdTopg;
    if (int retval = nc_def_var(fileId, "topg", NC_DOUBLE, 3, dimIds.data(), &varIdTopg)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    std::vector<PetscScalar> time(TIME);
    for (int i = 0; i < TIME; ++i) {
      time[i] = i;
    }

    // create 3d Matrix
    PetscScalar topgTime[TIME][GRID_SIZE_Y][GRID_SIZE_X];
    for (int t = 0; t < TIME; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        for (int x = 0; x < GRID_SIZE_X; ++x) {
          topgTime[t][y][x] = t;
        }
      }
    }
    if (int retval = nc_put_var_double(fileId, varIdTime, time.data())) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    if (int retval = nc_put_var_double(fileId, varIdTopg, &topgTime[0][0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  CUAS::ModelReader reader(timeForcingFile);
  // check that isTimeDependentField fails. The error message is arbitrary
  EXPECT_EXIT(reader.isTimeDependentField(timeForcingFile, "bmelt"), testing::ExitedWithCode(1), ".*");
}

TEST(TimeForcingDeathTest, tooManyDimensions) {
  std::string timeForcingFile = "timeForcingTooManyDimensions.nc";
  // create test file that is not correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds(4);
    nc_def_dim(fileId, "time", TIME, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[2]);
    nc_def_dim(fileId, "corrupted_dim", CORRUPTED_DIM, &dimIds[3]);

    int varIdTime;
    if (int retval = nc_def_var(fileId, "time", NC_DOUBLE, 1, &dimIds[0], &varIdTime)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineVector(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 4, dimIds.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    std::vector<PetscScalar> time(TIME);
    for (int i = 0; i < TIME; ++i) {
      time[i] = i;
    }

    // create 4d Matrix
    PetscScalar bmeltTime[CORRUPTED_TIME][GRID_SIZE_Y][GRID_SIZE_X][CORRUPTED_DIM];
    for (int t = 0; t < CORRUPTED_TIME; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        for (int x = 0; x < GRID_SIZE_X; ++x) {
          for (int dim = 0; dim < CORRUPTED_DIM; ++dim) {
            bmeltTime[t][y][x][dim] = t;
          }
        }
      }
    }
    if (int retval = nc_put_var_double(fileId, varIdTime, time.data())) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PetscVec: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0][0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  CUAS::ModelReader reader(timeForcingFile);
  // check that isTimeDependentField fails. The error message is arbitrary
  EXPECT_EXIT(reader.isTimeDependentField(timeForcingFile, "bmelt"), testing::ExitedWithCode(1), ".*");
}

TEST(TimeForcingDeathTest, twoDimsButOneIsTime) {
  std::string timeForcingFile = "timeForcingTwoDimsButOneIsTime.nc";
  // create test file that is not correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds(3);
    nc_def_dim(fileId, "time", TIME, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    // it is necessary to keep the dimension x in the file as it is required by our netcdf design
    nc_def_dim(fileId, "x", GRID_SIZE_Y, &dimIds[2]);

    int varIdTime;
    if (int retval = nc_def_var(fileId, "time", NC_DOUBLE, 1, &dimIds[0], &varIdTime)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineVector(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 2, dimIds.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    std::vector<PetscScalar> time(TIME);
    for (int i = 0; i < TIME; ++i) {
      time[i] = i;
    }

    // create 4d Matrix
    PetscScalar bmeltTime[TIME][GRID_SIZE_Y];
    for (int t = 0; t < TIME; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        bmeltTime[t][y] = t;
      }
    }
    if (int retval = nc_put_var_double(fileId, varIdTime, time.data())) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PetscVec: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  CUAS::ModelReader reader(timeForcingFile);
  // check that isTimeDependentField fails. The error message is arbitrary
  EXPECT_EXIT(reader.isTimeDependentField(timeForcingFile, "bmelt"), testing::ExitedWithCode(1), ".*");
}

TEST(TimeForcingDeathTest, threeDimsButNoDimTime) {
  std::string timeForcingFile = "timeForcingThreeDimsButNoDimTime.nc";
  // create test file that is not correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds(3);
    nc_def_dim(fileId, "z", CORRUPTED_DIM, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[2]);

    int varIdTime;
    if (int retval = nc_def_var(fileId, "time", NC_DOUBLE, 1, &dimIds[0], &varIdTime)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineVector(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 3, dimIds.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    std::vector<PetscScalar> time(TIME);
    for (int i = 0; i < TIME; ++i) {
      time[i] = i;
    }

    // create 3d Matrix
    PetscScalar bmeltTime[CORRUPTED_DIM][GRID_SIZE_Y][GRID_SIZE_X];
    for (int t = 0; t < CORRUPTED_DIM; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        for (int x = 0; x < GRID_SIZE_X; ++x) {
          bmeltTime[t][y][x] = t;
        }
      }
    }
    if (int retval = nc_put_var_double(fileId, varIdTime, time.data())) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  CUAS::ModelReader reader(timeForcingFile);
  // check that isTimeDependentField fails. The error message is arbitrary
  EXPECT_EXIT(reader.isTimeDependentField(timeForcingFile, "bmelt"), testing::ExitedWithCode(1), ".*");
}

TEST(TimeForcingDeathTest, threeDimsButNoFieldTime) {
  std::string timeForcingFile = "timeForcingThreeDimsButNoFieldTime.nc";
  // create test file that is not correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds(3);
    nc_def_dim(fileId, "time", TIME, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[2]);

    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 3, dimIds.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    // create 3d Matrix
    PetscScalar bmeltTime[CORRUPTED_DIM][GRID_SIZE_Y][GRID_SIZE_X];
    for (int time = 0; time < CORRUPTED_DIM; ++time) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        for (int x = 0; x < GRID_SIZE_X; ++x) {
          bmeltTime[time][y][x] = time;
        }
      }
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  CUAS::ModelReader reader(timeForcingFile);
  // check that isTimeDependentField fails. The error message is arbitrary
  EXPECT_EXIT(reader.isTimeDependentField(timeForcingFile, "bmelt"), testing::ExitedWithCode(1), ".*");
}

TEST(TimeForcingDeathTest, timeIsMultiDimensional) {
  std::string timeForcingFile = "timeForcingTimeIsMultiDimensional.nc";
  // create test file that is not correct
  {
    int fileId;
    if (int retval =
            nc_create_par(timeForcingFile.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    std::vector<int> dimIds(3);
    nc_def_dim(fileId, "time", TIME, &dimIds[0]);
    nc_def_dim(fileId, "y", GRID_SIZE_Y, &dimIds[1]);
    nc_def_dim(fileId, "x", GRID_SIZE_X, &dimIds[2]);

    int varIdTime;
    if (int retval = nc_def_var(fileId, "time", NC_DOUBLE, 2, &dimIds[0], &varIdTime)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineVector(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }
    int varIdBmelt;
    if (int retval = nc_def_var(fileId, "bmelt", NC_DOUBLE, 3, dimIds.data(), &varIdBmelt)) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_enddef(fileId);

    PetscScalar time[TIME][GRID_SIZE_Y];
    for (int t = 0; t < TIME; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        time[t][y] = t;
      }
    }

    // create 3d Matrix
    PetscScalar bmeltTime[TIME][GRID_SIZE_Y][GRID_SIZE_X];
    for (int t = 0; t < TIME; ++t) {
      for (int y = 0; y < GRID_SIZE_Y; ++y) {
        for (int x = 0; x < GRID_SIZE_X; ++x) {
          bmeltTime[t][y][x] = t;
        }
      }
    }
    if (int retval = nc_put_var_double(fileId, varIdTime, &time[0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    if (int retval = nc_put_var_double(fileId, varIdBmelt, &bmeltTime[0][0][0])) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError + "Exiting.")
      exit(1);
    }

    nc_close(fileId);
  }
  CUAS::CUASArgs args;
  CUAS::ModelReader reader(timeForcingFile);
  // check that isTimeDependentField fails. The error message is arbitrary
  EXPECT_EXIT(reader.isTimeDependentField(timeForcingFile, "bmelt"), testing::ExitedWithCode(1), ".*");
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
