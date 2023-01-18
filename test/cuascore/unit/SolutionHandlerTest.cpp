/**
 * File: SolutionHandlerTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "SolutionHandler.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define GRID_SIZE_X 9
#define GRID_SIZE_Y 5

TEST(SolutionHandlerTest, constructor) {
  CUAS::SolutionHandler handler("SolutionHandlerConstructor.nc", GRID_SIZE_X, GRID_SIZE_Y, "normal");
}

TEST(SolutionHandlerTest, storeInitialSetup) {
  PETScGrid u(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid melt(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid creep(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid cavity(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid temp(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid flux(GRID_SIZE_X, GRID_SIZE_Y);
  CUAS::CUASModel model(GRID_SIZE_X, GRID_SIZE_Y);
  CUAS::SolutionHandler handler("SolutionHandlerSaveInitalSetup.nc", GRID_SIZE_X, GRID_SIZE_Y, "normal");

  int argc = 2;
  char arg0[] = "test";
  char arg1[] = "--verbose";
  char *argv[] = {arg0, arg1};
  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

  handler.storeInitialSetup(u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp, args);
}

TEST(SolutionHandlerTest, storeSolution) {
  PETScGrid u(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid melt(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid creep(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid cavity(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid temp(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid flux(GRID_SIZE_X, GRID_SIZE_Y);
  CUAS::CUASModel model(GRID_SIZE_X, GRID_SIZE_Y);

  CUAS::SolutionHandler handler("SolutionHandlerSaveSolution.nc", GRID_SIZE_X, GRID_SIZE_Y, "normal");
  handler.storeSolution(10, u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp);
  handler.storeSolution(15, u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp);
  handler.storeSolution(20, u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp);
}

TEST(SolutionHandlerTest, storeSolutionOutputSize) {
  PETScGrid u(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid T(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid melt(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid creep(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid cavity(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid temp(GRID_SIZE_X, GRID_SIZE_Y);
  PETScGrid flux(GRID_SIZE_X, GRID_SIZE_Y);
  CUAS::CUASModel model(GRID_SIZE_X, GRID_SIZE_Y);

  {
    CUAS::SolutionHandler handler("SolutionHandlerSaveSolutionSmall.nc", GRID_SIZE_X, GRID_SIZE_Y, "small");
    handler.storeSolution(10, u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp);
  }
  {
    CUAS::SolutionHandler handler("SolutionHandlerSaveSolutionNormal.nc", GRID_SIZE_X, GRID_SIZE_Y, "normal");
    handler.storeSolution(10, u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp);
  }
  {
    CUAS::SolutionHandler handler("SolutionHandlerSaveSolutionLarge.nc", GRID_SIZE_X, GRID_SIZE_Y, "large");
    handler.storeSolution(10, u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp);
  }
  {
    CUAS::SolutionHandler handler("SolutionHandlerSaveSolutionXLarge.nc", GRID_SIZE_X, GRID_SIZE_Y, "xlarge");
    handler.storeSolution(10, u, T, model, flux, melt, creep, cavity, temp, temp, temp, temp);
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
