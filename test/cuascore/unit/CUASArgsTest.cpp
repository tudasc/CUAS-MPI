/**
 * File: CUASArgsTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASArgs.h"

#include "gtest/gtest.h"

#include "petsc.h"

int mpiRank;
int mpiSize;

TEST(CUASArgs, defaults) {
  char arg0[] = "test";
  char *argv[] = {arg0};
  int argc = sizeof(argv) / sizeof(argv[0]);

  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

  ASSERT_EQ(args.verbose, false);
  ASSERT_EQ(args.verboseSolver, false);
  ASSERT_EQ(args.directSolver, false);
  ASSERT_EQ(args.Tmax, 20.0);
  ASSERT_EQ(args.Tmin, 0.0000001);
  ASSERT_EQ(args.totaltime, "");
  ASSERT_EQ(args.dt, "");
  ASSERT_EQ(args.saveEvery, 0);
  ASSERT_EQ(args.conductivity, 10.0);
  ASSERT_EQ(args.doAllChannels, false);
  ASSERT_EQ(args.doAnyChannel, false);
  ASSERT_EQ(args.doCreep, false);
  ASSERT_EQ(args.doMelt, false);
  ASSERT_EQ(args.doCavity, false);
  ASSERT_EQ(args.disableUnconfined, false);
  ASSERT_EQ(args.flowConstant, 5e-25);
  ASSERT_EQ(args.roughnessFactor, 1.0);
  ASSERT_EQ(args.supplyMultiplier, 1.0);
  ASSERT_EQ(args.layerThickness, 0.1);
  ASSERT_EQ(args.unconfSmooth, 0.0);
  ASSERT_EQ(args.restart, "");
  ASSERT_EQ(args.restartNoneZeroInitialGuess, true);
  ASSERT_EQ(args.seaLevelForcing, "");
  ASSERT_EQ(args.output, "out.nc");
  ASSERT_EQ(args.input, "");
  ASSERT_EQ(args.timeStepFile, "");
  ASSERT_EQ(args.forcingFile, "");
  ASSERT_EQ(args.timeSteppingTheta, 1.0);
  ASSERT_EQ(args.sizeOfForcingBuffer, -1);
  ASSERT_EQ(args.starttime, "");
  ASSERT_EQ(args.endtime, "");
  ASSERT_EQ(args.enableUDS, false);
  ASSERT_EQ(args.thresholdThicknessUDS, 0.0);
  ASSERT_EQ(args.disableNonNegative, false);
  ASSERT_EQ(args.nonLinearIters, 0);
}

TEST(CUASArgs, allOpts) {
  char arg0[] = "test";
  char arg1[] = "--Tmax=7";
  char arg2[] = "--Tmin=3";
  char arg3[] = "--totaltime='4 weeks'";
  char arg4[] = "--dt='1 day'";
  char arg5[] = "--saveEvery=12";
  char arg6[] = "--conductivity=13.37";
  char arg7[] = "--doChannels";
  // disable unconfined = false
  char arg8[] = "--flowConstant=2.1";
  char arg9[] = "--roughnessFactor=2";
  char arg10[] = "--supplyMultiplier=1.2";
  char arg11[] = "--layerThickness=1.334";
  char arg12[] = "--unconfSmooth=4.3";
  char arg13[] = "--restart='restartFile.nc'";
  char arg14[] = "--seaLevelForcing='seaForcingFile.nc'";
  char arg15[] = "--output='output.nc'";
  char arg16[] = "--input='input.nc'";
  char arg17[] = "--selectedChannels=creep";
  char arg18[] = "--verbose";
  char arg19[] = "--timeStepFile='timesteps.nc'";
  char arg20[] = "--verboseSolver";
  char arg21[] = "--directSolver";
  char arg22[] = "--forcingFile='timeforcing.nc'";
  char arg23[] = "--timeSteppingTheta=0.5";
  char arg24[] = "--sizeOfForcingBuffer=5";
  char arg25[] = "--starttime='1 year'";
  char arg26[] = "--endtime='1 year 4 weeks'";
  char arg27[] = "--enableUDS";
  char arg28[] = "--thresholdThicknessUDS=0.1";
  char arg29[] = "--disableNonNegative";
  char arg30[] = "--nonLinearIters=1";

  // char arg25[] = "--restartNoneZeroInitialGuess=false";

  char *argv[] = {arg0,  arg1,  arg2,  arg3,  arg4,  arg5,  arg6,  arg7,  arg8,  arg9,  arg10,
                  arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20, arg21,
                  arg22, arg23, arg24, arg25, arg26, arg27, arg28, arg29, arg30};
  int argc = sizeof(argv) / sizeof(argv[0]);

  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

  ASSERT_EQ(args.verbose, true);
  ASSERT_EQ(args.verboseSolver, true);
  ASSERT_EQ(args.directSolver, true);
  ASSERT_EQ(args.Tmax, 7);
  ASSERT_EQ(args.Tmin, 3);
  ASSERT_EQ(args.totaltime, "'4 weeks'");
  ASSERT_EQ(args.dt, "'1 day'");
  ASSERT_EQ(args.saveEvery, 12);
  ASSERT_EQ(args.conductivity, 13.37);
  ASSERT_EQ(args.doAllChannels, false);
  ASSERT_EQ(args.doAnyChannel, true);
  ASSERT_EQ(args.doCreep, true);
  ASSERT_EQ(args.doMelt, false);
  ASSERT_EQ(args.doCavity, false);
  ASSERT_EQ(args.disableUnconfined, false);
  ASSERT_EQ(args.flowConstant, 2.1);
  ASSERT_EQ(args.roughnessFactor, 2);
  ASSERT_EQ(args.supplyMultiplier, 1.2);
  ASSERT_EQ(args.layerThickness, 1.334);
  ASSERT_EQ(args.unconfSmooth, 4.3);
  ASSERT_EQ(args.restart, "'restartFile.nc'");
  ASSERT_EQ(args.restartNoneZeroInitialGuess, true);
  ASSERT_EQ(args.seaLevelForcing, "'seaForcingFile.nc'");
  ASSERT_EQ(args.output, "'output.nc'");
  ASSERT_EQ(args.input, "'input.nc'");
  ASSERT_EQ(args.timeStepFile, "'timesteps.nc'");
  ASSERT_EQ(args.forcingFile, "'timeforcing.nc'");
  ASSERT_EQ(args.timeSteppingTheta, 0.5);
  ASSERT_EQ(args.sizeOfForcingBuffer, 5);
  ASSERT_EQ(args.starttime, "'1 year'");
  ASSERT_EQ(args.endtime, "'1 year 4 weeks'");
  ASSERT_EQ(args.enableUDS, true);
  ASSERT_EQ(args.thresholdThicknessUDS, 0.1);
  ASSERT_EQ(args.disableNonNegative, true);
  ASSERT_EQ(args.nonLinearIters, 1);
}

TEST(CUASArgs, selectiveChannels) {
  constexpr int argc = 3;
  char *argv[argc];
  char arg0[] = "test";
  argv[0] = arg0;

  // doChannels==true
  {
    char doChannels[] = "--doChannels";
    argv[1] = doChannels;
    {
      char selectedChannels[] = "";
      argv[2] = selectedChannels;

      CUAS::CUASArgs args;
      CUAS::parseArgs(argc, argv, args);

      ASSERT_EQ(args.doAllChannels, true);
      ASSERT_EQ(args.doAnyChannel, true);
      ASSERT_EQ(args.doCreep, true);
      ASSERT_EQ(args.doMelt, true);
      ASSERT_EQ(args.doCavity, true);
    }

    {
      char selectedChannels[] = "--selectedChannels=melt";
      argv[2] = selectedChannels;

      CUAS::CUASArgs args;
      CUAS::parseArgs(argc, argv, args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, true);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, true);
      ASSERT_EQ(args.doCavity, false);
    }

    {
      char selectedChannels[] = "--selectedChannels=";
      argv[2] = selectedChannels;

      CUAS::CUASArgs args;
      CUAS::parseArgs(argc, argv, args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, false);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, false);
      ASSERT_EQ(args.doCavity, false);
    }
  }

  // doChannels==false
  {
    char doChannels[] = "";
    argv[2] = doChannels;
    {
      char selectedChannels[] = "";
      argv[1] = selectedChannels;

      CUAS::CUASArgs args;
      CUAS::parseArgs(argc, argv, args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, false);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, false);
      ASSERT_EQ(args.doCavity, false);
    }

    {
      char selectedChannels[] = "--selectedChannels=melt,creep";
      argv[1] = selectedChannels;

      CUAS::CUASArgs args;
      CUAS::parseArgs(argc, argv, args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, true);
      ASSERT_EQ(args.doCreep, true);
      ASSERT_EQ(args.doMelt, true);
      ASSERT_EQ(args.doCavity, false);
    }

    {
      char selectedChannels[] = "--selectedChannels=";
      argv[1] = selectedChannels;

      CUAS::CUASArgs args;
      CUAS::parseArgs(argc, argv, args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, false);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, false);
      ASSERT_EQ(args.doCavity, false);
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
