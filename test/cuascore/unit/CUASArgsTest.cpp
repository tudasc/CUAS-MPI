#include "CUASArgs.h"

#include "gtest/gtest.h"

#include "petsc.h"

int mpiRank;
int mpiSize;

TEST(CUASArgs, allOpts) {
  char arg0[] = "test";
  char arg1[] = "--Tmax=7";
  char arg2[] = "--Tmin=3";
  char arg3[] = "--totaltime='4 Weeks'";
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
  char arg14[] = "--noSmoothMelt";
  char arg15[] = "--seaLevelForcing='seaForcingFile.nc'";
  char arg16[] = "--output='output.nc'";
  char arg17[] = "--input='input.nc'";
  char arg18[] = "--selectedChannels=creep";
  char arg19[] = "--verbose";
  char arg20[] = "--timeStepFile='timesteps.nc'";

  char *argv[] = {arg0,  arg1,  arg2,  arg3,  arg4,  arg5,  arg6,  arg7,  arg8,  arg9, arg10,
                  arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20};
  int argc = sizeof(argv) / sizeof(argv[0]);

  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

  ASSERT_EQ(args.verbose, true);
  ASSERT_EQ(args.Tmax, 7);
  ASSERT_EQ(args.Tmin, 3);
  ASSERT_EQ(args.totaltime, "'4 Weeks'");
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
  ASSERT_EQ(args.noSmoothMelt, true);
  ASSERT_EQ(args.seaLevelForcing, "'seaForcingFile.nc'");
  ASSERT_EQ(args.output, "'output.nc'");
  ASSERT_EQ(args.input, "'input.nc'");
  ASSERT_EQ(args.timeStepFile, "'timesteps.nc'");
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
