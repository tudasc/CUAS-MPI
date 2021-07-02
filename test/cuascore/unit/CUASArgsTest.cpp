#include "CUASArgs.h"

#include "gtest/gtest.h"

TEST(CUASArgs, allOpts) {
  int argc = 19;
  char arg0[] = "test";
  char arg05[] = "--verbose";
  char arg1[] = "--Tmax=7";
  char arg2[] = "--Tmin=3";
  char arg3[] = "--totaltime='4 Weeks'";
  char arg4[] = "--dt='1 day'";
  char arg5[] = "--saveEvery=12";
  char arg6[] = "--conductivity=13.37";
  char arg7[] = "--dochannels";
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
  char *argv[] = {arg0, arg05, arg1,  arg2,  arg3,  arg4,  arg5,  arg6,  arg7, arg8,
                  arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17};
  CUAS::CUASArgs args;
  CUAS::parseArgs(argc, argv, args);

  ASSERT_EQ(args.verbose, true);
  ASSERT_EQ(args.tMax, 7);
  ASSERT_EQ(args.tMin, 3);
  ASSERT_EQ(args.totaltime, "'4 Weeks'");
  ASSERT_EQ(args.dt, "'1 day'");
  ASSERT_EQ(args.saveEvery, 12);
  ASSERT_EQ(args.conductivity, 13.37);
  ASSERT_EQ(args.dochannels, true);
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
  ASSERT_EQ(args.netcdf, "'input.nc'");
}
