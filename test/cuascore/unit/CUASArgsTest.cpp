/**
 * File: CUASArgsTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASArgs.h"

#include "gtest/gtest.h"

#include "petsc.h"

TEST(CUASArgs, defaults) {
  std::vector input = {"CUASArgsTest.exe"};

  CUAS::CUASArgs args;
  parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

  // verbosity
  ASSERT_EQ(args.verbose, false);
  ASSERT_EQ(args.verboseSolver, false);

  // input and output
  ASSERT_EQ(args.input, "");
  ASSERT_EQ(args.output, "out.nc");
  ASSERT_EQ(args.coordinatesFile, "");
  ASSERT_EQ(args.restart, "");
  ASSERT_EQ(args.restartNoneZeroInitialGuess, true);

  // time stepping
  ASSERT_EQ(args.starttime, "");
  ASSERT_EQ(args.endtime, "");
  ASSERT_EQ(args.totaltime, "");
  ASSERT_EQ(args.dt, "");
  ASSERT_EQ(args.timeStepFile, "");

  // output behavior
  ASSERT_EQ(args.saveEvery, 0);
  ASSERT_EQ(args.saveInterval, "");
  ASSERT_EQ(args.outputSize, "normal");

  // forcing
  ASSERT_EQ(args.forcingFile, "");
  ASSERT_EQ(args.sizeOfForcingBuffer, -1);
  ASSERT_EQ(args.loopForcing, false);
  ASSERT_EQ(args.seaLevelForcing, "");

  // solver behavior
  ASSERT_EQ(args.directSolver, false);
  ASSERT_EQ(args.nonLinearIters, 0);
  ASSERT_EQ(args.timeSteppingTheta, 1.0);
  ASSERT_EQ(args.enableUDS, false);
  ASSERT_EQ(args.disableNonNegative, false);

  // channel configuration
  ASSERT_EQ(args.doChannels, false);
  ASSERT_EQ(args.selectedChannels, "noselected");
  ASSERT_EQ(args.doAllChannels, false);
  ASSERT_EQ(args.doAnyChannel, false);
  ASSERT_EQ(args.doCavity, false);
  ASSERT_EQ(args.doMelt, false);
  ASSERT_EQ(args.doCreep, false);

  // physics
  ASSERT_EQ(args.initialHead, "Nzero");
  ASSERT_EQ(args.Tmax, 20.0);
  ASSERT_EQ(args.Tmin, 0.0000001);
  ASSERT_EQ(args.Tinit, 0.2);
  ASSERT_EQ(args.disableUnconfined, false);
  ASSERT_EQ(args.conductivity, 10.0);
  ASSERT_EQ(args.flowConstant, 5e-25);
  ASSERT_EQ(args.roughnessFactor, 1.0);
  ASSERT_EQ(args.supplyMultiplier, 1.0);
  ASSERT_EQ(args.layerThickness, 0.1);
  ASSERT_EQ(args.unconfSmooth, 0.0);
  ASSERT_EQ(args.specificStorage, 0.0000982977696);
  ASSERT_EQ(args.specificYield, 0.4);
  ASSERT_EQ(args.thresholdThicknessUDS, 0.0);
  ASSERT_EQ(args.basalVelocityIce, 1e-6);
  ASSERT_EQ(args.cavityBeta, 5e-4);

  // outflow boundary conditions
  ASSERT_EQ(args.Twater, 100.0);
  ASSERT_EQ(args.dirichletBCWaterDepth, 1.0);
  ASSERT_EQ(args.blockInflow, 1);
  ASSERT_EQ(args.applyRestartChecks, false);
}

TEST(CUASArgs, allOpts) {
  std::vector input = {"CUASArgsTest.exe",
                       "--Tmax=7",
                       "--Tmin=3",
                       "--Tinit=0.4",
                       "--totaltime='4 weeks'",
                       "--dt='1 day'",
                       "--saveEvery=12",
                       "--conductivity=13.37",
                       "--doChannels",
                       "--flowConstant=2.1",
                       "--roughnessFactor=2",
                       "--supplyMultiplier=1.2",
                       "--layerThickness=1.334",
                       "--unconfSmooth=4.3",
                       "--restart='restartFile.nc'",
                       "--restartNoneZeroInitialGuess=false",
                       "--seaLevelForcing='seaForcingFile.nc'",
                       "--output='output.nc'",
                       "--input='input.nc'",
                       "--selectedChannels=creep",
                       "--verbose",
                       "--timeStepFile='timesteps.nc'",
                       "--verboseSolver",
                       "--directSolver",
                       "--forcingFile='timeforcing.nc'",
                       "--timeSteppingTheta=0.5",
                       "--sizeOfForcingBuffer=5",
                       "--loopForcing",
                       "--starttime='1 year'",
                       "--endtime='1 year 4 weeks'",
                       "--enableUDS",
                       "--thresholdThicknessUDS=0.1",
                       "--disableNonNegative",
                       "--nonLinearIters=1",
                       "--cavityBeta=3.14",
                       "--coordinatesFile='coords.nc'",
                       "--saveEvery=1",
                       "--saveInterval='1 week'",
                       "--outputSize=small",
                       "--selectedChannels=creep,melt",
                       "--initialHead=topg",
                       "--specificStorage=0.000977696",
                       "--specificYield=0.8",
                       "--basalVelocityIce=3e-6",
                       "--Twater=80.0",
                       "--dirichletBCWaterDepth=0.1",
                       "--blockInflow=2",
                       "--applyRestartChecks"};

  CUAS::CUASArgs args;
  parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

  // verbosity
  ASSERT_EQ(args.verbose, true);
  ASSERT_EQ(args.verboseSolver, true);

  // input and output
  ASSERT_EQ(args.input, "'input.nc'");
  ASSERT_EQ(args.output, "'output.nc'");
  ASSERT_EQ(args.coordinatesFile, "'coords.nc'");
  ASSERT_EQ(args.restart, "'restartFile.nc'");
  ASSERT_EQ(args.restartNoneZeroInitialGuess, false);

  // time stepping
  ASSERT_EQ(args.starttime, "'1 year'");
  ASSERT_EQ(args.endtime, "'1 year 4 weeks'");
  ASSERT_EQ(args.totaltime, "'4 weeks'");
  ASSERT_EQ(args.dt, "'1 day'");
  ASSERT_EQ(args.timeStepFile, "'timesteps.nc'");

  // output behavior
  ASSERT_EQ(args.saveEvery, 1);
  ASSERT_EQ(args.saveInterval, "'1 week'");
  ASSERT_EQ(args.outputSize, "small");

  // forcing
  ASSERT_EQ(args.forcingFile, "'timeforcing.nc'");
  ASSERT_EQ(args.sizeOfForcingBuffer, 5);
  ASSERT_EQ(args.loopForcing, true);
  ASSERT_EQ(args.seaLevelForcing, "'seaForcingFile.nc'");

  // solver behavior
  ASSERT_EQ(args.directSolver, true);
  ASSERT_EQ(args.nonLinearIters, 1);
  ASSERT_EQ(args.timeSteppingTheta, 0.5);
  ASSERT_EQ(args.enableUDS, true);
  ASSERT_EQ(args.disableNonNegative, true);

  // channel configuration
  ASSERT_EQ(args.doChannels, true);
  ASSERT_EQ(args.selectedChannels, "creep,melt");
  ASSERT_EQ(args.doAllChannels, false);
  ASSERT_EQ(args.doAnyChannel, true);
  ASSERT_EQ(args.doCavity, false);
  ASSERT_EQ(args.doMelt, true);
  ASSERT_EQ(args.doCreep, true);

  // physics
  ASSERT_EQ(args.initialHead, "topg");
  ASSERT_EQ(args.Tmax, 7);
  ASSERT_EQ(args.Tmin, 3);
  ASSERT_EQ(args.Tinit, 0.4);
  ASSERT_EQ(args.disableUnconfined, false);
  ASSERT_EQ(args.conductivity, 13.37);
  ASSERT_EQ(args.flowConstant, 2.1);
  ASSERT_EQ(args.roughnessFactor, 2);
  ASSERT_EQ(args.supplyMultiplier, 1.2);
  ASSERT_EQ(args.layerThickness, 1.334);
  ASSERT_EQ(args.unconfSmooth, 4.3);
  ASSERT_EQ(args.specificStorage, 0.000977696);
  ASSERT_EQ(args.specificYield, 0.8);
  ASSERT_EQ(args.thresholdThicknessUDS, 0.1);
  ASSERT_EQ(args.basalVelocityIce, 3e-6);
  ASSERT_EQ(args.cavityBeta, 3.14);

  // outflow boundary conditions
  ASSERT_EQ(args.Twater, 80.0);
  ASSERT_EQ(args.dirichletBCWaterDepth, 0.1);
  ASSERT_EQ(args.blockInflow, 2);
  ASSERT_EQ(args.applyRestartChecks, true);
}

TEST(CUASArgs, shorts) {
  std::vector input = {"CUASArgsTest.exe", "-x", "7", "-i", "3", "-v"};

  CUAS::CUASArgs args;
  parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

  // verbosity
  ASSERT_EQ(args.verbose, true);
  ASSERT_EQ(args.verboseSolver, false);

  // input and output
  ASSERT_EQ(args.input, "");
  ASSERT_EQ(args.output, "out.nc");
  ASSERT_EQ(args.coordinatesFile, "");
  ASSERT_EQ(args.restart, "");
  ASSERT_EQ(args.restartNoneZeroInitialGuess, true);

  // time stepping
  ASSERT_EQ(args.starttime, "");
  ASSERT_EQ(args.endtime, "");
  ASSERT_EQ(args.totaltime, "");
  ASSERT_EQ(args.dt, "");
  ASSERT_EQ(args.timeStepFile, "");

  // output behavior
  ASSERT_EQ(args.saveEvery, 0);
  ASSERT_EQ(args.saveInterval, "");
  ASSERT_EQ(args.outputSize, "normal");

  // forcing
  ASSERT_EQ(args.forcingFile, "");
  ASSERT_EQ(args.sizeOfForcingBuffer, -1);
  ASSERT_EQ(args.loopForcing, false);
  ASSERT_EQ(args.seaLevelForcing, "");

  // solver behavior
  ASSERT_EQ(args.directSolver, false);
  ASSERT_EQ(args.nonLinearIters, 0);
  ASSERT_EQ(args.timeSteppingTheta, 1.0);
  ASSERT_EQ(args.enableUDS, false);
  ASSERT_EQ(args.disableNonNegative, false);

  // channel configuration
  ASSERT_EQ(args.doChannels, false);
  ASSERT_EQ(args.selectedChannels, "noselected");
  ASSERT_EQ(args.doAllChannels, false);
  ASSERT_EQ(args.doAnyChannel, false);
  ASSERT_EQ(args.doCavity, false);
  ASSERT_EQ(args.doMelt, false);
  ASSERT_EQ(args.doCreep, false);

  // physics
  ASSERT_EQ(args.initialHead, "Nzero");
  ASSERT_EQ(args.Tmax, 7);
  ASSERT_EQ(args.Tmin, 3);
  ASSERT_EQ(args.Tinit, 0.2);
  ASSERT_EQ(args.disableUnconfined, false);
  ASSERT_EQ(args.conductivity, 10.0);
  ASSERT_EQ(args.flowConstant, 5e-25);
  ASSERT_EQ(args.roughnessFactor, 1.0);
  ASSERT_EQ(args.supplyMultiplier, 1.0);
  ASSERT_EQ(args.layerThickness, 0.1);
  ASSERT_EQ(args.unconfSmooth, 0.0);
  ASSERT_EQ(args.specificStorage, 0.0000982977696);
  ASSERT_EQ(args.specificYield, 0.4);
  ASSERT_EQ(args.thresholdThicknessUDS, 0.0);
  ASSERT_EQ(args.basalVelocityIce, 1e-6);
  ASSERT_EQ(args.cavityBeta, 5e-4);
}

TEST(CUASArgs, selectiveChannels) {
  std::vector input = {"CUASArgsTest.exe"};
  input.resize(3);

  // doChannels==true
  {
    input[1] = "--doChannels";
    {
      input[2] = "";

      CUAS::CUASArgs args;
      parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

      ASSERT_EQ(args.doAllChannels, true);
      ASSERT_EQ(args.doAnyChannel, true);
      ASSERT_EQ(args.doCreep, true);
      ASSERT_EQ(args.doMelt, true);
      ASSERT_EQ(args.doCavity, true);
    }

    {
      input[2] = "--selectedChannels=melt";

      CUAS::CUASArgs args;
      parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, true);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, true);
      ASSERT_EQ(args.doCavity, false);
    }

    {
      input[2] = "--selectedChannels=''";

      CUAS::CUASArgs args;
      parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, false);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, false);
      ASSERT_EQ(args.doCavity, false);
    }
  }

  // doChannels==false
  {
    input[1] = "";
    {
      input[2] = "";

      CUAS::CUASArgs args;
      parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, false);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, false);
      ASSERT_EQ(args.doCavity, false);
    }

    {
      input[2] = "--selectedChannels=melt,creep";

      CUAS::CUASArgs args;
      parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, true);
      ASSERT_EQ(args.doCreep, true);
      ASSERT_EQ(args.doMelt, true);
      ASSERT_EQ(args.doCavity, false);
    }

    {
      input[2] = "--selectedChannels=''";

      CUAS::CUASArgs args;
      parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

      ASSERT_EQ(args.doAllChannels, false);
      ASSERT_EQ(args.doAnyChannel, false);
      ASSERT_EQ(args.doCreep, false);
      ASSERT_EQ(args.doMelt, false);
      ASSERT_EQ(args.doCavity, false);
    }
  }
}

TEST(CUASArgs, positionalOpts) {
  std::vector input = {"CUASArgsTest.exe", "--doChannels", "--input=input.nc", "output.nc"};

  CUAS::CUASArgs args;
  parseArgs(static_cast<int>(input.size()), const_cast<char **>(input.data()), args);

  // verbosity
  ASSERT_EQ(args.verbose, false);
  ASSERT_EQ(args.verboseSolver, false);

  // input and output
  ASSERT_EQ(args.input, "input.nc");
  ASSERT_EQ(args.output, "output.nc");
  ASSERT_EQ(args.coordinatesFile, "");
  ASSERT_EQ(args.restart, "");
  ASSERT_EQ(args.restartNoneZeroInitialGuess, true);

  // time stepping
  ASSERT_EQ(args.starttime, "");
  ASSERT_EQ(args.endtime, "");
  ASSERT_EQ(args.totaltime, "");
  ASSERT_EQ(args.dt, "");
  ASSERT_EQ(args.timeStepFile, "");

  // output behavior
  ASSERT_EQ(args.saveEvery, 0);
  ASSERT_EQ(args.saveInterval, "");
  ASSERT_EQ(args.outputSize, "normal");

  // forcing
  ASSERT_EQ(args.forcingFile, "");
  ASSERT_EQ(args.sizeOfForcingBuffer, -1);
  ASSERT_EQ(args.loopForcing, false);
  ASSERT_EQ(args.seaLevelForcing, "");

  // solver behavior
  ASSERT_EQ(args.directSolver, false);
  ASSERT_EQ(args.nonLinearIters, 0);
  ASSERT_EQ(args.timeSteppingTheta, 1.0);
  ASSERT_EQ(args.enableUDS, false);
  ASSERT_EQ(args.disableNonNegative, false);

  // channel configuration
  ASSERT_EQ(args.doChannels, true);
  ASSERT_EQ(args.selectedChannels, "noselected");
  ASSERT_EQ(args.doAllChannels, true);
  ASSERT_EQ(args.doAnyChannel, true);
  ASSERT_EQ(args.doCavity, true);
  ASSERT_EQ(args.doMelt, true);
  ASSERT_EQ(args.doCreep, true);

  // physics
  ASSERT_EQ(args.initialHead, "Nzero");
  ASSERT_EQ(args.Tmax, 20.0);
  ASSERT_EQ(args.Tmin, 0.0000001);
  ASSERT_EQ(args.Tinit, 0.2);
  ASSERT_EQ(args.disableUnconfined, false);
  ASSERT_EQ(args.conductivity, 10.0);
  ASSERT_EQ(args.flowConstant, 5e-25);
  ASSERT_EQ(args.roughnessFactor, 1.0);
  ASSERT_EQ(args.supplyMultiplier, 1.0);
  ASSERT_EQ(args.layerThickness, 0.1);
  ASSERT_EQ(args.unconfSmooth, 0.0);
  ASSERT_EQ(args.specificStorage, 0.0000982977696);
  ASSERT_EQ(args.specificYield, 0.4);
  ASSERT_EQ(args.thresholdThicknessUDS, 0.0);
  ASSERT_EQ(args.basalVelocityIce, 1e-6);
  ASSERT_EQ(args.cavityBeta, 5e-4);
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  auto result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
