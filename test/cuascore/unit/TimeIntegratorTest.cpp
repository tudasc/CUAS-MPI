/**
 * File: TimeIntegratorTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASTimeIntegrator.h"

#include "gtest/gtest.h"

#include "petsc.h"

int mpiRank;
int mpiSize;

TEST(TimeIntegratorTest, initialization) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(11);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timesteps = integrator.getTimesteps();
  auto currentTime = integrator.getCurrentTime();
  auto dt = integrator.getCurrentDt();
  auto timeStepIndex = integrator.getTimestepIndex();

  for (int i = 0; i < inputTimeSteps.size(); ++i) {
    ASSERT_EQ(inputTimeSteps[i], timesteps[i]);
  }
  ASSERT_EQ(currentTime, 2);
  ASSERT_EQ(dt, 0);
  ASSERT_EQ(timeStepIndex, 0);
}

TEST(TimeIntegratorTest, getTimestepInformation) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(11);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(10);

  auto timesteps = integrator.getTimesteps();
  auto currentTime = integrator.getCurrentTime();
  auto dt = integrator.getCurrentDt();
  auto timeStepIndex = integrator.getTimestepIndex();

  for (int i = 0; i < inputTimeSteps.size(); ++i) {
    ASSERT_EQ(inputTimeSteps[i], timesteps[i]);
  }
  ASSERT_EQ(timestepInformation.first, 2);
  ASSERT_EQ(timestepInformation.second, 3);
  ASSERT_EQ(currentTime, 2);
  ASSERT_EQ(dt, 3);
  ASSERT_EQ(timeStepIndex, 0);
}

TEST(TimeIntegratorTest, finalizeTimestep) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(11);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(3);

  auto timesteps = integrator.getTimesteps();
  auto currentTime = integrator.getCurrentTime();
  auto dt = integrator.getCurrentDt();
  auto timeStepIndex = integrator.getTimestepIndex();

  for (int i = 0; i < inputTimeSteps.size(); ++i) {
    ASSERT_EQ(inputTimeSteps[i], timesteps[i]);
  }
  ASSERT_EQ(timestepInformation.first, 2);
  ASSERT_EQ(timestepInformation.second, 3);
  ASSERT_EQ(currentTime, 2);
  ASSERT_EQ(dt, 3);
  ASSERT_EQ(timeStepIndex, 0);

  integrator.finalizeTimestep(dt);

  currentTime = integrator.getCurrentTime();
  dt = integrator.getCurrentDt();
  timeStepIndex = integrator.getTimestepIndex();

  ASSERT_EQ(currentTime, 5);
  ASSERT_EQ(dt, 0);
  ASSERT_EQ(timeStepIndex, 1);
}

TEST(TimeIntegratorTest, finalizeIncompleteTimeStep) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(11);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(2);

  auto timesteps = integrator.getTimesteps();
  auto currentTime = integrator.getCurrentTime();
  auto dt = integrator.getCurrentDt();
  auto timeStepIndex = integrator.getTimestepIndex();

  for (int i = 0; i < inputTimeSteps.size(); ++i) {
    ASSERT_EQ(inputTimeSteps[i], timesteps[i]);
  }
  ASSERT_EQ(timestepInformation.first, 2);
  ASSERT_EQ(timestepInformation.second, 2);
  ASSERT_EQ(currentTime, 2);
  ASSERT_EQ(dt, 2);
  ASSERT_EQ(timeStepIndex, 0);

  integrator.finalizeTimestep(dt);

  currentTime = integrator.getCurrentTime();
  dt = integrator.getCurrentDt();
  timeStepIndex = integrator.getTimestepIndex();

  ASSERT_EQ(currentTime, 4);
  ASSERT_EQ(dt, 0);
  ASSERT_EQ(timeStepIndex, 0);

  timestepInformation = integrator.getTimestepInformation(2);

  currentTime = integrator.getCurrentTime();
  dt = integrator.getCurrentDt();
  timeStepIndex = integrator.getTimestepIndex();

  ASSERT_EQ(timestepInformation.first, 4);
  ASSERT_EQ(timestepInformation.second, 1);
  ASSERT_EQ(currentTime, 4);
  ASSERT_EQ(dt, 1);
  ASSERT_EQ(timeStepIndex, 0);

  integrator.finalizeTimestep(dt);

  currentTime = integrator.getCurrentTime();
  dt = integrator.getCurrentDt();
  timeStepIndex = integrator.getTimestepIndex();

  ASSERT_EQ(currentTime, 5);
  ASSERT_EQ(dt, 0);
  ASSERT_EQ(timeStepIndex, 1);
}

TEST(TimeIntegratorTest, emptyTimeStepArray) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(10);

  auto timesteps = integrator.getTimesteps();
  auto currentTime = integrator.getCurrentTime();
  auto dt = integrator.getCurrentDt();
  auto timeStepIndex = integrator.getTimestepIndex();

  for (int i = 0; i < inputTimeSteps.size(); ++i) {
    ASSERT_EQ(inputTimeSteps[i], timesteps[i]);
  }
  ASSERT_EQ(timestepInformation.first, 0);
  ASSERT_EQ(timestepInformation.second, 0);
  ASSERT_EQ(currentTime, 0);
  ASSERT_EQ(dt, 0);
  ASSERT_EQ(timeStepIndex, 0);
}

TEST(TimeIntegratorTest, oneEntryTimeStepArray) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(7);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(10);

  auto timesteps = integrator.getTimesteps();
  auto currentTime = integrator.getCurrentTime();
  auto dt = integrator.getCurrentDt();
  auto timeStepIndex = integrator.getTimestepIndex();

  for (int i = 0; i < inputTimeSteps.size(); ++i) {
    ASSERT_EQ(inputTimeSteps[i], timesteps[i]);
  }
  ASSERT_EQ(timestepInformation.first, 7);
  ASSERT_EQ(timestepInformation.second, 0);
  ASSERT_EQ(currentTime, 7);
  ASSERT_EQ(dt, 0);
  ASSERT_EQ(timeStepIndex, 0);
}

TEST(TimeIntegratorTest, twoEntriesTimeStepArray) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(18);
  inputTimeSteps.push_back(132);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(77);

  auto timesteps = integrator.getTimesteps();
  auto currentTime = integrator.getCurrentTime();
  auto dt = integrator.getCurrentDt();
  auto timeStepIndex = integrator.getTimestepIndex();

  for (int i = 0; i < inputTimeSteps.size(); ++i) {
    ASSERT_EQ(inputTimeSteps[i], timesteps[i]);
  }
  ASSERT_EQ(timestepInformation.first, 18);
  ASSERT_EQ(timestepInformation.second, 77);
  ASSERT_EQ(currentTime, 18);
  ASSERT_EQ(dt, 77);
  ASSERT_EQ(timeStepIndex, 0);
}

TEST(TimeIntegratorTest, DEATHnegativeTimeStep) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(-1);
  inputTimeSteps.push_back(3);
  inputTimeSteps.push_back(7);

  // execution
  ASSERT_DEATH(CUAS::CUASTimeIntegrator integrator(inputTimeSteps), "");
}

TEST(TimeIntegratorTest, DEATHnotContinuousSteps1) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(4);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(9);

  // execution
  ASSERT_DEATH(CUAS::CUASTimeIntegrator integrator(inputTimeSteps), "");
}

TEST(TimeIntegratorTest, DEATHnotContinuousSteps2) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(0);
  inputTimeSteps.push_back(3);
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(7);

  // execution
  ASSERT_DEATH(CUAS::CUASTimeIntegrator integrator(inputTimeSteps), "");
}

TEST(TimeIntegratorTest, DEATHnotContinuousSteps3) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(0);
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(4);

  // execution
  ASSERT_DEATH(CUAS::CUASTimeIntegrator integrator(inputTimeSteps), "");
}

TEST(TimeIntegratorTest, DEATHequalTimeStepsSteps) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(0);
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(4);

  // execution
  ASSERT_DEATH(CUAS::CUASTimeIntegrator integrator(inputTimeSteps), "");
}

TEST(TimeIntegratorTest, DEATHfinalizeTimestepSmallerZero) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(11);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(3);

  ASSERT_EQ(timestepInformation.first, 2);
  ASSERT_EQ(timestepInformation.second, 3);

  ASSERT_DEATH(integrator.finalizeTimestep(-1), "");
}

TEST(TimeIntegratorTest, DEATHfinalizeTimestepGreaterMaxNextDt) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(11);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(3);

  ASSERT_EQ(timestepInformation.first, 2);
  ASSERT_EQ(timestepInformation.second, 3);

  ASSERT_DEATH(integrator.finalizeTimestep(10), "");
}

TEST(TimeIntegratorTest, DEATHdoubleGetTimeStepInformation) {
  // setup
  std::vector<CUAS::timeSecs> inputTimeSteps;
  inputTimeSteps.push_back(2);
  inputTimeSteps.push_back(5);
  inputTimeSteps.push_back(8);
  inputTimeSteps.push_back(11);
  CUAS::CUASTimeIntegrator integrator(inputTimeSteps);

  // execution
  auto timestepInformation = integrator.getTimestepInformation(3);

  ASSERT_EQ(timestepInformation.first, 2);
  ASSERT_EQ(timestepInformation.second, 3);

  ASSERT_DEATH(timestepInformation = integrator.getTimestepInformation(3), "");
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