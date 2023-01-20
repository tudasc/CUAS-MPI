/**
 * File: timeparseTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "timeparse.h"

#include "Logger.h"

#include "gtest/gtest.h"

#include "petsc.h"

int mpiRank;
int mpiSize;

TEST(timeparseTest, parse10y1m2w15d5h) {
  std::string inputTime = "10 years 1 month 2 weeks 15 days 5 hours";
  CUAS::timeSecs secs = CUAS::parseTime(inputTime);
  ASSERT_EQ(secs, 320475600);
}

TEST(timeparseTest, parse50y10h) {
  std::string inputTime = "50 years 10 hours";
  CUAS::timeSecs secs = CUAS::parseTime(inputTime);
  ASSERT_EQ(secs, 1576836000);
}

TEST(timeparseTest, parse3h) {
  std::string inputTime = "3 hours";
  CUAS::timeSecs secs = CUAS::parseTime(inputTime);
  ASSERT_EQ(secs, 10800);
}

TEST(timeparseTest, parse1y2y1m) {
  std::string inputTime = "1 year 2 years 1 month";
  ASSERT_EXIT(CUAS::parseTime(inputTime), ::testing::ExitedWithCode(1),
              "timeparse.h: Wrong format! You either used a non-existing time-unit or used a time-unit multiple times "
              "in the input string! Exiting.");
}

TEST(timeparseTest, parseFractionalValue) {
  std::string inputTime = "0.5 hours";
  ASSERT_EXIT(CUAS::parseTime(inputTime), ::testing::ExitedWithCode(1),
              "timeparse.h: Wrong format for timeValue: '0.5'! Exiting.");
}

TEST(timeparseReverseTest, parse1y1mReverse) {
  std::string inputTime = "1 year 1 month";
  CUAS::timeSecs secs = 34128000;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, parse3hReverse) {
  std::string inputTime = "3 hours";
  CUAS::timeSecs secs = 10800;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, parse50y10hReverse) {
  std::string inputTime = "50 years 10 hours";
  CUAS::timeSecs secs = 1576836000;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, parse10y1m2w6d5hReverse) {
  std::string inputTime = "10 years 1 month 2 weeks 6 days 5 hours";
  CUAS::timeSecs secs = 319698000;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, parse5h28min16sReverse) {
  std::string inputTime = "5 hours 28 minutes 16 seconds";
  CUAS::timeSecs secs = 19696;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, wrongInput) {
  CUAS::timeSecs secs = -1;
  ASSERT_EXIT(CUAS::parseTime(secs), ::testing::ExitedWithCode(1),
              "timeparse.h: Invalid input. secs cannot be less than 0. Exiting.");
}

TEST(timeStepArray, getTimeStepArray) {
  auto timeStepArray = CUAS::getTimeStepArray(0, 23, 4);
  ASSERT_EQ(timeStepArray.size(), 7);
  for (int i = 0; i < timeStepArray.size(); ++i) {
    ASSERT_EQ(timeStepArray[i], 4 * i);
  }
}

TEST(timeparseReverseTest, lessThan1second) {
  CUAS::timeSecs secs = 0;
  std::string timeString = CUAS::parseTime(secs);
  ASSERT_EQ(timeString, "0 second");
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
