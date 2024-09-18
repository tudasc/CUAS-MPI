/**
 * File: forcingTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASConstants.h"
#include "Forcing/MultiForcing.h"
#include "Forcing/SteadyForcing.h"
#include "Forcing/TimeDependentForcing.h"
#include "timeparse.h"

#include "gtest/gtest.h"

int mpiRank;
int mpiSize;

#define MPI_SIZE 6

TEST(forcingTest, constant) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  PETScGrid bmelt(20, 10);
  bmelt.setConst(1);

  auto supplyMultiplier = 1.2;
  std::unique_ptr<CUAS::Forcing> forcing = std::make_unique<CUAS::SteadyForcing>(bmelt, supplyMultiplier / SPY);
  // constant forcing: 0, empty vector, false
  auto &read = forcing->getCurrent(0).getReadHandle();
  auto &readBmelt = bmelt.getReadHandle();
  for (int i = 0; i < bmelt.getLocalNumOfRows(); ++i) {
    for (int j = 0; j < bmelt.getLocalNumOfCols(); ++j) {
      ASSERT_DOUBLE_EQ(read(i, j), readBmelt(i, j) / SPY * supplyMultiplier);
    }
  }
}

TEST(forcingTest, timeForcing) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  std::vector<std::unique_ptr<PETScGrid>> qs;
  for (int i = 0; i < 4; ++i) {
    qs.push_back(std::make_unique<PETScGrid>(20, 10));
  }
  qs[0]->setConst(1);
  qs[1]->setConst(7);
  qs[2]->setConst(5);
  qs[3]->setConst(11);

  std::vector<CUAS::timeSecs> time_forcing = {10, 20, 100, 150};
  auto supplyMultiplier = 1.2;
  std::unique_ptr<CUAS::Forcing> forcing =
      std::make_unique<CUAS::TimeDependentForcing>(qs, time_forcing, supplyMultiplier / SPY, 0.0, false);

  // check interpolation
  {
    auto temp1 = 1.0 / SPY * supplyMultiplier;
    auto temp2 = 7.0 / SPY * supplyMultiplier;
    auto result = 0.5 * temp1 + 0.5 * temp2;
    auto &Q = forcing->getCurrent(15);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  {
    auto temp1 = 7.0 / SPY * supplyMultiplier;
    auto temp2 = 5.0 / SPY * supplyMultiplier;
    auto result = 0.8 * temp1 + 0.2 * temp2;
    auto &Q = forcing->getCurrent(36);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  {
    auto temp1 = 7.0 / SPY * supplyMultiplier;
    auto temp2 = 5.0 / SPY * supplyMultiplier;
    auto result = 0.4375 * temp1 + 0.5625 * temp2;
    auto &Q = forcing->getCurrent(65);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  // check exact hit
  {
    auto result = 5.0 / SPY * supplyMultiplier;
    auto &Q = forcing->getCurrent(100);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  // check higher currTime
  {
    auto result = 11.0 / SPY * supplyMultiplier;
    auto &Q = forcing->getCurrent(1000);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  // check lower currTime
  {
    auto result = 1.0 / SPY * supplyMultiplier;
    auto &Q = forcing->getCurrent(1);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
}

TEST(forcingTest, loopForcing) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  std::vector<std::unique_ptr<PETScGrid>> qs;
  for (int i = 0; i < 4; ++i) {
    qs.push_back(std::make_unique<PETScGrid>(20, 10));
  }
  qs[0]->setConst(1);
  qs[1]->setConst(7);
  qs[2]->setConst(5);
  qs[3]->setConst(11);

  std::vector<CUAS::timeSecs> time_forcing = {10, 20, 100, 150};
  auto supplyMultiplier = 1.2;
  std::unique_ptr<CUAS::Forcing> forcing =
      std::make_unique<CUAS::TimeDependentForcing>(qs, time_forcing, supplyMultiplier / SPY, 0.0, true);

  // check interpolation within the range
  {
    auto temp1 = 1.0 / SPY * supplyMultiplier;
    auto temp2 = 7.0 / SPY * supplyMultiplier;
    auto result = 0.5 * temp1 + 0.5 * temp2;
    auto &Q = forcing->getCurrent(15);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  {
    auto temp1 = 7.0 / SPY * supplyMultiplier;
    auto temp2 = 5.0 / SPY * supplyMultiplier;
    auto result = 0.8 * temp1 + 0.2 * temp2;
    auto &Q = forcing->getCurrent(36);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  {
    auto temp1 = 7.0 / SPY * supplyMultiplier;
    auto temp2 = 5.0 / SPY * supplyMultiplier;
    auto result = 0.4375 * temp1 + 0.5625 * temp2;
    auto &Q = forcing->getCurrent(65);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  // check exact hit
  {
    auto result = 5.0 / SPY * supplyMultiplier;
    auto &Q = forcing->getCurrent(100);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  // check higher currTime to test loop forcing
  {
    auto temp1 = 7.0 / SPY * supplyMultiplier;
    auto temp2 = 5.0 / SPY * supplyMultiplier;
    auto result = 0.5 * temp1 + 0.5 * temp2;
    auto &Q = forcing->getCurrent(210);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  {
    auto result = 7.0 / SPY * supplyMultiplier;
    auto &Q = forcing->getCurrent(470);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  // check lower currTime
  {
    auto result = 1.0 / SPY * supplyMultiplier;
    auto &Q = forcing->getCurrent(1);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
}

TEST(forcingTest, multiForcing) {
  ASSERT_EQ(mpiSize, MPI_SIZE);

  CUAS::MultiForcing multiForcing(20, 10);
  {
    std::vector<std::unique_ptr<PETScGrid>> forcingStack;
    std::vector<CUAS::timeSecs> time;
    for (int i = 0; i < 3; ++i) {
      forcingStack.push_back(std::make_unique<PETScGrid>(20, 10));
      forcingStack[i]->setConst(i + 1);
      time.push_back((i + 1) * 10);
    }
    std::unique_ptr<CUAS::Forcing> forcing = std::make_unique<CUAS::TimeDependentForcing>(forcingStack, time);
    multiForcing.registerNewForcing(forcing);
  }

  {
    PETScGrid bmelt(20, 10);
    bmelt.setConst(3.14);
    std::unique_ptr<CUAS::Forcing> forcing = std::make_unique<CUAS::SteadyForcing>(bmelt);
    multiForcing.registerNewForcing(forcing);
  }

  auto &t1 = multiForcing.getCurrent(15);
  auto &t1Read = t1.getReadHandle();
  ASSERT_DOUBLE_EQ(t1Read(0, 0), 3.14 + 1.5);

  auto &t2 = multiForcing.getCurrent(23);
  auto &t2Read = t1.getReadHandle();
  ASSERT_DOUBLE_EQ(t1Read(0, 0), 3.14 + 2.3);

  auto &t3 = multiForcing.getCurrent(33);
  auto &t3Read = t1.getReadHandle();
  ASSERT_DOUBLE_EQ(t1Read(0, 0), 3.14 + 3);
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