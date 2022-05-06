#include "CUASConstants.h"
#include "Forcing/ConstantForcing.h"
#include "Forcing/TimeForcing.h"
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
  std::unique_ptr<CUAS::Forcing> forcing = std::make_unique<CUAS::ConstantForcing>(bmelt, supplyMultiplier / SPY);
  // constant forcing: 0, empty vector, false
  auto &read = forcing->getCurrentQ().getReadHandle();
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

  std::vector<CUAS::timeSecs> time_forcing = {1, 2, 10, 15};
  auto supplyMultiplier = 1.2;
  std::unique_ptr<CUAS::Forcing> forcing =
      std::make_unique<CUAS::TimeForcing>(qs, time_forcing, supplyMultiplier / SPY, 0.0, false);

  // check interpolation
  {
    auto temp1 = 1.0 / SPY * supplyMultiplier;
    auto temp2 = 7.0 / SPY * supplyMultiplier;
    auto result = 0.5 * temp1 + 0.5 * temp2;
    auto &Q = forcing->getCurrentQ(1.5);
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
    auto &Q = forcing->getCurrentQ(3.6);
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
    auto &Q = forcing->getCurrentQ(6.5);
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
    auto &Q = forcing->getCurrentQ(10);
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
    auto &Q = forcing->getCurrentQ(100);
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
    auto &Q = forcing->getCurrentQ(0.1);
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

  std::vector<CUAS::timeSecs> time_forcing = {1, 2, 10, 15};
  auto supplyMultiplier = 1.2;
  std::unique_ptr<CUAS::Forcing> forcing =
      std::make_unique<CUAS::TimeForcing>(qs, time_forcing, supplyMultiplier / SPY, 0.0, true);

  // check interpolation
  {
    auto temp1 = 1.0 / SPY * supplyMultiplier;
    auto temp2 = 7.0 / SPY * supplyMultiplier;
    auto result = 0.5 * temp1 + 0.5 * temp2;
    auto &Q = forcing->getCurrentQ(1.5);
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
    auto &Q = forcing->getCurrentQ(3.6);
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
    auto &Q = forcing->getCurrentQ(6.5);
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
    auto &Q = forcing->getCurrentQ(10);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  // check higher currTime
  {
    auto temp1 = 7.0 / SPY * supplyMultiplier;
    auto temp2 = 5.0 / SPY * supplyMultiplier;
    auto result = 0.5 * temp1 + 0.5 * temp2;
    auto &Q = forcing->getCurrentQ(21);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
    }
  }
  {
    auto result = 7.0 / SPY * supplyMultiplier;
    auto &Q = forcing->getCurrentQ(47);
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
    auto &Q = forcing->getCurrentQ(0.1);
    auto &read1 = Q.getReadHandle();
    for (int i = 0; i < Q.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < Q.getLocalNumOfCols(); ++j) {
        ASSERT_DOUBLE_EQ(read1(i, j), result);
      }
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