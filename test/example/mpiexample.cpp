#include "gtest/gtest.h"

#include "mpi.h"

int Factorial(int n) {
  int r = 1;
  for (int i = 1; i <= n; ++i) {
    r*=i;
  }
  return r;
}

// Tests factorial of 0.
TEST(FactorialTest, HandlesZeroInput) {
  EXPECT_EQ(Factorial(0), 1);
}

// Tests factorial of positive numbers.
TEST(FactorialTest, HandlesPositiveInput) {
  EXPECT_EQ(Factorial(1), 1);
  EXPECT_EQ(Factorial(2), 2);
  EXPECT_EQ(Factorial(3), 6);
  EXPECT_EQ(Factorial(8), 40320);
}

int main(int argc, char* argv[]) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}
