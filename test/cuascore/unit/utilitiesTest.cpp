/**
 * File: utilitiesTest.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "utilities.h"

#include "gtest/gtest.h"

#include "petsc.h"

int mpiRank;
int mpiSize;

TEST(splitTest, standard) {
  constexpr auto delimiter = ';';
  auto input = "/abspath/tofile.nc;relpath/tofile.nc;file.nc";

  auto strings = CUAS::split(input, delimiter);

  ASSERT_EQ(strings.size(), 3);
  ASSERT_EQ(strings[0], "/abspath/tofile.nc");
  ASSERT_EQ(strings[1], "relpath/tofile.nc");
  ASSERT_EQ(strings[2], "file.nc");
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
