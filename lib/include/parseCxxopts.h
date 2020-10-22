#ifndef CUAS_PARSECXXOPTS_H
#define CUAS_PARSECXXOPTS_H

#include "cxxopts.hpp"

#include "petsc.h"

struct CUASArgs {
  PetscScalar tMax;
  PetscScalar tMin;
};

void parseArgs(int argc, char **argv, CUASArgs &args);

#endif
