#ifndef CUAS_PARSECXXOPTS_H
#define CUAS_PARSECXXOPTS_H

#include "cxxopts.hpp"

#include "petsc.h"

#include <string.h>

struct CUASArgs {
  PetscScalar tMax;
  PetscScalar tMin;
  std::string totaltime;
  std::string dt;
  int saveEvery;
  PetscScalar conductivity;
  bool dochannels;
  bool disableUnconfined;
  PetscScalar flowConstant;
  PetscScalar roughnessFactor;
  PetscScalar supplyMultiplier;
  PetscScalar layerThickness;
  PetscScalar unconfSmooth;
  std::string restart;
  PetscScalar Ssmulti;
  PetscScalar Sy;
  PetscScalar Texp;
  bool noSmoothMelt;
  bool loopForcing;
  PetscScalar basalVelocityIce;
  PetscScalar cavityBeta;
  std::string initialHead;
  std::string tempResults;
  bool version;
  std::string seaLevelForcing;
  bool verbose;
};

void parseArgs(int argc, char **argv, CUASArgs &args);

#endif