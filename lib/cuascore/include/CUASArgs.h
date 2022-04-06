#ifndef CUAS_PARSECXXOPTS_H
#define CUAS_PARSECXXOPTS_H

#include "petsc.h"

#include <string>

namespace CUAS {

struct CUASArgs {
  PetscScalar Tmax;
  PetscScalar Tmin;
  std::string totaltime;
  std::string dt;
  std::string timeStepFile;
  int saveEvery;
  PetscScalar conductivity;
  bool doAllChannels;
  bool doAnyChannel;
  bool doCavity;
  bool doMelt;
  bool doCreep;
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
  std::string seaLevelForcing;
  bool verbose;
  bool verboseSolver;
  bool directSolver;
  std::string input;
  std::string output;
  std::string outputSize;
};

void parseArgs(int argc, char **argv, CUASArgs &args);

}  // namespace CUAS

#endif
