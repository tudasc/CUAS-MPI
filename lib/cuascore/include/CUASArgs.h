/**
 * File: CUASArgs.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PARSECXXOPTS_H
#define CUAS_PARSECXXOPTS_H

#include "petsc.h"

#include <string>

namespace CUAS {

struct CUASArgs {
  PetscScalar Tmax;
  PetscScalar Tmin;
  PetscScalar Tinit;
  std::string starttime;
  std::string endtime;
  std::string totaltime;
  std::string dt;
  PetscScalar timeSteppingTheta;
  std::string timeStepFile;
  int saveEvery;
  std::string saveInterval;
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
  bool restartNoneZeroInitialGuess;
  PetscScalar specificStorage;  // Ss
  PetscScalar specificYield;    // Sy
  bool noSmoothMelt;
  std::string coordinatesFile;
  std::string forcingFile;
  int sizeOfForcingBuffer;
  bool loopForcing;
  PetscScalar basalVelocityIce;
  PetscScalar cavityBeta;
  std::string initialHead;
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
