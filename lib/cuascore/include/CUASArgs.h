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

class CUASArgs {
 public:
  // verbosity
  bool verbose;
  bool verboseSolver;

  // input and output
  std::string input;
  std::string output;
  std::string coordinatesFile;
  std::string restart;
  bool restartNoneZeroInitialGuess;

  // time stepping
  std::string starttime;
  std::string endtime;
  std::string totaltime;
  std::string dt;
  std::string timeStepFile;

  // output behavior
  int saveEvery;
  std::string saveInterval;
  std::string outputSize;

  // forcing
  std::string forcingFile;
  int sizeOfForcingBuffer;
  bool loopForcing;
  std::string seaLevelForcing;

  // solver behavior
  bool directSolver;
  int nonLinearIters;
  PetscScalar timeSteppingTheta;
  bool enableUDS;  // upwind scheme
  bool disableNonNegative;

  // channel configuration
  bool doChannels;
  std::string selectedChannels;
  bool doAllChannels;
  bool doAnyChannel;
  bool doCavity;
  bool doMelt;
  bool doCreep;

  // physics
  std::string initialHead;
  PetscScalar Tmax;
  PetscScalar Tmin;
  PetscScalar Tinit;
  bool disableUnconfined;
  PetscScalar conductivity;
  PetscScalar flowConstant;
  PetscScalar roughnessFactor;
  PetscScalar supplyMultiplier;
  PetscScalar layerThickness;
  PetscScalar unconfSmooth;
  PetscScalar specificStorage;  // Ss
  PetscScalar specificYield;    // Sy
  PetscScalar thresholdThicknessUDS;
  PetscScalar basalVelocityIce;
  PetscScalar cavityBeta;
};

void parseArgs(int argc, char **argv, CUASArgs &args);

}  // namespace CUAS

#endif
