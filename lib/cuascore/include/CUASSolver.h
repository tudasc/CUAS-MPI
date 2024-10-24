/**
 * File: CUASSolver.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_SOLVER_H
#define CUAS_SOLVER_H

#include "CUASArgs.h"
#include "CUASModel.h"
#include "SolutionHandler.h"
#include "timeparse.h"

#include "PETScGrid.h"
#include "PETScMatrix.h"
#include "fillgrid.h"

namespace CUAS {

class CUASSolver {
 public:
  explicit CUASSolver(CUASModel *model, CUASArgs const *args, SolutionHandler *solutionHandler = nullptr);
  CUASSolver(CUASSolver &) = delete;
  CUASSolver(CUASSolver &&) = delete;
  CUASSolver &operator=(CUASSolver const &) = delete;
  CUASSolver &operator=(CUASSolver const &&) = delete;
  ~CUASSolver() { DMDestroy(&dm); }

  // member functions
 public:
  void solve(std::vector<timeSecs> &timeSteps);
  void setup();

  // member
 public:
  std::unique_ptr<PETScGrid> nextHead;  // unknown head at new time level
  std::unique_ptr<PETScGrid> currHead;  // head at the current time level
  std::unique_ptr<PETScGrid> nextTransmissivity;
  std::unique_ptr<PETScGrid> currTransmissivity;

  std::unique_ptr<PETScGrid> melt;
  std::unique_ptr<PETScGrid> tmpMelt;
  std::unique_ptr<PETScGrid> creep;
  std::unique_ptr<PETScGrid> cavity;
  std::unique_ptr<PETScGrid> pEffective;
  std::unique_ptr<PETScGrid> fluxXDir;       //!< water flux in x-direction
  std::unique_ptr<PETScGrid> fluxYDir;       //!< water flux in y-direction
  std::unique_ptr<PETScGrid> fluxMagnitude;  //!< water flux magnitude

  std::unique_ptr<PETScGrid> rateFactorIce;
  std::unique_ptr<PETScGrid> basalVelocityIce;

  std::unique_ptr<PETScGrid> Seff;       //!< effective Storativity
  std::unique_ptr<PETScGrid> Teff;       //!< effective Transmissivity
  std::unique_ptr<PETScGrid> workSpace;  //!< a temporary grid to store e.g. waterPressure

  std::unique_ptr<PETScGrid> effTransEast;   //!< effective transmissivity at the eastern (i + 1/2) boundary
  std::unique_ptr<PETScGrid> effTransWest;   //!< effective transmissivity at the western (i - 1/2) boundary
  std::unique_ptr<PETScGrid> effTransNorth;  //!< effective transmissivity at the norther (j + 1/2) boundary
  std::unique_ptr<PETScGrid> effTransSouth;  //!< effective transmissivity at the southern (j - 1/2) boundary

  std::unique_ptr<PETScGrid> dirichletValues;  // public to be used in the coupler

  PetscScalar eps = 0.0;   //!<  \f$ max(|h^n - h^{n-1}|)/dt \f$
  PetscScalar Teps = 0.0;  //!<  \f$ max(|T^n - T^{n-1}|)/dt \f$

  // member
 private:
  CUASModel *const model;
  CUASArgs const *const args;
  SolutionHandler *const solutionHandler;

  DM dm;
  std::unique_ptr<PETScMatrix> matA;
  std::unique_ptr<PETScGrid> bGrid;
  std::unique_ptr<PETScGrid> solGrid;

  int const numOfCols;
  int const numOfRows;

  std::unique_ptr<PETScGrid> globalIndicesBlocked;

  std::unique_ptr<PETScGrid> waterSource;
  std::vector<WaterSource *> waterSources;

  // member functions
 private:
  void prepare();
  void updateWaterSource(timeSecs currTime);
  void preComputation();
  void storeData(CUASTimeIntegrator const &timeIntegrator);
  [[nodiscard]] bool updateHeadAndTransmissivity(CUASTimeIntegrator const &timeIntegrator);
  void addWaterSource(WaterSource *newWaterSource);
};

}  // namespace CUAS

#endif
