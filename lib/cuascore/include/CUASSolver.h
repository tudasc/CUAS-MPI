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
  explicit CUASSolver(CUASModel *model, CUASArgs const *args, CUAS::SolutionHandler *solutionHandler = nullptr);
  CUASSolver(CUASSolver &) = delete;
  CUASSolver(CUASSolver &&) = delete;
  CUASSolver &operator=(CUASSolver const &) = delete;
  CUASSolver &operator=(CUASSolver const &&) = delete;
  ~CUASSolver() { DMDestroy(&dm); }

  // member functions
 public:
  void solve(std::vector<CUAS::timeSecs> &timeSteps);
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
  std::unique_ptr<PETScGrid> pEffective;  // same as model->pIce
  std::unique_ptr<PETScGrid> gradHeadSquared;
  std::unique_ptr<PETScGrid> fluxMagnitude;

  std::unique_ptr<PETScGrid> Seff;  //!< effective Storativity
  std::unique_ptr<PETScGrid> Teff;  //!< effective Transmissivity

  PetscScalar eps = 0.0;   //!<  \f$ max(|h^n - h^{n-1}|)/dt \f$
  PetscScalar Teps = 0.0;  //!<  \f$ max(|T^n - T^{n-1}|)/dt \f$

  // member
 private:
  CUASModel *const model;
  CUASArgs const *const args;
  CUAS::SolutionHandler *const solutionHandler;

  DM dm;
  std::unique_ptr<PETScMatrix> matA;
  std::unique_ptr<PETScGrid> bGrid;
  std::unique_ptr<PETScGrid> solGrid;

  int const numOfCols;
  int const numOfRows;

  std::unique_ptr<PETScGrid> gradMask;
  std::unique_ptr<PETScGrid> dirichletValues;

  std::unique_ptr<PETScGrid> globalIndicesBlocked;

  std::unique_ptr<PETScGrid> rateFactorIce;
  std::unique_ptr<PETScGrid> basalVelocityIce;

  // member functions
 private:
  void prepare();
  void preComputation();
  void storeData(PETScGrid const &currentQ, CUASTimeIntegrator const &timeIntegrator);
  [[nodiscard]] bool updateHeadAndTransmissivity(PETScGrid const &currentQ, CUASTimeIntegrator const &timeIntegrator);
};

}  // namespace CUAS

#endif
