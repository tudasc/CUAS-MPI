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
  explicit CUASSolver(CUASModel *model, CUASArgs const *const args, CUAS::SolutionHandler *solutionHandler = nullptr)
      : model(model), args(args), solutionHandler(solutionHandler), numOfCols(model->Ncols), numOfRows(model->Nrows) {
    nextHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    currHead = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    nextTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    currTransmissivity = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    gradMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    dirichletValues = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    // rate factor from flow law (ice rheology)
    rateFactorIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    // basal velocity of ice
    basalVelocityIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    melt = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    tmpMelt = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    creep = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    cavity = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    pEffective = std::make_unique<PETScGrid>(numOfCols, numOfRows);  // same as model->pIce
    gradHeadSquared = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    fluxMagnitude = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    Seff = std::make_unique<PETScGrid>(numOfCols, numOfRows);  // effective Storativity
    Teff = std::make_unique<PETScGrid>(numOfCols, numOfRows);  // effective Transmissivity

    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, numOfCols, numOfRows,
                 PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr, &dm);
    DMSetFromOptions(dm);
    DMSetUp(dm);
    Mat petMat;
    DMCreateMatrix(dm, &petMat);
    matA = std::make_unique<PETScMatrix>(petMat);
    bGrid = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    solGrid = std::make_unique<PETScGrid>(numOfCols, numOfRows);

    globalIndicesBlocked = std::make_unique<PETScGrid>(numOfCols, numOfRows);
    fillGlobalIndicesBlocked(*globalIndicesBlocked);
  }
  CUASSolver(CUASSolver &) = delete;
  CUASSolver(CUASSolver &&) = delete;
  ~CUASSolver() { DMDestroy(&dm); }

  void solve(std::vector<CUAS::timeSecs> &timeSteps);

  void setup();

 public:
  std::unique_ptr<PETScGrid> nextHead;  // unknown head at new time level
  std::unique_ptr<PETScGrid> currHead;  // head at the current time level
  std::unique_ptr<PETScGrid> nextTransmissivity;
  std::unique_ptr<PETScGrid> currTransmissivity;

 private:
  int const numOfCols;
  int const numOfRows;

  std::unique_ptr<PETScGrid> gradMask;
  std::unique_ptr<PETScGrid> dirichletValues;

  std::unique_ptr<PETScGrid> rateFactorIce;
  std::unique_ptr<PETScGrid> basalVelocityIce;

  std::unique_ptr<PETScGrid> melt;
  std::unique_ptr<PETScGrid> tmpMelt;
  std::unique_ptr<PETScGrid> creep;
  std::unique_ptr<PETScGrid> cavity;
  std::unique_ptr<PETScGrid> pEffective;  // same as model->pIce
  std::unique_ptr<PETScGrid> gradHeadSquared;
  std::unique_ptr<PETScGrid> fluxMagnitude;

  std::unique_ptr<PETScGrid> Seff;  //!< effective Storativity
  std::unique_ptr<PETScGrid> Teff;  //!< effective Transmissivity

  CUASModel *const model;
  CUASArgs const *const args;
  CUAS::SolutionHandler *const solutionHandler;

  std::unique_ptr<PETScGrid> globalIndicesBlocked;

  DM dm;
  std::unique_ptr<PETScMatrix> matA;
  std::unique_ptr<PETScGrid> bGrid;
  std::unique_ptr<PETScGrid> solGrid;

  PetscScalar eps = 0.0;   //!<  \f$ max(|h^n - h^{n-1}|)/dt \f$
  PetscScalar Teps = 0.0;  //!<  \f$ max(|T^n - T^{n-1}|)/dt \f$

 private:
  static std::pair<timeSecs, timeSecs> getTimeStepInformation(std::vector<CUAS::timeSecs> const &timeSteps,
                                                              int timeStepIndex);
  void prepare();
  void preComputation();
  void storeData(PETScGrid const &currentQ, timeSecs dt, std::vector<CUAS::timeSecs> const &timeSteps,
                 int timeStepIndex);
  void updateHeadAndTransmissivity(timeSecs dt, PETScGrid const &currentQ, std::vector<CUAS::timeSecs> const &timeSteps,
                                   int timeStepIndex);
};

}  // namespace CUAS

#endif
