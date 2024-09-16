/**
 * File: SolutionHandler.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_SOLUTION_HANDLER_H
#define CUAS_SOLUTION_HANDLER_H

#include "CUASArgs.h"
#include "NetCDFFile.h"
#include "timeparse.h"

#include <memory>
#include <string>
#include <vector>

namespace CUAS {

class CUASModel;
class CUASSolver;
class CUASTimeIntegrator;

enum class OutputSize { SMALL, NORMAL, LARGE, XLARGE };
enum class OutputReason { NONE, INITIAL, NORMAL };
enum class SaveStrategy { DEFAULT, TIMEINTERVAL, INDEX };

class SolutionHandler {
 public:
  // use this constructor if you want to determine the shape of the solution on your own
  // TODO should we base this on CUASModel?
  SolutionHandler(std::string const &fileName, int dimX, int dimY, std::string const &outputSize,
                  bool storeMutable = false);
  SolutionHandler(SolutionHandler &) = delete;
  SolutionHandler(SolutionHandler &&) = delete;
  SolutionHandler &operator=(SolutionHandler const &) = delete;
  SolutionHandler &operator=(SolutionHandler const &&) = delete;
  ~SolutionHandler() = default;

  // member functions
 public:
  void storeData(CUASSolver const &solver, CUASModel const &model, CUASArgs const &args, PETScGrid const &waterSource,
                 CUASTimeIntegrator const &timeIntegrator);
  void storePETScOptions();

  void setTimeUnits(std::string const &s);
  void setCalendar(std::string const &s);
  void setSaveStrategy(SaveStrategy strategy, long saveInterval, int saveEvery);

 private:
  // write the values passed as parameters to the NetCDF file
  void storeInitialSetup(CUASSolver const &solver, CUASModel const &model, PETScGrid const &waterSource,
                         CUASArgs const &args, CUASTimeIntegrator const &timeIntegrator);
  void storeCUASArgs(CUASArgs const &args);
  void storeConstantModelInformation(CUASModel const &model);
  void storeMutableModelInformation(CUASModel const &model);
  // write the values passed as parameters to the NetCDF file
  void storeSolution(CUAS::timeSecs currTime, CUASSolver const &solver, PETScGrid const &waterSource,
                     PetscScalar eps_inf = NC_FILL_DOUBLE, PetscScalar Teps_inf = NC_FILL_DOUBLE);
  void finalizeSolution();

 private:
  // NetCDF file that the SolutionHandler uses to store results
  std::unique_ptr<NetCDFFile> file;
  // index of the next solution
  int nextSolution = 0;
  // output configuration
  OutputSize osize = OutputSize::SMALL;
  // save strategy
  SaveStrategy strategy = SaveStrategy::DEFAULT;
  timeSecs saveInterval = 1;
  int saveEvery = 1;
  // defines if all fields are written per time step (usually coupled mode): true
  // or only fields which change in pure CUAS-MPI: false (default)
  const bool storeMutable;

  // member functions
 private:
  // defines all grids that are written in storeSolution. saveEvery determines at which timeStep the values are written
  // to the NetCDF file. This is important to get the naming for the NetCDF variables right.
  void defineSolution();
  // reason why we need to save this time step
  OutputReason getOutputReason(CUASTimeIntegrator const &timeIntegrator) const;
};

}  // namespace CUAS

#endif
