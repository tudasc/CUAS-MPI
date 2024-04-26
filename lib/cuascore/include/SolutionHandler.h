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

enum class OutputSize { SMALL, NORMAL, LARGE, XLARGE };
enum class OutputReason { NONE, INITIAL, NORMAL };

class SolutionHandler {
 private:
  // netcdf file that the SolutionHandler handles
  std::unique_ptr<NetCDFFile> file;
  // mapping of the name of a PETScGrid and the corresponding values
  // std::unordered_map<std::string, const PETScGrid *> grids;
  // mapping og the name of a variable to its attributes
  // std::unordered_map<std::string, std::vector<std::vector<std::string>>> attributes;
  // mapping of the name of a PetscScalar and the corresponding value
  // std::unordered_map<std::string, PetscScalar> scalars;
  // vector of all timesteps that get executed by the SolutionHandler
  // std::vector<int> timeSteps;
  // total number of timeSteps
  // const int Nt;
  // save every
  // const int saveEvery;
  // index of the next solution
  int nextSolution = 0;
  //
  OutputSize osize = OutputSize::SMALL;

  // defines all grids that are written in storeSolution. saveEvery determines at which timeStep the values are written
  // to the netcdf file. This is important to get the naming for the netcdf variables right.
  void defineSolution();

 public:
  // use this constructor if you want to determine the shape of the solution on your own
  // TODO do we need Nt and saveEvery
  SolutionHandler(std::string const &fileName, int dimX, int dimY, std::string const &outputSize);
  ~SolutionHandler() = default;

  void setTimeUnits(std::string const &s);
  void setCalendar(std::string const &s);

  // write the values passed as parameters to the netcdf file
  void storeInitialSetup(PETScGrid const &hydraulicHead, PETScGrid const &hydraulicTransmissivity,
                         CUASModel const &model, PETScGrid const &fluxMagnitude, PETScGrid const &melt,
                         PETScGrid const &creep, PETScGrid const &cavity, PETScGrid const &effectivePressure,
                         PETScGrid const &effectiveStorativity, PETScGrid const &effectiveTransmissivity,
                         PETScGrid const &waterSource, CUASArgs const &args);
  // write the values passed as parameters to the netcdf file
  void storeSolution(CUAS::timeSecs currTime, PETScGrid const &hydraulicHead, PETScGrid const &hydraulicTransmissivity,
                     CUASModel const &model, PETScGrid const &fluxMagnitude, PETScGrid const &melt,
                     PETScGrid const &creep, PETScGrid const &cavity, PETScGrid const &effectivePressure,
                     PETScGrid const &effectiveStorativity, PETScGrid const &effectiveTransmissivity,
                     PETScGrid const &waterSource, PetscScalar eps_inf = NC_FILL_DOUBLE,
                     PetscScalar Teps_inf = NC_FILL_DOUBLE);

  void storePETScOptions();
  // reason why we need to save this time step
  static OutputReason getOutputReason(int timeStepIndex, int numberOfTimeSteps, int saveEvery);
};
}  // namespace CUAS

#endif
