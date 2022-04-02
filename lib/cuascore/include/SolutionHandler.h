#ifndef CUAS_SOLUTION_HANDLER_H
#define CUAS_SOLUTION_HANDLER_H

#include "CUASArgs.h"
#include "CUASModel.h"
#include "NetCDFFile.h"
#include "timeparse.h"

#include <memory>
#include <string>
#include <vector>

namespace CUAS {

enum class OutputSize { SMALL, NORMAL, LARGE };

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

  // write the values passed as parameters to the netcdf file
  void storeInitialSetup(PETScGrid const &hydraulicHead, PETScGrid const &hydraulicTransmissivity,
                         CUASModel const &model, PETScGrid const &melt, PETScGrid const &creep, PETScGrid const &cavity,
                         CUASArgs const &args);
  // write the values passed as parameters to the netcdf file
  void storeSolution(CUAS::timeSecs currTime, PETScGrid const &hydraulicHead, PETScGrid const &hydraulicTransmissivity,
                     CUASModel const &model, PETScGrid const &melt, PETScGrid const &creep, PETScGrid const &cavity);

  void setTimeUnits(std::string const &s);
  void setCalendar(std::string const &s);

  ~SolutionHandler();
};
}  // namespace CUAS

#endif
