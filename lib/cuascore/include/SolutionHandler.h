#ifndef CUAS_SOLUTION_HANDLER_H
#define CUAS_SOLUTION_HANDLER_H

#include "CUASFile.h"

#include <map>

namespace CUAS {

class SolutionHandler {
 private:
  // netcdf file that the SolutionHandler handles
  std::unique_ptr<CUASFile> file;
  // mapping of the name of a PETScGrid and the corresponding values
  std::map<std::string, const PETScGrid *> grids;
  // mapping of the name of a PetscScalar and the corresponding value
  std::map<std::string, PetscScalar> scalars;
  // vector of all timesteps that get executed by the SolutionHandler
  std::vector<int> timeSteps;
  // total number of timeSteps
  const int Nt;

 public:
  // create a new CUASFile and netcdf file with the name fileName to write the solution to. Nt is the total number of
  // time steps and saveEvery determines at what time steps the solution should be saved.
  SolutionHandler(std::string const &fileName, const int Nt, const int saveEvery, int dimX, int dimY, int mpiRank);
  // defines all grids that are written in saveSolution. saveEvery determines at which timeStep the values are written
  // to the netcdf file. This is important to get the naming for the netcdf variables right.
  void defineSolution(int saveEvery);
  // write the values passed as parameters to the netcdf file
  void saveSolution(int timeStep, CUASArgs const &args, int rank, PETScGrid const &u, PETScGrid const &u_n,
                    CUASModel const &model, PETScGrid const &melt, PetscScalar cavity_opening);
};
}  // namespace CUAS

#endif
