#ifndef CUAS_FILE_H
#define CUAS_FILE_H

#include "CUASArgs.h"
#include "CUASModel.h"
#include "PETScGrid.h"
#include "PETScVector.h"

#include <netcdf.h>
#include <netcdf_par.h>

#include <string>
#include <vector>

namespace CUAS {

class CUASFile {
 private:
  // file name of the netcdf file
  std::string const fileName;
  // id of the netcdf file. This is important for the netcdf file
  int fileId;
  // dimension id for the grids. This is important for the netcdf file
  std::vector<int> dimIdsGrid;
  // gets the variable id of the specifed variable name in the netcdf file
  int getVarId(std::string const &varName) const;
  // checks if the input PETScGrid has the same number of cols/rows like the grid in the netcdf file
  bool checkDimensions(std::string const &varName, PETScGrid const &input);
  // checks if the input PetscVector has the same dimensions as the corresponding vector in the netcdf file
  bool checkDimensions(std::string const &varName, PETScVector const &input);
  // checks if the input vector has the same dimensions as the corresponding vector in the netcdf file
  bool checkDimensions(std::string const &varName, std::vector<PetscScalar> const &input);

 public:
  // open file
  CUASFile(std::string const &fileName, char const mode);
  // create file
  CUASFile(std::string const &fileName, int dimX, int dimY);
  // get the fileName of the netcdf file
  int getFileId() { return fileId; };
  // defines a grid with the name varName in a netcdf file
  void defineGrid(std::string const &varName);
  // defines a scalar with the name varName in a netcdf file
  void defineScalar(std::string const &varName);
  // writes the grid with the name varName from a netcdf file in the input PETScGrid
  void read(std::string const &varName, PETScGrid &input);
  // writes the vector with the name varName from a netcdf file in the input PETScVector
  void read(std::string const &varName, PETScVector &input);
  // writes the vector with the name varName from a netcdf file in the input std::vector
  // the vector read by this function is not shared. For shared vectors use the read method with the PetscVec
  void read(std::string const &varName, std::vector<PetscScalar> &input);
  // writes the scalar with the name varName from a netcdf file in the input PetscScalar
  void read(std::string const &varName, PetscScalar &input);
  // writes the input PETScGrid to the variable in the netcdf file with the name varName
  void write(std::string const &varName, PETScGrid const &input);
  // writes the input PetscScalar to the variable in the netcdf file with the name varName
  void write(std::string const &varName, PetscScalar const &input);
  // writes a vector with all timesteps to the netcdf file
  void writeTimeSteps(std::vector<int>);
  ~CUASFile();
};
}  // namespace CUAS

#endif  // CUAS_FILE_H
