/**
 * File: NetCDFFile.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_NETCDF_FILE_H
#define CUAS_NETCDF_FILE_H

#include "PETScGrid.h"
#include "PETScVector.h"

#include <netcdf.h>
#include <netcdf_par.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#define UNLIMITED true
#define LIMITED false

namespace CUAS {

// forward declaration of ModelReader to make it a friend class
class ModelReader;

class NetCDFFile {
  struct NetCDFVar {
    std::string name;
    int varId;
    bool isUnlimited;
  };
  friend class ModelReader;

 private:
  // map from varname to netcdfvar
  std::unordered_map<std::string, NetCDFVar> netcdfVars;
  // file name of the netcdf file
  std::string const fileName;
  // id of the netcdf file. This is important for the netcdf file
  int fileId;

  int dimX;

  int dimY;
  // dimension id for all structures. This is important for the netcdf file
  std::vector<int> dimIds;
  // gets the dimensions id the the specified dimension in the netcdf file
  int getDimId(std::string const &dimName) const;
  // gets the variable id of the specifed variable name in the netcdf file
  int getVarId(std::string const &varName) const;
  // checks if the input PETScGrid has the same number of cols/rows like the grid in the netcdf file
  bool checkDimensions(std::string const &varName, PETScGrid const &input) const;
  // checks dimensions in case of a file with unlimited dimension
  bool checkDimensionsUnlimited(std::string const &varName, PETScGrid const &input) const;
  // checks if the input PetscVector has the same dimensions as the corresponding vector in the netcdf file
  bool checkDimensions(std::string const &varName, PETScVector const &input) const;
  // checks dimension in case of a file with unlimited dimension
  bool checkDimensionsUnlimited(std::string const &varName, PETScVector const &input) const;
  // checks dimensions in when dealing with a forcing file
  bool checkDimensionsTimeForcing(std::string const &varName, std::vector<std::unique_ptr<PETScGrid>> &forcing);
  // checks if the input vector has the same dimensions as the corresponding vector in the netcdf file
  bool checkDimensions(std::string const &varName, std::vector<PetscScalar> const &input) const;
  // checks if the input vector has the same dimensions as the corresponding vector in the netcdf file
  bool checkDimensions(std::string const &varName, std::vector<long> const &input) const;
  // defines either a vector of size x or y
  void defineVectorWithOffset(std::string const &varName, int offset);
  // returns the position of the current timestep
  int getCurrentHead() const;
  // gets the dimension ids for a variable given by varname
  std::vector<int> getDimIdsForVariable(std::string const &varName) const;

 public:
  // open file
  NetCDFFile(std::string const &fileName, char mode);
  // create file
  NetCDFFile(std::string const &fileName, int dimX, int dimY);
  // gets the fileId
  int getFileId() const { return fileId; };
  // gets the number of dimensions in the netcdf file
  int getNumberOfDimensions() const;
  // checks if the file has the dimensions specified by name
  bool hasDimension(std::string const &name) const;
  // gets the x dimension of the file
  int getDimX() const { return dimX; };
  // gets the y dimension of the file
  int getDimY() const { return dimY; };
  // a generalized method to get the size of a specific dimension by the dimName
  int getDimLength(std::string const &name) const;
  // a generalized method to get the size of a specific dimension by the dimId
  int getDimLength(int dimId) const;
  // checks if the file has the variable specified by name
  bool hasVariable(std::string const &name) const;
  // returns the number of dimensions the specified variable has
  int getDimensionsForVariable(std::string const &name) const;
  // checks if the variable specified by varName has the dimensions specified by dimName
  bool variableHasDimensionByName(std::string const &varName, std::string const &dimName) const;
  // defines a PETSCGrid with the name varName in a netcdf file
  void defineGrid(std::string const &varName, bool isUnlimted = LIMITED);
  // defines a PetscVector with the name varName and size of dimension x in a netcdf file
  void defineVectorX(std::string const &varName);
  // defines a PetscVector with the name varName and size of dimension y in a netcdf file
  void defineVectorY(std::string const &varName);
  // defines a PetscScalar with the name varName in a netcdf file
  void defineScalar(std::string const &varName, bool isUnlimited = LIMITED);
  // adds an attribute to a variable
  void addAttributeToVariable(std::string const &varName, std::string const &attributeName,
                              std::string const &attributeText);

  // adds global attribute to file
  void addGlobalAttribute(std::string const &attributeName, std::string const &attributeText);
  void addGlobalAttribute(std::string const &attributeName, PetscScalar attribute);
  void addGlobalAttribute(std::string const &attributeName, bool attribute);
  void addGlobalAttribute(std::string const &attributeName, int attribute);

  // writes the grid with the name varName from a netcdf file in the input PETScGrid
  void read(std::string const &varName, PETScGrid &dest) const;
  // writes the vector with the name varName from a netcdf file in the input PETScVector
  void read(std::string const &varName, PETScVector &dest) const;
  // writes the vector with the name varName from a netcdf file in the input std::vector
  // the vector read by this function is not shared. For shared vectors use the read method with the PetscVec
  void read(std::string const &varName, std::vector<PetscScalar> &dest) const;
  //! \brief  Read in a 1d variable of type long
  void read(std::string const &varName, std::vector<long> &dest) const;
  // writes the scalar with the name varName from a netcdf file in the input PetscScalar
  void read(std::string const &varName, PetscScalar &dest) const;
  // writes the 3 dimensional forcing grid into a vector of grids
  void read(std::string const &forcingName, std::vector<std::unique_ptr<PETScGrid>> &forcing);
  // flushes cached data to file
  void sync() const;
  // writes the input PETScGrid to the variable in the netcdf file with the name varName
  void write(std::string const &varName, PETScGrid const &input, int currentTimeStep = -1);
  // writes the vector to the variable in the netcdf file with the name varName
  void write(std::string const &varName, std::vector<PetscScalar> const &input, int currentTimeStep = -1);
  // writes the input PETScVector to the variable in the netcdf file with the name varName
  // TODO const ref to input should be possible
  void write(std::string const &varName, PETScVector &input, int currentTimeStep = -1);
  // writes the input PetscScalar to the variable in the netcdf file with the name varName
  void write(std::string const &varName, PetscScalar input, int currentTimeStep = -1);
  // writes a vector with all timesteps to the netcdf file
  // TODO deprecated?
  // void writeTimeSteps(std::vector<int>);

  //! \brief Get a text attribute.
  std::string readTextAttribute(std::string const &varName, std::string const &attName);

  ~NetCDFFile();

  // Return a string containing the name of the type defined in netcdf.h
  static std::string getTypeName(nc_type varType);
  static bool isNetCDFFloatingPoint(nc_type varType);
  static bool isNetCDFIntegral(nc_type varType);
};
}  // namespace CUAS

#endif  // CUAS_NETCDF_FILE_H
