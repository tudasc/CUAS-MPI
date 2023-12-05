/**
 * File: NetCDFFile.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "NetCDFFile.h"

#include <array>
#include <cmath>

// GRID_ROWWISE causes a deadlock in nc_put_vara_double
// probably caused by NC_COLLECTIVE
// we currently use GRID_BLOCKWISE presupposing PETScGrid data stored en bloc
#define GRID_ROWWISE 0
#define GRID_BLOCKWISE 1
#ifndef GRID_IMPLEMENTATION
#define GRID_IMPLEMENTATION GRID_BLOCKWISE
#endif

namespace CUAS {

#define SECURED_NETCDF_EXECUTION(val) check((val), #val, __func__, __FILE__, __LINE__)
template <typename T>
void check(T ncerr, const char *callee, const char *func, const char *file, int line) {
  if (ncerr) {
    std::string netcdfError = nc_strerror(ncerr);
    CUAS_ERROR("{}:{} {}\ncalling {}:\nA netcdf error occurred:\n{}\nExiting.",  //
               file, line, func, callee, netcdfError)
    std::exit(EXIT_FAILURE);
  }
}

int NetCDFFile::getVarId(std::string const &varName) const {
  int varId;
  SECURED_NETCDF_EXECUTION(nc_inq_varid(fileId, varName.c_str(), &varId));
  return varId;
}

int NetCDFFile::getDimId(std::string const &dimName) const {
  int dimId;
  SECURED_NETCDF_EXECUTION(nc_inq_dimid(fileId, dimName.c_str(), &dimId));
  return dimId;
}

int NetCDFFile::getNumberOfDimensions() const {
  int numberOfDimensions = 0;
  SECURED_NETCDF_EXECUTION(nc_inq_ndims(fileId, &numberOfDimensions));
  return numberOfDimensions;
}

bool NetCDFFile::hasDimension(std::string const &name) const {
  int dimId;
  if (int retval = nc_inq_dimid(fileId, name.c_str(), &dimId)) {
    return false;
  }
  return true;
}

int NetCDFFile::getDimLength(std::string const &name) const {
  // get id of dimension
  int dimId = getDimId(name);
  // get length of dimension
  size_t lengthOfDim;
  nc_inq_dimlen(fileId, dimId, &lengthOfDim);
  return lengthOfDim;
}

int NetCDFFile::getDimLength(int dimId) const {
  // get length of dimension
  size_t lengthOfDim;
  nc_inq_dimlen(fileId, dimId, &lengthOfDim);
  return lengthOfDim;
}

std::string NetCDFFile::getDimName(int dimId) const {
  char name[NC_MAX_NAME + 1];
  SECURED_NETCDF_EXECUTION(nc_inq_dimname(fileId, dimId, &name[0]));
  return std::string(name);
}

bool NetCDFFile::checkDimensions(std::string const &varName, PETScGrid const &input) const {
  int varId = getVarId(varName);
  size_t numOfCols;
  size_t numOfRows;
  std::array<int, 2> dimIds;
  nc_inq_vardimid(fileId, varId, dimIds.data());
  nc_inq_dimlen(fileId, dimIds[0], &numOfRows);
  nc_inq_dimlen(fileId, dimIds[1], &numOfCols);

  return numOfCols == input.getTotalNumOfCols() && numOfRows == input.getTotalNumOfRows();
}

bool NetCDFFile::checkDimensionsUnlimited(std::string const &varName, PETScGrid const &input) const {
  int varId = getVarId(varName);
  size_t numOfCols;
  size_t numOfRows;
  std::array<int, 3> dimIds;
  nc_inq_vardimid(fileId, varId, dimIds.data());
  nc_inq_dimlen(fileId, dimIds[1], &numOfRows);
  nc_inq_dimlen(fileId, dimIds[2], &numOfCols);
  size_t unlimited;
  nc_inq_dimlen(fileId, dimIds[0], &unlimited);

  return numOfCols == input.getTotalNumOfCols() && numOfRows == input.getTotalNumOfRows();
}

bool NetCDFFile::checkDimensionsTimeForcing(std::string const &varName,
                                            std::vector<std::unique_ptr<PETScGrid>> &forcing) {
  int varId = getVarId(varName);
  size_t time;
  size_t numOfCols;
  size_t numOfRows;
  std::array<int, 3> dimIds;

  nc_inq_vardimid(fileId, varId, dimIds.data());
  nc_inq_dimlen(fileId, dimIds[0], &time);
  nc_inq_dimlen(fileId, dimIds[1], &numOfRows);
  nc_inq_dimlen(fileId, dimIds[2], &numOfCols);

  if (time == forcing.size()) {
    for (const std::unique_ptr<PETScGrid> &grid : forcing) {
      if (numOfCols != grid->getTotalNumOfCols() && numOfRows != grid->getTotalNumOfRows()) {
        return false;
      }
      return true;
    }
  }
  return false;
}

bool NetCDFFile::checkDimensions(std::string const &varName, PETScVector const &input) const {
  int varId = getVarId(varName);
  size_t size;
  int dimId;
  nc_inq_vardimid(fileId, varId, &dimId);
  nc_inq_dimlen(fileId, dimId, &size);

  return size == input.getSize();
}

bool NetCDFFile::checkDimensionsUnlimited(std::string const &varName, PETScVector const &input) const {
  int varId = getVarId(varName);
  size_t size;
  std::vector<int> dimIds(2);
  nc_inq_vardimid(fileId, varId, dimIds.data());
  nc_inq_dimlen(fileId, dimIds[1], &size);

  return size == input.getSize();
}

bool NetCDFFile::checkDimensions(std::string const &varName, std::vector<PetscScalar> const &input) const {
  int varId = getVarId(varName);
  size_t size;
  int dimId;
  nc_inq_vardimid(fileId, varId, &dimId);
  nc_inq_dimlen(fileId, dimId, &size);

  return size == input.size();
}

bool NetCDFFile::checkDimensions(std::string const &varName, std::vector<long> const &input) const {
  int varId = getVarId(varName);
  size_t size;
  int dimId;
  nc_inq_vardimid(fileId, varId, &dimId);
  nc_inq_dimlen(fileId, dimId, &size);

  return size == input.size();
}

void NetCDFFile::defineVectorWithOffset(std::string const &varName, int offset) {
  // maybe we need a different dimIdGrids array for unlimited dimensions
  int varId;
  SECURED_NETCDF_EXECUTION(nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 1, &dimIds[offset], &varId));

  netcdfVars[varName] = NetCDFVar{varName, varId, LIMITED};
}

int NetCDFFile::getCurrentHead() const {
  int unlimId;
  nc_inq_unlimdim(fileId, &unlimId);
  size_t unlimitedLength;
  nc_inq_dimlen(fileId, unlimId, &unlimitedLength);
  return unlimitedLength;
}

std::vector<int> NetCDFFile::getDimIdsForVariable(std::string const &varName) const {
  int varId = getVarId(varName);
  int numberOfDimensions = getDimensionsForVariable(varName);
  std::vector<int> dimIdsForVar(numberOfDimensions);
  SECURED_NETCDF_EXECUTION(nc_inq_vardimid(fileId, varId, dimIdsForVar.data()));
  return dimIdsForVar;
}

NetCDFFile::NetCDFFile(std::string const &fileName, char mode) : fileName(fileName) {
  switch (mode) {
    // openfile for reading
    case 'r': {
      SECURED_NETCDF_EXECUTION(nc_open_par(fileName.c_str(), NC_NOWRITE, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId));
      break;
    }
    // openfile for writing
    case 'w': {
      SECURED_NETCDF_EXECUTION(nc_open_par(fileName.c_str(), NC_WRITE, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId));
      break;
    }
    default: {
      CUAS_ERROR(
          "NetCDFFile.cpp: NetCDFFile() in write-mode: Please enter 'r' or 'w' as a mode to read or write the netcdf "
          "file. Exiting.")
      exit(1);
    }
  }
  nc_enddef(fileId);

  // get dimId of unlimtedDimension
  int dimIdUnlimited;
  size_t sizeOfUnlimitedDim;
  // TODO check
  if (int retval = nc_inq_dimid(fileId, "time", &dimIdUnlimited)) {
    std::string netcdfError = nc_strerror(retval);
    CUAS_INFO_RANK0("NetCDFFile.cpp: The netcdf file doesn't have unlimited dimensions.")
  }
  // get lenght of dimension timeSteps
  else if (int retval = nc_inq_dimlen(fileId, dimIdUnlimited, &sizeOfUnlimitedDim)) {
    std::string netcdfError = nc_strerror(retval);
    CUAS_ERROR("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError + "Exiting.")
    exit(1);
  }

  // fixme: Here is the problem: all files are expected to have "x" and "y"! This is not the case
  //        for, e.g. a file that just contains time steps.
  dimX = getDimLength("x");
  dimY = getDimLength("y");
}

NetCDFFile::NetCDFFile(std::string const &fileName, int dimX, int dimY) : fileName(fileName), dimX(dimX), dimY(dimY) {
  SECURED_NETCDF_EXECUTION(
      nc_create_par(fileName.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId));

  // define dimensions for grids
  dimIds.reserve(3);
  nc_def_dim(fileId, "time", NC_UNLIMITED, &dimIds[0]);
  nc_def_dim(fileId, "y", dimY, &dimIds[1]);
  nc_def_dim(fileId, "x", dimX, &dimIds[2]);

  nc_enddef(fileId);
}

bool NetCDFFile::hasVariable(std::string const &name) const {
  int varId;
  if (int retval = nc_inq_varid(fileId, name.c_str(), &varId)) {
    return false;
  }
  return true;
}

int NetCDFFile::getDimensionsForVariable(std::string const &name) const {
  int varId = getVarId(name);
  int numberOfDimensions;
  SECURED_NETCDF_EXECUTION(nc_inq_varndims(fileId, varId, &numberOfDimensions));
  return numberOfDimensions;
}

bool NetCDFFile::variableHasDimensionByName(std::string const &varName, std::string const &dimName) const {
  int dimId;

  // Return early, if the requested dimension does not exist at all
  if (NC_NOERR != nc_inq_dimid(fileId, dimName.c_str(), &dimId)) {
    return false;
  }

  int const numberOfDimensions = getDimensionsForVariable(varName);
  std::vector<int> dimIdsForVar(numberOfDimensions);
  int varId = getVarId(varName);
  SECURED_NETCDF_EXECUTION(nc_inq_vardimid(fileId, varId, dimIdsForVar.data()));
  // checks if dimIdsForVar contains dimId
  return std::find(std::begin(dimIdsForVar), std::end(dimIdsForVar), dimId) != std::end(dimIdsForVar);
}

void NetCDFFile::defineGrid(std::string const &varName, bool isUnlimited) {
  int varId;
  if (isUnlimited) {
    SECURED_NETCDF_EXECUTION(nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 3, dimIds.data(), &varId));
  } else {
    // maybe we need a different dimIdGrids array for unlimited dimensions
    SECURED_NETCDF_EXECUTION(nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 2, &dimIds[1], &varId));
  }

  netcdfVars[varName] = NetCDFVar{varName, varId, isUnlimited};
}

void NetCDFFile::defineVectorX(std::string const &varName) { defineVectorWithOffset(varName, 2); }

void NetCDFFile::defineVectorY(std::string const &varName) { defineVectorWithOffset(varName, 1); }

void NetCDFFile::defineScalar(std::string const &varName, bool isUnlimited) {
  int varId;
  if (isUnlimited) {
    SECURED_NETCDF_EXECUTION(nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 1, &dimIds[0], &varId));
  } else {
    SECURED_NETCDF_EXECUTION(nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 0, nullptr, &varId));
  }

  netcdfVars[varName] = NetCDFVar{varName, varId, isUnlimited};
}

void NetCDFFile::addAttributeToVariable(std::string const &varName, std::string const &attributeName,
                                        std::string const &attributeText) {
  auto varId = getVarId(varName);
  SECURED_NETCDF_EXECUTION(
      nc_put_att_text(fileId, varId, attributeName.c_str(), strlen(attributeText.c_str()), attributeText.c_str()));
}

void NetCDFFile::addAttributeToVariable(std::string const &varName, std::string const &attributeName,
                                        PetscScalar attribute) {
  auto varId = getVarId(varName);
  double petscToDouble(attribute);  // PetscScalar is double, but this could change
  SECURED_NETCDF_EXECUTION(nc_put_att_double(fileId, varId, attributeName.c_str(), NC_DOUBLE, 1, &petscToDouble));
}

void NetCDFFile::addAttributeToVariable(std::string const &varName, std::string const &attributeName, int attribute) {
  auto varId = getVarId(varName);
  SECURED_NETCDF_EXECUTION(nc_put_att_int(fileId, varId, attributeName.c_str(), NC_INT, 1, &attribute));
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, std::string const &attributeText) {
  SECURED_NETCDF_EXECUTION(
      nc_put_att_text(fileId, NC_GLOBAL, attributeName.c_str(), strlen(attributeText.c_str()), attributeText.c_str()));
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, PetscScalar attribute) {
  SECURED_NETCDF_EXECUTION(nc_put_att_double(fileId, NC_GLOBAL, attributeName.c_str(), NC_DOUBLE, 1, &attribute));
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, bool attribute) {
  int boolToInt(attribute);
  SECURED_NETCDF_EXECUTION(nc_put_att_int(fileId, NC_GLOBAL, attributeName.c_str(), NC_INT, 1, &boolToInt));
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, int attribute) {
  SECURED_NETCDF_EXECUTION(nc_put_att_int(fileId, NC_GLOBAL, attributeName.c_str(), NC_INT, 1, &attribute));
}

void NetCDFFile::read(std::string const &varName, PETScGrid &dest) const {
  auto varId = getVarId(varName);
  int numOfDims;
  nc_inq_varndims(fileId, varId, &numOfDims);
  bool isUnlimited;
  if (numOfDims == 2) {
    isUnlimited = LIMITED;
  } else if (numOfDims == 3) {
    isUnlimited = UNLIMITED;
  } else {
    CUAS_ERROR("NetCDFFile.cpp: read('{}') with PETScGrid: The grid is not limited or unlimited! Exiting.", varName)
    exit(1);
  }
  // determine start and count of values
  if (isUnlimited) {
    if (!checkDimensionsUnlimited(varName, dest)) {
      CUAS_ERROR("NetCDFFile.cpp: read('{}') with PETScGrid: sizes do not fit! Exiting.", varName)
      exit(1);
    }
    int headIndex = getCurrentHead() - 1;

    std::array<size_t, 3> start;
    std::array<size_t, 3> count;

    auto grid = dest.getWriteHandle();
    auto gridFromNetcdf = grid.getRaw();

#if GRID_IMPLEMENTATION == GRID_ROWWISE
    start[0] = headIndex;
    start[2] = dest.getCornerX();
    count[0] = 1;
    count[1] = 1;
    count[2] = dest.getLocalNumOfCols();

    for (int i = 0; i < dest.getLocalNumOfRows(); ++i) {
      start[1] = dest.getCornerY() + i;
      SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0]));
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = headIndex;
    start[1] = dest.getCornerY();
    start[2] = dest.getCornerX();
    count[0] = 1;
    count[1] = dest.getLocalNumOfRows();
    count[2] = dest.getLocalNumOfCols();

    SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0]));
#endif
  } else {
    if (!checkDimensions(varName, dest)) {
      CUAS_ERROR("NetCDFFile.cpp: read('{}') with PETScGrid: sizes do not fit! Exiting.", varName)
      exit(1);
    }

    std::array<size_t, 2> start;
    std::array<size_t, 2> count;

    auto grid = dest.getWriteHandle();
    auto gridFromNetcdf = grid.getRaw();

#if GRID_IMPLEMENTATION == GRID_ROWWISE
    start[1] = dest.getCornerX();
    count[0] = 1;
    count[1] = dest.getLocalNumOfCols();

    for (int i = 0; i < dest.getLocalNumOfRows(); ++i) {
      start[0] = dest.getCornerY() + i;
      SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0]));
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = dest.getCornerY();
    start[1] = dest.getCornerX();
    count[0] = dest.getLocalNumOfRows();
    count[1] = dest.getLocalNumOfCols();

    SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0]));
#endif
  }
}

void NetCDFFile::read(std::string const &varName, PETScVector &dest) const {
  auto varId = getVarId(varName);
  int numOfDims;
  nc_inq_varndims(fileId, varId, &numOfDims);
  bool isUnlimited;
  if (numOfDims == 1) {
    isUnlimited = LIMITED;
  } else if (numOfDims == 2) {
    isUnlimited = UNLIMITED;
  } else {
    CUAS_ERROR("NetCDFFile.cpp: read('{}') with PETScVector: The grid is not limited or unlimited! Exiting.", varName)
    exit(1);
  }

  if (isUnlimited) {
    if (!checkDimensionsUnlimited(varName, dest)) {
      CUAS_ERROR("NetCDFFile.cpp: read('{}') with PETScVector: size does not fit! Exiting.", varName)
      exit(1);
    }

    int headIndex = getCurrentHead() - 1;

    int startIndex;
    int endIndex;
    dest.getOwnershipRange(startIndex, endIndex);

    std::array<size_t, 2> start;
    std::array<size_t, 2> count;

    start[0] = headIndex;
    start[1] = startIndex;
    count[0] = 1;
    count[1] = endIndex - startIndex;

    int localSizeOfInput = count[1];
    std::vector<PetscScalar> vecFromNetcdf;
    vecFromNetcdf.resize(localSizeOfInput);

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, start.data(), count.data(), vecFromNetcdf.data()));

    for (int i = startIndex; i < endIndex; ++i) {
      dest.setValue(i, vecFromNetcdf[i - startIndex]);
    }

    dest.assemble();
  } else {
    if (!checkDimensions(varName, dest)) {
      CUAS_ERROR("NetCDFFile.cpp: read('{}') with PETScVector: size does not fit! Exiting.", varName)
      exit(1);
    }

    int startIndex;
    int endIndex;
    dest.getOwnershipRange(startIndex, endIndex);

    size_t count = endIndex - startIndex;
    size_t start = startIndex;

    int localSizeOfInput = count;
    std::vector<PetscScalar> vecFromNetcdf;
    vecFromNetcdf.resize(localSizeOfInput);

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, &start, &count, vecFromNetcdf.data()));

    for (int i = startIndex; i < endIndex; ++i) {
      dest.setValue(i, vecFromNetcdf[i - startIndex]);
    }

    dest.assemble();
  }
}

void NetCDFFile::read(std::string const &varName, std::vector<PetscScalar> &dest) const {
  int varId = getVarId(varName);

  if (!checkDimensions(varName, dest)) {
    CUAS_ERROR("NetCDFFile.cpp: read('{}') with std::vector<PetscScalar>: size does not fit! Exiting.", varName)
    exit(1);
  }

  // the std::vector is not distributed. Therefore, the sequential get method from netcdf is called
  SECURED_NETCDF_EXECUTION(nc_get_var_double(fileId, varId, dest.data()));
}

void NetCDFFile::read(std::string const &varName, std::vector<long> &dest) const {
  int varId = getVarId(varName);
  int numOfDims;
  int dimId;
  size_t dimLen;
  nc_type varType;

  nc_inq_varndims(fileId, varId, &numOfDims);
  if (numOfDims != 1) {
    // todo:   "variable '%s' in '%s' should to have 1 dimension (got %d)"
    CUAS_ERROR("NetCDFFile.cpp: read('{}') with std::vector<long>: Wrong number of dimensions! Exiting.", varName)
    exit(1);
  }

  nc_inq_vardimid(fileId, varId, &dimId);
  nc_inq_dimlen(fileId, dimId, &dimLen);

  dest.resize(dimLen);  // memory allocation happens here

  if (!checkDimensions(varName, dest)) {
    CUAS_ERROR("NetCDFFile.cpp: read with std::vector<long>: size does not fit! Exiting.")
    exit(1);
  }

  SECURED_NETCDF_EXECUTION(nc_inq_vartype(fileId, varId, &varType));

  // The std::vector is not distributed. Therefore, the sequential get method from netcdf is called.

  if (isNetCDFFloatingPoint(varType)) {
    std::vector<double> tmp(dest.size());
    double fractpart, intpart;

    // read as double and check for truncation error during cast to long
    SECURED_NETCDF_EXECUTION(nc_get_var_double(fileId, varId, tmp.data()));
    for (std::vector<long>::size_type i = 0; i != dest.size(); ++i) {
      fractpart = modf(tmp[i], &intpart);  // fractional part with sign
      dest[i] = (long)intpart;
      if (fractpart != 0.0) {
        // We don't want the code to behave in an unexpected way with truncated values, thus error instead of warn.
        CUAS_ERROR("NetCDFFile.cpp: read('{}') with std::vector<long>: Loss of decimal places reading as long!",
                   varName)
        exit(1);
      }
    }
  } else {
    // NC_LONG would be the correct type, but NC_INT, NC_SHORT and probably other types work as well
    if (!isNetCDFIntegral(varType)) {
      CUAS_WARN(
          "NetCDFFile.cpp: read('{}') with std::vector<long>: Unsupported variable type {}. Try reading as NC_LONG",
          varName, getTypeName(varType))
    }
    SECURED_NETCDF_EXECUTION(nc_get_var_long(fileId, varId, dest.data()));
  }
}

void NetCDFFile::read(std::string const &varName, PetscScalar &dest) const {
  int varId = getVarId(varName);
  int numOfDims;
  nc_inq_varndims(fileId, varId, &numOfDims);
  bool isUnlimited;
  if (numOfDims == 0) {
    isUnlimited = LIMITED;
  } else if (numOfDims == 1) {
    isUnlimited = UNLIMITED;
  } else {
    CUAS_ERROR("NetCDFFile.cpp: read('{}') with PETScScalar: The grid is not limited or unlimited! Exiting.", varName)
    exit(1);
  }

  if (isUnlimited) {
    int headIndex = getCurrentHead() - 1;

    size_t start[1];
    size_t count[1];

    start[0] = headIndex;
    count[0] = 1;

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, start, count, &dest));
  }

  else {
    SECURED_NETCDF_EXECUTION(nc_get_var_double(fileId, varId, &dest));
  }
}

void NetCDFFile::read(std::string const &forcingName, std::vector<std::unique_ptr<PETScGrid>> &forcing) {
  auto varId = getVarId(forcingName);
  int numOfDims;
  nc_inq_varndims(fileId, varId, &numOfDims);
  if (numOfDims != 3) {
    CUAS_ERROR("NetCDFFile.cpp: read() with std::vector of PETScGrid: The timeforcing is not 3 dimensional! Exiting.")
    exit(1);
  }

  // determine start and count of values
  if (!checkDimensionsTimeForcing(forcingName, forcing)) {
    CUAS_ERROR("NetCDFFile.cpp: read() with std::vector of PETScGrid: sizes do not fit! Exiting.")
    exit(1);
  }
  // int headIndex = getCurrentHead() - 1;

  for (int i = 0; i < forcing.size(); ++i) {
    PETScGrid &currentGrid = *forcing[i];

    std::array<size_t, 3> start;
    std::array<size_t, 3> count;

    auto grid = currentGrid.getWriteHandle();
    auto gridFromNetcdf = grid.getRaw();

    start[0] = i;
    start[1] = currentGrid.getCornerY();
    start[2] = currentGrid.getCornerX();
    count[0] = 1;
    count[1] = currentGrid.getLocalNumOfRows();
    count[2] = currentGrid.getLocalNumOfCols();

    SECURED_NETCDF_EXECUTION(nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0]));
  }
}

void NetCDFFile::write(std::string const &varName, PETScGrid const &input, int currentTimeStep) {
  auto currentVar = netcdfVars[varName];
  int varId = currentVar.varId;
  bool isUnlimited = currentVar.isUnlimited;

  if (isUnlimited) {
    if (!checkDimensionsUnlimited(varName, input)) {
      CUAS_ERROR("NetCDFFile.cpp: write with PETScGrid: sizes do not fit! Exiting.")
      exit(1);
    }
    // determine start and count of values
    std::array<size_t, 3> start;
    std::array<size_t, 3> count;

    auto &grid = input.getReadHandle();
    auto gridFromNetcdf = grid.getRaw();

#if GRID_IMPLEMENTATION == GRID_ROWWISE
    start[0] = currentTimeStep;
    start[2] = input.getCornerX();
    count[0] = 1;
    count[1] = 1;
    count[2] = input.getLocalNumOfCols();

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    for (int i = 0; i < input.getLocalNumOfRows(); ++i) {
      start[1] = input.getCornerY() + i;
      SECURED_NETCDF_EXECUTION(nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0]));
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = currentTimeStep;
    start[1] = input.getCornerY();
    start[2] = input.getCornerX();
    count[0] = 1;
    count[1] = input.getLocalNumOfRows();
    count[2] = input.getLocalNumOfCols();

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0]));
#endif
  } else {
    if (!checkDimensions(varName, input)) {
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScGrid: sizes do not fit! Exiting.")
      exit(1);
    }

    std::array<size_t, 2> start;
    std::array<size_t, 2> count;

    auto &grid = input.getReadHandle();
    auto gridFromNetcdf = grid.getRaw();

#if GRID_IMPLEMENTATION == GRID_ROWWISE
    start[1] = input.getCornerX();
    count[0] = 1;
    count[1] = input.getLocalNumOfCols();

    for (int i = 0; i < input.getLocalNumOfRows(); ++i) {
      start[0] = input.getCornerY() + i;
      SECURED_NETCDF_EXECUTION(nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0]));
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = input.getCornerY();
    start[1] = input.getCornerX();
    count[0] = input.getLocalNumOfRows();
    count[1] = input.getLocalNumOfCols();

    SECURED_NETCDF_EXECUTION(nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0]));
#endif
  }
}

void NetCDFFile::write(std::string const &varName, PETScVector &input, const int currentTimeStep) {
  auto currentVar = netcdfVars[varName];
  int varId = currentVar.varId;
  bool isUnlimited = currentVar.isUnlimited;
  if (isUnlimited) {
    CUAS_ERROR("writing unlimited PETScVectors is not supported")
    exit(1);
    /*
    if (!checkDimensionsUnlimited(varName, input)) {
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScVector: size does not fit! Exiting.");
      exit(1);
    }

    int startIndex;
    int endIndex;
    input.getOwnershipRange(startIndex, endIndex);

    std::array<size_t, 2> start;
    std::array<size_t, 2> count;

    start[0] = currentTimeStep;
    start[1] = startIndex;
    count[0] = 1;
    count[1] = endIndex - startIndex;

    int localSizeOfInput = count[1];

    // TODO this seams to be an unnecessary copy
    std::vector<PetscScalar> vecFromNetcdf;
    {
      vecFromNetcdf.resize(localSizeOfInput);

      std::vector<int> indicies;
      indicies.resize(localSizeOfInput);

      for (int i = 0; i < localSizeOfInput; ++i) {
        indicies[i] = startIndex + i;
      }

      VecGetValues(input.getRaw(), localSizeOfInput, indicies.data(), vecFromNetcdf.data());
    }

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_put_vara_double(fileId, varId, start.data(), count.data(), vecFromNetcdf.data()));*/
  } else {
    if (!checkDimensions(varName, input)) {
      CUAS_ERROR("NetCDFFile.cpp: write() with PETScVector: size does not fit! Exiting.")
      exit(1);
    }

    int startIndex;
    int endIndex;
    input.getOwnershipRange(startIndex, endIndex);

    size_t count = endIndex - startIndex;
    size_t start = startIndex;

    // TODO this seams to be an unnecessary copy
    std::vector<PetscScalar> vecFromNetcdf;
    {
      vecFromNetcdf.resize(count);

      std::vector<int> indicies;
      indicies.resize(count);

      for (int i = 0; i < count; ++i) {
        indicies[i] = startIndex + i;
      }

      VecGetValues(input.getRaw(), count, indicies.data(), vecFromNetcdf.data());
    }

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_put_vara_double(fileId, varId, &start, &count, vecFromNetcdf.data()));
  }
}

void NetCDFFile::write(std::string const &varName, std::vector<PetscScalar> const &input, const int currentTimeStep) {
  // Y.F. this is always limited as PetscVec
  // Y.F. btw. we never store PetscVec :-P
  // Y.F. i guess we are able to recapture this from the old method writeTimeStep (sse below)
  // Y.F. after writeTimeStep was recaptured we can delete the old code below
  // Y.F. if we need writeTimeStep as we have used it before again, we will use this new method
  auto currentVar = netcdfVars[varName];
  int varId = currentVar.varId;
  bool isUnlimited = currentVar.isUnlimited;
  if (isUnlimited) {
    CUAS_ERROR("writing unlimited PETScVectors is not supported")
    exit(1);
  } else {
    if (!checkDimensions(varName, input)) {
      CUAS_ERROR("NetCDFFile.cpp: write() with std::vector: size does not fit! Exiting.")
      exit(1);
    }

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_put_var_double(fileId, varId, input.data()));
  }
}

// change this to scalars that are written as an unlimited var
/*void NetCDFFile::writeTimeSteps(std::vector<int> timeSteps) {
  nc_redef(fileId);
  std::vector<int> dimTimeSteps(1);
  SECURED_NETCDF_EXECUTION(nc_def_dim(fileId, "numOfTimeSteps", timeSteps.size(), &dimTimeSteps[0]));
  int varId;
  SECURED_NETCDF_EXECUTION(nc_def_var(fileId, "time", NC_INT, 1, dimTimeSteps.data(), &varId));
  nc_enddef(fileId);
  SECURED_NETCDF_EXECUTION(nc_put_var_int(fileId, varId, timeSteps.data()));
}*/

void NetCDFFile::write(std::string const &varName, PetscScalar input, const int currentTimeStep) {
  auto currentVar = netcdfVars[varName];
  auto varId = currentVar.varId;
  auto isUnlimited = currentVar.isUnlimited;
  if (isUnlimited) {
    std::array<size_t, 1> start;
    std::array<size_t, 1> count;

    start[0] = currentTimeStep;
    count[0] = 1;

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    SECURED_NETCDF_EXECUTION(nc_put_vara_double(fileId, varId, start.data(), count.data(), &input));
  } else {
    SECURED_NETCDF_EXECUTION(nc_put_var_double(fileId, varId, &input));
  }
}

void NetCDFFile::sync() const { SECURED_NETCDF_EXECUTION(nc_sync(fileId)); }

NetCDFFile::~NetCDFFile() { SECURED_NETCDF_EXECUTION(nc_close(fileId)); }

std::string NetCDFFile::readTextAttribute(const std::string &varName, const std::string &attName) {
  int varId = getVarId(varName);
  nc_type ncType;
  size_t attLen = 0;
  std::string result;

  int retval = nc_inq_atttype(fileId, varId, attName.c_str(), &ncType);
  if (retval == NC_NOERR) {
    //
    // attribute found by name
    //
    if (ncType == NC_CHAR) {
      int stat = nc_inq_attlen(fileId, varId, attName.c_str(), &attLen);
      if (stat != NC_NOERR) {
        result = "";
      } else {
        // we have a text attribute with a given length yeah!
        std::vector<char> buffer(attLen + 1, 0);
        stat = nc_get_att_text(fileId, varId, attName.c_str(), buffer.data());
        result = (stat == NC_NOERR) ? buffer.data() : "";
      }
    } else {
      //
      CUAS_ERROR("NetCDFFile.cpp: readTextAttribute(): Invalid attribute type. Exiting.")
      exit(1);
    }
  } else if (retval == NC_ENOTATT) {
    //
    // attribute not found
    //
    CUAS_WARN("NetCDFFile.cpp: readTextAttribute(): attribute '" + attName + "' not found", attName)
    result = "";
  } else {
    //
    // something went wrong
    //
    std::string netcdfError = nc_strerror(retval);
    CUAS_ERROR("NetCDFFile.cpp: readTextAttribute(): A netcdf error occurred: " + netcdfError + "Exiting.")
    exit(1);
  }

  return result;
}

// modified after ncType.cpp|h from netcdf4 c++ API
std::string NetCDFFile::getTypeName(nc_type varType) {
  switch (varType) {
    case NC_BYTE:
      return "NC_BYTE";
    case NC_UBYTE:
      return "NC_UBYTE";
    case NC_CHAR:
      return "NC_CHAR";
    case NC_SHORT:
      return "NC_SHORT";
    case NC_USHORT:
      return "NC_USHORT";
    case NC_INT:
      return "NC_INT";
    case NC_UINT:
      return "NC_UINT";
    case NC_INT64:
      return "NC_INT64";
    case NC_UINT64:
      return "NC_UINT64";
    case NC_FLOAT:
      return "NC_FLOAT";
    case NC_DOUBLE:
      return "NC_DOUBLE";
    case NC_STRING:
      return "NC_STRING";
    case NC_VLEN:
      return "NC_VLEN";
    case NC_OPAQUE:
      return "NC_OPAQUE";
    case NC_ENUM:
      return "NC_ENUM";
    case NC_COMPOUND:
      return "NC_COMPOUND";
    default:
      // we should never get here!
      return "NOT_A_KNOWN_TYPE";
  }
}

bool NetCDFFile::isNetCDFFloatingPoint(nc_type varType) { return varType == NC_FLOAT || varType == NC_DOUBLE; }
bool NetCDFFile::isNetCDFIntegral(nc_type varType) {
  return (varType == NC_LONG) || (varType == NC_INT) || (varType == NC_SHORT);
}

// Use nc_copy_var (int ncid_in, int varid_in, int ncid_out)
// This will copy a variable that is an array of primitive type and its
// attributes from one file to another, assuming dimensions in the output
// file are already defined and have same dimension IDs and length.
void NetCDFFile::copyCoordinatesFrom(const std::string &srcFileName) {
  constexpr static std::string_view prefix{"NetCDFFile::copyCoordinatesFrom()"};

  if (srcFileName.empty()) {
    CUAS_ERROR_RANK0("{}: called with empty fileName. Exiting.", prefix)
    exit(1);
  }

  auto ncin = std::make_unique<NetCDFFile>(srcFileName, 'r');

  // copy lat(y,x), lon(y,x), lat_bnds(y, x, nv4), lon_bnds(y, x, nv4),  ;
  // The name of the third dimension name 'nv4' may vary.
  const std::array<std::string, 4> varNames{"lat", "lon", "lat_bnds", "lon_bnds"};
  for (auto &varName : varNames) {
    if (!ncin->hasVariable(varName)) {
      if (varName == "lat" || varName == "lon") {
        CUAS_ERROR_RANK0("{}: Variabe '{}' not found in file '{}'. Exiting.", prefix, varName, srcFileName)
        exit(1);
      } else {
        // lat_bnds and lon_bnds are only needed for later remapping using CDO
        CUAS_WARN_RANK0("{}: Variabe '{}' not found in file '{}'. Continue.", prefix, varName, srcFileName)
        continue;  // terminate current iteration for varName
      }
    }
    // check the dimensions and create new dimensions if missing
    auto srcDimIds = ncin->getDimIdsForVariable(varName);
    for (auto dimId : srcDimIds) {
      auto dimName = ncin->getDimName(dimId);
      if (!hasDimension(dimName)) {
        if (dimName == "x" || dimName == "y") {
          CUAS_ERROR_RANK0("{}: Dimension '{}' should already exist in output file. Exiting.", prefix, dimName)
          exit(1);
        }
        auto dimLen = ncin->getDimLength(dimName);
        int newDimId;
        SECURED_NETCDF_EXECUTION(nc_redef(fileId));
        SECURED_NETCDF_EXECUTION(nc_def_dim(fileId, dimName.c_str(), dimLen, &newDimId));
        SECURED_NETCDF_EXECUTION(nc_enddef(fileId));
      }
    }
    // now all dimensions exist and we can start copy all over
    auto srcVarId = ncin->getVarId(varName);
    auto srcFileId = ncin->getFileId();
    int retval;
    if ((retval = nc_copy_var(srcFileId, srcVarId, fileId)) != NC_NOERR) {
      std::string netcdfError = nc_strerror(retval);
      CUAS_ERROR("{}: nc_copy_var() failed for '{}' with error {}. Exiting.", prefix, varName, netcdfError)
      exit(1);
    }
  }
}

void NetCDFFile::setCoordinatesAttribute() {
  const char attText[] = "lat lon";
  for (auto &ncVar : netcdfVars) {
    auto varName = ncVar.first;
    auto varId = ncVar.second.varId;
    // only variable in the map plane are qualified
    if (variableHasDimensionByName(varName, "x") && variableHasDimensionByName(varName, "y")) {
      SECURED_NETCDF_EXECUTION(nc_put_att_text(fileId, varId, "coordinates", strlen(attText), attText));
    }
  }
}

}  // namespace CUAS
