#include "NetCDFFile.h"

#include <array>

#include <iostream>

// GRID_ROWWISE causes a deadlock in nc_put_vara_double
// probably caused by NC_COLLECTIVE
// we currently use GRID_BLOCKWISE presupposing PETScGrid data stored en bloc
#define GRID_ROWWISE 0
#define GRID_BLOCKWISE 1
#define GRID_IMPLEMENTATION GRID_BLOCKWISE

namespace CUAS {

int NetCDFFile::getVarId(std::string const &varName) const {
  int varId;
  // get varId
  if (int retval = nc_inq_varid(fileId, varName.c_str(), &varId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: getVarId(): A netcdf error occurred: " + netcdfError + "Exiting.");
    exit(1);
  }
  return varId;
};

bool NetCDFFile::checkDimensions(std::string const &varName, PETScGrid const &input) const {
  int varId = getVarId(varName);
  size_t numOfCols;
  size_t numOfRows;
  std::array<int, 2> dimIds;
  nc_inq_vardimid(fileId, varId, dimIds.data());
  nc_inq_dimlen(fileId, dimIds[0], &numOfRows);
  nc_inq_dimlen(fileId, dimIds[1], &numOfCols);
  if (numOfCols == input.getTotalNumOfCols() && numOfRows == input.getTotalNumOfRows()) {
    return true;
  } else {
    return false;
  }
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
  if (numOfCols == input.getTotalNumOfCols() && numOfRows == input.getTotalNumOfRows()) {
    return true;
  } else {
    return false;
  }
}

bool NetCDFFile::checkDimensions(std::string const &varName, PETScVector const &input) const {
  int varId = getVarId(varName);
  size_t size;
  int dimId;
  nc_inq_vardimid(fileId, varId, &dimId);
  nc_inq_dimlen(fileId, dimId, &size);
  if (size == input.getSize()) {
    return true;
  } else {
    return false;
  }
}

bool NetCDFFile::checkDimensionsUnlimited(std::string const &varName, PETScVector const &input) const {
  int varId = getVarId(varName);
  size_t size;
  std::vector<int> dimIds(2);
  nc_inq_vardimid(fileId, varId, dimIds.data());
  nc_inq_dimlen(fileId, dimIds[1], &size);
  if (size == input.getSize()) {
    return true;
  } else {
    return false;
  }
}

bool NetCDFFile::checkDimensions(std::string const &varName, std::vector<PetscScalar> const &input) const {
  int varId = getVarId(varName);
  size_t size;
  int dimId;
  nc_inq_vardimid(fileId, varId, &dimId);
  nc_inq_dimlen(fileId, dimId, &size);
  if (size == input.size()) {
    return true;
  } else {
    return false;
  }
}

void NetCDFFile::defineVectorWithOffset(std::string const &varName, int offset) {
  // maybe we need a different dimIdGrids array for unlimited dimensions
  int varId;
  if (int retval = nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 1, &dimIds[offset], &varId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.");
    exit(1);
  }

  netcdfVars[varName] = NetCDFVar{varName, varId, LIMITED};
}

int NetCDFFile::getCurrentHead() const {
  int unlimId;
  nc_inq_unlimdim(fileId, &unlimId);
  size_t unlimitedLength;
  nc_inq_dimlen(fileId, unlimId, &unlimitedLength);
  return unlimitedLength;
}

NetCDFFile::NetCDFFile(std::string const &fileName, char mode) : fileName(fileName) {
  switch (mode) {
    // openfile for reading
    case 'r': {
      if (int retval = nc_open_par(fileName.c_str(), NC_NOWRITE, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("NetCDFFile.cpp: NetCDFFile() in read-mode: A netcdf error occurred: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
      break;
    }
    // openfile for writing
    case 'w': {
      if (int retval = nc_open_par(fileName.c_str(), NC_WRITE, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("NetCDFFile.cpp: NetCDFFile() in write-mode: A netcdf error occurred: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
      break;
    }
    default: {
      Logger::instance().error(
          "NetCDFFile.cpp: NetCDFFile() in write-mode: Please enter 'r' or 'w' as a mode to read or write the netcdf "
          "file. "
          "Exiting.");
      exit(1);
    }
  }
  nc_enddef(fileId);

  // get dimId of unlimtedDimension
  int dimIdUnlimited;
  size_t sizeOfUnlimitedDim;
  if (int retval = nc_inq_dimid(fileId, "time", &dimIdUnlimited)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().info("NetCDFFile.cpp: The netcdf file doesn't have unlimited dimensions.");
  }
  // get lenght of dimension timeSteps
  else if (int retval = nc_inq_dimlen(fileId, dimIdUnlimited, &sizeOfUnlimitedDim)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }

  // get dimId of x
  int dimIdX;
  nc_inq_dimid(fileId, "x", &dimIdX);
  // get lenght of dimension x
  size_t modelDimX;
  nc_inq_dimlen(fileId, dimIdX, &modelDimX);
  // get dimId of y
  int dimIdY;
  nc_inq_dimid(fileId, "y", &dimIdY);
  // get lenght of dimension x
  size_t modelDimY;
  nc_inq_dimlen(fileId, dimIdY, &modelDimY);

  dimX = modelDimX;
  dimY = modelDimY;
}

NetCDFFile::NetCDFFile(std::string const &fileName, int dimX, int dimY) : fileName(fileName), dimX(dimX), dimY(dimY) {
  if (int retval = nc_create_par(fileName.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: NetCDFFile() in create-mode: A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }

  // define dimensions for grids
  dimIds.reserve(3);
  nc_def_dim(fileId, "time", NC_UNLIMITED, &dimIds[0]);
  nc_def_dim(fileId, "y", dimY, &dimIds[1]);
  nc_def_dim(fileId, "x", dimX, &dimIds[2]);

  nc_enddef(fileId);
}

void NetCDFFile::defineGrid(std::string const &varName, bool isUnlimited) {
  int varId;
  if (isUnlimited) {
    if (int retval = nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 3, dimIds.data(), &varId)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.");
      exit(1);
    }
  } else {
    // maybe we need a different dimIdGrids array for unlimited dimensions
    if (int retval = nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 2, &dimIds[1], &varId)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: defineGrid(): A netcdf error occurred: " + netcdfError + "Exiting.");
      exit(1);
    }
  }

  netcdfVars[varName] = NetCDFVar{varName, varId, isUnlimited};
}

void NetCDFFile::defineVectorX(std::string const &varName) { defineVectorWithOffset(varName, 2); }

void NetCDFFile::defineVectorY(std::string const &varName) { defineVectorWithOffset(varName, 1); }

void NetCDFFile::defineScalar(std::string const &varName, bool isUnlimited) {
  int varId;
  if (isUnlimited) {
    if (int retval = nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 1, &dimIds[0], &varId)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: defineScalar(): A netcdf error occurred: " + netcdfError + "Exiting.");
      exit(1);
    }
  } else {
    if (int retval = nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 0, nullptr, &varId)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: defineScalar(): A netcdf error occurred: " + netcdfError + "Exiting.");
      exit(1);
    }
  }

  netcdfVars[varName] = NetCDFVar{varName, varId, isUnlimited};
}

void NetCDFFile::addAttributeToVariable(std::string const &varName, std::string const &attributeName,
                                        std::string const &attributeText) {
  auto varId = getVarId(varName);
  if (int retval =
          nc_put_att_text(fileId, varId, attributeName.c_str(), strlen(attributeText.c_str()), attributeText.c_str())) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: addAttributeToVariable(): A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, std::string const &attributeText) {
  if (int retval = nc_put_att_text(fileId, NC_GLOBAL, attributeName.c_str(), strlen(attributeText.c_str()),
                                   attributeText.c_str())) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: addGlobalAttribute(): A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, PetscScalar attribute) {
  if (int retval = nc_put_att_double(fileId, NC_GLOBAL, attributeName.c_str(), NC_DOUBLE, 1, &attribute)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: addGlobalAttribute(): A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, bool attribute) {
  int boolToInt(attribute);
  if (int retval = nc_put_att_int(fileId, NC_GLOBAL, attributeName.c_str(), NC_INT, 1, &boolToInt)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: addGlobalAttribute(): A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
}

void NetCDFFile::addGlobalAttribute(std::string const &attributeName, int attribute) {
  if (int retval = nc_put_att_int(fileId, NC_GLOBAL, attributeName.c_str(), NC_INT, 1, &attribute)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: addGlobalAttribute(): A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
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
    Logger::instance().error("NetCDFFile.cpp: read() with PETScGrid: The grid is not limited or unlimited! Exiting.");
    exit(1);
  }
  // determine start and count of values
  if (isUnlimited) {
    if (!checkDimensionsUnlimited(varName, dest)) {
      Logger::instance().error("NetCDFFile.cpp: read() with PETScGrid: sizes do not fit! Exiting.");
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
      if (int retval = nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0])) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("NetCDFFile.cpp: read() with PETScGrid: A netcdf error occurred: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = headIndex;
    start[1] = dest.getCornerY();
    start[2] = dest.getCornerX();
    count[0] = 1;
    count[1] = dest.getLocalNumOfRows();
    count[2] = dest.getLocalNumOfCols();

    if (int retval = nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0])) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: read() with PETScGrid: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
#endif
  } else {
    if (!checkDimensions(varName, dest)) {
      Logger::instance().error("NetCDFFile.cpp: read() with PETScGrid: sizes do not fit! Exiting.");
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
      if (int retval = nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0])) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("NetCDFFile.cpp: read() with PETScGrid: A netcdf error occurred: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = dest.getCornerY();
    start[1] = dest.getCornerX();
    count[0] = dest.getLocalNumOfRows();
    count[1] = dest.getLocalNumOfCols();

    if (int retval = nc_get_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0])) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: read() with PETScGrid: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
#endif
  }
};

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
    Logger::instance().error("NetCDFFile.cpp: read() with PETScVector: The grid is not limited or unlimited! Exiting.");
    exit(1);
  }

  if (isUnlimited) {
    if (!checkDimensionsUnlimited(varName, dest)) {
      Logger::instance().error("NetCDFFile.cpp: read with PETScVector: size does not fit! Exiting.");
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
    if (int retval = nc_get_vara_double(fileId, varId, start.data(), count.data(), vecFromNetcdf.data())) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp:read() with PETScVector: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }

    for (int i = startIndex; i < endIndex; ++i) {
      dest.setValue(i, vecFromNetcdf[i - startIndex]);
    }

    dest.assemble();
  } else {
    if (!checkDimensions(varName, dest)) {
      Logger::instance().error("NetCDFFile.cpp: read with PETScVector: size does not fit! Exiting.");
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
    if (int retval = nc_get_vara_double(fileId, varId, &start, &count, vecFromNetcdf.data())) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp:read() with PETScVector: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }

    for (int i = startIndex; i < endIndex; ++i) {
      dest.setValue(i, vecFromNetcdf[i - startIndex]);
    }

    dest.assemble();
  }
};

void NetCDFFile::read(std::string const &varName, std::vector<PetscScalar> &dest) const {
  int varId = getVarId(varName);

  if (!checkDimensions(varName, dest)) {
    Logger::instance().error("NetCDFFile.cpp: read with std::vector: size does not fit! Exiting.");
    exit(1);
  }

  // the std::vector is not distributed. Therefore the sequential get method from netcdf is called
  if (int retval = nc_get_var_double(fileId, varId, dest.data())) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: read() with std::vector: A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
}

void NetCDFFile::read(std::string const &varName, std::vector<long> &dest) const {
  int varId = getVarId(varName);
  int numOfDims;
  int dimId;
  size_t dimLen;

  nc_inq_varndims(fileId, varId, &numOfDims);
  if (numOfDims != 1) {
    // todo:   "variable '%s' in '%s' should to have 1 dimension (got %d)"
    Logger::instance().error("NetCDFFile.cpp: read() with std::vector<long>: Wrong number of dimensions! Exiting.");
    exit(1);
  }

  nc_inq_vardimid(fileId, varId, &dimId);
  nc_inq_dimlen(fileId, dimId, &dimLen);

  dest.resize(dimLen);  // memory allocation happens here

  if (int retval = nc_get_var_long(fileId, varId, dest.data())) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: read() with std::vector<long>: A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
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
    Logger::instance().error("NetCDFFile.cpp: read() with PETScScalar: The grid is not limited or unlimited! Exiting.");
    exit(1);
  }

  if (isUnlimited) {
    int headIndex = getCurrentHead() - 1;

    size_t start[1];
    size_t count[1];

    start[0] = headIndex;
    count[0] = 1;

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    if (int retval = nc_get_vara_double(fileId, varId, start, count, &dest)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: read() with PetscScalar: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
  }

  else {
    if (int retval = nc_get_var_double(fileId, varId, &dest)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: read() with PetscScalar: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
  }
}

void NetCDFFile::write(std::string const &varName, PETScGrid const &input, int currentTimeStep) {
  auto currentVar = netcdfVars[varName];
  int varId = currentVar.varId;
  bool isUnlimited = currentVar.isUnlimited;

  if (isUnlimited) {
    if (!checkDimensionsUnlimited(varName, input)) {
      Logger::instance().error("NetCDFFile.cpp: write with PETScGrid: sizes do not fit! Exiting.");
      exit(1);
    }
    // determine start and count of values
    std::array<size_t, 3> start;
    std::array<size_t, 3> count;

    auto grid = input.getReadHandle();
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
      if (int retval = nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0])) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = currentTimeStep;
    start[1] = input.getCornerY();
    start[2] = input.getCornerX();
    count[0] = 1;
    count[1] = input.getLocalNumOfRows();
    count[2] = input.getLocalNumOfCols();

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    if (int retval = nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0])) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
#endif
  } else {
    if (!checkDimensions(varName, input)) {
      Logger::instance().error("NetCDFFile.cpp: write() with PETScGrid: sizes do not fit! Exiting.");
      exit(1);
    }

    std::array<size_t, 2> start;
    std::array<size_t, 2> count;

    auto grid = input.getReadHandle();
    auto gridFromNetcdf = grid.getRaw();

#if GRID_IMPLEMENTATION == GRID_ROWWISE
    start[1] = input.getCornerX();
    count[0] = 1;
    count[1] = input.getLocalNumOfCols();

    for (int i = 0; i < input.getLocalNumOfRows(); ++i) {
      start[0] = input.getCornerY() + i;
      if (int retval = nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[i][0])) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
    }
#elif GRID_IMPLEMENTATION == GRID_BLOCKWISE
    start[0] = input.getCornerY();
    start[1] = input.getCornerX();
    count[0] = input.getLocalNumOfRows();
    count[1] = input.getLocalNumOfCols();

    if (int retval = nc_put_vara_double(fileId, varId, start.data(), count.data(), &gridFromNetcdf[0][0])) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: write() with PETScGrid: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
#endif
  }
}

void NetCDFFile::write(std::string const &varName, PETScVector &input, const int currentTimeStep) {
  auto currentVar = netcdfVars[varName];
  int varId = currentVar.varId;
  bool isUnlimited = currentVar.isUnlimited;
  if (isUnlimited) {
    Logger::instance().error("writing unlimited PETScVectors is not supported");
    exit(1);
    /*
    if (!checkDimensionsUnlimited(varName, input)) {
      Logger::instance().error("NetCDFFile.cpp: write() with PETScVector: size does not fit! Exiting.");
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
    if (int retval = nc_put_vara_double(fileId, varId, start.data(), count.data(), vecFromNetcdf.data())) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: write() with PETScVector: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }*/
  } else {
    if (!checkDimensions(varName, input)) {
      Logger::instance().error("NetCDFFile.cpp: write() with PETScVector: size does not fit! Exiting.");
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
    if (int retval = nc_put_vara_double(fileId, varId, &start, &count, vecFromNetcdf.data())) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: write() with PETScVector: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
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
    Logger::instance().error("writing unlimited PETScVectors is not supported");
    exit(1);
  } else {
    if (!checkDimensions(varName, input)) {
      Logger::instance().error("NetCDFFile.cpp: write() with std::vector: size does not fit! Exiting.");
      exit(1);
    }

    nc_var_par_access(fileId, varId, NC_COLLECTIVE);
    if (int retval = nc_put_var_double(fileId, varId, input.data())) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: write() with std::vector: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
  }
}

// change this to scalars that are written as an unlimited var
/*void NetCDFFile::writeTimeSteps(std::vector<int> timeSteps) {
  nc_redef(fileId);
  std::vector<int> dimTimeSteps(1);
  if (int retval = nc_def_dim(fileId, "numOfTimeSteps", timeSteps.size(), &dimTimeSteps[0])) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: writeTimeSteps() when defining dimensions: A netcdf error occurred: " +
                             netcdfError + "Exiting.");
    exit(1);
  }
  int varId;
  if (int retval = nc_def_var(fileId, "time", NC_INT, 1, dimTimeSteps.data(), &varId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: writeTimeSteps() when defining a variable: A netcdf error occurred: " +
                             netcdfError + "Exiting.");
    exit(1);
  }
  nc_enddef(fileId);
  if (int retval = nc_put_var_int(fileId, varId, timeSteps.data())) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error(
        "NetCDFFile.cpp: writeTimeSteps() when defining writing timeSteps: A netcdf error occurred: " + netcdfError +
        "Exiting.");
    exit(1);
  }
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
    if (auto retval = nc_put_vara_double(fileId, varId, start.data(), count.data(), &input)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: write() with PetscScalar: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
  } else {
    if (auto retval = nc_put_var_double(fileId, varId, &input)) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("NetCDFFile.cpp: read() with PetscScalar: A netcdf error occurred: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
  }
}

NetCDFFile::~NetCDFFile() {
  if (int retval = nc_close(fileId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: NetCDFFilele(): A netcdf error occurred: " + netcdfError + "Exiting.");
    exit(1);
  }
}

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
      Logger::instance().error("NetCDFFile.cpp: readTextAttribute(): Invalid attribute type. Exiting.");
      exit(1);
    }
  } else if (retval == NC_ENOTATT) {
    //
    // attribute not found
    //
    Logger::instance().warn("NetCDFFile.cpp: readTextAttribute(): attribute '" + attName + "' not found", attName);
    result = "";
  } else {
    //
    // something went wrong
    //
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("NetCDFFile.cpp: readTextAttribute(): A netcdf error occurred: " + netcdfError +
                             "Exiting.");
    exit(1);
  }

  return result;
}

}  // namespace CUAS
