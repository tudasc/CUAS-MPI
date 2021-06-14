#include "CUASFile.h"

namespace CUAS {

int CUASFile::getVarId(std::string const &varName) const {
  int varId;
  // get varId
  if (int retval = nc_inq_varid(fileId, varName.c_str(), &varId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: getVarId(): A netcdf error occured: " + netcdfError + "Exiting.");
    exit(1);
  }
  return varId;
};

bool CUASFile::checkDimensions(std::string const &varName, PETScGrid const &input) {
  int varId = getVarId(varName);
  size_t numOfCols;
  size_t numOfRows;
  std::vector<int> dimIds(2);
  nc_inq_vardimid(fileId, varId, dimIds.data());
  nc_inq_dimlen(fileId, dimIds[0], &numOfRows);
  nc_inq_dimlen(fileId, dimIds[1], &numOfCols);
  if (numOfCols == input.getTotalNumOfCols() && numOfRows == input.getTotalNumOfRows()) {
    return true;
  } else {
    return false;
  }
}

bool CUASFile::checkDimensions(std::string const &varName, PETScVec const &input) {
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

bool CUASFile::checkDimensions(std::string const &varName, std::vector<PetscScalar> const &input) {
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

void CUASFile::defineGrid(std::string const &varName) {
  int varId;
  if (int retval = nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 2, dimIdsGrid.data(), &varId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: defineGrid(): A netcdf error occured: " + netcdfError + "Exiting.");
    exit(1);
  }
}

void CUASFile::defineScalar(std::string const &varName) {
  int varId;
  if (int retval = nc_def_var(fileId, varName.c_str(), NC_DOUBLE, 0, nullptr, &varId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: defineScalar(): A netcdf error occured: " + netcdfError + "Exiting.");
    exit(1);
  }
}

CUASFile::CUASFile(std::string const &fileName, char const mode) : fileName(fileName) {
  switch (mode) {
    // openfile for reading
    case 'r': {
      if (int retval = nc_open_par(fileName.c_str(), NC_NOWRITE, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("CUASFile.cpp: CUASFile() in read-mode: A netcdf error occured: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
      break;
    }
    // openfile for writing
    case 'w': {
      if (int retval = nc_open_par(fileName.c_str(), NC_WRITE, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
        std::string netcdfError = nc_strerror(retval);
        Logger::instance().error("CUASFile.cpp: CUASFile() in write-mode: A netcdf error occured: " + netcdfError +
                                 "Exiting.");
        exit(1);
      }
      break;
    }
    default: {
      Logger::instance().error(
          "CUASFile.cpp: CUASFile() in write-mode: Please enter 'r' or 'w' as a mode to read or write the netcdf file. "
          "Exiting.");
      exit(1);
    }
  }
  nc_enddef(fileId);
}

CUASFile::CUASFile(std::string const &fileName, int dimX, int dimY, int mpiRank) : fileName(fileName) {
  if (int retval = nc_create_par(fileName.c_str(), NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &fileId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: CUASFile() in create-mode: A netcdf error occured: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
  dimIdsGrid.reserve(2);
  nc_def_dim(fileId, "x", dimX, &dimIdsGrid[1]);
  nc_def_dim(fileId, "y", dimY, &dimIdsGrid[0]);
  nc_enddef(fileId);
}

void CUASFile::read(std::string const &varName, PETScGrid &input) {
  int varId = getVarId(varName);

  if (!checkDimensions(varName, input)) {
    Logger::instance().error("CUASFile.cpp: read() with PETScGrid: sizes do not fit! Exiting.");
    exit(1);
  }

  // determine start and count of values
  size_t start[2];
  size_t count[2];

  auto grid = input.getWriteHandle();
  auto gridFromNetcdf = grid.getRaw();

  start[1] = input.getCornerX();
  count[0] = 1;
  count[1] = input.getLocalNumOfCols();

  for (int i = 0; i < input.getLocalNumOfRows(); ++i) {
    start[0] = input.getCornerY() + i;
    if (int retval = nc_get_vara_double(fileId, varId, start, count, &gridFromNetcdf[i][0])) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("CUASFile.cpp: write() with PETScGrid: A netcdf error occured: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
  }
};

void CUASFile::read(std::string const &varName, PETScVec &input) {
  int varId = getVarId(varName);

  if (!checkDimensions(varName, input)) {
    Logger::instance().error("CUASFile.cpp: read with PETScVec: size does not fit! Exiting.");
    exit(1);
  }

  int startIndex;
  int endIndex;
  input.getOwnershipRange(startIndex, endIndex);

  size_t count = endIndex - startIndex;
  size_t start = startIndex;

  int localSizeOfInput = count;
  double vecFromNetcdf[localSizeOfInput];

  if (int retval = nc_get_vara_double(fileId, varId, &start, &count, &vecFromNetcdf[0])) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp:read() with PETScVec: A netcdf error occured: " + netcdfError + "Exiting.");
    exit(1);
  }

  for (int i = startIndex; i < endIndex; ++i) {
    input.setValue(i, vecFromNetcdf[i - startIndex]);
  }

  input.assemble();
};

void CUASFile::read(std::string const &varName, std::vector<PetscScalar> &input) {
  int varId = getVarId(varName);

  if (!checkDimensions(varName, input)) {
    Logger::instance().error("CUASFile.cpp: read with std::vector: size does not fit! Exiting.");
    exit(1);
  }

  int globalSizeOfInput = input.size();

  double vecFromNetcdf[globalSizeOfInput];

  // the std::vector is not distributed. Therefore the sequential get method from netcdf is called
  if (int retval = nc_get_var_double(fileId, varId, &vecFromNetcdf[0])) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: read() with std::vector: A netcdf error occured: " + netcdfError +
                             "Exiting.");
    exit(1);
  }

  for (int i = 0; i < globalSizeOfInput; ++i) {
    input[i] = vecFromNetcdf[i];
  }
}

void CUASFile::read(std::string const &varName, PetscScalar &input) {
  int varId = getVarId(varName);

  if (int retval = nc_get_var_double(fileId, varId, &input)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: read() with PetscScalar: A netcdf error occured: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
}

void CUASFile::write(std::string const &varName, PETScGrid const &input) {
  int varId = getVarId(varName);

  if (!checkDimensions(varName, input)) {
    Logger::instance().error("CUASFile.cpp: write with PETScGrid: sizes do not fit! Exiting.");
    exit(1);
  }
  // determine start and count of values
  size_t start[2];
  size_t count[2];

  auto grid = input.getReadHandle();
  auto gridFromNetcdf = grid.getRaw();

  start[1] = input.getCornerX();
  count[0] = 1;
  count[1] = input.getLocalNumOfCols();

  for (int i = 0; i < input.getLocalNumOfRows(); ++i) {
    start[0] = input.getCornerY() + i;
    if (int retval = nc_put_vara_double(fileId, varId, start, count, &gridFromNetcdf[i][0])) {
      std::string netcdfError = nc_strerror(retval);
      Logger::instance().error("CUASFile.cpp: write() with PETScGrid: A netcdf error occured: " + netcdfError +
                               "Exiting.");
      exit(1);
    }
  }
}

void CUASFile::write(std::string const &varName, PetscScalar const &input) {
  int varId = getVarId(varName);

  if (int retval = nc_put_var_double(fileId, varId, &input)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: write() with PetscScalar: A netcdf error occured: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
}

void CUASFile::writeTimeSteps(std::vector<int> timeSteps) {
  nc_redef(fileId);
  std::vector<int> dimTimeSteps(1);
  if (int retval = nc_def_dim(fileId, "numOfTimeSteps", timeSteps.size(), &dimTimeSteps[0])) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error(
        "CUASFile.cpp: writeTimeSteps() when defining dimensions: A netcdf error occured: " + netcdfError + "Exiting.");
    exit(1);
  }
  int varId;
  if (int retval = nc_def_var(fileId, "timeSteps", NC_INT, 1, dimTimeSteps.data(), &varId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error(
        "CUASFile.cpp: writeTimeSteps() when defining a variable: A netcdf error occured: " + netcdfError + "Exiting.");
    exit(1);
  }
  nc_enddef(fileId);
  if (int retval = nc_put_var_int(fileId, varId, timeSteps.data())) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error(
        "CUASFile.cpp: writeTimeSteps() when defining writing timeSteps: A netcdf error occured: " + netcdfError +
        "Exiting.");
    exit(1);
  }
}

CUASFile::~CUASFile() {
  if (int retval = nc_close(fileId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: ~CUASFile(): A netcdf error occured: " + netcdfError + "Exiting.");
    exit(1);
  }
};
}  // namespace CUAS
