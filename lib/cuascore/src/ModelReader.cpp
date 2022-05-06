#include "ModelReader.h"

#include "timeparse.h"

namespace CUAS {

ModelReader::ModelReader(std::string const &fileName) { file = std::make_unique<NetCDFFile>(fileName, 'r'); }

std::unique_ptr<CUAS::CUASModel> ModelReader::fillModelFromNetcdf() {
  auto pmodel = std::make_unique<CUAS::CUASModel>(file->getDimX(), file->getDimY());
  auto &model = *pmodel;

  file->read("x", model.xAxis);
  file->read("y", model.yAxis);

  file->read("topg", *model.topg);
  file->read("thk", *model.thk);
  file->read("bnd_mask", *model.bndMask);

  model.init();

  return pmodel;
}

void ModelReader::restartFromFile(CUAS::CUASSolver &solver, std::string const &restartFile,
                                  bool restartNoneZeroInitialGuess) {
  auto restartNetcdfFile = std::make_unique<NetCDFFile>(restartFile, 'r');
  restartNetcdfFile->read("head", *solver.currHead);
  restartNetcdfFile->read("transmissivity", *solver.currTransmissivity);

  // TODO: remove after merge request !112 feat/updateCUASKernels
  solver.nextHead->copy(*solver.currHead);                      // obsolete
  solver.nextTransmissivity->copy(*solver.nextTransmissivity);  // obsolete

  if (restartNoneZeroInitialGuess) {
    CUAS_WARN("restartNoneZeroInitialGuess not implemented yet")
  }
}

std::unique_ptr<CUAS::TimeForcing> ModelReader::getTimeForcing(std::string const &timeForcingFileName,
                                                               std::string const &fieldName,
                                                               PetscScalar const multiplier, PetscScalar const offset,
                                                               bool const loopForcing) {
  auto timeForcingFile = std::make_unique<NetCDFFile>(timeForcingFileName, 'r');
  auto timeDim = timeForcingFile->getDimLength("time");
  auto x = timeForcingFile->getDimX();
  auto y = timeForcingFile->getDimY();
  std::vector<timeSecs> time(timeDim);
  std::vector<std::unique_ptr<PETScGrid>> forcing(timeDim);
  for (std::unique_ptr<PETScGrid> &grid : forcing) {
    grid = std::make_unique<PETScGrid>(x, y);
  }
  timeForcingFile->read("time", time);
  timeForcingFile->read(fieldName, forcing);
  return std::make_unique<CUAS::TimeForcing>(forcing, time, multiplier, offset, loopForcing);
}

std::unique_ptr<CUAS::ConstantForcing> ModelReader::getConstantForcing(std::string const &timeForcingFileName,
                                                                       std::string const &fieldName,
                                                                       PetscScalar const multiplier,
                                                                       PetscScalar const offset) {
  auto timeForcingFile = std::make_unique<NetCDFFile>(timeForcingFileName, 'r');
  auto x = timeForcingFile->getDimX();
  auto y = timeForcingFile->getDimY();
  PETScGrid forcing(x, y);
  timeForcingFile->read(fieldName, forcing);
  return std::make_unique<CUAS::ConstantForcing>(forcing, multiplier, offset);
}

bool ModelReader::isTimeDependentField(std::string const &timeForcingFileName, std::string const &fieldName) {
  auto timeForcingFile = std::make_unique<NetCDFFile>(timeForcingFileName, 'r');

  if (!timeForcingFile->hasVariable(fieldName)) {
    CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): The timeForcingFile does not have the variable " + fieldName +
               " that was passed to the function. Exiting.")
    exit(1);
  }

  auto dims = timeForcingFile->getDimensionsForVariable(fieldName);
  auto hasTime = timeForcingFile->variableHasDimensionByName(fieldName, "time");

  if (dims == 2) {
    if (hasTime) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): The variable " + fieldName +
                 " is 2 dimensional and has time dimension. Exiting.")
      exit(1);
    }

    return false;
  } else if (dims == 3) {
    if (!hasTime) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): The variable " + fieldName +
                 " is 3 dimensional but has not a time dimension. Exiting.")
      exit(1);
    }
    if (!timeForcingFile->hasVariable("time")) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): " + fieldName +
                 " is 3 dimensional but the timeForcingFile does not have the required variable 'time'. "
                 "Exiting.")
      exit(1);
    }
    if (!timeForcingFile->hasDimension("time")) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): " + fieldName +
                 " is 3 dimensional but The timeForcingFile does not have the required dimension 'time'. "
                 "Exiting.")
      exit(1);
    }

    // check if the dimension time and the size of the variable "time" are the same
    auto lengthOfDimTime = timeForcingFile->getDimLength("time");
    auto numberOfDimensionsTime = timeForcingFile->getDimensionsForVariable("time");
    if (numberOfDimensionsTime != 1) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): the variable time is multidimensional. Exiting.")
      exit(1);
    }

    // before there was fieledname.c_str() but it should be the variable "time" right?
    auto dimIdsForVariableTime = timeForcingFile->getDimIdsForVariable("time");
    // compare the dimension time and the variable "time". The time dimension is the 0th dimId of "time"
    auto lengthOfVarTime = timeForcingFile->getDimLength(dimIdsForVariableTime[0]);
    if (lengthOfDimTime != lengthOfVarTime) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): The size of time in " + fieldName +
                 " is not equal to the size of the variable time in the netcdf file. Exiting.")
      exit(1);
    }
    return true;
  } else {
    CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): The variable " + fieldName +
               " is not 2 or 3 dimensional. Exiting.")
    exit(1);
  }
}

}  // namespace CUAS
