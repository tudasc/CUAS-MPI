/**
 * File: ModelReader.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

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

  if (restartNoneZeroInitialGuess) {
    CUAS_WARN("restartNoneZeroInitialGuess not implemented yet")
  }
}

std::unique_ptr<CUAS::TimeDependentForcing> ModelReader::getTimeDependentForcing(
    std::string const &ncFileName, std::string const &variableName, std::vector<PetscScalar> const &xAxis,
    std::vector<PetscScalar> const &yAxis, PetscScalar multiplier, PetscScalar offset, bool loopForcing) {
  auto ncFile = std::make_unique<NetCDFFile>(ncFileName, 'r');
  auto nx = xAxis.size();
  auto ny = yAxis.size();
  auto mx = ncFile->getDimX();
  auto my = ncFile->getDimY();

  if (nx != mx || ny != my) {
    CUAS_ERROR("ModelReader.cpp: getTimeDependentForcing(): nx: {} != {} or ny: {} != {} in file '{}'. Exiting.", nx,
               mx, ny, my, ncFileName)
    exit(1);
  }

  // TODO: check if x- and y- axis match

  auto nt = ncFile->getDimLength("time");
  std::vector<timeSecs> time(nt);
  std::vector<std::unique_ptr<PETScGrid>> forcing(nt);
  for (std::unique_ptr<PETScGrid> &grid : forcing) {
    grid = std::make_unique<PETScGrid>(nx, ny);
  }

  ncFile->read("time", time);
  ncFile->read(variableName, forcing);

  return std::make_unique<CUAS::TimeDependentForcing>(forcing, time, multiplier, offset, loopForcing);
}

std::unique_ptr<CUAS::BufferedForcing> ModelReader::getBufferedForcing(
    std::string const &ncFileName, std::string const &variableName, std::vector<PetscScalar> const &xAxis,
    std::vector<PetscScalar> const &yAxis, int numberOfSlicesPerLoad, PetscScalar multiplier, PetscScalar offset,
    bool loopForcing) {
  auto ncFile = std::make_unique<NetCDFFile>(ncFileName, 'r');

  auto nx = xAxis.size();
  auto ny = yAxis.size();
  auto mx = ncFile->getDimX();
  auto my = ncFile->getDimY();

  if (nx != mx || ny != my) {
    CUAS_ERROR("{}::{} {}: nx: {} != {} or ny: {} != {} in file '{}'. Exiting.", __FILE__, __LINE__,
               __PRETTY_FUNCTION__, nx, mx, ny, my, ncFileName)
    exit(1);
  }

  return std::make_unique<CUAS::BufferedForcing>(ncFile, variableName, nx, ny, numberOfSlicesPerLoad, multiplier,
                                                 offset, loopForcing);
}

std::unique_ptr<CUAS::ScalarTimeDependentForcing> ModelReader::getScalarTimeDependentForcing(
    std::string const &ncFileName, std::string const &variableName, std::vector<PetscScalar> const &xAxis,
    std::vector<PetscScalar> const &yAxis, PetscScalar multiplier, PetscScalar offset, bool loopForcing) {
  // expected dimension names and dimension order. Note, len(x) = 1 and len(y) = 1
  const std::vector<std::string> expectedDimNames{"time", "y", "x"};

  auto ncFile = std::make_unique<NetCDFFile>(ncFileName, 'r');
  if (!ncFile->hasVariable(variableName)) {
    CUAS_ERROR("ModelReader.cpp: getScalarTimeDependentForcing(): Variable '{}' not found in file '{}'. Exiting.",
               variableName, ncFileName)
    exit(1);
  }

  // size of the model, not size from the forcing file
  auto nx = xAxis.size();
  auto ny = yAxis.size();

  // exit with error if dimensions have the wrong order or if something is missing
  auto dimIds = ncFile->getDimIdsForVariable(variableName);
  for (int i = 0; i < expectedDimNames.size(); ++i) {
    auto dimName = ncFile->getDimName(dimIds[i]);
    if (dimName != expectedDimNames[i]) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): Unexpected dim name ('{}' != '{}'). Exiting.", dimName,
                 expectedDimNames[i])
      exit(1);
    }
  }

  if (!ncFile->variableHasDimensionByName(variableName, "time")) {
    // has been checked already in tools/main.cpp
    CUAS_ERROR(
        "ModelReader.cpp: getScalarTimeDependentForcing(): Time dimension not found for variable '{}' in file '{}'. "
        "Exiting.",
        variableName, ncFileName)
    exit(1);
  }

  auto mx = ncFile->getDimLength("x");
  auto my = ncFile->getDimLength("y");
  if (mx != 1 || my != 1) {
    CUAS_ERROR(
        "ModelReader.cpp: getScalarTimeDependentForcing(): Incorrect number of points in x- or y-dimension (!= 1) for "
        "variable "
        "'{}' in file "
        "'{}'. Exiting.",
        variableName, ncFileName)
    exit(1);
  }

  auto nt = ncFile->getDimLength("time");
  std::vector<timeSecs> time(nt);
  std::vector<PetscScalar> timeSeries(nt);
  ncFile->read("time", time);
  ncFile->read(variableName, timeSeries);

  std::vector<PetscScalar> xPos({NC_FILL_DOUBLE});
  std::vector<PetscScalar> yPos({NC_FILL_DOUBLE});
  ncFile->read("x", xPos);
  ncFile->read("y", yPos);

  // find index in grid, if any
  int rowIndex, colIndex;
  auto allowedError = 1.0e-3 * (xAxis[1] - xAxis[0]);  // 0.1% of dx, todo: update after MR181

  // search for rowIndex
  {
    auto p = xPos.front();
    auto upperBound = std::upper_bound(xAxis.begin(), xAxis.end(), p) - xAxis.begin();
    auto distUpper = xAxis[upperBound] - p;
    auto distLower = p - xAxis[upperBound - 1];
    rowIndex = (distUpper > distLower) ? upperBound - 1 : upperBound;
    if (std::fabs(p - xAxis[rowIndex]) > allowedError) {
      CUAS_ERROR("ModelReader.cpp: getScalarTimeDependentForcing(): |xAxis[index] - x| > {}. Exiting.", allowedError)
      exit(1);
    }
    CUAS_INFO_RANK0("ModelReader.cpp: getScalarTimeDependentForcing(): Found xAxis[{}] = {} for x = {}", rowIndex,
                    xAxis[rowIndex], p)
  }

  // search for colIndex
  {
    auto p = yPos.front();
    auto upperBound = std::upper_bound(yAxis.begin(), yAxis.end(), p) - yAxis.begin();
    auto distUpper = yAxis[upperBound] - p;
    auto distLower = p - yAxis[upperBound - 1];
    colIndex = (distUpper > distLower) ? upperBound - 1 : upperBound;
    if (std::fabs(p - yAxis[rowIndex]) > allowedError) {
      CUAS_ERROR("ModelReader.cpp: getScalarTimeDependentForcing(): |yAxis[index] - y| > {}. Exiting.", allowedError)
      exit(1);
    }
    CUAS_INFO_RANK0("ModelReader.cpp: getScalarTimeDependentForcing(): Found yAxis[{}] = {} for y = {}", colIndex,
                    yAxis[colIndex], p)
  }

  return std::make_unique<CUAS::ScalarTimeDependentForcing>(nx, ny, rowIndex, colIndex, timeSeries, time, multiplier,
                                                            offset, loopForcing);
}

std::unique_ptr<CUAS::SteadyForcing> ModelReader::getSteadyForcing(std::string const &ncFileName,
                                                                   std::string const &variableName,
                                                                   std::vector<PetscScalar> const &xAxis,
                                                                   std::vector<PetscScalar> const &yAxis,
                                                                   PetscScalar multiplier, PetscScalar offset) {
  auto ncFile = std::make_unique<NetCDFFile>(ncFileName, 'r');
  // size of the model, not size from the forcing file
  auto nx = (int)xAxis.size();
  auto ny = (int)yAxis.size();
  auto mx = ncFile->getDimX();
  auto my = ncFile->getDimY();
  if (nx != mx || ny != my) {
    CUAS_ERROR("ModelReader.cpp: getSteadyForcing(): nx: {} != {} or ny: {} != {} in file '{}'. Exiting.", nx, mx, ny,
               my, ncFileName)
    exit(1);
  }

  // TODO: check if x- and y- axis match

  PETScGrid forcing(nx, ny);
  ncFile->read(variableName, forcing);
  return std::make_unique<CUAS::SteadyForcing>(forcing, multiplier, offset);
}

bool ModelReader::isTimeDependent(std::string const &ncFileName, std::string const &variableName) {
  auto ncFile = std::make_unique<NetCDFFile>(ncFileName, 'r');
  if (!ncFile->hasVariable(variableName)) {
    CUAS_ERROR("ModelReader.cpp: isTimeDependent(): Variable '{}' not found in file '{}'. Exiting.", variableName,
               ncFileName)
    exit(1);
  }
  return ncFile->variableHasDimensionByName(variableName, "time");
}

bool ModelReader::isTimeDependentField(std::string const &ncFileName, std::string const &variableName) {
  // expected dimension names and dimension order
  const std::vector<std::string> expectedDimNames{"time", "y", "x"};

  auto ncFile = std::make_unique<NetCDFFile>(ncFileName, 'r');
  if (!ncFile->hasVariable(variableName)) {
    CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): Variable '{}' not found in file '{}'. Exiting.", variableName,
               ncFileName)
    exit(1);
  }

  // exit early if wrong number of dimensions
  auto dims = ncFile->getDimensionsForVariable(variableName);
  if (dims != expectedDimNames.size()) {
    return false;
  }

  // exit early if the wrong dimensions
  auto hasTime = ncFile->variableHasDimensionByName(variableName, "time");
  auto hasX = ncFile->variableHasDimensionByName(variableName, "y");
  auto hasY = ncFile->variableHasDimensionByName(variableName, "x");
  if (!(hasTime && hasX && hasY)) {
    return false;
  }

  // exit early, if we have degenerated (len==1) spatial dimensions, e.g. for a point
  if (2 > ncFile->getDimLength("x") || 2 > ncFile->getDimLength("y")) {
    return false;
  }

  // exit with error if dimensions have the wrong order
  auto dimIds = ncFile->getDimIdsForVariable(variableName);
  for (int i = 0; i < expectedDimNames.size(); ++i) {
    auto dimName = ncFile->getDimName(dimIds[i]);
    if (dimName != expectedDimNames[i]) {
      CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): Unexpected dim name ('{}' != '{}'). Exiting.", dimName,
                 expectedDimNames[i])
      exit(1);
    }
  }

  // check if the dimension time and the size of the variable "time" are the same
  auto numberOfDimensionsTime = ncFile->getDimensionsForVariable("time");
  if (numberOfDimensionsTime != 1) {
    CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): the variable time is multidimensional. Exiting.")
    exit(1);
  }

  // compare the dimension time and the variable "time". The time dimension is the 0th dimId of "time"
  auto lengthOfDimTime = ncFile->getDimLength("time");
  if (lengthOfDimTime < 1) {
    CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): Invalid dimLen = {} for 'time'. Exiting.", lengthOfDimTime)
    exit(1);
  }

  auto dimIdsForVariableTime = ncFile->getDimIdsForVariable("time");
  auto lengthOfVarTime = ncFile->getDimLength(dimIdsForVariableTime[0]);
  if (lengthOfDimTime != lengthOfVarTime) {
    CUAS_ERROR("ModelReader.cpp: isTimeDependentField(): The size of time in " + variableName +
               " is not equal to the size of the variable time in the netcdf file. Exiting.")
    exit(1);
  }

  // finally
  return lengthOfVarTime > 0;
}

}  // namespace CUAS
