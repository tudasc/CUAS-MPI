/**
 * File: ModelReader.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_MODEL_READER_H
#define CUAS_MODEL_READER_H

#include "CUASModel.h"
#include "CUASSolver.h"
#include "Forcing/BufferedForcing.h"
#include "Forcing/ScalarTimeDependentForcing.h"
#include "Forcing/SteadyForcing.h"
#include "Forcing/TimeDependentForcing.h"
#include "NetCDFFile.h"

#include <memory>
#include <string>
#include <vector>

// example usage:
// create a model reader
//  CUAS::ModelReader reader("NEGIS1200m_for_CUAS_SICOc_IDBMG4v3N_nLakes.nc");
// create a model with the reader
//  auto model = reader.fillModelFromNetcdf(GRID_SIZE_X, GRID_SIZE_Y);
// use the model for calculations
//  dump(*model.get()->topg);

namespace CUAS {

class ModelReader {
 private:
  std::unique_ptr<NetCDFFile> file;

 public:
  explicit ModelReader(std::string const &fileName);
  ModelReader(ModelReader const &) = delete;
  ModelReader &operator=(ModelReader const &) = delete;
  ModelReader(ModelReader &&) = delete;
  ModelReader &operator=(ModelReader &&) = delete;
  ~ModelReader() = default;

  std::unique_ptr<CUASModel> fillModelFromNetcdf();
  static void restartFromFile(CUASSolver &solver, std::string const &restartFile, bool restartNoneZeroInitialGuess);

  static std::unique_ptr<TimeDependentForcing> getTimeDependentForcing(
      std::string const &ncFileName, std::string const &variableName, std::vector<PetscScalar> const &xAxis,
      std::vector<PetscScalar> const &yAxis, PetscScalar multiplier = 1.0, PetscScalar offset = 0.0,
      bool loopForcing = false);

  static std::unique_ptr<BufferedForcing> getBufferedForcing(std::string const &ncFileName,
                                                             std::string const &variableName,
                                                             std::vector<PetscScalar> const &xAxis,
                                                             std::vector<PetscScalar> const &yAxis,
                                                             int numberOfSlicesPerLoad, PetscScalar multiplier = 1.0,
                                                             PetscScalar offset = 0.0, bool loopForcing = false);

  static std::unique_ptr<ScalarTimeDependentForcing> getScalarTimeDependentForcing(
      std::string const &ncFileName, std::string const &variableName, std::vector<PetscScalar> const &xAxis,
      std::vector<PetscScalar> const &yAxis, PetscScalar multiplier = 1.0, PetscScalar offset = 0.0,
      bool loopForcing = false);

  static std::unique_ptr<SteadyForcing> getSteadyForcing(std::string const &ncFileName, std::string const &variableName,
                                                         std::vector<PetscScalar> const &xAxis,
                                                         std::vector<PetscScalar> const &yAxis,
                                                         PetscScalar multiplier = 1.0, PetscScalar offset = 0.0);

  static bool isTimeDependent(std::string const &ncFileName, std::string const &variableName);

  static bool isTimeDependentField(std::string const &ncFileName, std::string const &variableName);
};

}  // namespace CUAS

#endif
