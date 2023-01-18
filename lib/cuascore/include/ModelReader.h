/**
 * File: ModelReader.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_MODEL_READER_H
#define CUAS_MODEL_READER_H

#include "CUASModel.h"
#include "CUASSolver.h"
#include "Forcing/ConstantForcing.h"
#include "Forcing/TimeForcing.h"
#include "NetCDFFile.h"

#include <memory>

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
  std::unique_ptr<CUAS::CUASModel> fillModelFromNetcdf();
  static void restartFromFile(CUAS::CUASSolver &solver, std::string const &restartFile,
                              bool restartNoneZeroInitialGuess);
  static std::unique_ptr<CUAS::TimeForcing> getTimeForcing(std::string const &timeForcingFileName,
                                                           std::string const &fieldName, PetscScalar multiplier = 1.0,
                                                           PetscScalar offset = 0.0, bool loopForcing = false);
  static std::unique_ptr<CUAS::ConstantForcing> getConstantForcing(std::string const &timeForcingFileName,
                                                                   std::string const &fieldName,
                                                                   PetscScalar multiplier = 1.0,
                                                                   PetscScalar offset = 0.0);
  static bool isTimeDependentField(std::string const &timeForcingFileName, std::string const &fieldName);
};

}  // namespace CUAS

#endif
