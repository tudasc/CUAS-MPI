#ifndef CUAS_MODEL_READER_H
#define CUAS_MODEL_READER_H

#include "CUASModel.h"
#include "CUASSolver.h"
#include "NetCDFFile.h"

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
};
}  // namespace CUAS

#endif
