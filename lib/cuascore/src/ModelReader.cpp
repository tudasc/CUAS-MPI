#include "ModelReader.h"

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
  file->read("bmelt", *model.bmelt);

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
}  // namespace CUAS
