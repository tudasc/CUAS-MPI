#include "ModelReader.h"

namespace CUAS {

ModelReader::ModelReader(std::string const &fileName) { file = std::make_unique<CUASFile>(fileName, 'r'); }

std::unique_ptr<CUAS::CUASModel> ModelReader::fillModelFromNetcdf() {
  // get fileId
  int fileId = file->getFileId();
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

  auto pmodel = std::make_unique<CUAS::CUASModel>(modelDimX, modelDimY);
  auto &model = *pmodel;

  file->read("x", model.cols);
  file->read("y", model.rows);

  file->read("usurf", *model.usurf.get());
  file->read("topg", *model.topg.get());
  file->read("thk", *model.thk.get());
  file->read("bnd_mask", *model.bndMask.get());
  file->read("bmelt", *model.bmelt.get());

  model.init();

  return pmodel;
}
}  // namespace CUAS
