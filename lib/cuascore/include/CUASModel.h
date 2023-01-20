/**
 * File: CUASModel.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_MODEL_H
#define CUAS_MODEL_H

#include "Forcing/Forcing.h"

#include "PETScGrid.h"

#include <memory>
#include <vector>

namespace CUAS {

class CUASModel {
 public:
  explicit CUASModel(int numOfCols, int numOfRows);
  CUASModel(CUASModel &) = delete;
  CUASModel(CUASModel &&) = delete;
  //~CUASModel();

  // x = cols, y = rows
  int const Ncols, Nrows;  // TODO equivalent xAxis.size() and yAxis.size()
  PetscScalar dx, dy;
  std::vector<PetscScalar> xAxis, yAxis;
  std::unique_ptr<PETScGrid> topg;
  std::unique_ptr<PETScGrid> thk;
  std::unique_ptr<PETScGrid> bndMask;
  std::unique_ptr<PETScGrid> pIce;
  std::unique_ptr<Forcing> Q;

  void init();
};

}  // namespace CUAS

#endif
