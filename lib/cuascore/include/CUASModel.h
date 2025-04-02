/**
 * File: CUASModel.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_MODEL_H
#define CUAS_MODEL_H

#include "Forcing/Forcing.h"
#include "WaterSource.h"

#include "PETScGrid.h"

#include <memory>
#include <string>
#include <vector>

namespace CUAS {

class CUASModel : public WaterSource {
 public:
  explicit CUASModel(int numOfCols, int numOfRows);
  CUASModel(CUASModel const &) = delete;
  CUASModel &operator=(CUASModel const &) = delete;
  CUASModel(CUASModel &&) = delete;
  CUASModel &operator=(CUASModel &&) = delete;
  ~CUASModel() = default;

  // member functions
 public:
  void init();

  [[nodiscard]] bool providesWaterSource() const override;
  PETScGrid const &getCurrentWaterSource(timeSecs currTime) override;
  void setWaterSource(std::unique_ptr<Forcing> waterSource);

  // member
 public:
  // x = cols, y = rows
  int const Ncols, Nrows;  // TODO equivalent xAxis.size() and yAxis.size()
  PetscScalar dx, dy;
  std::vector<PetscScalar> xAxis, yAxis;
  std::unique_ptr<PETScGrid> topg;
  std::unique_ptr<PETScGrid> thk;
  std::unique_ptr<PETScGrid> bndMask;
  std::unique_ptr<PETScGrid> pIce;

  // member
 private:
  std::unique_ptr<Forcing> waterSource;
  // member functions
 private:
  static PetscScalar getAxisSpacing(std::vector<PetscScalar> const &axis, const std::string &attName);
};

}  // namespace CUAS

#endif
