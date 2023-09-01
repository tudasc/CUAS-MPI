/**
 * File: ScalarTimeForcing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_SCALARTIMEFORCING_H
#define CUAS_SCALARTIMEFORCING_H

#include "Forcing.h"

#include "timeparse.h"

#include "Logger.h"
#include "utilities.h"

#include <algorithm>
#include <cmath>
#include <functional>  // vector times scalar or vector plus scalar
#include <memory>
#include <vector>

namespace CUAS {

class ScalarTimeForcing : public Forcing {
 public:
  explicit ScalarTimeForcing(int const nx, int const ny, int const rowIndex, int const colIndex,
                             std::vector<PetscScalar> &timeSeries, std::vector<timeSecs> const &time,
                             const PetscScalar multiplier = 1.0, const PetscScalar offset = 0.0,
                             const bool loopForcing = false)
      : rowIndex(rowIndex), colIndex(colIndex), time(time), timeSeries(timeSeries), loopForcing(loopForcing) {
    if (time.size() != timeSeries.size()) {
      CUAS_ERROR("ScalarTimeForcing.h: time and time series sizes are not compatible. Exiting.")
      exit(1);
    }
    if (time.empty()) {
      CUAS_ERROR("ScalarTimeForcing.h: time dimension length is less than 1. Exiting.")
      exit(1);
    }
    if (!isIncreasing(time)) {
      CUAS_ERROR("ScalarTimeForcing.h: time is not strictly increasing. Exiting.")
      exit(1);
    }

    currQ = std::make_unique<PETScGrid>(nx, ny);
    currQ->setZero();

    if (multiplier != 1.0) {
      std::transform(this->timeSeries.begin(), this->timeSeries.end(), this->timeSeries.begin(),
                     [&multiplier](auto &el) { return el * multiplier; });
    }
    if (offset != 0.0) {
      std::transform(this->timeSeries.begin(), this->timeSeries.end(), this->timeSeries.begin(),
                     [&offset](auto &el) { return el + offset; });
    }
  }
  ScalarTimeForcing(ScalarTimeForcing &) = delete;
  ScalarTimeForcing(ScalarTimeForcing &&) = delete;

  PETScGrid const &getCurrentQ(timeSecs currTime = 0) override {
    if (currTime < 0) {
      CUAS_ERROR("ScalarTimeForcing.h: getCurrentQ was called with currTime < 0. Exiting.")
      exit(1);
    }

    if (loopForcing) {
      currTime = currTime % time.back();
    } else if (currTime > time.back()) {
      CUAS_WARN(
          "ScalarTimeForcing.h: getCurrentQ was called with currTime > time.back(). Using last Q of forcingStack. "
          "Consider "
          "using --loopForcing argument.")
      auto value = timeSeries.back();
      return setValueAtPos(value);
    }

    if (currTime < time.front()) {
      CUAS_WARN(
          "ScalarTimeForcing.h: getCurrentQ was called with currTime < time.front(). Using first Q of forcingStack.")
      auto value = timeSeries.front();
      return setValueAtPos(value);
    }

    // this requires time to be sorted
    auto upperBound = std::upper_bound(time.begin(), time.end(), currTime) - time.begin();
    auto lowerBound = upperBound - 1;

    if (time[lowerBound] == currTime) {
      auto value = timeSeries[lowerBound];
      return setValueAtPos(value);
    }

    // compute weights for linear interpolation
    auto diff = time[upperBound] - time[lowerBound];
    auto wUpper = (PetscScalar)(currTime - time[lowerBound]) / (PetscScalar)diff;
    auto wLower = (PetscScalar)(time[upperBound] - currTime) / (PetscScalar)diff;
    auto value = timeSeries[upperBound] * wUpper + timeSeries[lowerBound] * wLower;
    return setValueAtPos(value);
  }

 private:
  std::unique_ptr<PETScGrid> currQ;
  std::vector<timeSecs> const time;
  std::vector<PetscScalar> timeSeries;
  bool const loopForcing;
  int rowIndex, colIndex;

  // Note, we apply the multiplier directly on the scalar time series input.
  void applyMultiplier(PetscScalar multiplier) override {}

  // Note, we apply the offset directly on the scalar time series input.
  void applyOffset(PetscScalar offset) override {}

  PETScGrid const &setValueAtPos(PetscScalar value) {
    auto rows = currQ->getLocalNumOfRows();
    auto cols = currQ->getLocalNumOfCols();
    auto cornerX = currQ->getCornerX();
    auto cornerY = currQ->getCornerY();
    auto currQWrite = currQ->getWriteHandle();

    // target point belongs to one processor only
    auto row = rowIndex - cornerY;
    auto col = colIndex - cornerX;
    if (row >= 0 && row < rows && col >= 0 && col < cols) {
      currQWrite(row, col) = value;
    }
    currQWrite.setValues();

    return *currQ;
  }
};
}  // namespace CUAS

#endif
