/**
 * File: ScalarTimeDependentForcing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_SCALARTIMEDEPENDENTFORCING_H
#define CUAS_SCALARTIMEDEPENDENTFORCING_H

#include "Forcing.h"

#include "Logger.h"
#include "PETScGrid.h"

#include "timeparse.h"
#include "utilities.h"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <vector>

namespace CUAS {

class ScalarTimeDependentForcing : public Forcing {
 public:
  explicit ScalarTimeDependentForcing(int nx, int ny, int rowIndex, int colIndex, std::vector<PetscScalar> &timeSeries,
                                      std::vector<timeSecs> const &time, PetscScalar multiplier = 1.0,
                                      PetscScalar offset = 0.0, bool loopForcing = false);
  ScalarTimeDependentForcing(const ScalarTimeDependentForcing &) = delete;
  ScalarTimeDependentForcing &operator=(ScalarTimeDependentForcing const &) = delete;
  ScalarTimeDependentForcing(const ScalarTimeDependentForcing &&) = delete;
  ScalarTimeDependentForcing &operator=(ScalarTimeDependentForcing const &&) = delete;
  ~ScalarTimeDependentForcing() override = default;

  // member functions
 public:
  PETScGrid const &getCurrentQ(timeSecs currTime) override;

  // member
 public:
  // member
 private:
  std::vector<timeSecs> const time;
  std::vector<PetscScalar> timeSeries;

  bool const loopForcing;
  int const rowIndex;
  int const colIndex;

  // member functions
 private:
  PETScGrid const &setValueAtPos(PetscScalar value);

  // Note, we apply the multiplier directly on the scalar time series input.
  void applyMultiplier(PetscScalar multiplier) override {
    if (multiplier == 1.0) {
      return;
    }
    std::transform(this->timeSeries.begin(), this->timeSeries.end(), this->timeSeries.begin(),
                   [&multiplier](auto &el) { return el * multiplier; });
  }

  // Note, we apply the offset directly on the scalar time series input.
  void applyOffset(PetscScalar offset) override {
    if (offset == 0.0) {
      return;
    }
    std::transform(this->timeSeries.begin(), this->timeSeries.end(), this->timeSeries.begin(),
                   [&offset](auto &el) { return el + offset; });
  }
};

inline ScalarTimeDependentForcing::ScalarTimeDependentForcing(int nx, int const ny, int rowIndex, int colIndex,
                                                              std::vector<PetscScalar> &timeSeries,
                                                              std::vector<timeSecs> const &time, PetscScalar multiplier,
                                                              PetscScalar offset, bool loopForcing)
    : rowIndex(rowIndex), colIndex(colIndex), time(time), timeSeries(timeSeries), loopForcing(loopForcing) {
  if (time.size() != timeSeries.size()) {
    CUAS_ERROR("{}: time and time series sizes are not compatible. Exiting.", __PRETTY_FUNCTION__)
    exit(1);
  }
  if (time.empty()) {
    CUAS_ERROR("{}: time dimension length is less than 1. Exiting.", __PRETTY_FUNCTION__)
    exit(1);
  }
  if (!isIncreasing(time)) {
    CUAS_ERROR("{}: time is not strictly increasing. Exiting.", __PRETTY_FUNCTION__)
    exit(1);
  }

  currQ = std::make_unique<PETScGrid>(nx, ny);
  currQ->setZero();

  ScalarTimeDependentForcing::applyMultiplier(multiplier);
  ScalarTimeDependentForcing::applyOffset(offset);
}

inline PETScGrid const &ScalarTimeDependentForcing::getCurrentQ(timeSecs currTime) {
  if (currTime < 0) {
    CUAS_ERROR("{} was called with currTime < 0. Exiting.", __PRETTY_FUNCTION__)
    exit(1);
  }

  if (loopForcing) {
    // this condition would change the behaviour of loop forcing,
    // keep consistent with BufferedForcing(loopForcing)
    // if (currTime > allTime.back()) {
    currTime = currTime % time.back();
    // }
  } else if (currTime > time.back()) {
    CUAS_WARN_RANK0(
        "{} was called with currTime > time.back(). Using last Q of forcingStack. "
        "Consider "
        "using --loopForcing argument.",
        __PRETTY_FUNCTION__)
    auto value = timeSeries.back();
    return setValueAtPos(value);
  }

  if (currTime < time.front()) {
    CUAS_WARN_RANK0("{} was called with currTime < time.front(). Using first Q of forcingStack.", __PRETTY_FUNCTION__)
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

inline PETScGrid const &ScalarTimeDependentForcing::setValueAtPos(PetscScalar value) {
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

}  // namespace CUAS

#endif
