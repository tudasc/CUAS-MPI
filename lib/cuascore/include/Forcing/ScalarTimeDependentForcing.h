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
  ScalarTimeDependentForcing(ScalarTimeDependentForcing const &) = delete;
  ScalarTimeDependentForcing &operator=(ScalarTimeDependentForcing const &) = delete;
  ScalarTimeDependentForcing(ScalarTimeDependentForcing &&) = delete;
  ScalarTimeDependentForcing &operator=(ScalarTimeDependentForcing &&) = delete;
  ~ScalarTimeDependentForcing() override = default;

  // member functions
 public:
  PETScGrid const &getCurrent(timeSecs currentTime) override;

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

  current = std::make_unique<PETScGrid>(nx, ny);
  current->setZero();

  ScalarTimeDependentForcing::applyMultiplier(multiplier);
  ScalarTimeDependentForcing::applyOffset(offset);
}

inline PETScGrid const &ScalarTimeDependentForcing::getCurrent(timeSecs currentTime) {
  if (currentTime < 0) {
    CUAS_ERROR("{} was called with currentTime < 0. Exiting.", __PRETTY_FUNCTION__)
    exit(1);
  }

  if (loopForcing) {
    // this condition would change the behaviour of loop forcing,
    // keep consistent with BufferedForcing(loopForcing)
    // if (currentTime > allTime.back()) {
    currentTime = currentTime % time.back();
    // }
  } else if (currentTime > time.back()) {
    CUAS_WARN_RANK0(
        "{} was called with currentTime > time.back(). Using last slice of forcingStack. "
        "Consider using --loopForcing argument.",
        __PRETTY_FUNCTION__)
    auto value = timeSeries.back();
    return setValueAtPos(value);
  }

  if (currentTime < time.front()) {
    CUAS_WARN_RANK0("{} was called with currentTime < time.front(). Using first slice of forcingStack.",
                    __PRETTY_FUNCTION__)
    auto value = timeSeries.front();
    return setValueAtPos(value);
  }

  // this requires time to be sorted
  auto upperBound = std::upper_bound(time.begin(), time.end(), currentTime) - time.begin();
  auto lowerBound = upperBound - 1;

  if (time[lowerBound] == currentTime) {
    auto value = timeSeries[lowerBound];
    return setValueAtPos(value);
  }

  // compute weights for linear interpolation
  auto diff = time[upperBound] - time[lowerBound];
  auto wUpper = (PetscScalar)(currentTime - time[lowerBound]) / (PetscScalar)diff;
  auto wLower = (PetscScalar)(time[upperBound] - currentTime) / (PetscScalar)diff;
  auto value = timeSeries[upperBound] * wUpper + timeSeries[lowerBound] * wLower;

  return setValueAtPos(value);
}

inline PETScGrid const &ScalarTimeDependentForcing::setValueAtPos(PetscScalar value) {
  auto rows = current->getLocalNumOfRows();
  auto cols = current->getLocalNumOfCols();
  auto cornerX = current->getCornerX();
  auto cornerY = current->getCornerY();

  {
    auto currentWrite = current->getWriteHandle();
    // target point belongs to one processor only
    auto row = rowIndex - cornerY;
    auto col = colIndex - cornerX;
    if (row >= 0 && row < rows && col >= 0 && col < cols) {
      currentWrite(row, col) = value;
    }
  }

  return *current;
}

}  // namespace CUAS

#endif
