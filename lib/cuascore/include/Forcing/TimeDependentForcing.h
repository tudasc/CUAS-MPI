/**
 * File: TimeDependentForcing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_TIMEDEPENDENTFORCING_H
#define CUAS_TIMEDEPENDENTFORCING_H

#include "Forcing.h"

#include "Logger.h"
#include "PETScGrid.h"

#include "timeparse.h"
#include "utilities.h"

#include <algorithm>
#include <cstdlib>   // exit
#include <iterator>  // for std::back_inserter
#include <memory>
#include <vector>

namespace CUAS {

class TimeDependentForcing : public Forcing {
 public:
  explicit TimeDependentForcing(std::vector<std::unique_ptr<PETScGrid>> &forcing, std::vector<timeSecs> const &time,
                                PetscScalar const multiplier = 1.0, PetscScalar const offset = 0.0,
                                bool const loopForcing = false)
      : time(time), loopForcing(loopForcing) {
    if (time.size() != forcing.size()) {
      CUAS_ERROR("{}: time and forcing sizes are not compatible. Exiting.", __PRETTY_FUNCTION__)
      exit(1);
    }
    if (time.size() < 2) {
      CUAS_ERROR("{}: time dimension length is less than 2. Did you want to use ConstantForcing? Exiting.",
                 __PRETTY_FUNCTION__)
      exit(1);
    }
    if (!isIncreasing(time)) {
      CUAS_ERROR("{}: time is not strictly increasing. Exiting.", __PRETTY_FUNCTION__)
      exit(1);
    }

    currQ = std::make_unique<PETScGrid>(forcing[0]->getTotalNumOfCols(), forcing[0]->getTotalNumOfRows());
    std::move(begin(forcing), end(forcing), std::back_inserter(forcingStack));

    if (multiplier != 1.0) {
      TimeDependentForcing::applyMultiplier(multiplier);
    }
    if (offset != 0.0) {
      TimeDependentForcing::applyOffset(offset);
    }
  }
  TimeDependentForcing(TimeDependentForcing &) = delete;
  TimeDependentForcing(TimeDependentForcing &&) = delete;

  PETScGrid const &getCurrentQ(timeSecs currTime = 0) override {
    if (currTime < 0) {
      CUAS_ERROR("{}: getCurrentQ was called with currTime < 0. Exiting.", __PRETTY_FUNCTION__)
      exit(1);
    }

    if (loopForcing) {
      currTime = currTime % time.back();
    } else if (currTime > time.back()) {
      CUAS_WARN_RANK0(
          "{} was called with currTime > time.back(). Using last Q of forcingStack. Consider "
          "using --loopForcing argument.",
          __PRETTY_FUNCTION__)
      return *forcingStack.back();
    }

    if (currTime < time.front()) {
      CUAS_WARN_RANK0("{} was called with currTime < time.front(). Using first Q of forcingStack.", __PRETTY_FUNCTION__)
      return *forcingStack.front();
    }

    // this requires time to be sorted
    auto upperBound = std::upper_bound(time.begin(), time.end(), currTime) - time.begin();
    auto lowerBound = upperBound - 1;

    if (time[lowerBound] == currTime) {
      return *forcingStack[lowerBound];
    }

    // compute weights for linear interpolation
    // diff must a real type (double|float); with "auto" diff would be of type timeSecs aka long
    // and division later would fail.
    auto diff = time[upperBound] - time[lowerBound];
    auto wUpper = (PetscScalar)(currTime - time[lowerBound]) / (PetscScalar)diff;
    auto wLower = (PetscScalar)(time[upperBound] - currTime) / (PetscScalar)diff;

    {
      auto &fLower = forcingStack[lowerBound]->getReadHandle();
      auto &fUpper = forcingStack[upperBound]->getReadHandle();
      auto rows = forcingStack[lowerBound]->getLocalNumOfRows();
      auto cols = forcingStack[lowerBound]->getLocalNumOfCols();
      auto currQWrite = currQ->getWriteHandle();
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
          currQWrite(i, j) = fUpper(i, j) * wUpper + fLower(i, j) * wLower;
        }
      }
    }
    return *currQ;
  }

 private:
  std::vector<timeSecs> const time;
  std::vector<std::unique_ptr<PETScGrid>> forcingStack;
  bool const loopForcing;
  std::unique_ptr<PETScGrid> currQ;

  void applyMultiplier(PetscScalar multiplier) override {
    for (auto &forcing : forcingStack) {
      forcing->applyMultiplier(multiplier);
    }
  }

  void applyOffset(PetscScalar offset) override {
    for (auto &forcing : forcingStack) {
      forcing->applyOffset(offset);
    }
  }
};
}  // namespace CUAS

#endif
