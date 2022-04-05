#ifndef CUAS_TIMEFORCING_H
#define CUAS_TIMEFORCING_H

#include "Forcing.h"

#include "Logger.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace CUAS {

class TimeForcing : public Forcing {
 public:
  explicit TimeForcing(std::vector<std::unique_ptr<PETScGrid>> &forcing, PetscScalar const supplyMultiplier,
                       std::vector<int> const &time_forcing, bool const loopForcing)
      : time_forcing(time_forcing), loopForcing(loopForcing) {
    if (time_forcing.size() != forcing.size()) {
      CUAS_ERROR("TimeForcing.h: time_forcing and forcing sizes are not compatible. Exiting.");
      exit(1);
    }
    if (time_forcing.size() < 2) {
      CUAS_ERROR("TimeForcing.h: time_forcing is smaller than 2. Did you want to use ConstantForcing? Exiting.");
      exit(1);
    }
    currQ = std::make_unique<PETScGrid>(forcing[0]->getTotalNumOfCols(), forcing[0]->getTotalNumOfRows());
    std::move(begin(forcing), end(forcing), std::back_inserter(forcingStack));
    setup(supplyMultiplier);
  }
  TimeForcing(TimeForcing &) = delete;
  TimeForcing(TimeForcing &&) = delete;

  virtual PETScGrid const &getCurrentQ(PetscScalar currTime = 0.0) override {
    if (currTime < 0) {
      CUAS_ERROR("TimeForcing.h: getCurrentQ was called with currTime < 0. Exiting.");
      exit(1);
    }

    if (loopForcing) {
      currTime = std::fmod(currTime, time_forcing.back());
    } else if (currTime >= time_forcing.back()) {
      CUAS_WARN(
          "TimeForcing.h: getCurrentQ was called with currTime >= time_forcing.back(). Using last Q of forcingStack. "
          "Consider using --loopForcing argument.");
      return *forcingStack.back();
    }

    if (currTime <= time_forcing.front()) {
      CUAS_WARN(
          "TimeForcing.h: getCurrentQ was called with currTime <= time_forcing.front(). Using first Q of "
          "forcingStack.");
      return *forcingStack.front();
    }

    // this requires currTime_forcing to be sorted
    auto upperBound = std::upper_bound(time_forcing.begin(), time_forcing.end(), currTime) - time_forcing.begin();
    auto lowerBound = upperBound - 1;

    if (time_forcing[lowerBound] == currTime) {
      return *forcingStack[lowerBound];
    }

    auto div = time_forcing[upperBound] - time_forcing[lowerBound];
    auto fac1 = (currTime - time_forcing[lowerBound]) / div;
    auto fac2 = (time_forcing[upperBound] - currTime) / div;

    {
      auto fLower = forcingStack[lowerBound]->getReadHandle();
      auto fUpper = forcingStack[upperBound]->getReadHandle();
      auto rows = forcingStack[lowerBound]->getLocalNumOfRows();
      auto cols = forcingStack[lowerBound]->getLocalNumOfCols();
      auto currQWrite = currQ->getWriteHandle();
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
          currQWrite(i, j) = fUpper(i, j) * fac1 + fLower(i, j) * fac2;
        }
      }
    }
    return *currQ;
  }

 private:
  std::vector<int> const time_forcing;
  std::vector<std::unique_ptr<PETScGrid>> forcingStack;
  bool const loopForcing;
  std::unique_ptr<PETScGrid> currQ;

  void setup(PetscScalar supplyMultiplier) {
    for (auto &forcing : forcingStack) {
      auto fWrite = forcing->getWriteHandle();
      auto rows = forcing->getLocalNumOfRows();
      auto cols = forcing->getLocalNumOfCols();
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
          fWrite(i, j) = fWrite(i, j) / SPY * supplyMultiplier;
        }
      }
    }
  }
};
}  // namespace CUAS

#endif
