/**
 * File: CUASTimeIntegrator.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_TIMEINTEGRATOR_H
#define CUAS_TIMEINTEGRATOR_H

#include "timeparse.h"

#include <utility>
#include <vector>

namespace CUAS {

class CUASTimeIntegrator {
 public:
  explicit CUASTimeIntegrator(std::vector<timeSecs> const &timeSteps)
      : timeSteps(timeSteps), timeStepIndex(0), currTime(0), currDt(0), maxNextDt(0) {
    if (!timeSteps.empty() && timeSteps[0] < 0) {
      CUAS_ERROR("{}::{}: first time step is smaller than 0: {}.", __PRETTY_FUNCTION__, __LINE__, timeSteps[0])
      exit(1);
    }
    for (int i = 1; i < timeSteps.size(); ++i) {
      if (timeSteps[i] <= timeSteps[i - 1]) {
        CUAS_ERROR("{}::{}: time steps are not increasing: timeSteps[{}]={}, timeSteps[{}]={}.", __PRETTY_FUNCTION__,
                   __LINE__, i - 1, timeSteps[i - 1], i, timeSteps[i])
        exit(1);
      }
    }

    if (!timeSteps.empty()) {
      currTime = timeSteps[0];
    }

    if (timeSteps.size() >= 2) {
      maxNextDt = timeSteps[1] - timeSteps[0];
    }
  }
  CUASTimeIntegrator(CUASTimeIntegrator &) = delete;
  CUASTimeIntegrator &operator=(CUASTimeIntegrator const &) = delete;
  CUASTimeIntegrator(CUASTimeIntegrator &&) = delete;
  CUASTimeIntegrator &operator=(CUASTimeIntegrator const &&) = delete;
  ~CUASTimeIntegrator() = default;

  // member
 private:
  const std::vector<timeSecs> timeSteps;
  decltype(timeSteps)::size_type timeStepIndex;
  timeSecs currTime;
  timeSecs currDt;
  timeSecs maxNextDt;

  // member functions
 public:
  void finalizeTimestep(timeSecs sizeOfComputedTimeStep) {
    if (sizeOfComputedTimeStep < 0) {
      CUAS_ERROR("{}::{}: sizeOfComputedTimeStep {} is smaller than 0.", __PRETTY_FUNCTION__, __LINE__,
                 sizeOfComputedTimeStep)
      exit(1);
    }
    if (sizeOfComputedTimeStep > maxNextDt) {
      // ERROR a larger time step has been calculated
      CUAS_ERROR("{}::{}: sizeOfComputedTimeStep {} is larger than the time step the integrator calculated {}.",
                 __PRETTY_FUNCTION__, __LINE__, sizeOfComputedTimeStep, currDt)
      exit(1);
    }

    if (sizeOfComputedTimeStep == maxNextDt) {
      // the time step has been fully calculated
      ++timeStepIndex;
      if (timeStepIndex < timeSteps.size() - 1) {
        // more time steps to be computed
        currTime = timeSteps[timeStepIndex];
        auto nextTime = timeSteps[timeStepIndex + 1];
        maxNextDt = nextTime - currTime;
      } else {
        // reached the last time step
        currTime = timeSteps.back();
        maxNextDt = 0;
      }
    } else if (sizeOfComputedTimeStep < maxNextDt) {
      // the time step has not been fully calculated
      // we do not increase timeStepIndex, because we are still in the same intervall
      currTime += sizeOfComputedTimeStep;
      maxNextDt -= sizeOfComputedTimeStep;
    }

    currDt = 0;
  }

  std::pair<timeSecs, timeSecs> getTimestepInformation(timeSecs input) {
    if (currDt != 0) {
      CUAS_ERROR("{}::{}: new time step information generated, while currDt is not 0. Call finalizeTimestep first.",
                 __PRETTY_FUNCTION__, __LINE__)
      exit(1);
    }

    currDt = std::min(input, maxNextDt);
    return std::make_pair(currTime, currDt);
  }

  [[nodiscard]] timeSecs getCurrentTime() const { return currTime; }
  [[nodiscard]] timeSecs getCurrentDt() const { return currDt; }
  [[nodiscard]] decltype(timeSteps) const &getTimesteps() const { return timeSteps; }
  [[nodiscard]] decltype(timeSteps)::size_type getTimestepIndex() const { return timeStepIndex; }
};

// deprecated
/*std::pair<timeSecs, timeSecs> getTimestepInformation(std::vector<timeSecs> const &timeSteps,
                                                                 int timeStepIndex) {
  timeSecs dt = 0;
  auto currTime = timeSteps[timeStepIndex];
  if (timeSteps.size() > 1) {
    if (timeStepIndex < timeSteps.size() - 1) {
      auto nextTime = timeSteps[timeStepIndex + 1];
      dt = nextTime - currTime;
    } else {
      // This is the last iteration only used to recompute the diagnostic fields for saving
      dt = -9999;  // Something stupid to let the solver crash if used
    }
  }

  return std::make_pair(currTime, dt);
}*/

}  // namespace CUAS

#endif
