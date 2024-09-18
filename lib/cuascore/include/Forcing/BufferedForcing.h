/**
 * File: BufferedForcing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_BUFFEREDFORCING_H
#define CUAS_BUFFEREDFORCING_H

#include "Forcing.h"

#include "NetCDFFile.h"
#include "timeparse.h"

#include "Logger.h"
#include "utilities.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace CUAS {

class BufferedForcing : public Forcing {
 public:
  explicit BufferedForcing(std::unique_ptr<NetCDFFile> &ncFile, std::string const &fieldName, int nx, int ny,
                           int numberOfSlicesPerLoad = 10, PetscScalar multiplier = 1.0, PetscScalar offset = 0.0,
                           bool loopForcing = false);
  BufferedForcing(const BufferedForcing &) = delete;
  BufferedForcing &operator=(BufferedForcing const &) = delete;
  BufferedForcing(const BufferedForcing &&) = delete;
  BufferedForcing &operator=(BufferedForcing const &&) = delete;
  ~BufferedForcing() override = default;

  // member functions
 public:
  PETScGrid const &getCurrent(timeSecs currentTime) override;

  // member
 public:
  // member
 private:
  std::vector<timeSecs> allTime;
  std::vector<timeSecs> time;
  std::vector<std::unique_ptr<PETScGrid>> forcingStack;

  int indexOfFirstSliceInBuffer;
  int indexOfLastSliceInBuffer;
  int endOfBuffer;
  bool reachedEndOfFile;

  std::unique_ptr<NetCDFFile> const ncFile;
  std::string const fieldName;
  int const numberOfSlicesPerLoad;

  PetscScalar const multiplier;
  PetscScalar const offset;
  bool const loopForcing;

  // member functions
 private:
  void loadSlices(int first);

  void applyMultiplier(PetscScalar multiplier) override {
    if (multiplier == 1.0) {
      return;
    }
    for (auto &forcing : forcingStack) {
      forcing->applyMultiplier(multiplier);
    }
  }

  void applyOffset(PetscScalar offset) override {
    if (offset == 0.0) {
      return;
    }
    for (auto &forcing : forcingStack) {
      forcing->applyOffset(offset);
    }
  }
};

inline BufferedForcing::BufferedForcing(std::unique_ptr<NetCDFFile> &ncFile, std::string const &fieldName, int nx,
                                        int ny, int numberOfSlicesPerLoad, PetscScalar multiplier, PetscScalar offset,
                                        bool loopForcing)
    : ncFile(std::move(ncFile)),
      fieldName(fieldName),
      numberOfSlicesPerLoad(numberOfSlicesPerLoad),
      multiplier(multiplier),
      offset(offset),
      loopForcing(loopForcing) {
  if (numberOfSlicesPerLoad < 2) {
    CUAS_ERROR("{}: buffer size ({}) has to be at least 2. Exiting.", __FILE__, numberOfSlicesPerLoad)
    exit(1);
  }

  // get number of discrete time steps
  auto nt = this->ncFile->getDimLength("time");
  if (nt < 2) {
    CUAS_ERROR("{}::{} {} time dimension length ({}) is less than 2. Did you want to use ConstantForcing? Exiting.",
               __FILE__, __LINE__, __PRETTY_FUNCTION__, nt)
    exit(1);
  }
  // make sure time dimension is also available for the requested field
  if (!this->ncFile->variableHasDimensionByName(fieldName, "time")) {
    CUAS_ERROR(
        "{}::{} {} The field <{}> has no time dimension. Please restart with "
        "matching time and forcing field. Exiting.",
        __FILE__, __LINE__, __PRETTY_FUNCTION__, fieldName)
    exit(1);
  }

  // read all discrete time steps, only forcing fields are partially loaded
  allTime.resize(nt);
  this->ncFile->read("time", allTime);

  if (loopForcing && allTime.front() != 0) {
    CUAS_ERROR("{}::{} {} loop forcing is enabled and time does not start at 0. Exiting.", __FILE__, __LINE__,
               __PRETTY_FUNCTION__)
    exit(1);
  }
  if (allTime.front() < 0) {
    CUAS_ERROR("{}::{} {} time is smaller than 0. Exiting.", __FILE__, __LINE__, __PRETTY_FUNCTION__)
    exit(1);
  }
  if (!isIncreasing(allTime)) {
    CUAS_ERROR("{}::{} {} time is not strictly increasing. Exiting.", __FILE__, __LINE__, __PRETTY_FUNCTION__)
    exit(1);
  }
  if (numberOfSlicesPerLoad > nt) {
    CUAS_ERROR("{}::{} {} number of slices per load is larger than the number of slices in the forcing file. Exiting.",
               __FILE__, __LINE__, __PRETTY_FUNCTION__)
  }

  // resize time and forcingStack as a buffer
  time.resize(numberOfSlicesPerLoad, 0);
  forcingStack.resize(numberOfSlicesPerLoad);

  for (std::unique_ptr<PETScGrid> &grid : forcingStack) {
    grid = std::make_unique<PETScGrid>(nx, ny);
  }

  loadSlices(0);

  current = std::make_unique<PETScGrid>(forcingStack[0]->getTotalNumOfCols(), forcingStack[0]->getTotalNumOfRows());
}

inline PETScGrid const &BufferedForcing::getCurrent(timeSecs currentTime) {
  if (currentTime < 0) {
    CUAS_ERROR("{}::{} {} getCurrent was called with currentTime < 0. Exiting.", __FILE__, __LINE__,
               __PRETTY_FUNCTION__)
    exit(1);
  }

  if (currentTime < time.front()) {
    CUAS_WARN(
        "BufferedForcing.h: getCurrent was called with currentTime < time.front(). Using first slice of forcingStack.")
    return *forcingStack.front();
  }

  if (loopForcing) {
    // this condition would change the behaviour of loop forcing,
    // keep consistent with SclarTimeDependentForcing(loopForcing))
    // if (currentTime > allTime.back() {
    currentTime = currentTime % allTime.back();
    // }
    if (currentTime < time.front() || currentTime > time[endOfBuffer]) {
      CUAS_INFO("BufferedForcing.h: getCurrent out of current scope, triggering load of slices (loop forcing).")
      // this requires time to be sorted
      auto upperBound = std::upper_bound(allTime.begin(), allTime.end(), currentTime) - allTime.begin();
      auto lowerBound = upperBound - 1;
      loadSlices(lowerBound);
    }
  } else {
    if (currentTime > time[endOfBuffer]) {
      CUAS_INFO("BufferedForcing.h: getCurrent was called with currentTime > time.back().")
      if (reachedEndOfFile) {
        CUAS_WARN("Return last slice, because loop forcing is disabled.")
        return *forcingStack[endOfBuffer];
      } else {
        CUAS_INFO("Triggering load of slices.")
        loadSlices(indexOfLastSliceInBuffer);
        return getCurrent(currentTime);
      }
    }
  }

  // this requires time to be sorted
  auto upperBound = std::upper_bound(time.begin(), time.begin() + endOfBuffer + 1, currentTime) - time.begin();
  auto lowerBound = upperBound - 1;

  if (time[lowerBound] == currentTime) {
    return *forcingStack[lowerBound];
  }

  // compute weights for linear interpolation
  // diff must a real type (double|float); with "auto" diff would be of type timeSecs aka long
  // and division later would fail.
  auto diff = time[upperBound] - time[lowerBound];
  auto wUpper = static_cast<PetscScalar>(currentTime - time[lowerBound]) / static_cast<PetscScalar>(diff);
  auto wLower = static_cast<PetscScalar>(time[upperBound] - currentTime) / static_cast<PetscScalar>(diff);

  {
    auto &fLower = forcingStack[lowerBound]->getReadHandle();
    auto &fUpper = forcingStack[upperBound]->getReadHandle();
    auto rows = forcingStack[lowerBound]->getLocalNumOfRows();
    auto cols = forcingStack[lowerBound]->getLocalNumOfCols();
    auto currentWrite = current->getWriteHandle();
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        currentWrite(i, j) = fUpper(i, j) * wUpper + fLower(i, j) * wLower;
      }
    }
  }
  return *current;
}

inline void BufferedForcing::loadSlices(int first) {
  indexOfFirstSliceInBuffer = first;
  // check if enough slices remain in the file to fill the entire buffer
  indexOfLastSliceInBuffer = std::min(first + numberOfSlicesPerLoad - 1, static_cast<int>(allTime.size() - 1));

  endOfBuffer = indexOfLastSliceInBuffer - indexOfFirstSliceInBuffer;
  reachedEndOfFile = indexOfLastSliceInBuffer == allTime.size() - 1;

  // copy currently buffered time steps from the list of all time steps
  for (int i = 0; i <= endOfBuffer; ++i) {
    time[i] = allTime[indexOfFirstSliceInBuffer + i];
  }

  auto numberOfLoadedDiscreteTimesteps = endOfBuffer + 1;
  ncFile->read(fieldName, forcingStack, indexOfFirstSliceInBuffer, numberOfLoadedDiscreteTimesteps);

  if (numberOfLoadedDiscreteTimesteps < numberOfSlicesPerLoad) {
    CUAS_INFO(
        "BufferedForcing.h: {} slices loaded, instead of {}. It seams, that we reached the end of the forcing "
        "series.",
        numberOfLoadedDiscreteTimesteps, numberOfSlicesPerLoad)
  }

  BufferedForcing::applyMultiplier(multiplier);
  BufferedForcing::applyOffset(offset);
}

}  // namespace CUAS

#endif
