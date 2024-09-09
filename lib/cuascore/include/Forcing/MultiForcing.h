/**
 * File: MultiForcing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_MULTIFORCING_H
#define CUAS_MULTIFORCING_H

#include "Forcing.h"

#include "PETScGrid.h"

#include "timeparse.h"

namespace CUAS {

class MultiForcing : public Forcing {
 public:
  explicit MultiForcing(int nCols, int nRows) : currQ(std::make_unique<PETScGrid>(nCols, nRows)) {}
  MultiForcing(const MultiForcing &) = delete;
  MultiForcing &operator=(MultiForcing const &) = delete;
  MultiForcing(const MultiForcing &&) = delete;
  MultiForcing &operator=(MultiForcing const &&) = delete;
  ~MultiForcing() override = default;

  // member functions
 public:
  PETScGrid const &getCurrentQ(timeSecs currTime) override;

  void registerNewForcing(std::unique_ptr<Forcing> &forcing);

  // member
 public:
  // member
 private:
  std::unique_ptr<PETScGrid> currQ;
  std::vector<std::unique_ptr<Forcing>> forcings;

  // member functions
 private:
  void applyMultiplier(PetscScalar multiplier) override {}

  void applyOffset(PetscScalar offset) override {}
};

inline PETScGrid const &MultiForcing::getCurrentQ(timeSecs currTime) {
  if (forcings.size() == 1) {
    CUAS_WARN("MultiForcing contains only one forcing.")
    return forcings[0]->getCurrentQ(currTime);
  }

  currQ->setZero();

  if (forcings.empty()) {
    CUAS_WARN("MultiForcing does not contain any forcings.")
    return *currQ;
  }

  auto writeHandle = currQ->getWriteHandle();

  for (auto &f : forcings) {
    auto &data = f->getCurrentQ(currTime);
    auto &readHandle = data.getReadHandle();
    // currQ += data;
    for (int j = 0; j < currQ->getLocalNumOfRows(); ++j) {
      for (int i = 0; i < currQ->getLocalNumOfCols(); ++i) {
        writeHandle(j, i) = writeHandle(j, i) + readHandle(j, i);
      }
    }
  }

  return *currQ;
}

inline void MultiForcing::registerNewForcing(std::unique_ptr<Forcing> &forcing) {
  if (!currQ->isCompatible(forcing->getCurrentQ(0))) {
    CUAS_ERROR("{}::{} {} Forcing is incompatible! Exiting.", __FILE__, __LINE__, __PRETTY_FUNCTION__)
    exit(1);
  }

  forcings.emplace_back(std::move(forcing));
}

}  // namespace CUAS

#endif
