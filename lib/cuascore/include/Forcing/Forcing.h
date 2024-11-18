/**
 * File: Forcing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_FORCING_H
#define CUAS_FORCING_H

#include "PETScGrid.h"

#include "timeparse.h"

namespace CUAS {

class Forcing {
 public:
  explicit Forcing() = default;
  Forcing(const Forcing &) = delete;
  Forcing &operator=(Forcing const &) = delete;
  Forcing(const Forcing &&) = delete;
  Forcing &operator=(Forcing const &&) = delete;
  virtual ~Forcing() = default;

  // member functions
 public:
  virtual PETScGrid const &getCurrent(timeSecs currentTime) = 0;

  // member
 public:
  // member
 protected:
  std::unique_ptr<PETScGrid> current;
  // member functions
 private:
  virtual void applyMultiplier(PetscScalar multiplier) = 0;
  virtual void applyOffset(PetscScalar offset) = 0;
};

}  // namespace CUAS

#endif
