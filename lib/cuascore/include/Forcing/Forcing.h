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
  virtual PETScGrid const &getCurrentQ(timeSecs currTime = 0) = 0;
  virtual ~Forcing() = default;

 private:
  virtual void applyMultiplier(PetscScalar multiplier) = 0;
  virtual void applyOffset(PetscScalar offset) = 0;
};

}  // namespace CUAS

#endif
