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
