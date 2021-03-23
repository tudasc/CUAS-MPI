#ifndef CUAS_FORCING_H
#define CUAS_FORCING_H

#include "PETScGrid.h"

namespace CUAS {

class Forcing {
 public:
  virtual PETScGrid const &getCurrentQ(PetscScalar currTime = 0.0) = 0;
};

}  // namespace CUAS

#endif
