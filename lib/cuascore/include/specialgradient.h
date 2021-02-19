#ifndef CUAS_SPECIALGRADIENT_H
#define CUAS_SPECIALGRADIENT_H

#include "PETScGrid.h"

namespace CUAS {

// PetscGrid h and PetscGrid gradient need to be created beforehand. The gradient-grid needs to be of the same size as
// the grid h. The same needs to be respected for the calculation of gradient2_central.
void gradient2(PetscGrid &gradient, PetscGrid const &h, PetscScalar const dx);

// not used, not tested
// void gradient2_central(PetscGrid &gradient, PetscGrid const &h, PetscScalar const dx);

}  // namespace CUAS

#endif
