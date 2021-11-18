#ifndef CUAS_FILL_MATRIX_COO_H
#define CUAS_FILL_MATRIX_COO_H

#include "PETScGrid.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

namespace CUAS {

void systemmatrix(PETScMatrix &A, PETScVector &b, int const Nx, int const Ny, PETScGrid const &hydraulicStorativity,
                  PETScGrid const &hydraulicTransmissivity, PetscScalar const dx, PetscScalar const dt,
                  PetscScalar const theta, PETScGrid const &hydraulicHead, PETScGrid const &Q,
                  PETScGrid const &dirichletValues, PETScGrid const &bndMask);
}  // namespace CUAS

#endif
