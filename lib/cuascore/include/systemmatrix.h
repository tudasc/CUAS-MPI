#ifndef CUAS_FILL_MATRIX_COO_H
#define CUAS_FILL_MATRIX_COO_H

#include "PETScGrid.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

namespace CUAS {

void systemmatrix(PETScMatrix &A, PETScVector &b, int const Nx, int const Ny, PETScGrid const &S, PETScGrid const &T,
                  PetscScalar const dx, PetscScalar const dt, PetscScalar const theta, PETScGrid const &u_n,
                  PETScGrid const &Q, PETScGrid const &dirichlet_values, PETScGrid const &dirichlet_mask);
}  // namespace CUAS

#endif
