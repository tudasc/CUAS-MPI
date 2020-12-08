#ifndef CUAS_FILL_MATRIX_COO_H
#define CUAS_FILL_MATRIX_COO_H

#include "PetscGrid.h"
#include "PetscMat.h"
#include "PetscVec.h"

namespace CUAS {

void fill_matrix_coo(PetscMat &A, PetscVec &b, int const Nx, int const Ny, PetscGrid &S, PetscGrid &T,
                     PetscScalar const dx, PetscScalar const dt, PetscScalar const theta, PetscGrid &u_n, PetscGrid &Q,
                     PetscGrid &dirichlet_values, PetscGrid &dirichlet_mask);
}  // namespace CUAS

#endif
