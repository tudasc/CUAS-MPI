#ifndef CUAS_FILL_MATRIX_COO_H
#define CUAS_FILL_MATRIX_COO_H

#include "PETScGrid.h"
#include "PETScMat.h"
#include "PETScVec.h"

namespace CUAS {

void fill_matrix_coo(PetscMat &A, PetscVec &b, int const Nx, int const Ny, PetscGrid const &S, PetscGrid const &T,
                     PetscScalar const dx, PetscScalar const dt, PetscScalar const theta, PetscGrid const &u_n,
                     PetscGrid const &Q, PetscGrid const &dirichlet_values, PetscGrid const &dirichlet_mask);
}  // namespace CUAS

#endif
