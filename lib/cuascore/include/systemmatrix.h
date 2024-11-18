/**
 * File: systemmatrix.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_FILL_MATRIX_COO_H
#define CUAS_FILL_MATRIX_COO_H

#include "PETScGrid.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

namespace CUAS {

void systemmatrixDeprecated(PETScMatrix &A, PETScGrid &b, PETScGrid const &hydraulicStorativity,
                            PETScGrid const &hydraulicTransmissivity, PetscScalar dx, PetscScalar dt, PetscScalar theta,
                            PETScGrid const &hydraulicHead, PETScGrid const &waterSource,
                            PETScGrid const &dirichletValues, PETScGrid const &bndMask, PETScGrid const &globalIndices);

/**
 * Implements the system matrix using a first order upwind scheme (UDS)
 * or second order central difference scheme (CDS) depending on the interface transmissivity
 */
void systemmatrix(PETScMatrix &A, PETScGrid &b, PETScGrid const &hydraulicStorativity, PETScGrid const &effTransEast,
                  PETScGrid const &effTransWest, PETScGrid const &effTransNorth, PETScGrid const &effTransSouth,
                  PetscScalar dx, PetscScalar dt, PetscScalar theta, PETScGrid const &hydraulicHead,
                  PETScGrid const &waterSource, PETScGrid const &dirichletValues, PETScGrid const &bndMask,
                  PETScGrid const &globalIndices);

}  // namespace CUAS

#endif
