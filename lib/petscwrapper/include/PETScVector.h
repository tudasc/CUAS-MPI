/**
 * File: PETScVector.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PETSCVEC_H
#define CUAS_PETSCVEC_H

#include "petsc.h"

class PETScVector {
 public:
  explicit PETScVector(int size) : size(size) {
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, PETSC_DECIDE, size);
    VecSetFromOptions(vec);
  }
  PETScVector(PETScVector &) = delete;
  PETScVector(PETScVector &&) = delete;
  ~PETScVector() { VecDestroy(&vec); }

  void setValue(int position, PetscScalar value) { VecSetValue(vec, position, value, INSERT_VALUES); }
  void setConst(PetscScalar value) { VecSet(vec, value); }
  void setZero() { VecZeroEntries(vec); };
  void assemble() {
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
  }

  void getOwnershipRange(int &firstIndex, int &lastIndex) const { VecGetOwnershipRange(vec, &firstIndex, &lastIndex); }

  int getSize() const { return size; }

  Vec getRaw() { return vec; }

 private:
  Vec vec;
  const int size;

  friend class PETScSolver;
};

#endif
