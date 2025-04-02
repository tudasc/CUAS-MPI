/**
 * File: PETScVector.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PETSCVEC_H
#define CUAS_PETSCVEC_H

#include "petsc.h"
#include "petscwrapperutils.h"

class PETScVector {
  friend class PETScSolver;

 public:
  explicit PETScVector(int size) : size(size) {
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, PETSC_DECIDE, size);
    VecSetFromOptions(vec);
  }
  PETScVector(PETScVector const &) = delete;
  PETScVector &operator=(PETScVector const &) = delete;
  PETScVector(PETScVector &&) = delete;
  PETScVector &operator=(PETScVector &&) = delete;
  ~PETScVector() { VecDestroy(&vec); }

  // member functions
 public:
  void setValue(int position, PetscScalar value) { VecSetValue(vec, position, value, INSERT_VALUES); }
  void setConst(PetscScalar value) { VecSet(vec, value); }
  void setZero() { VecZeroEntries(vec); };

  void getOwnershipRange(int &firstIndex, int &lastIndex) const { VecGetOwnershipRange(vec, &firstIndex, &lastIndex); }
  int getSize() const { return size; }
  Vec getRaw() { return vec; }

  void assemble() {
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
  }

  // member
 public:
  // member
 private:
  Vec vec;
  const int size;

  // member functions
 private:
};

#endif
