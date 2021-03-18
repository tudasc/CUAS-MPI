#ifndef CUAS_PETSCVEC_H
#define CUAS_PETSCVEC_H

#include "petsc.h"

class PETScVec {
 public:
  explicit PETScVec(int size) : size(size) {
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, PETSC_DECIDE, size);
    VecSetFromOptions(vec);
  }
  ~PETScVec() { VecDestroy(&vec); }

  void setValue(int position, PetscScalar value) { VecSetValue(vec, position, value, INSERT_VALUES); }
  void setConst(PetscScalar value) { VecSet(vec, value); }
  void setZero() { VecZeroEntries(vec); };
  void assemble() {
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
  }

  int getSize() const { return size; }

  Vec getRaw() { return vec; }

 private:
  Vec vec;
  const int size;

  friend class PETScSolver;
};

#endif