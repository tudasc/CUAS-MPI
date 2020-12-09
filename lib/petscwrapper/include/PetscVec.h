#ifndef CUAS_PETSCVEC_H
#define CUAS_PETSCVEC_H

#include "petsc.h"

class PetscVec {
 public:
  PetscVec(int size);
  ~PetscVec();
  void setValue(int position, PetscScalar value);
  void setConst(PetscScalar value);
  void setZero() { setConst(0); };
  void assemble();
  void zero();

  Vec getPetscRaw() { return vec; }
  int getSize() const { return size; }

 private:
  Vec vec;
  const int size;

  friend class PetscSolver;
};

#endif
