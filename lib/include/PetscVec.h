#ifndef CUAS_PETSCVEC_H
#define CUAS_PETSCVEC_H

#include "petsc.h"

class PetscVec {
 public:
  PetscVec(int size);
  ~PetscVec();
  void setValue(int position, PetscScalar value);
  void assemble();
  void view();
  void zero();

 private:
  Vec vec;
};

#endif