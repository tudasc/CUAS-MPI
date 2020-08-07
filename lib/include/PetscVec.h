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