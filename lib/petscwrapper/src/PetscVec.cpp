#include "PetscVec.h"

PetscVec::PetscVec(int size) : size(size) {
  VecCreate(PETSC_COMM_WORLD, &vec);
  VecSetSizes(vec, PETSC_DECIDE, size);
  VecSetFromOptions(vec);
}

void PetscVec::setValue(int position, PetscScalar value) { VecSetValue(vec, position, value, INSERT_VALUES); }

void PetscVec::assemble() {
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);
}

void PetscVec::setConst(PetscScalar value) { VecSet(vec, value); }

void PetscVec::view() { VecView(vec, PETSC_VIEWER_STDOUT_WORLD); }

void PetscVec::zero() { VecZeroEntries(vec); }

PetscVec::~PetscVec() { VecDestroy(&vec); }
