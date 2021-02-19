#ifndef CUAS_SAVETIMESTEPS_H
#define CUAS_SAVETIMESTEPS_H

#include "CUASArgs.h"
#include "CUASModel.h"
#include "specialgradient.h"

#include "PETScGrid.h"

namespace CUAS {

inline void saveSolution(int timeStep, CUASArgs const &args, int rank, PetscGrid const &u, PetscGrid const &u_n,
                         CUASModel &model, PetscGrid &melt, PetscScalar cavity_opening) {
  if (args.verbose && rank == 0) {
    // TODO: implement time_display_units from helper to get time_scaling
    // std::cout << "time(" << i << "/" << Nt << ") = " << time_current/time_scaling << " (" << time_units << ")"
    // << std::endl;
  }
  // eps_head: get highest differing value from u and u_n
  // eps_T: get highest differing value from T and T_n
  PetscScalar eps_head;
  PetscScalar eps_T;
  auto &uGlobal = u.getReadHandle();
  auto &u_nGlobal = u_n.getReadHandle();
  auto &TGlobal = model.T->getReadHandle();
  auto &T_nGlobal = model.T_n->getReadHandle();
  for (int j = 0; j < u.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < u.getLocalNumOfCols(); ++i) {
      PetscScalar diffHead = abs(uGlobal(j, i) - u_nGlobal(j, i));
      PetscScalar diffT = abs(TGlobal(j, i) - T_nGlobal(j, i));
      if (diffHead > eps_head) {
        eps_head = diffHead;
      }
      if (diffT > eps_T) {
        eps_T = diffT;
      }
    }
  }

  auto meltGlobal = melt.getWriteHandle();
  for (int j = 0; j < melt.getLocalNumOfRows(); ++j) {
    for (int i = 0; i < melt.getLocalNumOfCols(); ++i) {
      meltGlobal(j, i) = meltGlobal(j, i) + cavity_opening;
    }
  }
  // should be out of scope here anyways
  meltGlobal.setValues();
  // solution.saveTimestep();
}

}  // namespace CUAS

#endif
