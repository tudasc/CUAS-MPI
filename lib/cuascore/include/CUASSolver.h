#ifndef CUAS_SOLVER_H
#define CUAS_SOLVER_H

#include "CUASArgs.h"
#include "CUASModel.h"

#include "PETScGrid.h"

namespace CUAS {

void solve(std::unique_ptr<PETScGrid> &u, std::unique_ptr<PETScGrid> &u_n, CUASModel &model, int const Nt,
           CUASArgs const &args, PetscScalar const totaltime_secs, PetscScalar const dt_secs);
}  // namespace CUAS

#endif