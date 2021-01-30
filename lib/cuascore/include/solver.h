#ifndef CUAS_SOLVER_H
#define CUAS_SOLVER_H

#include "CUASModel.h"
#include "PetscGrid.h"
#include "parseCxxopts.h"

void solve(std::unique_ptr<PetscGrid> &u, std::unique_ptr<PetscGrid> &u_n, int const Nt, CUAS::CUASModel &model,
           CUAS::CUASArgs const &args, PetscScalar const totaltime_secs, PetscScalar const dt_secs);
#endif
