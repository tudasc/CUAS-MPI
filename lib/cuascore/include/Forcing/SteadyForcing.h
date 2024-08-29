/**
 * File: SteadyForcing.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_STEADYFORCING_H
#define CUAS_STEADYFORCING_H

#include "Forcing.h"

#include "PETScGrid.h"

#include "timeparse.h"

namespace CUAS {

class SteadyForcing : public Forcing {
 public:
  explicit SteadyForcing(PETScGrid const &m_forcing, PetscScalar const multiplier = 1.0, PetscScalar const offset = 0.0)
      : forcing(m_forcing.getTotalNumOfCols(), m_forcing.getTotalNumOfRows()) {
    forcing.copy(m_forcing);

    SteadyForcing::applyMultiplier(multiplier);
    SteadyForcing::applyOffset(offset);
  };
  SteadyForcing(const SteadyForcing &) = delete;
  SteadyForcing &operator=(SteadyForcing const &) = delete;
  SteadyForcing(const SteadyForcing &&) = delete;
  SteadyForcing &operator=(SteadyForcing const &&) = delete;
  ~SteadyForcing() = default;

  PETScGrid const &getCurrentQ(timeSecs /*currTime*/) override { return forcing; }

 private:
  PETScGrid forcing;

  void applyMultiplier(PetscScalar multiplier) override {
    if (multiplier == 1.0) {
      return;
    }
    forcing.applyMultiplier(multiplier);
  }

  void applyOffset(PetscScalar offset) override {
    if (offset == 0.0) {
      return;
    }
    forcing.applyOffset(offset);
  }
};

}  // namespace CUAS

#endif
