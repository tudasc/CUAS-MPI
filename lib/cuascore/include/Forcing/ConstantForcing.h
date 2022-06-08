#ifndef CUAS_CONSTANTFORCING_H
#define CUAS_CONSTANTFORCING_H

#include "Forcing.h"

#include "CUASConstants.h"

namespace CUAS {

class ConstantForcing : public Forcing {
 public:
  explicit ConstantForcing(PETScGrid const &m_forcing, PetscScalar const multiplier = 1.0,
                           PetscScalar const offset = 0.0)
      : forcing(m_forcing.getTotalNumOfCols(), m_forcing.getTotalNumOfRows()) {
    forcing.copy(m_forcing);

    ConstantForcing::applyMultiplier(multiplier);
    ConstantForcing::applyOffset(offset);
  };
  ConstantForcing(ConstantForcing &) = delete;
  ConstantForcing(ConstantForcing &&) = delete;

  PETScGrid const &getCurrentQ(timeSecs currTime = 0) override { return forcing; }

 private:
  PETScGrid forcing;

  void applyMultiplier(PetscScalar multiplier) override { forcing.applyMultiplier(multiplier); }

  void applyOffset(PetscScalar offset) override { forcing.applyOffset(offset); }
};

}  // namespace CUAS

#endif
