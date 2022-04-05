#ifndef CUAS_CONSTANTFORCING_H
#define CUAS_CONSTANTFORCING_H

#include "CUASConstants.h"
#include "Forcing.h"

namespace CUAS {

class ConstantForcing : public Forcing {
 public:
  explicit ConstantForcing(PETScGrid const &m_forcing, PetscScalar const supplyMultiplier)
      : forcing(m_forcing.getTotalNumOfCols(), m_forcing.getTotalNumOfRows()) {
    // TODO this copy is not necessary we can read the data of m_forcing in the loop
    forcing.copy(m_forcing);

    auto fWrite = forcing.getWriteHandle();
    for (int i = 0; i < forcing.getLocalNumOfRows(); ++i) {
      for (int j = 0; j < forcing.getLocalNumOfCols(); ++j) {
        fWrite(i, j) = fWrite(i, j) / SPY * supplyMultiplier;
      }
    }
  };
  ConstantForcing(ConstantForcing &) = delete;
  ConstantForcing(ConstantForcing &&) = delete;

  virtual PETScGrid const &getCurrentQ(PetscScalar currTime = 0.0) override { return forcing; }

 private:
  PETScGrid forcing;
};

}  // namespace CUAS

#endif