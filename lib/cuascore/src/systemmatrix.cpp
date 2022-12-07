#include "systemmatrix.h"

#include "CUASConstants.h"

namespace CUAS {

inline PetscScalar hmean(PetscScalar x1, PetscScalar x2) { return 2.0 * x1 * x2 / (x1 + x2 + TINY); }

void systemmatrix(PETScMatrix &A, PETScGrid &b, PETScGrid const &hydraulicStorativity,
                  PETScGrid const &hydraulicTransmissivity, PetscScalar dx, PetscScalar dt, PetscScalar theta,
                  PETScGrid const &hydraulicHead, PETScGrid const &Q, PETScGrid const &dirichletValues,
                  PETScGrid const &bndMask, PETScGrid const &globalIndices) {
  auto &mask = bndMask.getReadHandle();
  auto &values = dirichletValues.getReadHandle();
  auto &storativity = hydraulicStorativity.getReadHandle();
  auto &transmissivity = hydraulicTransmissivity.getReadHandle();
  auto &head = hydraulicHead.getReadHandle();
  auto &source = Q.getReadHandle();
  auto &globalIndicesHandle = globalIndices.getReadHandle();

  auto bHandle = b.getWriteHandle();

  // not sure if necessary
  b.setConst(0);
  A.setZero();

  const auto one_dx2 = 1.0 / (dx * dx);

  auto cornerX = hydraulicStorativity.getCornerX();
  auto cornerY = hydraulicStorativity.getCornerY();
  auto cornerXGhost = hydraulicStorativity.getCornerXGhost();
  auto cornerYGhost = hydraulicStorativity.getCornerYGhost();

  // we take localnumOfCols, so that we can iterate over the normal boundaries
  for (int row = 0; row < hydraulicStorativity.getLocalNumOfRows(); ++row) {
    for (int col = 0; col < hydraulicStorativity.getLocalNumOfCols(); ++col) {
      auto p = globalIndicesHandle(row, col);

      if (mask(row, col) != (PetscScalar)COMPUTE_FLAG) {
        A.setValue(p, p, 1.0);
        bHandle(row, col) = values(row, col);
      } else {
        // needed for transmissivity and head in order to use LocalArray rather than GlobalArray
        auto ghostLocalRow = row + (cornerY - cornerYGhost);
        auto ghostLocalCol = col + (cornerX - cornerXGhost);

        auto sP = storativity(row, col);
        auto tP = transmissivity(ghostLocalRow, ghostLocalCol, GHOSTED);

        // fill A
        auto tE = transmissivity(ghostLocalRow, ghostLocalCol + 1, GHOSTED);
        auto dE = hmean(tP, tE) * one_dx2;
        auto aE = -theta * dt / sP * dE;
        auto iE = globalIndicesHandle(ghostLocalRow, ghostLocalCol + 1, GHOSTED);
        A.setValue(p, iE, aE);

        auto tW = transmissivity(ghostLocalRow, ghostLocalCol - 1, GHOSTED);
        auto dW = hmean(tP, tW) * one_dx2;
        auto aW = -theta * dt / sP * dW;
        auto iW = globalIndicesHandle(ghostLocalRow, ghostLocalCol - 1, GHOSTED);
        A.setValue(p, iW, aW);

        auto tN = transmissivity(ghostLocalRow - 1, ghostLocalCol, GHOSTED);
        auto dN = hmean(tP, tN) * one_dx2;
        auto aN = -theta * dt / sP * dN;
        auto iN = globalIndicesHandle(ghostLocalRow - 1, ghostLocalCol, GHOSTED);
        A.setValue(p, iN, aN);

        auto tS = transmissivity(ghostLocalRow + 1, ghostLocalCol, GHOSTED);
        auto dS = hmean(tP, tS) * one_dx2;
        auto aS = -theta * dt / sP * dS;
        auto iS = globalIndicesHandle(ghostLocalRow + 1, ghostLocalCol, GHOSTED);
        A.setValue(p, iS, aS);

        auto dP = -(dE + dW + dN + dS);
        auto aP = 1.0 - theta * dt / sP * dP;
        A.setValue(p, p, aP);

        // fill b
        PetscScalar value = head(ghostLocalRow, ghostLocalCol, GHOSTED) +
                            (1.0 - theta) * dt / sP *
                                (dE * head(ghostLocalRow, ghostLocalCol + 1, GHOSTED) +
                                 dW * head(ghostLocalRow, ghostLocalCol - 1, GHOSTED) +
                                 dP * head(ghostLocalRow, ghostLocalCol, GHOSTED) +
                                 dN * head(ghostLocalRow - 1, ghostLocalCol, GHOSTED) +
                                 dS * head(ghostLocalRow + 1, ghostLocalCol, GHOSTED)) +
                            dt / sP * source(row, col);
        bHandle(row, col) = value;
      }
    }
  }
  A.assemble();
}

}  // namespace CUAS
