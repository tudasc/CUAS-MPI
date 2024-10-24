/**
 * File: systemmatrix.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "systemmatrix.h"

#include "CUASConstants.h"

namespace CUAS {

inline PetscScalar hmean(PetscScalar x1, PetscScalar x2) { return 2.0 * x1 * x2 / (x1 + x2 + TINY); }

/*
 * We use compass notation to simplify indexing, thus "E" means to the East i+1 while "e" means
 * at the eastern interface at i+1/2 and in a similar fashion for N:North, S:South, W:West.
 * East = col + 1, West = col - 1
 * North = row - 1, South = row + 1  (old)
 * North = row + 1, South = row - 1  (new)
 *
 * P: just means the current location at (i,j) aka (row,col)
 *
 * S dh/dt + (q_e - q_w)/dx + (q_n - q_s)/dy = Q
 * S dh/dt + div_q_x + div_q_y = Q
 *
 * Flux divergence in the upwind scheme (Jarosh et al. 2013, eq. 36, MUSCL 1D)
 *
 * div_q_x = (T_up_x * (h[i+1] - h[i])/dx - T_down_x * (h[i] - h[i-1])/dx)/dx
 * ->      = 1/dx^2 * ( T_up_x * h[i+1] - T_up_x * h[i] - T_down_x * h[i] + T_down_x * h[i-1] )
 * ->      = 1/dx^2 * ( T_up_x * h[i+1] - (T_up_x + T_down_x) * h[i] + T_down_x * h[i-1])
 *
 * similar for y-direction
 * div_q_y = (T_up_y * (h[j+1] - h[j])/dy - T_dn_y * (h[j] - h[j-1])/dy)/dy
 * ->      = 1/dy^2 * ( T_up_y * h[j+1] - (T_up_y + T_dn_y) * h[j] + T_dn_y * h[j-1])
 *
 */
void systemmatrix(PETScMatrix &A, PETScGrid &b, PETScGrid const &hydraulicStorativity, PETScGrid const &effTransEast,
                  PETScGrid const &effTransWest, PETScGrid const &effTransNorth, PETScGrid const &effTransSouth,
                  PetscScalar dx, PetscScalar dt, PetscScalar theta, PETScGrid const &hydraulicHead,
                  PETScGrid const &waterSource, PETScGrid const &dirichletValues, PETScGrid const &bndMask,
                  PETScGrid const &globalIndices) {
  auto &mask = bndMask.getReadHandle();
  auto &values = dirichletValues.getReadHandle();
  auto &storativity = hydraulicStorativity.getReadHandle();

  // up: upwind direction, dn: downwind
  auto &T_up_x = effTransEast.getReadHandle();   // interface (eff.) transmissivity in x-dir
  auto &T_dn_x = effTransWest.getReadHandle();   // interface (eff.) transmissivity in x-dir
  auto &T_up_y = effTransNorth.getReadHandle();  // interface (eff.) transmissivity in y-dir
  auto &T_dn_y = effTransSouth.getReadHandle();  // interface (eff.) transmissivity in y-dir

  auto &head = hydraulicHead.getReadHandle();
  auto &source = waterSource.getReadHandle();
  // Indices are of integer type but stored as PetsScalar aka double for convenience.
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
      auto iP = static_cast<int>(globalIndicesHandle(row, col));

      if (mask(row, col) != (PetscScalar)COMPUTE_FLAG) {
        A.setValue(iP, iP, 1.0);
        bHandle(row, col) = values(row, col);
      } else {
        // needed for transmissivity and head in order to use LocalArray rather than GlobalArray
        auto j = row + (cornerY - cornerYGhost);
        auto i = col + (cornerX - cornerXGhost);

        auto sP = storativity(row, col);

        // some abbreviations
        auto head_E = head(j, i + 1, GHOSTED);
        auto head_W = head(j, i - 1, GHOSTED);
        auto head_N = head(j + 1, i, GHOSTED);
        auto head_S = head(j - 1, i, GHOSTED);
        auto head_P = head(j, i, GHOSTED);

        //
        // fill A
        //
        // auto dE = harmonicmean(T_P, T_E) * one_dx2;
        auto dE = T_up_x(row, col) * one_dx2;
        auto aE = -theta * dt / sP * dE;
        auto iE = static_cast<int>(globalIndicesHandle(j, i + 1, GHOSTED));
        A.setValue(iP, iE, aE);

        // auto dW = harmonicmean(T_P, T_W) * one_dx2;
        auto dW = T_dn_x(row, col) * one_dx2;
        auto aW = -theta * dt / sP * dW;
        auto iW = static_cast<int>(globalIndicesHandle(j, i - 1, GHOSTED));
        A.setValue(iP, iW, aW);

        // auto dN = harmonicmean(T_P, T_N) * one_dx2;
        auto dN = T_up_y(row, col) * one_dx2;
        auto aN = -theta * dt / sP * dN;
        auto iN = static_cast<int>(globalIndicesHandle(j + 1, i, GHOSTED));
        A.setValue(iP, iN, aN);

        // auto dS = harmonicmean(T_P, T_S) * one_dx2;
        auto dS = T_dn_y(row, col) * one_dx2;
        auto aS = -theta * dt / sP * dS;
        auto iS = static_cast<int>(globalIndicesHandle(j - 1, i, GHOSTED));
        A.setValue(iP, iS, aS);

        auto dP = -(dE + dW + dN + dS);  // this could be zero and thus aP = 1.0 (sources only)
        auto aP = 1.0 - theta * dt / sP * dP;
        A.setValue(iP, iP, aP);

        // fill b
        PetscScalar value =
            head_P + (1.0 - theta) * dt / sP * (dE * head_E + dW * head_W + dP * head_P + dN * head_N + dS * head_S) +
            dt / sP * source(row, col);
        bHandle(row, col) = value;
      }
    }
  }
  A.assemble();
}

void systemmatrixDeprecated(PETScMatrix &A, PETScGrid &b, PETScGrid const &hydraulicStorativity,
                            PETScGrid const &hydraulicTransmissivity, PetscScalar dx, PetscScalar dt, PetscScalar theta,
                            PETScGrid const &hydraulicHead, PETScGrid const &waterSource,
                            PETScGrid const &dirichletValues, PETScGrid const &bndMask,
                            PETScGrid const &globalIndices) {
  auto &mask = bndMask.getReadHandle();
  auto &values = dirichletValues.getReadHandle();
  auto &storativity = hydraulicStorativity.getReadHandle();
  auto &transmissivity = hydraulicTransmissivity.getReadHandle();
  auto &head = hydraulicHead.getReadHandle();
  auto &source = waterSource.getReadHandle();
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
