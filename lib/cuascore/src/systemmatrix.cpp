#include "systemmatrix.h"

#include "CUASConstants.h"

namespace CUAS {

inline PetscScalar hmean(PetscScalar x1, PetscScalar x2) { return 2.0 * x1 * x2 / (x1 + x2 + TINY); }

inline int m(int i, int j, int Nx) { return j * Nx + i; }

void systemmatrix(PETScMatrix &A, PETScVector &b, int const Nx, int const Ny, PETScGrid const &hydraulicStorativity,
                  PETScGrid const &hydraulicTransmissivity, PetscScalar const dx, PetscScalar const dt,
                  PetscScalar const theta, PETScGrid const &hydraulicHead, PETScGrid const &Q,
                  PETScGrid const &dirichletValues, PETScGrid const &bndMask) {
  // we take localnumOfCols, so that we can iterate over the normal boundaries
  auto numOfCols = hydraulicStorativity.getLocalNumOfCols();
  auto numOfRows = hydraulicStorativity.getLocalNumOfRows();
  auto cornerX = hydraulicStorativity.getCornerX();
  auto cornerY = hydraulicStorativity.getCornerY();
  auto cornerXGhost = hydraulicStorativity.getCornerXGhost();
  auto cornerYGhost = hydraulicStorativity.getCornerYGhost();
  // note: T and u_n are taken local, because we potentially need to access ghost-cells
  auto &mask = bndMask.getReadHandle();
  auto &values = dirichletValues.getReadHandle();
  auto &storativity = hydraulicStorativity.getReadHandle();
  auto &transmissivity = hydraulicTransmissivity.getReadHandle();  // local
  auto &head = hydraulicHead.getReadHandle();                      // local
  auto &source = Q.getReadHandle();                                // local
  PetscScalar S_P, d_N, d_S, d_W, d_E, d_P;
  PetscScalar A_N, A_S, A_W, A_E, A_P;
  // not sure if needed
  b.setZero();
  int p;
  int currJ;
  int currI;
  int ghostI;
  int ghostJ;

  const auto one_dx2 = 1.0 / (dx * dx);

  for (int j = 0; j < numOfCols; ++j) {
    currJ = j + cornerX;
    for (int i = 0; i < numOfRows; ++i) {
      currI = i + cornerY;
      // needed for T and u_n in order to use LocalArray rather than GlobalArray
      ghostI = i + (cornerY - cornerYGhost);
      ghostJ = j + (cornerX - cornerXGhost);
      if (mask(i, j) != (PetscScalar)COMPUTE_FLAG) {
        // is filling the Mat with setValue performant?
        p = m(currI, currJ, Nx);
        A.setValue(p, p, 1);
        b.setValue(p, values(i, j));
      } else {
        S_P = storativity(i, j);
        d_N = hmean(transmissivity(ghostI, ghostJ, GHOSTED), transmissivity(ghostI, ghostJ + 1, GHOSTED)) * one_dx2;
        d_S = hmean(transmissivity(ghostI, ghostJ, GHOSTED), transmissivity(ghostI, ghostJ - 1, GHOSTED)) * one_dx2;
        d_W = hmean(transmissivity(ghostI, ghostJ, GHOSTED), transmissivity(ghostI - 1, ghostJ, GHOSTED)) * one_dx2;
        d_E = hmean(transmissivity(ghostI, ghostJ, GHOSTED), transmissivity(ghostI + 1, ghostJ, GHOSTED)) * one_dx2;
        d_P = -(d_N + d_S + d_W + d_E);
        A_N = -theta * dt / S_P * d_N;
        A_S = -theta * dt / S_P * d_S;
        A_W = -theta * dt / S_P * d_W;
        A_E = -theta * dt / S_P * d_E;
        A_P = 1 - theta * dt / S_P * d_P;
        p = m(currI, currJ, Nx);
        // fill A matrix
        // might be beneficial to set those values together using MatSetValues,
        // which is faster.
        // A[p][m(i,j-1)] = A_S
        A.setValue(p, m(currI, currJ - 1, Nx), A_S);
        // A[p, m(i-1,j)] = A_W
        A.setValue(p, m(currI - 1, currJ, Nx), A_W);
        // A[p, p] = A_P
        A.setValue(p, p, A_P);
        // A[p, m(i+1,j)] = A_E
        A.setValue(p, m(currI + 1, currJ, Nx), A_E);
        // A[p, m(i,j+1)] = A_N
        A.setValue(p, m(currI, currJ + 1, Nx), A_N);
        // fill b
        PetscScalar value = head(ghostI, ghostJ, GHOSTED) +
                            (1 - theta) * dt / S_P *
                                (d_N * head(ghostI, ghostJ + 1, GHOSTED) + d_S * head(ghostI, ghostJ - 1, GHOSTED) +
                                 d_P * head(ghostI, ghostJ, GHOSTED) + d_W * head(ghostI - 1, ghostJ, GHOSTED) +
                                 d_E * head(ghostI + 1, ghostJ, GHOSTED)) +
                            dt / S_P * source(i, j);
        b.setValue(p, value);
      }
    }
  }
  A.assemble();
  b.assemble();
}

}  // namespace CUAS
