#include "fill_matrix_coo.h"

#include "CUASConstants.h"

namespace CUAS {

inline PetscScalar hmean(PetscScalar x1, PetscScalar x2) { return 2.0 * x1 * x2 / (x1 + x2 + TINY); }

inline int m(int i, int j, int Nx) { return j * Nx + i; }

void fill_matrix_coo(PETScMat &A, PETScVec &b, int const Nx, int const Ny, PETScGrid const &S, PETScGrid const &T,
                     PetscScalar const dx, PetscScalar const dt, PetscScalar const theta, PETScGrid const &u_n,
                     PETScGrid const &Q, PETScGrid const &dirichlet_values, PETScGrid const &dirichlet_mask) {
  // we take localnumOfCols, so that we can iterate over the normal boundaries
  auto numOfCols = S.getLocalNumOfCols();
  auto numOfRows = S.getLocalNumOfRows();
  auto cornerX = S.getCornerX();
  auto cornerY = S.getCornerY();
  auto cornerXGhost = S.getCornerXGhost();
  auto cornerYGhost = S.getCornerYGhost();
  // note: T and u_n are taken local, because we potentially need to access ghost-cells
  auto &diri_mask2d = dirichlet_mask.getReadHandle();
  auto &diri_values2d = dirichlet_values.getReadHandle();
  auto &S_2d = S.getReadHandle();
  auto &T_2d = T.getReadHandle();     // local
  auto &u_n2d = u_n.getReadHandle();  // local
  auto &Q_2d = Q.getReadHandle();
  PetscScalar S_P, d_N, d_S, d_W, d_E, d_P;
  PetscScalar A_N, A_S, A_W, A_E, A_P;
  // not sure if needed
  b.setZero();
  int p;
  int currJ;
  int currI;
  int GhostI;
  int GhostJ;
  for (int j = 0; j < numOfCols; ++j) {
    currJ = j + cornerX;
    for (int i = 0; i < numOfRows; ++i) {
      currI = i + cornerY;
      // needed for T and u_n in order to use LocalArray rather than GlobalArray
      GhostI = i + (cornerY - cornerYGhost);
      GhostJ = j + (cornerX - cornerXGhost);
      if (diri_mask2d(i, j)) {
        // is filling the Mat with setValue performant?
        p = m(currI, currJ, Nx);
        A.setValue(p, p, 1);
        b.setValue(p, diri_values2d(i, j));
      } else {
        S_P = S_2d(i, j);
        d_N = hmean(T_2d(GhostI, GhostJ, GHOSTED), T_2d(GhostI, GhostJ + 1, GHOSTED)) / (dx * dx);
        d_S = hmean(T_2d(GhostI, GhostJ, GHOSTED), T_2d(GhostI, GhostJ - 1, GHOSTED)) / (dx * dx);
        d_W = hmean(T_2d(GhostI, GhostJ, GHOSTED), T_2d(GhostI - 1, GhostJ, GHOSTED)) / (dx * dx);
        d_E = hmean(T_2d(GhostI, GhostJ, GHOSTED), T_2d(GhostI + 1, GhostJ, GHOSTED)) / (dx * dx);
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
        PetscScalar value = u_n2d(GhostI, GhostJ, GHOSTED) +
                            (1 - theta) * dt / S_P *
                                (d_N * u_n2d(GhostI, GhostJ + 1, GHOSTED) + d_S * u_n2d(GhostI, GhostJ - 1, true) +
                                 d_P * u_n2d(GhostI, GhostJ, GHOSTED) + d_W * u_n2d(GhostI - 1, GhostJ, true) +
                                 d_E * u_n2d(GhostI + 1, GhostJ, GHOSTED)) +
                            dt / S_P * Q_2d(i, j);
        b.setValue(p, value);
      }
    }
  }
  A.assemble();
  b.assemble();
}

}  // namespace CUAS
