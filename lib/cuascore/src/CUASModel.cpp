/**
 * File: CUASModel.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASModel.h"

#include "CUASConstants.h"
#include "Forcing/SteadyForcing.h"

#include "Logger.h"

namespace CUAS {

CUASModel::CUASModel(int numOfCols, int numOfRows) : Ncols(numOfCols), Nrows(numOfRows) {
  xAxis.resize(Ncols);
  yAxis.resize(Nrows);
  dx = 0.0;
  dy = 0.0;

  topg = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  thk = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  bndMask = std::make_unique<PETScGrid>(numOfCols, numOfRows);
  PETScGrid zeroGrid(numOfCols, numOfRows);
  zeroGrid.setZero();
  waterSource = std::make_unique<SteadyForcing>(zeroGrid);
  pIce = std::make_unique<PETScGrid>(numOfCols, numOfRows);
}

// call init before using Model but after thk has been set
void CUASModel::init() {
  dx = getAxisSpacing(xAxis, "xAxis");
  dy = getAxisSpacing(yAxis, "yAxis");
  auto err = fabs(dx - dy);
  if (err > dx * RELATIVE_GRID_SPACING_ERROR || err > dy * RELATIVE_GRID_SPACING_ERROR) {
    CUAS_ERROR("CUASModel.cpp: init(): dx = {} and dy = {} are not equal. Exiting.", dx, dy)
    exit(1);
  }

  // pIce = thk * RHO_ICE * GRAVITY (python)
  // ensure that
  auto &thk2d = thk->getReadHandle();
  auto pIce2d = pIce->getWriteHandleGhost();
  for (int j = 0; j < thk->getLocalGhostNumOfRows(); ++j) {
    for (int i = 0; i < thk->getLocalGhostNumOfCols(); ++i) {
      pIce2d(j, i) = thk2d(j, i, GHOSTED) * RHO_ICE * GRAVITY;
    }
  }

  // outer (ghost) boundary is no-flow (Neumann BC)
  bndMask->setGhostBoundary((PetscScalar)NOFLOW_FLAG);

  // cuas grid nodes are not allowed to be active (compute) at the margin
  // This could be removed later, if cuas-python is obsolete, or we only
  // require that the outer (ghost) nodes are no-flow in cuas-mpi.
  bndMask->findAndReplaceRealBoundary((PetscScalar)COMPUTE_FLAG, (PetscScalar)NOFLOW_FLAG);
}

PETScGrid const &CUASModel::getCurrentWaterSource(timeSecs currTime) { return waterSource->getCurrentQ(currTime); }

void CUASModel::setWaterSource(std::unique_ptr<Forcing> waterSource) { this->waterSource = std::move(waterSource); }

PetscScalar CUASModel::getAxisSpacing(std::vector<PetscScalar> const &axis, const std::string &axisName) {
  auto n = axis.size();
  if (n < 3) {
    CUAS_ERROR("CUASModel.cpp: getAxisSpacing(): The size n = {} of {} is invalid. Exiting.", n, axisName)
  }

  auto spacing = (axis.back() - axis.front()) / (PetscScalar)(n - 1);
  if (spacing <= 0.0) {
    CUAS_ERROR("CUASModel.cpp: getAxisSpacing(): Invalid spacing {} for {}. Exiting.", spacing, axisName)
  }

  for (int i = 0; i < n - 1; ++i) {
    auto err = fabs(axis[i + 1] - axis[i] - spacing);
    if (err > spacing * RELATIVE_GRID_SPACING_ERROR) {
      CUAS_ERROR(
          "CUASModel.cpp: getAxisSpacing(): The values of {} are not evenly spaced (err={}, index={i}). Exiting.",
          axisName, err, i)
      exit(1);
    }
  }

  return spacing;
}

}  // namespace CUAS
