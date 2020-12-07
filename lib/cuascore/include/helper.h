#ifndef CUAS_HELPER_H
#define CUAS_HELPER_H

#include "PetscGrid.h"
#include "physicalConstants.h"

namespace CUAS {

void pressure2head(PetscGrid *result, PetscGrid *pressure, PetscGrid *bed_elevation, PetscScalar sea_level) {
  PetscScalar **result2d = result->getAsGlobal2dArr();
  PetscScalar **pressure2d = pressure->getAsGlobal2dArr();
  PetscScalar **bed_elevation2d = bed_elevation->getAsGlobal2dArr();
  for (int j = 0; j < result->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result->getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d[j][i] - sea_level;
      result2d[j][i] = pressure2d[j][i] / (RHO_WATER * GRAVITY) + effective_bed_elevation;
    }
  }
  result->setAsGlobal2dArr(result2d);
  pressure->restoreGlobal2dArr(pressure2d);
  bed_elevation->restoreGlobal2dArr(bed_elevation2d);
}

void head2pressure(PetscGrid *result, PetscGrid *head, PetscGrid *bed_elevation, PetscScalar sea_level) {
  PetscScalar **result2d = result->getAsGlobal2dArr();
  PetscScalar **head2d = head->getAsGlobal2dArr();
  PetscScalar **bed_elevation2d = bed_elevation->getAsGlobal2dArr();
  for (int j = 0; j < result->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result->getLocalNumOfCols(); ++i) {
      double effective_bed_elevation = bed_elevation2d[j][i] - sea_level;
      result2d[j][i] = RHO_WATER * GRAVITY * (head2d[j][i] - effective_bed_elevation);
    }
  }
  result->setAsGlobal2dArr(result2d);
  head->restoreGlobal2dArr(head2d);
  bed_elevation->restoreGlobal2dArr(bed_elevation2d);
}

void overburdenPressure(PetscGrid *result, PetscGrid *thk) {
  PetscScalar **result2d = result->getAsGlobal2dArr();
  PetscScalar **thk2d = thk->getAsGlobal2dArr();
  for (int j = 0; j < result->getLocalNumOfRows(); ++j) {
    for (int i = 0; i < result->getLocalNumOfCols(); ++i) {
      result2d[j][i] = thk2d[j][i] * RHO_ICE * GRAVITY;
    }
  }
  result->setAsGlobal2dArr(result2d);
  thk->restoreGlobal2dArr(thk2d);
}

}  // namespace CUAS

#endif
