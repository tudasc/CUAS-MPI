/**
 * File: WaterSource.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_WATERSOURCE_H
#define CUAS_WATERSOURCE_H

#include "timeparse.h"

#include "PETScGrid.h"

namespace CUAS {

/**
 * Interface class defining the interface of a water source in CUAS-MPI
 *    providesWaterSource
 *    getCurrentWaterSource
 */
class WaterSource {
 public:
  /**
   * Checks if a water source exists.
   * Should be checked before getCurrentWaterSource is called to prevent overhead.
   */
  [[nodiscard]] virtual bool providesWaterSource() const = 0;
  /**
   * Calculates and returns the water source at the current time stamp.
   * Check providesWaterSource before to prevent overhead of generating a zero grid.
   *
   * @param currTime current time stamp (s)
   *
   * @return cons reference to the grid representing the water source at the current time stamp
   */
  virtual PETScGrid const &getCurrentWaterSource(timeSecs currTime) = 0;
};

}  // namespace CUAS

#endif
