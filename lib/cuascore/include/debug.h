#ifndef CUAS_DEBUG_H
#define CUAS_DEBUG_H

#include "NetCDFFile.h"

#include "PETScGrid.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace CUAS {
struct namedGrid {
  std::string name;
  PETScGrid const *grid;
};

inline void dumpGridsToNetCDF(std::string function_name, std::string postfix, std::vector<namedGrid> const &grids) {
  for (auto &g : grids) {
    if (!g.grid->isCompatible(*grids[0].grid)) {
      CUAS_ERROR("{}:{}:{} was called with incompatible PETScGrids. Exiting.", __FILE__, __LINE__, __func__);
      exit(1);
    }
  }

  static std::unordered_map<std::string, int> numbering;
  int mpiSize;
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  std::string filename = function_name + "_ncpus" + std::to_string(mpiSize) + "_" + postfix;
  int number = numbering[filename];
  numbering[filename]++;
  filename += std::to_string(number);
  filename += ".nc";

  NetCDFFile file(filename, grids[0].grid->getTotalNumOfCols(), grids[0].grid->getTotalNumOfRows());
  for (auto &g : grids) {
    file.defineGrid(g.name);
  }
  for (auto &g : grids) {
    file.write(g.name, *g.grid);
  }
}
}  // namespace CUAS

#endif
