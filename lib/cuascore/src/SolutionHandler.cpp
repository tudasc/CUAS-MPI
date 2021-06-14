#include "SolutionHandler.h"

namespace CUAS {

SolutionHandler::SolutionHandler(std::string const &fileName, const int Nt, const int saveEvery, int dimX, int dimY,
                                 int mpiRank)
    : Nt(Nt) {
  file = std::make_unique<CUASFile>(fileName, dimX, dimY, mpiRank);
  SolutionHandler::defineSolution(saveEvery);
}

void SolutionHandler::defineSolution(int saveEvery) {
  // define grids and scalars
  grids["u"];
  grids["u_n"];

  grids["usurf"];
  grids["topg"];
  grids["thk"];
  grids["bndMask"];
  grids["pIce"];

  grids["melt"];

  scalars["cavity_opening"];

  int fileId = file->getFileId();
  nc_redef(fileId);
  for (int i = saveEvery; i < Nt + 1; i = i + saveEvery) {
    for (auto const &grid : grids) {
      std::string gridNameTimeStep = grid.first + std::to_string(i);
      file->defineGrid(gridNameTimeStep);
    }
    for (auto const &scalar : scalars) {
      std::string scalarNameTimeStep = scalar.first + std::to_string(i);
      file->defineScalar(scalarNameTimeStep);
    }
  }
  nc_enddef(fileId);
}

void SolutionHandler::saveSolution(int timeStep, CUASArgs const &args, int rank, PETScGrid const &u,
                                   PETScGrid const &u_n, CUASModel const &model, PETScGrid const &melt,
                                   PetscScalar cavity_opening) {
  grids["u"] = &u;
  grids["u_n"] = &u_n;

  grids["usurf"] = model.usurf.get();
  grids["topg"] = model.topg.get();
  grids["thk"] = model.thk.get();
  grids["bndMask"] = model.bndMask.get();
  grids["pIce"] = model.pIce.get();

  grids["melt"] = &melt;

  scalars["cavity_opening"] = cavity_opening;

  // write grids and scalars
  for (auto const &grid : grids) {
    std::string gridNameTimeStep = grid.first + std::to_string(timeStep);
    file->write(gridNameTimeStep, *grid.second);
  }
  for (auto const &scalar : scalars) {
    std::string scalarNameTimeStep = scalar.first + std::to_string(timeStep);
    file->write(scalarNameTimeStep, scalar.second);
  }

  timeSteps.push_back(timeStep);
  if (timeStep >= Nt) {
    file->writeTimeSteps(timeSteps);
  }
}
};  // namespace CUAS
