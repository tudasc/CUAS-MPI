#include "SolutionHandler.h"

namespace CUAS {

SolutionHandler::SolutionHandler(std::string const &fileName, const int Nt, const int saveEvery,
                                 std::string const &inputFileName)
    : Nt(Nt) {
  // read from input
  int inputFileId;
  if (int retval = nc_open_par(inputFileName.c_str(), NC_NOWRITE, PETSC_COMM_WORLD, MPI_INFO_NULL, &inputFileId)) {
    std::string netcdfError = nc_strerror(retval);
    Logger::instance().error("CUASFile.cpp: CUASFile() in read-mode: A netcdf error occured: " + netcdfError +
                             "Exiting.");
    exit(1);
  }
  // get dimId of x
  int dimIdX;
  nc_inq_dimid(inputFileId, "x", &dimIdX);
  // get lenght of dimension x
  size_t dimXAsSize;
  nc_inq_dimlen(inputFileId, dimIdX, &dimXAsSize);
  // get dimId of y
  int dimIdY;
  nc_inq_dimid(inputFileId, "y", &dimIdY);
  // get lenght of dimension y
  size_t dimYAsSize;
  nc_inq_dimlen(inputFileId, dimIdY, &dimYAsSize);
  nc_close(inputFileId);
  int dimX = dimXAsSize;
  int dimY = dimYAsSize;
  file = std::make_unique<CUASFile>(fileName, dimX, dimY);
  SolutionHandler::defineSolution(saveEvery);
}

SolutionHandler::SolutionHandler(std::string const &fileName, const int Nt, const int saveEvery, int dimX, int dimY)
    : Nt(Nt) {
  file = std::make_unique<CUASFile>(fileName, dimX, dimY);
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
  grids["cavity_opening"];

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
                                   PETScGrid const &cavityOpening) {
  grids["u"] = &u;
  grids["u_n"] = &u_n;

  grids["usurf"] = model.usurf.get();
  grids["topg"] = model.topg.get();
  grids["thk"] = model.thk.get();
  grids["bndMask"] = model.bndMask.get();
  grids["pIce"] = model.pIce.get();

  grids["melt"] = &melt;
  grids["cavity_opening"] = &cavityOpening;

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
}

SolutionHandler::~SolutionHandler() { file->writeTimeSteps(timeSteps); }
};  // namespace CUAS
