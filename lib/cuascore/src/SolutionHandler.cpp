/**
 * File: SolutionHandler.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "SolutionHandler.h"

#include "CUASConstants.h"
#include "CUASKernels.h"
#include "CUASModel.h"
#include "CUASSolver.h"
#include "CUASTimeIntegrator.h"
#include "utilities.h"  // for version

#include "petscoptionswrapper.h"

namespace CUAS {

SolutionHandler::SolutionHandler(std::string const &fileName, int dimX, int dimY, std::string const &outputSize,
                                 bool storeMutable)
    : storeMutable(storeMutable) {
  file = std::make_unique<NetCDFFile>(fileName, dimX, dimY);

  if (outputSize == "small") {
    osize = OutputSize::SMALL;
  } else if (outputSize == "normal") {
    osize = OutputSize::NORMAL;
  } else if (outputSize == "large") {
    osize = OutputSize::LARGE;
  } else if (outputSize == "xlarge" || outputSize == "x-large") {
    osize = OutputSize::XLARGE;
  } else {
    CUAS_ERROR("SolutionHandler.cpp: SolutionHandler(...): unknown output size keyword : " + outputSize + "Exiting.")
    exit(1);
  }

  SolutionHandler::defineSolution();
}

void SolutionHandler::defineSolution() {
  /*
   * grid description
   */
  file->defineScalar("time", UNLIMITED);
  file->addAttributeToVariable("time", "units", "seconds since 01-01-01 00:00:00");
  file->addAttributeToVariable("time", "standard_name", "time");
  file->addAttributeToVariable("time", "calendar", "365_day");
  file->addAttributeToVariable("time", "axis", "T");

  file->defineVectorX("x");
  file->addAttributeToVariable("x", "units", "m");
  file->addAttributeToVariable("x", "standard_name", "projection_x_coordinate");
  file->addAttributeToVariable("x", "long_name", "X-coordinate in Cartesian system");
  file->addAttributeToVariable("x", "axis", "X");

  file->defineVectorY("y");
  file->addAttributeToVariable("y", "units", "m");
  file->addAttributeToVariable("y", "standard_name", "projection_y_coordinate");
  file->addAttributeToVariable("y", "long_name", "Y-coordinate in Cartesian system");
  file->addAttributeToVariable("y", "axis", "Y");

  // in the python version the noflow_mask was stored, here we store the bnd_mask
  file->defineGrid("bnd_mask", storeMutable);
  file->addAttributeToVariable("bnd_mask", "units", "1");
  file->addAttributeToVariable("bnd_mask", "flag_meanings",
                               "COMPUTE_FLAG DIRICHLET_FLAG NOFLOW_FLAG DIRICHLET_OCEAN_FLAG DIRICHLET_LAKE_FLAG");
  file->addAttributeToVariable("bnd_mask", "standard_name", "bnd_mask");
  // TODO use c++20 std::format, we only allow for c++17 in CMakeLists.txt
  // clang-format off
  std::string s = "bnd_mask (" +
                  std::to_string(COMPUTE_FLAG) + " = cuas active, " +
                  std::to_string(DIRICHLET_FLAG) + " = domain boundary (Dirichlet), " +
                  std::to_string(NOFLOW_FLAG) + " = no-flow (Neumann), " +
                  std::to_string(DIRICHLET_OCEAN_FLAG) + " = ocean (Dirichlet), " +
                  std::to_string(DIRICHLET_LAKE_FLAG) + " = periglacial lake or river (Dirichlet))";
  // clang-format on
  file->addAttributeToVariable("bnd_mask", "long_name", s);

  /*
   * state variables
   */
  file->defineGrid("head", UNLIMITED);
  file->addAttributeToVariable("head", "units", "m");
  file->addAttributeToVariable("head", "standard_name", "hydraulic_head");
  file->addAttributeToVariable("head", "long_name", "hydraulic head");

  file->defineGrid("transmissivity", UNLIMITED);
  file->addAttributeToVariable("transmissivity", "units", "m2 s-1");
  file->addAttributeToVariable("transmissivity", "standard_name", "hydraulic_transmissivity");
  file->addAttributeToVariable("transmissivity", "long_name", "hydraulic transmissivity");

  /*
   * solver and convergence
   */
  file->defineScalar("eps_inf", UNLIMITED);
  file->addAttributeToVariable("eps_inf", "units", "m s-1");
  file->addAttributeToVariable("eps_inf", "standard_name", "head_change_rate");
  file->addAttributeToVariable("eps_inf", "long_name",
                               "rate of change in head inf-norm: eps_inf = max(|h^n-h^(n-1)|)/dt");

  file->defineScalar("Teps_inf", UNLIMITED);
  file->addAttributeToVariable("Teps_inf", "units", "m2 s-2");
  file->addAttributeToVariable("Teps_inf", "standard_name", "transmissivity_change_rate");
  file->addAttributeToVariable("Teps_inf", "long_name",
                               "rate of change in transmissivity inf-norm: eps_inf = max(|T^n-T^(n-1)|)/dt");

  if (osize >= OutputSize::NORMAL) {
    file->defineGrid("topg", storeMutable);  // sometimes called bedrock
    file->addAttributeToVariable("topg", "units", "m");
    file->addAttributeToVariable("topg", "standard_name", "land_ice_bed_elevation");
    file->addAttributeToVariable("topg", "long_name", "land_ice_bed_elevation");

    file->defineGrid("watersource", UNLIMITED);
    file->addAttributeToVariable("watersource", "units", "m s-1");
    file->addAttributeToVariable("watersource", "standard_name", "water_source");
    file->addAttributeToVariable("watersource", "long_name", "Water input into model");

    /*
     * channel evolution
     */
    file->defineGrid("a_melt", UNLIMITED);
    file->addAttributeToVariable("a_melt", "units", "m2 s-2");
    file->addAttributeToVariable("a_melt", "standard_name", "a_melt");
    file->addAttributeToVariable("a_melt", "long_name", "dT/dt due to channel wall melt");
    file->addAttributeToVariable("a_melt", "doc", "dT/dt = a_melt + a_creep + a_cavity");

    file->defineGrid("a_creep", UNLIMITED);
    file->addAttributeToVariable("a_creep", "units", "m2 s-2");
    file->addAttributeToVariable("a_creep", "standard_name", "a_creep");
    file->addAttributeToVariable("a_creep", "long_name", "dT/dt due to creep opening");
    file->addAttributeToVariable("a_creep", "doc", "dT/dt = a_melt + a_creep + a_cavity");

    file->defineGrid("a_cavity", UNLIMITED);
    file->addAttributeToVariable("a_cavity", "units", "m2 s-2");
    file->addAttributeToVariable("a_cavity", "standard_name", "a_creep");
    file->addAttributeToVariable("a_cavity", "long_name", "dT/dt due to cavity opening");
    file->addAttributeToVariable("a_cavity", "doc", "dT/dt = a_melt + a_creep + a_cavity");

    /*
     * diagnostics
     */
    file->defineGrid("peffective", UNLIMITED);  // could be derived from head and geometry in post-processing
    file->addAttributeToVariable("peffective", "units", "Pa");
    file->addAttributeToVariable("peffective", "standard_name", "effective_pressure");
    file->addAttributeToVariable("peffective", "long_name", "effective pressure");
    file->addAttributeToVariable("peffective", "doc", "effective pressure at the aquifer surface");

    file->defineGrid("flux", UNLIMITED);  // could be derived from head and geometry in post-processing
    file->addAttributeToVariable("flux", "units", "m2 s-1");
    file->addAttributeToVariable("flux", "standard_name", "flux");
    file->addAttributeToVariable("flux", "long_name", "flux");
    file->addAttributeToVariable("flux", "doc", "flux = T * |grad(head)|");  // fixme: use Texp
  }

  if (osize >= OutputSize::LARGE) {
    file->defineGrid("thk", storeMutable);
    file->addAttributeToVariable("thk", "units", "m");
    file->addAttributeToVariable("thk", "standard_name", "land_ice_thickness");
    file->addAttributeToVariable("thk", "long_name", "land ice thickness");

    file->defineGrid("effective_transmissivity", UNLIMITED);
    file->addAttributeToVariable("effective_transmissivity", "units", "m2 s-1");
    file->addAttributeToVariable("effective_transmissivity", "standard_name", "effective_transmissivity");
    file->addAttributeToVariable("effective_transmissivity", "long_name", "Effective layer transmissivity");

    file->defineGrid("effective_storativity", UNLIMITED);
    file->addAttributeToVariable("effective_storativity", "units", "1");
    file->addAttributeToVariable("effective_storativity", "standard_name", "effective_storativity");
    file->addAttributeToVariable("effective_storativity", "long_name", "Effective layer storativity");

    file->defineGrid("fluxXDir", UNLIMITED);  // could be derived from head and geometry in post-processing
    file->addAttributeToVariable("fluxXDir", "units", "m2 s-1");
    file->addAttributeToVariable("fluxXDir", "standard_name", "flux_x_dir");
    file->addAttributeToVariable("fluxXDir", "long_name", "x-component of the flux");

    file->defineGrid("fluxYDir", UNLIMITED);  // could be derived from head and geometry in post-processing
    file->addAttributeToVariable("fluxYDir", "units", "m2 s-1");
    file->addAttributeToVariable("fluxYDir", "standard_name", "flux_y_dir");
    file->addAttributeToVariable("fluxYDir", "long_name", "y-component of the flux");
  }

  if (osize >= OutputSize::XLARGE) {
    file->defineGrid("pice", storeMutable);  // direct input from ice sheet model or computed from geometry
    file->addAttributeToVariable("pice", "units", "Pa");
    file->addAttributeToVariable("pice", "standard_name", "land_ice_pressure");
    file->addAttributeToVariable("pice", "long_name", "land ice pressure");
  }
}

void SolutionHandler::storeData(CUASSolver const &solver, CUASModel const &model, CUASArgs const &args,
                                PETScGrid const &waterSource, CUASTimeIntegrator const &timeIntegrator) {
  //
  // STORE DATA, IF NEEDED
  //
  OutputReason reason = getOutputReason(timeIntegrator);

  if (reason != OutputReason::NONE) {
    if (!args.verboseSolver && args.verbose) {
      // show only if verboseSolver is off
      CUAS_INFO_RANK0("time({}/{}) = {} s, dt = {} s", timeIntegrator.getTimestepIndex(),
                      timeIntegrator.getTimesteps().size() - 1, timeIntegrator.getCurrentTime(),
                      timeIntegrator.getCurrentDt())
    }

    // enhancement: Process diagnostic variables that are for output only here
    // instead of in CUASSolver. We don't need them for every time step

    if (reason == OutputReason::INITIAL) {
      // storeInitialSetup() calls storeSolution() to store initial values for time dependent fields
      storeInitialSetup(solver, model, waterSource, args, timeIntegrator);
    } else {
      if (storeMutable) {
        storeMutableModelInformation(model);
      }
      storeSolution(timeIntegrator.getCurrentTime(), solver, waterSource, solver.eps, solver.Teps);
    }

    finalizeSolution();
  }
}

void SolutionHandler::storeCUASArgs(CUASArgs const &args) {
  // TODO: The code below will break soon or later if somebody changes CUASArgs
  //       This should be done in an automated way, e.g. by using 'reflection'.
  file->addGlobalAttribute("Tmax", args.Tmax);
  file->addGlobalAttribute("Tmin", args.Tmin);
  file->addGlobalAttribute("Tinit", args.Tinit);
  file->addGlobalAttribute("totaltime", args.totaltime);
  file->addGlobalAttribute("starttime", args.starttime);
  file->addGlobalAttribute("endtime", args.endtime);
  file->addGlobalAttribute("dt", args.dt);
  file->addGlobalAttribute("timeSteppingTheta", args.timeSteppingTheta);
  file->addGlobalAttribute("timeStepFile", args.timeStepFile);
  file->addGlobalAttribute("saveEvery", args.saveEvery);
  file->addGlobalAttribute("conductivity", args.conductivity);
  file->addGlobalAttribute("doAllChannels", args.doAllChannels);
  file->addGlobalAttribute("doAnyChannels", args.doAnyChannel);
  file->addGlobalAttribute("doCavity", args.doCavity);
  file->addGlobalAttribute("doMelt", args.doMelt);
  file->addGlobalAttribute("doCreep", args.doCreep);
  file->addGlobalAttribute("disableUnconfined", args.disableUnconfined);
  file->addGlobalAttribute("flowConstant", args.flowConstant);
  file->addGlobalAttribute("roughnessFactor", args.roughnessFactor);
  file->addGlobalAttribute("supplyMultiplier", args.supplyMultiplier);
  file->addGlobalAttribute("layerThickness", args.layerThickness);
  file->addGlobalAttribute("unconfSmooth", args.unconfSmooth);
  file->addGlobalAttribute("restart", args.restart);
  file->addGlobalAttribute("restartNoneZeroInitialGuess", args.restartNoneZeroInitialGuess);
  file->addGlobalAttribute("specificStorage", args.specificStorage);
  file->addGlobalAttribute("specificYield", args.specificYield);
  file->addGlobalAttribute("loopForcing", args.loopForcing);
  file->addGlobalAttribute("coordinatesFile", args.coordinatesFile);
  file->addGlobalAttribute("forcingFile", args.forcingFile);
  file->addGlobalAttribute("basalVelocityIce", args.basalVelocityIce);
  file->addGlobalAttribute("cavityBeta", args.cavityBeta);
  file->addGlobalAttribute("initialHead", args.initialHead);
  // ignore (std::string) tempResults,
  // ignore (bool) version,
  file->addGlobalAttribute("seaLevelForcing", args.seaLevelForcing);
  // ignore (bool) verbose and (bool) verboseSolver,
  file->addGlobalAttribute("directSolver", args.directSolver);
  file->addGlobalAttribute("input", args.input);
  file->addGlobalAttribute("output", args.output);
  file->addGlobalAttribute("outputSize", args.outputSize);

  file->addGlobalAttribute("enableUDS", args.enableUDS);
  file->addGlobalAttribute("thresholdThicknessUDS", args.thresholdThicknessUDS);

  // compile time CUASConstants
  file->addGlobalAttribute("version", version());
  file->addGlobalAttribute("TINY", TINY);
  file->addGlobalAttribute("NOFLOW_VALUE", NOFLOW_VALUE);
  file->addGlobalAttribute("RHO_ICE", RHO_ICE);
  file->addGlobalAttribute("SPY", SPY);

  // copy lat/lon lat_bnds/lon_bnds from input file if needed
  if (!args.coordinatesFile.empty()) {
    CUAS_INFO_RANK0("Copy coordinates from {}", args.coordinatesFile)
    file->copyCoordinatesFrom(args.coordinatesFile);
    file->setCoordinatesAttribute();
  }
}

void SolutionHandler::storeConstantModelInformation(const CUASModel &model) {
  file->write("x", model.xAxis);
  file->write("y", model.yAxis);
}

void SolutionHandler::storeMutableModelInformation(const CUASModel &model) {
  file->write("bnd_mask", *model.bndMask, nextSolution);

  if (osize >= OutputSize::NORMAL) {
    file->write("topg", *model.topg, nextSolution);
  }

  if (osize >= OutputSize::LARGE) {
    file->write("thk", *model.thk, nextSolution);
  }

  if (osize >= OutputSize::XLARGE) {
    file->write("pice", *model.pIce, nextSolution);
  }
}

void SolutionHandler::storeInitialSetup(CUASSolver const &solver, CUASModel const &model, PETScGrid const &waterSource,
                                        CUASArgs const &args, CUASTimeIntegrator const &timeIntegrator) {
  storeConstantModelInformation(model);

  storeMutableModelInformation(model);

  storeCUASArgs(args);

  // store initial conditions if needed
  storeSolution(timeIntegrator.getCurrentTime(), solver, waterSource);
}

void SolutionHandler::storeSolution(timeSecs currTime, CUASSolver const &solver, PETScGrid const &waterSource,
                                    PetscScalar eps_inf, PetscScalar Teps_inf) {
  // write scalars
  file->write("time", currTime, nextSolution);

  // TODO make eps_inf and Teps_inf available for storage
  file->write("eps_inf", eps_inf, nextSolution);
  file->write("Teps_inf", Teps_inf, nextSolution);

  // write grids
  file->write("head", *solver.currHead, nextSolution);
  file->write("transmissivity", *solver.currTransmissivity, nextSolution);

  if (osize >= OutputSize::NORMAL) {
    file->write("a_melt", *solver.melt, nextSolution);
    file->write("a_creep", *solver.creep, nextSolution);
    file->write("a_cavity", *solver.cavity, nextSolution);
    file->write("peffective", *solver.pEffective, nextSolution);
    file->write("watersource", waterSource, nextSolution);
    file->write("flux", *solver.fluxMagnitude, nextSolution);
  }

  if (osize >= OutputSize::LARGE) {
    file->write("effective_transmissivity", *solver.Teff, nextSolution);
    file->write("effective_storativity", *solver.Seff, nextSolution);
    file->write("fluxXDir", *solver.fluxXDir, nextSolution);
    file->write("fluxYDir", *solver.fluxYDir, nextSolution);
  }

  if (osize >= OutputSize::XLARGE) {
    // we would need #include "CUASKernels.h", args.layerThickness and headToPressure() for water pressure
    // output
    // todo: add other fields
  }
}

void SolutionHandler::finalizeSolution() {
  // Write everything to the file
  file->sync();

  // Increment time-slice index for next call
  ++nextSolution;
}

void SolutionHandler::setTimeUnits(std::string const &s) {
  // TODO: do some checks
  file->addAttributeToVariable("time", "units", s);
}

void SolutionHandler::setCalendar(std::string const &s) {
  // TODO: do some checks
  file->addAttributeToVariable("time", "calendar", s);
}

void SolutionHandler::setSaveStrategy(SaveStrategy strategy, long saveInterval, int saveEvery) {
  // TODO add sanity checks for input

  this->strategy = strategy;
  this->saveInterval = saveInterval;
  this->saveEvery = saveEvery;
}

void SolutionHandler::storePETScOptions() {
  file->addGlobalAttribute("PETSC_OPTIONS", getPETScOptionsAll());
  file->addGlobalAttribute("PETSC_OPTIONS_UNUSED", getPETScOptionsUnused());
  file->addGlobalAttribute("PETSC_OPTIONS_USED", getPETScOptionsUsed());
}

OutputReason SolutionHandler::getOutputReason(CUASTimeIntegrator const &timeIntegrator) const {
  // TODO clean up get output reason

  auto timeStepIndex = timeIntegrator.getTimestepIndex();

  if (timeStepIndex == 0) {
    return OutputReason::INITIAL;
  }

  if (strategy == SaveStrategy::INDEX) {
    if (saveEvery > 0) {  // avoid division by zero in modulo operation
      if ((timeStepIndex % saveEvery == 0) || (timeStepIndex == timeIntegrator.getTimesteps().size() - 1)) {
        return OutputReason::NORMAL;
      } else {
        return OutputReason::NONE;
      }
    } else {
      return OutputReason::NONE;
    }
  } else if (strategy == SaveStrategy::TIMEINTERVAL) {
    if (saveInterval > 0) {
      if ((timeIntegrator.getCurrentTime() % saveInterval == 0) ||
          (timeIntegrator.getCurrentTime() == timeIntegrator.getTimesteps().back())) {
        return OutputReason::NORMAL;
      } else {
        return OutputReason::NONE;
      }
    } else {
      return OutputReason::NONE;
    }
  } else if (strategy == SaveStrategy::DEFAULT) {
    return OutputReason::NORMAL;
  }

  return OutputReason::NONE;
}

}  // namespace CUAS
