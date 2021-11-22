#include "SolutionHandler.h"

namespace CUAS {

SolutionHandler::SolutionHandler(std::string const &fileName, int dimX, int dimY, std::string const &outputSize) {
  file = std::make_unique<NetCDFFile>(fileName, dimX, dimY);

  if (outputSize == "small") {
    osize = OutputSize::SMALL;
  } else if (outputSize == "normal") {
    osize = OutputSize::NORMAL;
  } else if (outputSize == "large") {
    osize = OutputSize::LARGE;
  } else {
    Logger::instance().error("SolutionHandler.cpp: SolutionHandler(...): unknown output size keyword : " + outputSize +
                             "Exiting.");
    exit(1);
  }

  SolutionHandler::defineSolution();
}

void SolutionHandler::defineSolution() {
  /*
   * grid description
   */

  file->defineScalar("time", true);
  file->addAttributeToVariable("time", "unit", "seconds since 01-01-01 00:00:00");
  file->addAttributeToVariable("time", "standard_name", "time");
  file->addAttributeToVariable("time", "calendar", "365_day");
  file->addAttributeToVariable("time", "axis", "T");

  file->defineVectorX("x");
  file->addAttributeToVariable("x", "unit", "m");
  file->addAttributeToVariable("x", "standard_name", "projection_x_coordinate");
  file->addAttributeToVariable("x", "long_name", "X-coordinate in Cartesian system");
  file->addAttributeToVariable("x", "axis", "X");

  file->defineVectorY("y");
  file->addAttributeToVariable("y", "unit", "m");
  file->addAttributeToVariable("y", "standard_name", "projection_y_coordinate");
  file->addAttributeToVariable("y", "long_name", "Y-coordinate in Cartesian system");
  file->addAttributeToVariable("y", "axis", "Y");

  file->defineGrid("bnd_mask");  // in the python version the noflow_mask was stored, here we store the bnd
  file->addAttributeToVariable("bnd_mask", "unit", "1");
  file->addAttributeToVariable("bnd_mask", "flag_meanings",
                               "COMPUTE_FLAG DIRICHLET_FLAG NOFLOW_FLAG DIRICHLET_LAKE_FLAG");
  // todo add list of flag values as: bnd_mask:flag_values = 0, 1, 2, 3 ;
  file->addAttributeToVariable("bnd_mask", "standard_name", "bnd_mask");
  file->addAttributeToVariable("bnd_mask", "long_name",
                               "bnd_mask (0 = cuas, 1 = floating ice or ocean, 2 = no-flow, 3 = periglacial lake)");

  /*
   * state variables
   */
  file->defineGrid("head", UNLIMITED);
  file->addAttributeToVariable("head", "unit", "m");
  file->addAttributeToVariable("head", "standard_name", "hydraulic_head");
  file->addAttributeToVariable("head", "long_name", "hydraulic head");

  file->defineGrid("transmissivity", UNLIMITED);
  file->addAttributeToVariable("transmissivity", "unit", "m2 s-1");
  file->addAttributeToVariable("transmissivity", "standard_name", "hydraulic_transmissivity");
  file->addAttributeToVariable("transmissivity", "long_name", "hydraulic transmissivity");

  /*
   * solver and convergence
   */
  file->defineScalar("eps_inf", UNLIMITED);
  file->defineScalar("Teps_inf", UNLIMITED);

  if (osize >= OutputSize::NORMAL) {
    /*
     * channel evolution
     */
    file->defineGrid("a_melt", UNLIMITED);
    file->addAttributeToVariable("a_melt", "unit", "m2 s-2");
    file->addAttributeToVariable("a_melt", "standard_name", "a_melt");
    file->addAttributeToVariable("a_melt", "long_name", "dT/dt due to channel wall melt");
    file->addAttributeToVariable("a_melt", "doc", "dT/dt = a_melt + a_creep + a_cavity");

    file->defineGrid("a_creep", UNLIMITED);
    file->addAttributeToVariable("a_creep", "unit", "m2 s-2");
    file->addAttributeToVariable("a_creep", "standard_name", "a_creep");
    file->addAttributeToVariable("a_creep", "long_name", "dT/dt due to creep opening");
    file->addAttributeToVariable("a_creep", "doc", "dT/dt = a_melt + a_creep + a_cavity");

    file->defineGrid("a_cavity", UNLIMITED);
    file->addAttributeToVariable("a_cavity", "unit", "m2 s-2");
    file->addAttributeToVariable("a_cavity", "standard_name", "a_creep");
    file->addAttributeToVariable("a_cavity", "long_name", "dT/dt due to cavity opening");
    file->addAttributeToVariable("a_cavity", "doc", "dT/dt = a_melt + a_creep + a_cavity");

    /*
     * diagnostics
     */
    file->defineGrid("peffective", UNLIMITED);  // could be derived from head and geometry in post-processing
    file->addAttributeToVariable("peffective", "unit", "Pa");
    file->addAttributeToVariable("peffective", "standard_name", "effective_pressure");
    file->addAttributeToVariable("peffective", "long_name", "effective pressure");
    file->addAttributeToVariable("peffective", "doc", "effective pressure at the aquifer surface");

    file->defineGrid("flux", UNLIMITED);  // could be derived from head and geometry in post-processing
    file->addAttributeToVariable("flux", "unit", "m2 s-1");
    file->addAttributeToVariable("flux", "standard_name", "flux");
    file->addAttributeToVariable("flux", "long_name", "flux");
    file->addAttributeToVariable("flux", "doc", "flux = T * |grad(head)|");
  }

  if (osize >= OutputSize::LARGE) {
    /*
     * geometry
     */
    file->defineGrid("usurf");  // upper surface
    file->addAttributeToVariable("usurf", "unit", "m");
    file->addAttributeToVariable("usurf", "standard_name", "land_ice_surface_elevation");
    file->addAttributeToVariable("usurf", "long_name", "land ice surface elevation");

    file->defineGrid("topg");  // sometimes called bedrock
    file->addAttributeToVariable("topg", "unit", "m");
    file->addAttributeToVariable("topg", "standard_name", "land_ice_bed_elevation");
    file->addAttributeToVariable("topg", "long_name", "land_ice_bed_elevation");

    file->defineGrid("thk");
    file->addAttributeToVariable("thk", "unit", "m");
    file->addAttributeToVariable("thk", "standard_name", "land_ice_thickness");
    file->addAttributeToVariable("thk", "long_name", "land ice thickness");

    file->defineGrid("pice");  // direct input from ice sheet model or computed from geometry
    file->addAttributeToVariable("pice", "unit", "Pa");
    file->addAttributeToVariable("pice", "standard_name", "land_ice_pressure");
    file->addAttributeToVariable("pice", "long_name", "land ice pressure");

    /*
     * diganostics
     */
    file->defineGrid("pwater");  // direct input from ice sheet model or computed from geometry
    file->addAttributeToVariable("pwater", "unit", "Pa");
    file->addAttributeToVariable("pwater", "standard_name", "water_pressure");
    file->addAttributeToVariable("pwater", "long_name", "water pressure");
  }
}

void SolutionHandler::storeInitialSetup(PETScGrid const &hydraulicHead, PETScGrid const &hydraulicTransmissivity,
                                        CUASModel const &model, PETScGrid const &melt, PETScGrid const &creep,
                                        PETScGrid const &cavity, CUASArgs const &args) {
  file->write("x", model.xAxis);
  file->write("y", model.yAxis);

  //
  file->write("bnd_mask", *model.bndMask, nextSolution);

  if (osize >= OutputSize::LARGE) {
    file->write("usurf", *model.usurf, nextSolution);
    file->write("topg", *model.topg, nextSolution);
    file->write("thk", *model.thk, nextSolution);
    file->write("pice", *model.pIce, nextSolution);
  }

  // Todo: this should be done in an automated way using 'reflection'.
  //  The code below will break soon or later if somebody changes CUASArgs
  file->addGlobalAttribute("Tmax", args.Tmax);
  file->addGlobalAttribute("Tmin", args.Tmin);
  file->addGlobalAttribute("totaltime", args.totaltime);
  file->addGlobalAttribute("dt", args.dt);
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
  file->addGlobalAttribute("Ssmulti", args.Ssmulti);
  file->addGlobalAttribute("Sy", args.Sy);
  file->addGlobalAttribute("Texp", args.Texp);
  file->addGlobalAttribute("noSmoothMelt", args.noSmoothMelt);
  file->addGlobalAttribute("loopForcing", args.loopForcing);
  file->addGlobalAttribute("basalVelocityIce", args.basalVelocityIce);
  file->addGlobalAttribute("cavityBeta", args.cavityBeta);
  file->addGlobalAttribute("initialHead", args.initialHead);
  // ignore (std::string) tempResults,
  // ignore (bool) version,
  file->addGlobalAttribute("seaLevelForcing", args.seaLevelForcing);
  // ignore (bool) verbose,
  file->addGlobalAttribute("input", args.input);
  file->addGlobalAttribute("output", args.output);
  file->addGlobalAttribute("outputSize", args.outputSize);

  // store initial conditions if needed
  storeSolution(0, hydraulicHead, hydraulicTransmissivity, model, melt, creep, cavity);
}

void SolutionHandler::storeSolution(CUAS::timeSecs currTime, PETScGrid const &hydraulicHead,
                                    PETScGrid const &hydraulicTransmissivity, CUASModel const &model,
                                    PETScGrid const &melt, PETScGrid const &creep, PETScGrid const &cavity) {
  // write scalars
  file->write("time", currTime, nextSolution);

  // TODO make eps_inf and Teps_inf available for storage
  PetscScalar dummy(-9999.0);
  file->write("eps_inf", dummy, nextSolution);
  file->write("Teps_inf", dummy, nextSolution);

  // write grids
  file->write("head", hydraulicHead, nextSolution);
  file->write("transmissivity", hydraulicTransmissivity, nextSolution);

  PETScGrid pwater(hydraulicHead.getTotalNumOfCols(), hydraulicHead.getTotalNumOfRows());
  pwater.setZero();

  if (osize >= OutputSize::NORMAL) {
    file->write("a_melt", melt, nextSolution);
    file->write("a_creep", creep, nextSolution);
    file->write("a_cavity", cavity, nextSolution);

    // todo update effective pressure based on current hydraulicHead
    PETScGrid peff(hydraulicHead.getTotalNumOfCols(), hydraulicHead.getTotalNumOfRows());
    peff.setZero();
    file->write("peffective", peff, nextSolution);

    // todo update flux based on current hydraulicHead and hydraulicTransmissivity
    PETScGrid flux(hydraulicHead.getTotalNumOfCols(), hydraulicHead.getTotalNumOfRows());
    flux.setZero();
    file->write("flux", flux, nextSolution);
  }

  if (osize >= OutputSize::LARGE) {
    file->write("pwater", pwater, nextSolution);
  }

  ++nextSolution;
}

SolutionHandler::~SolutionHandler() {}
};  // namespace CUAS
