/**
 * File: CUASArgs.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASArgs.h"

#include "utilities.h"

#include "Logger.h"

#include <iostream>

namespace CUAS {

inline void defineArgs(CUASArgs &args, cxxopts::Options &options);

inline void handleHelpAndVersion(cxxopts::Options const &options, cxxopts::ParseResult const &result);

inline void parseCUASArgs(CUASArgs &args, cxxopts::ParseResult const &result);

inline void evaluateDoChannels(CUASArgs &args, cxxopts::ParseResult const &result);

inline void sanityChecks(CUASArgs const &args);

void parseArgs(int argc, char **argv, CUASArgs &args) {
  cxxopts::Options options("CUAS", "MPI parallel version of CUAS");

  defineArgs(args, options);

  auto result = options.parse(argc, argv);

  handleHelpAndVersion(options, result);

  parseCUASArgs(args, result);

  // check if all channels or any channels should be applied
  evaluateDoChannels(args, result);

  sanityChecks(args);

  if (args.verbose) {
    CUAS_INFO_RANK0("CUASArgs.cpp: parseArgs:\n\tinput: {}\n\toutput: {}.", result["input"].as<std::string>(),
                    result["output"].as<std::string>())
  }
}

inline void defineArgs(CUASArgs &args, cxxopts::Options &options) {
  options.positional_help("INPUT [OUTPUT]").show_positional_help();

  options.add_options()("h,help", "Print help")("version", "Show version information");

  for (auto &cuasoption : args.cuasOptions) {
    cuasoption->init(options);
  }

  options.add_options()("positional", "Do not use.", cxxopts::value<std::vector<std::string>>());

  options.parse_positional({"input", "output", "positional"});
}

inline void handleHelpAndVersion(cxxopts::Options const &options, cxxopts::ParseResult const &result) {
  if (result.count("help")) {
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

  if (result.count("version")) {
    std::cout << version() << std::endl;
    exit(0);
  }
}

inline void parseCUASArgs(CUASArgs &args, cxxopts::ParseResult const &result) {
  // by using explicit positional arguments for input and output this should not happen
  if (result.count("positional")) {
    CUAS_ERROR("CUASArgs.cpp: parseArgs(): Only two positional arguments allowed.")
    CUAS_ERROR("CUASArgs.cpp: parseArgs(): input = <{}>", result["input"].as<std::string>(),
               result["output"].as<std::string>())
    CUAS_ERROR("CUASArgs.cpp: parseArgs(): output = <{}>", result["output"].as<std::string>())
    auto &positional = result["positional"].as<std::vector<std::string>>();
    for (int i = 0; i < result.count("positional"); ++i) {
      CUAS_ERROR("CUASArgs.cpp: parseArgs(): positional[{}] = <{}>", i, positional[i])
    }
    CUAS_ERROR("Exiting.")
    exit(1);
  }

  for (auto &cuasoption : args.cuasOptions) {
    cuasoption->parse(result);
  }
}

inline void evaluateDoChannels(CUASArgs &args, cxxopts::ParseResult const &result) {
  // selectedChannels is not given, doChannels enables or disables all
  if (args.selectedChannels == std::string("noselected")) {
    args.doAllChannels = args.doAnyChannel = args.doMelt = args.doCreep = args.doCavity = args.doChannels;
  }
  // selectedChannels is given but its empty. We assume that no channels will be applied.
  else if (args.selectedChannels.empty()) {
    args.doAllChannels = args.doAnyChannel = false;
    args.doMelt = args.doCreep = args.doCavity = false;
  }
  // selectedChannels is given and contains some input.
  else {
    // check if melt is applied
    args.doMelt = args.selectedChannels.find("melt") != std::string::npos;
    // check if creep is applied
    args.doCreep = args.selectedChannels.find("creep") != std::string::npos;
    // check if cavity is applied
    args.doCavity = args.selectedChannels.find("cavity") != std::string::npos;
    // if channels were applied set doAllChannels and doAnyChannel
    args.doAllChannels = (args.doMelt && args.doCreep && args.doCavity);
    args.doAnyChannel = (args.doMelt || args.doCreep || args.doCavity);
  }
}

void sanityChecks(CUASArgs const &args) {
  if (args.nonLinearIters < 0) {
    CUAS_ERROR("CUASArgs.cpp: parseArgs(): args.nonLinearIters = <{}> < 0. Exiting.", args.nonLinearIters);
    exit(1);
  }
}

void CUASArgs::setup() {
  std::string description;

  // verbosity
  description = "Verbose output. Disables Progressbar";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<bool>>("v,verbose", description, &verbose));
  description = "Verbose Solver output.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<bool>>("verboseSolver", description, &verboseSolver));

  // input and output
  description = "NetCDF input file.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("input", description, &input));
  description = "NetCDF output file.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("output", description, &output));
  description = "file containing lat/lon grids to copied into the output file (NetCDF)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<std::string>>("coordinatesFile", description, &coordinatesFile));
  description = "Restart from this file.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("restart", description, &restart));
  description = "sets solution vector during restart";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<bool>>("restartNoneZeroInitialGuess", description,
                                                                     &restartNoneZeroInitialGuess));

  // time stepping
  description = "The time of the first point in time. Example: --starttime '3 years 1 week'";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("starttime", description, &starttime));
  description =
      "The last point in time (starttime + totaltime == endtime). Use either totaltime or endtime. Example: --endtime "
      "'3 years 1 week'";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("endtime", description, &endtime));
  description =
      "The total time, which is simulated (starttime + totaltime == endtime). Use either totaltime or endtime. "
      "Example: --totaltime '3 years 1 week'";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("totaltime", description, &totaltime));
  description = "Time step length. Example: --dt '12 hours' or --dt 1 day";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("dt", description, &dt));
  description = "NetCDF input file to read a time step array";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<std::string>>("timeStepFile", description, &timeStepFile));

  // output behavior
  description = "Save to NetCDF every nth timestep (deprecated).";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<int>>("saveEvery", description, &saveEvery));
  description = "Save to NetCDF whenever the interval is finished. Example: --saveInterval '5 days'";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<std::string>>("saveInterval", description, &saveInterval));
  description = "Netcdf output file size. ('small', 'normal', 'large', 'xlarge')";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("outputSize", description, &outputSize));

  // forcing
  description =
      "Defines an input file (NetCDF) to create a forcing. A list of multiple files is possible (Separator ;). Each "
      "file has to provide a field named 'bmelt'. Example: --forcingFile 'filepath;filepath'";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("forcingFile", description, &forcingFile));
  description =
      "Time slices to use in Forcing Buffer (TimeForcing --> BufferedForcing). Value has to be -1 or >= 2 (default -1 "
      "--> load all).";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<int>>("sizeOfForcingBuffer", description, &sizeOfForcingBuffer));
  description =
      "Loop the forcing when total time is longer than forcing. Otherwise the last step of the forcing is used.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<bool>>("loopForcing", description, &loopForcing));
  description = "Apply sea level forcing from NetCDF scalar time series file.";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<std::string>>("seaLevelForcing", description, &seaLevelForcing));

  // solver behavior
  description = "Set PETSc options for MUMPS+PARMETIS.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<bool>>("directSolver", description, &directSolver));
  description = "Number of non-linear sub-iterations.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<int>>("nonLinearIters", description, &nonLinearIters));
  description = "Time stepping family, e.g. theta=1 -> backward Euler, theta=0.5 -> Crank-Nicolson (0 <= theta <= 1)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("timeSteppingTheta", description, &timeSteppingTheta));
  description = "Enable upwind difference scheme (UDS). The default is the central difference scheme (CDS).";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<bool>>("enableUDS", description, &enableUDS));
  description = "disable maintain non-negativity (psi >= 0)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<bool>>("disableNonNegative", description, &disableNonNegative));

  // channel configuration
  description = "Evolve channels?";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<bool>>("doChannels", description, &doChannels));
  description = "select an individual channel configuration";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<std::string>>("selectedChannels", description, &selectedChannels));

  // physics
  description =
      "Initial value for head. Argument must be 'Nzero', 'Nopc', 'low', 'mid', 'high', 'topg' or a valid number";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<std::string>>("initialHead", description, &initialHead));
  description = "Maximum T to be allowed in the evolution.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<PetscScalar>>("x,Tmax", description, &Tmax));
  description = "Minimum T to be allowed in the evolution.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<PetscScalar>>("i,Tmin", description, &Tmin));
  description = "Initial T to be used in the evolution.";
  cuasOptions.emplace_back(std::make_unique<CUASOptionGeneric<PetscScalar>>("Tinit", description, &Tinit));
  description = "Disable unconfined aquifer case.";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<bool>>("disableUnconfined", description, &disableUnconfined));
  description = "Conductivity of layer.";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("conductivity", description, &conductivity));
  description = "Ice Flow Constant A.";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("flowConstant", description, &flowConstant, "5e-25"));
  description = "Roughness factor for opening term.";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("roughnessFactor", description, &roughnessFactor));
  description = "Multiplier for supply.";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("supplyMultiplier", description, &supplyMultiplier));
  description = "Water layer thickness (m)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("layerThickness", description, &layerThickness));
  description = "Unconfined confined transition (m)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("unconfSmooth", description, &unconfSmooth));
  description = "Specific storage, Ss (unit: m-1)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("specificStorage", description, &specificStorage));
  description = "Specific yield, Sy (unit: 1)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("specificYield", description, &specificYield));
  description = "Threshold for UDS scheme (m). Common choices are zero or layer thickness.";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("thresholdThicknessUDS", description, &thresholdThicknessUDS));
  description = "Basal velocity of the ice (m/s)";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("basalVelocityIce", description, &basalVelocityIce, "1e-6"));
  description = "cavity opening parameter";
  cuasOptions.emplace_back(
      std::make_unique<CUASOptionGeneric<PetscScalar>>("cavityBeta", description, &cavityBeta, "5e-4"));
}

}  // namespace CUAS
