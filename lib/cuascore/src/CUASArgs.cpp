/**
 * File: CUASArgs.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASArgs.h"

#include "utilities.h"

#include "Logger.h"

#include "cxxopts.hpp"

#include <iostream>

namespace CUAS {

inline void defineArgs(cxxopts::Options &options);

inline void handleHelpAndVersion(cxxopts::Options const &options, cxxopts::ParseResult const &result);

inline void parseCUASArgs(CUASArgs &args, cxxopts::ParseResult const &result);

inline void evaluateDoChannels(CUASArgs &args, cxxopts::ParseResult const &result);

void parseArgs(int argc, char **argv, CUASArgs &args) {
  cxxopts::Options options("CUAS", "MPI parallel version of CUAS");

  defineArgs(options);

  auto result = options.parse(argc, argv);

  handleHelpAndVersion(options, result);

  parseCUASArgs(args, result);

  // check if all channels or any channels should be applied
  evaluateDoChannels(args, result);

  if (args.verbose) {
    CUAS_INFO_RANK0("CUASArgs.cpp: parseArgs:\n\tinput: {}\n\toutput: {}.", result["input"].as<std::string>(),
                    result["output"].as<std::string>())
  }
}

inline void defineArgs(cxxopts::Options &options) {
  // clang-format off
  options
    .positional_help("INPUT [OUTPUT]")
    .show_positional_help()
    .add_options()
      ("h,help",
       "Print help")
      ("v,verbose",
       "Verbose output. Disables Progressbar")
      ("verboseSolver",
       "Verbose Solver output.")
      ("directSolver",
       "Set PETSc options for MUMPS+PARMETIS.")
      ("version",
       "Show version information")
      ("input",
       "Netcdf input file.",
       cxxopts::value<std::string>()->default_value(""))
      ("output",
       "Netcdf output file.",
       cxxopts::value<std::string>()->default_value("out.nc"))
      ("outputSize",
       "Netcdf output file size. ('small', 'normal', 'large')",
        cxxopts::value<std::string>()->default_value("normal"))
      ("totaltime",
       "The total time, which is simulated (starttime + totaltime == endtime). Use either totaltime or endtime. Example: --totaltime '3 years 1 week'",
       cxxopts::value<std::string>()->default_value(""))
      ("starttime",
       "The time of the first point in time. Example: --starttime '3 years 1 week'",
       cxxopts::value<std::string>()->default_value(""))
      ("endtime",
       "The last point in time (starttime + totaltime == endtime). Use either totaltime or endtime. Example: --endtime '3 years 1 week'",
       cxxopts::value<std::string>()->default_value(""))
      ("dt",
       "Time step length. Example: --dt '12 hours' or --dt 1 day",
       cxxopts::value<std::string>()->default_value(""))
      ("timeSteppingTheta",
       "Time stepping family, e.g. theta=1 -> backward Euler, theta=0.5 -> Crank-Nicolson (0 <= theta <= 1)",
       cxxopts::value<PetscScalar>()->default_value("1.0"))
      ("timeStepFile",
       "NetCDF input file to read a time step array",
       cxxopts::value<std::string>()->default_value(""))
      ("saveEvery",
       "Save to NetCDF every nth timestep (deprecated).",
       cxxopts::value<int>()->default_value("0"))
      ("saveInterval",
       "Save to NetCDF whenever the interval is finished. Example: --saveInterval '5 days'",
       cxxopts::value<std::string>()->default_value(""))
      ("conductivity",
       "Conductivity of layer.",
       cxxopts::value<PetscScalar>()->default_value("10"))
      ("doChannels",
       "Evolve channels?")
      ("selectedChannels",
       "select an individual channel configuration",
       cxxopts::value<std::string>()->default_value("noselected"))
      ("disableUnconfined",
       "Disable unconfined aquifer case.")
      ("x,Tmax",
       "Maximum T to be allowed in the evolution.",
       cxxopts::value<PetscScalar>()->default_value("20.0"))
      ("i,Tmin",
       "Minimum T to be allowed in the evolution.",
       cxxopts::value<PetscScalar>()->default_value("0.0000001"))
      ("Tinit",
       "Inital T to be used in the evolution.",
       cxxopts::value<PetscScalar>()->default_value("0.2"))
      ("flowConstant",
       "Ice Flow Constant A.",
       cxxopts::value<PetscScalar>()->default_value("5e-25"))
      ("roughnessFactor",
       "Roughness factor for opening term.",
       cxxopts::value<PetscScalar>()->default_value("1.0"))
      ("supplyMultiplier",
       "Multiplier for supply.",
       cxxopts::value<PetscScalar>()->default_value("1.0"))
      ("layerThickness",
       "Water layer thickness (m)",
       cxxopts::value<PetscScalar>()->default_value("0.1"))
      ("unconfSmooth",
       "Unconfined confined transition (m)",
       cxxopts::value<PetscScalar>()->default_value("0.0"))
      ("restart",
       "Restart from this file.",
       cxxopts::value<std::string>()->default_value(""))
      ("restartNoneZeroInitialGuess",
       "sets solution vector during restart",
       cxxopts::value<bool>()->default_value("true"))
      ("specificStorage",
       "Specific storage, Ss (unit: m-1)",
       cxxopts::value<PetscScalar>()->default_value("0.0000982977696"))
      ("specificYield",
       "Specific yield, Sy (unit: 1)",
       cxxopts::value<PetscScalar>()->default_value("0.4"))
      ("enableUDS",
       "Enable upwind difference scheme (UDS). The default is the central difference scheme (CDS).")
      ("thresholdThicknessUDS",
       "Threshold for UDS scheme (m). Common choices are zero or layer thickness.",
        cxxopts::value<PetscScalar>()->default_value("0.0"))
      ("basalVelocityIce",
       "Basal velocity of the ice (m/s)",
       cxxopts::value<PetscScalar>()->default_value("1e-6"))
      ("cavityBeta",
       "cavity opening parameter",
       cxxopts::value<PetscScalar>()->default_value("5e-4"))
      ("initialHead",
       "Initial value for head. Argument must be 'Nzero', 'Nopc', 'low', 'mid', 'high', 'topg' or a valid number",
       cxxopts::value<std::string>()->default_value("Nzero"))
      ("sizeOfForcingBuffer",
       "Time slices to use in Forcing Buffer (TimeForcing --> BufferedForcing). Value has to be -1 or >= 2 (default -1 --> load all).",
       cxxopts::value<int>()->default_value("-1"))
      ("loopForcing",
       "Loop the forcing when total time is longer than forcing. Otherwise the last step of the forcing is used.")
      ("forcingFile",
       "Defines an input file (NetCDF) to create a forcing. "
       "A list of multiple files is possible (Separator ;). "
       "Each file has to provide a field named 'bmelt'. "
       "Example: --forcingFile 'filepath;filepath'",
       cxxopts::value<std::string>()->default_value(""))
      ("coordinatesFile",
       "file containing lat/lon grids to copied into the output file (NetCDF)",
       cxxopts::value<std::string>()->default_value(""))
      ("seaLevelForcing",
       "Apply sea level forcing from NetCDF scalar time series file.",
       cxxopts::value<std::string>()->default_value(""))
      ("positional",
       "Positional arguments: these are the arguments that are entered without an option.",
       cxxopts::value<std::vector<std::string>>());
  // clang-format on

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

  args.Tmax = result["Tmax"].as<PetscScalar>();
  args.Tmin = result["Tmin"].as<PetscScalar>();
  args.Tinit = result["Tinit"].as<PetscScalar>();

  // need to be parsed
  args.totaltime = result["totaltime"].as<std::string>();
  args.starttime = result["starttime"].as<std::string>();
  args.endtime = result["endtime"].as<std::string>();
  args.dt = result["dt"].as<std::string>();
  args.timeStepFile = result["timeStepFile"].as<std::string>();
  args.saveEvery = result["saveEvery"].as<int>();
  args.saveInterval = result["saveInterval"].as<std::string>();
  args.conductivity = result["conductivity"].as<PetscScalar>();
  args.disableUnconfined = result["disableUnconfined"].as<bool>();
  args.flowConstant = result["flowConstant"].as<PetscScalar>();
  args.roughnessFactor = result["roughnessFactor"].as<PetscScalar>();
  args.supplyMultiplier = result["supplyMultiplier"].as<PetscScalar>();
  args.layerThickness = result["layerThickness"].as<PetscScalar>();
  args.unconfSmooth = result["unconfSmooth"].as<PetscScalar>();
  args.restart = result["restart"].as<std::string>();
  args.restartNoneZeroInitialGuess = result["restartNoneZeroInitialGuess"].as<bool>();
  args.specificStorage = result["specificStorage"].as<PetscScalar>();
  args.specificYield = result["specificYield"].as<PetscScalar>();
  args.sizeOfForcingBuffer = result["sizeOfForcingBuffer"].as<int>();
  args.loopForcing = result["loopForcing"].as<bool>();
  args.coordinatesFile = result["coordinatesFile"].as<std::string>();
  args.forcingFile = result["forcingFile"].as<std::string>();
  args.basalVelocityIce = result["basalVelocityIce"].as<PetscScalar>();
  args.cavityBeta = result["cavityBeta"].as<PetscScalar>();
  args.initialHead = result["initialHead"].as<std::string>();
  args.seaLevelForcing = result["seaLevelForcing"].as<std::string>();
  args.verbose = result["verbose"].as<bool>();
  args.verboseSolver = result["verboseSolver"].as<bool>();
  args.directSolver = result["directSolver"].as<bool>();
  args.input = result["input"].as<std::string>();
  args.output = result["output"].as<std::string>();
  args.outputSize = result["outputSize"].as<std::string>();  // todo: check valid keywords ('small', 'normal', 'large')
  args.timeSteppingTheta = result["timeSteppingTheta"].as<PetscScalar>();
  args.enableUDS = result["enableUDS"].as<bool>();
  args.thresholdThicknessUDS = result["thresholdThicknessUDS"].as<PetscScalar>();
}

inline void evaluateDoChannels(CUASArgs &args, cxxopts::ParseResult const &result) {
  auto doChannels = result["doChannels"].as<bool>();
  auto selectedChannels = result["selectedChannels"].as<std::string>();

  // doChannels should not be executed as neither doChannels nor selectedChannels is given

  // selectedChannels is not given, doChannels enables or disables all
  if (selectedChannels == std::string("noselected")) {
    args.doAllChannels = args.doAnyChannel = args.doMelt = args.doCreep = args.doCavity = doChannels;
  }
  // selectedChannels is given but its empty. We assume that no channels will be applied.
  else if (selectedChannels.empty()) {
    args.doAllChannels = args.doAnyChannel = false;
    args.doMelt = args.doCreep = args.doCavity = false;
  }
  // selectedChannels is given and contains some input.
  else {
    // check if melt is applied
    args.doMelt = selectedChannels.find("melt") != std::string::npos;
    // check if creep is applied
    args.doCreep = selectedChannels.find("creep") != std::string::npos;
    // check if cavity is applied
    args.doCavity = selectedChannels.find("cavity") != std::string::npos;
    // if channels were applied set doAllChannels and doAnyChannel
    args.doAllChannels = (args.doMelt && args.doCreep && args.doCavity);
    args.doAnyChannel = (args.doMelt || args.doCreep || args.doCavity);
  }
}

}  // namespace CUAS
