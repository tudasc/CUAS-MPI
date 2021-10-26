#include "CUASArgs.h"

#include "cxxopts.hpp"

namespace CUAS {

void parseArgs(int argc, char **argv, CUASArgs &args) {
  cxxopts::Options options("CUAS", "MPI parallel version of CUAS");

  // clang-format off
  options
    .add_options()
      ("h,help", "Print help")
      ("x,Tmax",
       "Maximum T to be allowed in the evolution.",
       cxxopts::value<PetscScalar>()->default_value("20.0"))
      ("i,Tmin",
       "Minimum T to be allowed in the evolution.",
       cxxopts::value<PetscScalar>()->default_value("0.0000001"))
      ("input", 
       "Netcdf input file.",
       cxxopts::value<std::string>()->default_value(""))
      ("output",
       "Netcdf output file.",
       cxxopts::value<std::string>()->default_value("out.nc"))
      ("totaltime",
       "Total time to run model. Example: --totaltime '4 weeks', --totaltime '3 years 6 months' or --totaltime "
       "'50 years'",
       cxxopts::value<std::string>()->default_value("10 years"))
      ("dt",
       "Time step length. Example: --dt '12 hours', --dt 1day", 
       cxxopts::value<std::string>()->default_value("12 hours"))
      ("saveEvery",
       "Save every nth timestep to netcdf.",
       cxxopts::value<int>()->default_value("0"))
      ("conductivity",
       "Conductivity of layer.",
       cxxopts::value<PetscScalar>()->default_value("10"))
      ("dochannels",
       "Evolve channels?")
      ("disableUnconfined",
       "Disable unconfined aquifer case.")
      ("flowConstant",
       "Ice Flow Constant A.",
       cxxopts::value<PetscScalar>()->default_value("5e-25"))
      ("roughnessFactor", "Roughness factor for opening term.",
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
      ("Ssmulti",
       "Multiplier for specific storage Ss.",
       cxxopts::value<PetscScalar>()->default_value("1.0"))
      ("Sy",
       "Specific yield for unconfined.",
       cxxopts::value<PetscScalar>()->default_value("0.4"))
      ("Texp",
       "Exponent of T.",
       cxxopts::value<PetscScalar>()->default_value("1"))
      ("noSmoothMelt",
       "Smooth melt term before computing change in T?")
      ("loopForcing",
       "Loop the forcing when total time is longer than forcing. Otherwise the last step of the forcing is used.")
      ("basalVelocityIce",
       "Basal velocity of the ice (m/s)",
       cxxopts::value<PetscScalar>()->default_value("1e-6"))
      ("cavityBeta",
       "cavity opening parameter",
       cxxopts::value<PetscScalar>()->default_value("5e-4"))
      ("initialHead",
       "Initial value for head. Nzero means that head is set such that effective pressure N is zero.",
       cxxopts::value<std::string>()->default_value("Nzero"))
      ("tempResults",
       "Save temporary results to this file(s) to later restart from them.",
       cxxopts::value<std::string>()->default_value(""))
      ("version",
       "Show version information")
      ("seaLevelForcing",
       "Apply sea level forcing from netcdf scalar time series file.",
       cxxopts::value<std::string>()->default_value(""))
      ("v,verbose",
       "Verbose output. Disables Progressbar");
  // clang-format on

  cxxopts::ParseResult result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

  args.tMax = result["Tmax"].as<PetscScalar>();
  args.tMin = result["Tmin"].as<PetscScalar>();

  // need to be parsed
  args.totaltime = result["totaltime"].as<std::string>();
  args.dt = result["dt"].as<std::string>();
  args.saveEvery = result["saveEvery"].as<int>();
  args.conductivity = result["conductivity"].as<PetscScalar>();
  args.dochannels = result["dochannels"].as<bool>();
  args.disableUnconfined = result["disableUnconfined"].as<bool>();
  args.flowConstant = result["flowConstant"].as<PetscScalar>();
  args.roughnessFactor = result["roughnessFactor"].as<PetscScalar>();
  args.supplyMultiplier = result["supplyMultiplier"].as<PetscScalar>();
  args.layerThickness = result["layerThickness"].as<PetscScalar>();
  args.unconfSmooth = result["unconfSmooth"].as<PetscScalar>();
  args.restart = result["restart"].as<std::string>();
  args.Ssmulti = result["Ssmulti"].as<PetscScalar>();
  args.Sy = result["Sy"].as<PetscScalar>();
  args.Texp = result["Texp"].as<PetscScalar>();
  args.noSmoothMelt = result["noSmoothMelt"].as<bool>();
  args.loopForcing = result["loopForcing"].as<bool>();
  args.basalVelocityIce = result["basalVelocityIce"].as<PetscScalar>();
  args.cavityBeta = result["cavityBeta"].as<PetscScalar>();
  args.initialHead = result["initialHead"].as<std::string>();
  args.tempResults = result["tempResults"].as<std::string>();
  args.version = result["version"].as<bool>();
  args.seaLevelForcing = result["seaLevelForcing"].as<std::string>();
  args.verbose = result["verbose"].as<bool>();
  args.output = result["output"].as<std::string>();
  args.netcdf = result["input"].as<std::string>();
}

}  // namespace CUAS
