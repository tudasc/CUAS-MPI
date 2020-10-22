#include "parseCxxopts.h"

void parseArgs(int argc, char **argv, CUASArgs &args) {
  cxxopts::Options options("CUAS", "MPI parallel version of CUAS");

  options.add_options()("x,Tmax", "Maximum T to be allowed in the evolution.",
                        cxxopts::value<PetscScalar>()->default_value("20.0"))(
      "i,Tmin", "Minimum T to be allowed in the evolution.", cxxopts::value<PetscScalar>()->default_value("0.0000001"))(
      "h,help", "Print help");

  cxxopts::ParseResult result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

  args.tMax = result["Tmax"].as<PetscScalar>();
  args.tMin = result["Tmin"].as<PetscScalar>();
}
