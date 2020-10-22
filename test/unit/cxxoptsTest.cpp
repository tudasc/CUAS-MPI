#include "parseCxxopts.h"

// run with --verbose or -v
int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  CUASArgs cuArgs;
  parseArgs(argc, argv, cuArgs);
  std::cout << "verbose? " << cuArgs.verbose << std::endl;
  PetscFinalize();
  return 0;
}
