#include "CUASSolver.h"
#include "Forcing/ConstantForcing.h"
#include "Forcing/TimeForcing.h"
#include "ModelReader.h"
#include "SolutionHandler.h"
#include "timeparse.h"

#include <math.h>

int main(int argc, char *argv[]) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  {
    CUAS::CUASArgs args;
    CUAS::parseArgs(argc, argv, args);
    std::string outfile = args.output;

    PetscScalar totaltime_secs = CUAS::parseTime(args.totaltime);
    PetscScalar dt_secs = CUAS::parseTime(args.dt);
    // rounds division to the nearest integer
    double relationTotalDt = (double)totaltime_secs / (double)dt_secs;
    int Nt = (int)std::rint(relationTotalDt);

    CUAS::ModelReader reader(args.input);
    auto model = reader.fillModelFromNetcdf();

    model->Q = std::make_unique<CUAS::ConstantForcing>(*model->bmelt, args.supplyMultiplier);

    std::unique_ptr<CUAS::SolutionHandler> solutionHandler;
    if (args.saveEvery > 0) {
      solutionHandler = std::make_unique<CUAS::SolutionHandler>(args.output, Nt, args.saveEvery, args.input);
    } else {
      solutionHandler = nullptr;
    }

    auto solver = std::make_unique<CUAS::CUASSolver>(model.get(), &args, solutionHandler.get());

    solver->setup();

    solver->solve(Nt, totaltime_secs, dt_secs);
  }
  PetscFinalize();
  return 0;
}
