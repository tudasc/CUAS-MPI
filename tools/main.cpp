#include "CUASSolver.h"
#include "Forcing/ConstantForcing.h"
#include "Forcing/TimeForcing.h"
#include "ModelReader.h"
#include "SolutionHandler.h"
#include "timeparse.h"

int main(int argc, char *argv[]) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  {
    CUAS::CUASArgs args;
    CUAS::parseArgs(argc, argv, args);
    std::string outfile = args.output;

    CUAS::ModelReader reader(args.input);
    auto model = reader.fillModelFromNetcdf();

    model->Q = std::make_unique<CUAS::ConstantForcing>(*model->bmelt, args.supplyMultiplier);

    std::unique_ptr<CUAS::SolutionHandler> solutionHandler;
    if (args.saveEvery > 0) {
      solutionHandler =
          std::make_unique<CUAS::SolutionHandler>(args.output, model->Ncols, model->Nrows, args.outputSize);
    } else {
      solutionHandler = nullptr;
    }

    auto solver = std::make_unique<CUAS::CUASSolver>(model.get(), &args, solutionHandler.get());

    solver->setup();

    std::vector<CUAS::timeSecs> timeSteps;
    timeSteps = CUAS::getTimeStepArray(0, CUAS::parseTime(args.totaltime), CUAS::parseTime(args.dt));

    solver->solve(timeSteps);
  }
  PetscFinalize();
  return 0;
}
