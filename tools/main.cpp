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
    if (!args.timeStepFile.empty()) {
      if (args.verbose) {
        Logger::instance().info("read time step array from file: " + args.timeStepFile);
      }
      CUAS::NetCDFFile file(args.timeStepFile, 'r');
      file.read("time", timeSteps);
      auto units = file.readTextAttribute("time", "units");        // --> "seconds since 01-01-01 00:00:00"
      auto calendar = file.readTextAttribute("time", "calendar");  // --> "365_day"
      if (args.verbose) {
        Logger::instance().info("time units '" + units + "'");
        Logger::instance().info("  calendar '" + calendar + "'");
      }
      if (solutionHandler != nullptr) {
        if (!units.empty()) {
          solutionHandler->setTimeUnits(units);
        }
        if (!calendar.empty()) {
          solutionHandler->setCalendar(calendar);
        }
      }
    } else {
      if (args.verbose) {
        Logger::instance().info("generates time step array using command line paramters");
      }
      timeSteps = CUAS::getTimeStepArray(0, CUAS::parseTime(args.totaltime), CUAS::parseTime(args.dt));
    }

    solver->solve(timeSteps);
  }
  PetscFinalize();
  return 0;
}
