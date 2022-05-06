#include "CUASSolver.h"
#include "Forcing/ConstantForcing.h"
#include "Forcing/TimeForcing.h"
#include "ModelReader.h"
#include "SolutionHandler.h"
#include "timeparse.h"

void setupTime(CUAS::Time &time, CUAS::CUASArgs const &args) {
  if (!args.timeStepFile.empty()) {
    if (args.verbose) {
      CUAS_INFO_RANK0("read time step array from file: " + args.timeStepFile);
    }
    CUAS::NetCDFFile file(args.timeStepFile, 'r');
    file.read("time", time.timeSteps);
    time.units = file.readTextAttribute("time", "units");        // --> "seconds since 01-01-01 00:00:00"
    time.calendar = file.readTextAttribute("time", "calendar");  // --> "365_day"
    if (args.verbose) {
      CUAS_INFO_RANK0("time units '" + time.units + "'");
      CUAS_INFO_RANK0("  calendar '" + time.calendar + "'");
    }
  } else {
    if (args.verbose) {
      CUAS_INFO_RANK0("generates time step array using command line paramters");
    }
    time.timeSteps = CUAS::getTimeStepArray(0, CUAS::parseTime(args.totaltime), CUAS::parseTime(args.dt));
  }

  if (args.verbose) {
    CUAS_INFO_RANK0("Number of time steps = {}", time.timeSteps.size() - 1);
  }
}

// store at least initial conditions and the final results
// unless saveEvery is negative --> no output
std::unique_ptr<CUAS::SolutionHandler> setupSolutionHandler(CUAS::CUASArgs &args, CUAS::Time const &time,
                                                            CUAS::CUASModel const &model) {
  if (args.saveEvery < 0)
    return nullptr;

  if (args.saveEvery == 0) {
    CUAS_WARN("Option --saveEvery is == 0, reset to {}", time.timeSteps.size());
    args.saveEvery = time.timeSteps.size();
  }

  return std::make_unique<CUAS::SolutionHandler>(args.output, model.Ncols, model.Nrows, args.outputSize);
}

void restartFromNetCDF(CUAS::CUASSolver &solver, std::string const &restartFile, bool restartNoneZeroInitialGuess) {
  CUAS::ModelReader::restartFromFile(solver, restartFile, restartNoneZeroInitialGuess);
}

int main(int argc, char *argv[]) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  {
    CUAS::CUASArgs args;
    CUAS::parseArgs(argc, argv, args);
    std::string outfile = args.output;

    CUAS::ModelReader reader(args.input);
    auto model = reader.fillModelFromNetcdf();

    model->Q = std::make_unique<CUAS::ConstantForcing>(*model->bmelt, args.supplyMultiplier);

    CUAS::Time time;
    setupTime(time, args);

    std::unique_ptr<CUAS::SolutionHandler> solutionHandler = setupSolutionHandler(args, time, *model);

    if (solutionHandler != nullptr) {
      if (!time.units.empty()) {
        solutionHandler->setTimeUnits(time.units);
      }
      if (!time.calendar.empty()) {
        solutionHandler->setCalendar(time.calendar);
      }
    }

    auto solver = std::make_unique<CUAS::CUASSolver>(model.get(), &args, solutionHandler.get());
    solver->setup();

    if (!args.restart.empty()) {
      restartFromNetCDF(*solver, args.restart, args.restartNoneZeroInitialGuess);
    }

    solver->solve(time.timeSteps);
  }
  PetscFinalize();
  return 0;
}
