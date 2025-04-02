/**
 * File: main.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASArgs.h"
#include "CUASModel.h"
#include "CUASSolver.h"
#include "Forcing/setup.h"
#include "ModelReader.h"
#include "SolutionHandler.h"
#include "timeparse.h"

/**
 * setup of CUAS-MPI main application
 */
namespace CUASSetup {

/**
 * differs different ways of defining the time in CUAS-MPI and sets the variables
 *
 * @param time
 * @param args
 */
void setupTime(CUAS::Time &time, CUAS::CUASArgs const &args) {
  if (!args.timeStepFile.empty()) {
    if (args.verbose) {
      CUAS_INFO_RANK0("read time step array from file: " + args.timeStepFile)
    }
    CUAS::NetCDFFile file(args.timeStepFile, 'r');
    file.read("time", time.timeSteps);
    time.units = file.readTextAttribute("time", "units");        // --> "seconds since 01-01-01 00:00:00"
    time.calendar = file.readTextAttribute("time", "calendar");  // --> "365_day"
    if (args.verbose) {
      CUAS_INFO_RANK0("time units '" + time.units + "'")
      CUAS_INFO_RANK0("  calendar '" + time.calendar + "'")
    }
  } else {
    if (args.verbose) {
      CUAS_INFO_RANK0("generates time step array using command line parameters")
    }
    auto starttime = CUAS::parseTime(args.starttime);
    auto endtime = CUAS::parseTime(args.endtime);
    auto totaltime = CUAS::parseTime(args.totaltime);
    auto dt = CUAS::parseTime(args.dt);

    if (totaltime > 0 && endtime > 0) {
      CUAS_WARN_RANK0("endtime and totaltime are defined, totaltime '{}' is ignored, endtime '{}' is used.",
                      args.totaltime, args.endtime)
    }
    if (totaltime <= 0 && endtime <= 0) {
      CUAS_WARN_RANK0("neither totaltime nor endtime are defined.")
    }

    if (endtime <= 0 && totaltime > 0) {
      endtime = starttime + totaltime;
    }

    CUAS_INFO_RANK0("generating time step array starttime {}, endtime {}, totaltime {}, dt {}.",
                    CUAS::parseTime(starttime), CUAS::parseTime(endtime), CUAS::parseTime(totaltime),
                    CUAS::parseTime(dt))
    time.timeSteps = CUAS::getTimeStepArray(starttime, endtime, dt);
  }

  if (args.verbose) {
    if (time.timeSteps.size() == 1) {
      CUAS_INFO_RANK0("No time steps requested. This run is diagnostic only!")
    } else {
      CUAS_INFO_RANK0("Number of time steps = {}", time.timeSteps.size() - 1)
    }
  }
}

/**
 * setup of the CUAS SolutionHandler
 *
 * store at least initial conditions and the final results
 * unless saveEvery is negative --> no output
 *
 * @param args
 * @param time
 * @param model
 */
std::unique_ptr<CUAS::SolutionHandler> setupSolutionHandler(CUAS::CUASArgs &args, CUAS::Time const &time,
                                                            CUAS::CUASModel const &model) {
  if (args.saveEvery < 0 && args.saveInterval.empty()) {
    return nullptr;
  }

  // TODO we currently do not differentiate whether we run with dt, timesteparray or coupled

  std::unique_ptr<CUAS::SolutionHandler> solutionHandler;

  if (!args.saveInterval.empty()) {
    auto saveInterval = CUAS::parseTime(args.saveInterval);
    if (args.saveEvery > 0) {
      CUAS_WARN("Both --saveInterval and --saveEvery are used: ignoring --saveEvery.")
    }
    solutionHandler = std::make_unique<CUAS::SolutionHandler>(args.output, model.Ncols, model.Nrows, args.outputSize);
    solutionHandler->setSaveStrategy(CUAS::SaveStrategy::TIMEINTERVAL, saveInterval, -1);
  } else {
    auto saveEvery = args.saveEvery;
    if (saveEvery == 0) {
      CUAS_WARN("Option --saveEvery is == 0, reset to {}", time.timeSteps.size())
      saveEvery = static_cast<int>(time.timeSteps.size());
    }
    solutionHandler = std::make_unique<CUAS::SolutionHandler>(args.output, model.Ncols, model.Nrows, args.outputSize);
    solutionHandler->setSaveStrategy(CUAS::SaveStrategy::INDEX, -1, saveEvery);
  }

  // TODO move to solution handler constructor?
  if (!time.units.empty()) {
    solutionHandler->setTimeUnits(time.units);
  }
  if (!time.calendar.empty()) {
    solutionHandler->setCalendar(time.calendar);
  }

  return solutionHandler;
}

/**
 * setup to restart CUAS from an existing solution in NetCDF
 *
 * @param solver
 * @param restartFile
 * @param restartNoneZeroInitialGuess
 */
void restartFromNetCDF(CUAS::CUASSolver &solver, std::string const &restartFile, bool restartNoneZeroInitialGuess) {
  CUAS::ModelReader::restartFromFile(solver, restartFile, restartNoneZeroInitialGuess);
}

}  // namespace CUASSetup

/**
 * CUAS-MPI main function
 */
int main(int argc, char *argv[]) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  {
    CUAS::CUASArgs args;
    parseArgs(argc, argv, args);

    CUAS::ModelReader reader(args.input);
    auto model = reader.fillModelFromNetcdf();

    setupForcing(*model, args);

    CUAS::Time time;
    CUASSetup::setupTime(time, args);

    std::unique_ptr<CUAS::SolutionHandler> solutionHandler = CUASSetup::setupSolutionHandler(args, time, *model);

    auto solver = std::make_unique<CUAS::CUASSolver>(model.get(), &args, solutionHandler.get());
    solver->setup();

    if (!args.restart.empty()) {
      CUASSetup::restartFromNetCDF(*solver, args.restart, args.restartNoneZeroInitialGuess);
    }

    solver->solve(time.timeSteps);

    // todo: store model run-time
    if (solutionHandler) {
      solutionHandler->storePETScOptions();
    }
  }
  PetscFinalize();
  return 0;
}
