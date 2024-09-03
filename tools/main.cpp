/**
 * File: main.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include "CUASArgs.h"
#include "CUASConstants.h"
#include "CUASModel.h"
#include "CUASSolver.h"
#include "ModelReader.h"
#include "SolutionHandler.h"
#include "timeparse.h"

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

// store at least initial conditions and the final results
// unless saveEvery is negative --> no output
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

void restartFromNetCDF(CUAS::CUASSolver &solver, std::string const &restartFile, bool restartNoneZeroInitialGuess) {
  CUAS::ModelReader::restartFromFile(solver, restartFile, restartNoneZeroInitialGuess);
}

void setupForcing(CUAS::CUASModel &model, CUAS::CUASArgs &args) {
  if (args.forcingFile.empty()) {
    CUAS_WARN_RANK0("no forcing file given --> try to use input file as forcing file")
    args.forcingFile = args.input;
  }

  // Note, bmelt is still in units m.a-1, but needs to be m.s-1
  if (CUAS::ModelReader::isTimeDependent(args.forcingFile, "bmelt")) {
    if (CUAS::ModelReader::isTimeDependentField(args.forcingFile, "bmelt")) {
      CUAS_INFO_RANK0("Using time forcing with file: " + args.forcingFile)
      if (args.sizeOfForcingBuffer < 2) {
        model.Q = CUAS::ModelReader::getTimeDependentForcing(args.forcingFile, "bmelt", model.xAxis, model.yAxis,
                                                             args.supplyMultiplier / SPY, 0.0, args.loopForcing);
      } else {
        model.Q = CUAS::ModelReader::getBufferedForcing(args.forcingFile, "bmelt", model.xAxis, model.yAxis,
                                                        args.sizeOfForcingBuffer, args.supplyMultiplier / SPY, 0.0,
                                                        args.loopForcing);
      }
    } else {
      CUAS_INFO_RANK0("Using scalar time series forcing with file: " + args.forcingFile)
      model.Q = CUAS::ModelReader::getScalarTimeDependentForcing(args.forcingFile, "bmelt", model.xAxis, model.yAxis,
                                                                 args.supplyMultiplier / SPY, 0.0, args.loopForcing);
    }
  } else {
    CUAS_INFO_RANK0("Using constant forcing with file: " + args.forcingFile)
    model.Q = CUAS::ModelReader::getSteadyForcing(args.forcingFile, "bmelt", model.xAxis, model.yAxis,
                                                  args.supplyMultiplier / SPY, 0.0);
  }
}

int main(int argc, char *argv[]) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  {
    CUAS::CUASArgs args;
    CUAS::parseArgs(argc, argv, args);
    std::string outfile = args.output;

    CUAS::ModelReader reader(args.input);
    auto model = reader.fillModelFromNetcdf();

    setupForcing(*model, args);

    CUAS::Time time;
    setupTime(time, args);

    std::unique_ptr<CUAS::SolutionHandler> solutionHandler = setupSolutionHandler(args, time, *model);

    auto solver = std::make_unique<CUAS::CUASSolver>(model.get(), &args, solutionHandler.get());
    solver->setup();

    if (!args.restart.empty()) {
      restartFromNetCDF(*solver, args.restart, args.restartNoneZeroInitialGuess);
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
