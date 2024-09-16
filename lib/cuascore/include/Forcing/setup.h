/**
 * File: setup.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_FORCING_SETUP_H
#define CUAS_FORCING_SETUP_H

#include "CUASArgs.h"
#include "CUASConstants.h"
#include "CUASModel.h"
#include "Forcing/MultiForcing.h"
#include "ModelReader.h"
#include "utilities.h"

namespace CUAS {

// FIXME: bmelt is still in units m.a-1, but needs to be m.s-1

std::unique_ptr<Forcing> createForcing(std::string const &fileName, CUASModel const &model, CUASArgs const &args,
                                       std::string const &fieldName);

std::unique_ptr<MultiForcing> createMultiForcing(std::vector<std::string> const &fileNames, CUASModel const &model,
                                                 CUASArgs const &args, std::string const &fieldName);

void setupForcing(CUASModel &model, CUASArgs &args, std::string const &fieldName = "bmelt") {
  constexpr auto delimiter = ';';
  if (!args.forcingFile.empty()) {
    if (args.forcingFile.find(delimiter) != std::string::npos) {
      CUAS_INFO_RANK0("Using multi forcing: {}.", args.forcingFile)
      auto fileNames = split(args.forcingFile, delimiter);
      model.setWaterSource(createMultiForcing(fileNames, model, args, fieldName));
    } else {
      CUAS_INFO_RANK0("Using forcing: {}.", args.forcingFile)
      model.setWaterSource(createForcing(args.forcingFile, model, args, fieldName));
    }
  } else {
    CUAS_WARN_RANK0("no forcing file given --> try to use input file as forcing file")
    args.forcingFile = args.input;

    model.setWaterSource(createForcing(args.forcingFile, model, args, fieldName));
  }
}

std::unique_ptr<Forcing> createForcing(std::string const &fileName, CUASModel const &model, CUASArgs const &args,
                                       std::string const &fieldName) {
  if (ModelReader::isTimeDependent(fileName, fieldName)) {
    if (ModelReader::isTimeDependentField(fileName, fieldName)) {
      CUAS_INFO_RANK0("Using time forcing with file: " + fileName)
      if (args.sizeOfForcingBuffer < 2) {
        return ModelReader::getTimeDependentForcing(fileName, fieldName, model.xAxis, model.yAxis,
                                                    args.supplyMultiplier / SPY, 0.0, args.loopForcing);
      } else {
        return ModelReader::getBufferedForcing(fileName, fieldName, model.xAxis, model.yAxis, args.sizeOfForcingBuffer,
                                               args.supplyMultiplier / SPY, 0.0, args.loopForcing);
      }
    } else {
      CUAS_INFO_RANK0("Using scalar time series forcing with file: " + fileName)
      return ModelReader::getScalarTimeDependentForcing(fileName, fieldName, model.xAxis, model.yAxis,
                                                        args.supplyMultiplier / SPY, 0.0, args.loopForcing);
    }
  } else {
    CUAS_INFO_RANK0("Using constant forcing with file: " + fileName)
    return ModelReader::getSteadyForcing(fileName, fieldName, model.xAxis, model.yAxis, args.supplyMultiplier / SPY,
                                         0.0);
  }
}

std::unique_ptr<MultiForcing> createMultiForcing(std::vector<std::string> const &fileNames, CUASModel const &model,
                                                 CUASArgs const &args, std::string const &fieldName) {
  auto multiForcing = std::make_unique<MultiForcing>(model.xAxis.size(), model.yAxis.size());
  for (auto &fileName : fileNames) {
    auto forcing = createForcing(fileName, model, args, fieldName);
    multiForcing->registerNewForcing(forcing);
  }
  return multiForcing;
}

}  // namespace CUAS

#endif
