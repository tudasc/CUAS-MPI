/**
 * File: CUASArgs.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PARSECXXOPTS_H
#define CUAS_PARSECXXOPTS_H

#include "petsc.h"

#include "Logger.h"
#include "utilities.h"

#include "cxxopts.hpp"
#include "yaml-cpp/yaml.h"

#include <string>
#include <vector>

namespace CUAS {

class CUASArgs {
 public:
  CUASArgs() { setup(); }
  ~CUASArgs() = default;
  // CUASArgs cannot be copied because cuasOptions is not copyable
  CUASArgs(CUASArgs const &) = delete;
  CUASArgs &operator=(CUASArgs const &) = delete;
  CUASArgs(CUASArgs &&) = delete;
  CUASArgs &operator=(CUASArgs &&) = delete;

 private:
  void setup();

 public:
  std::string configFile = "";

 public:
  // verbosity
  bool verbose = false;
  bool verboseSolver = false;

  // input and output
  std::string input = "";
  std::string output = "out.nc";
  std::string coordinatesFile = "";
  std::string restart = "";
  bool restartNoneZeroInitialGuess = true;

  // time stepping
  std::string starttime = "";
  std::string endtime = "";
  std::string totaltime = "";
  std::string dt = "";
  std::string timeStepFile = "";

  // output behavior
  int saveEvery = 0;
  std::string saveInterval = "";
  std::string outputSize = "normal";

  // forcing
  std::string forcingFile = "";
  int sizeOfForcingBuffer = -1;
  bool loopForcing = false;
  std::string seaLevelForcing = "";

  // solver behavior
  bool directSolver = false;
  int nonLinearIters = 0;
  PetscScalar timeSteppingTheta = 1.0;
  bool enableUDS = false;  // upwind scheme
  bool disableNonNegative = false;

  // channel configuration
  bool doChannels = false;
  std::string selectedChannels = "noselected";
  bool doAllChannels = false;
  bool doAnyChannel = false;
  bool doCavity = false;
  bool doMelt = false;
  bool doCreep = false;

  // physics
  std::string initialHead = "Nzero";
  PetscScalar Tmax = 20.0;
  PetscScalar Tmin = 0.0000001;
  PetscScalar Tinit = 0.2;
  bool disableUnconfined = false;
  PetscScalar conductivity = 10.0;
  PetscScalar flowConstant = 5e-25;
  PetscScalar roughnessFactor = 1.0;
  PetscScalar supplyMultiplier = 1.0;
  PetscScalar layerThickness = 0.1;
  PetscScalar unconfSmooth = 0.0;
  PetscScalar specificStorage = 0.0000982977696;  // Ss
  PetscScalar specificYield = 0.4;                // Sy
  PetscScalar thresholdThicknessUDS = 0.0;
  PetscScalar basalVelocityIce = 1e-6;
  PetscScalar cavityBeta = 5e-4;

  // outflow boundary conditions
  PetscScalar Twater = 100.0;
  PetscScalar dirichletBCWaterDepth = 1.0;
  int blockInflow = 1;
  bool applyRestartChecks = false;

 public:
  class CUASOption {
   public:
    CUASOption(std::string const &optionName, std::string const &description, std::string const &defaultValue)
        : optionIdentifier(optionName), description(description), defaultValue(defaultValue) {
      auto optionTokens = split(optionName, ',');
      if (optionTokens.size() == 1) {
        this->optionName = optionTokens[0];
      } else if (optionTokens.size() == 2) {
        this->optionName = optionTokens[1];
      } else {
        CUAS_ERROR("{}::{}, {} Option is not in the correct format (e.g. 'h,help'): {}", __FILE__, __LINE__, __func__,
                   optionIdentifier);
        exit(1);
      }
    };
    virtual ~CUASOption() = default;
    CUASOption(CUASOption const &) = delete;
    CUASOption &operator=(CUASOption const &) = delete;
    CUASOption(CUASOption &&) = delete;
    CUASOption &operator=(CUASOption &&) = delete;

   public:
    std::string optionIdentifier;
    std::string optionName;
    std::string description;
    std::string defaultValue;

   public:
    virtual void parse(cxxopts::ParseResult const &result, YAML::Node *configFile) = 0;
    virtual void init(cxxopts::Options &options) = 0;
  };

  template <typename ValueType>
  class CUASOptionGeneric : public CUASOption {
   public:
    CUASOptionGeneric(std::string const &optionName, std::string const &description, ValueType *destination,
                      std::string const &defaultValue)
        : CUASOption(optionName, description, defaultValue), destination(destination) {
      if (!destination) {
        CUAS_ERROR("{}::{}, {} nullptr", __FILE__, __LINE__, __func__)
        exit(1);
      }
    }
    CUASOptionGeneric(std::string const &optionName, std::string const &description, ValueType *destination)
        : CUASOptionGeneric(optionName, description, destination, toString(*destination)) {}
    ~CUASOptionGeneric() override = default;
    CUASOptionGeneric(CUASOptionGeneric const &) = delete;
    CUASOptionGeneric &operator=(CUASOptionGeneric const &) = delete;
    CUASOptionGeneric(CUASOptionGeneric &&) = delete;
    CUASOptionGeneric &operator=(CUASOptionGeneric &&) = delete;

    // ValueType defaultValueGeneric;
    ValueType *destination;

   public:
    void parse(cxxopts::ParseResult const &result, YAML::Node *configFile) override {
      // prefer cxxopts, if the option is not explicitly defined, check definition in option file, else use default
      if (result.count(optionName)) {
        *destination = result[optionName].template as<ValueType>();
      } else if (configFile) {
        auto &config = *configFile;
        if (config[optionName].IsDefined()) {
          *destination = config[optionName].template as<ValueType>();
        }
      }
      // else do nothing and keep default value
    }
    void init(cxxopts::Options &options) override {
      options.add_options()(optionIdentifier, description, cxxopts::value<ValueType>()->default_value(defaultValue));
    }
  };

 public:
  std::vector<std::unique_ptr<CUASOption>> cuasOptions;
};

void parseArgs(int argc, char **argv, CUASArgs &args);

}  // namespace CUAS

#endif
