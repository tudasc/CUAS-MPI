/**
 * File: Logger.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef PETSC_LOGGER_H
#define PETSC_LOGGER_H

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include "mpi.h"

#include <memory>

#define CUAS_INFO(...) Logger::instance().info(__VA_ARGS__);
#define CUAS_WARN(...) Logger::instance().warn(__VA_ARGS__);
#define CUAS_ERROR(...) Logger::instance().error(__VA_ARGS__);
#define CUAS_INFO_RANK0(...) Logger::instance().info0(__VA_ARGS__);
#define CUAS_WARN_RANK0(...) Logger::instance().warn0(__VA_ARGS__);
#define CUAS_ERROR_RANK0(...) Logger::instance().error0(__VA_ARGS__);

class Logger {
 private:
  std::shared_ptr<spdlog::logger> console;
  std::shared_ptr<spdlog::logger> err_logger;
  Logger() {
    console = spdlog::stdout_color_mt("console");
    err_logger = spdlog::stderr_color_mt("stderr");
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  int rank;

 public:
  Logger(const Logger &) = delete;
  Logger &operator=(const Logger &) = delete;
  Logger(Logger &&) = delete;
  Logger &operator=(Logger &&) = delete;

  static auto &instance() {
    static Logger logger;
    return logger;
  }

  template <typename... ArgsT>
  void error(const ArgsT... args) const {
    err_logger->error(args...);
  }
  template <typename... ArgsT>
  void warn(const ArgsT... args) const {
    console->warn(args...);
  }
  template <typename... ArgsT>
  void info(const ArgsT... args) const {
    console->info(args...);
  }
  template <typename... ArgsT>
  void error0(const ArgsT... args) const {
    if (rank == 0) {
      err_logger->error(args...);
    }
  }
  template <typename... ArgsT>
  void warn0(const ArgsT... args) const {
    if (rank == 0) {
      console->warn(args...);
    }
  }
  template <typename... ArgsT>
  void info0(const ArgsT... args) const {
    if (rank == 0) {
      console->info(args...);
    }
  }
};

#endif
