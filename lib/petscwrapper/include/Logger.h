#ifndef PETSC_LOGGER_H
#define PETSC_LOGGER_H

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <memory>

class Logger {
 private:
  std::shared_ptr<spdlog::logger> console;
  std::shared_ptr<spdlog::logger> err_logger;
  Logger() {
    console = spdlog::stdout_color_mt("console");
    err_logger = spdlog::stderr_color_mt("stderr");
  }

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
};

#endif
