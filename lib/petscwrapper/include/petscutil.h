#include "Logger.h"

#include "petsc.h"

#include <string>

std::string getPETScOptionsAll() {
  char *copts;
  PetscOptionsGetAll(PETSC_NULL, &copts);
  auto result = std::string(copts);
  PetscFree(copts);
  return result;
}

void getPETScOptionsAll(std::vector<std::string> &names, std::vector<std::string> &values) {
  // TODO getPETScOptionsAll --> split string
  CUAS_WARN("{} not implemented yet", __PRETTY_FUNCTION__)
  return;
}

void getPETScOptionsUnused(std::vector<std::string> &names, std::vector<std::string> &values) {
  PetscInt N;
  char **tnames;
  char **tvalues;
  PetscOptionsLeftGet(PETSC_NULL, &N, &tnames, &tvalues);

  names.resize(N);
  values.resize(N);

  for (int i = 0; i < N; ++i) {
    auto name = std::string(tnames[i]);
    if (!name.empty()) {
      names[i] = "-" + name;
      auto value = std::string(tvalues[i]);
      if (!value.empty()) {
        values[i] = value;
      }
    } else {
      CUAS_WARN("option with empty name")
    }
  }

  PetscOptionsLeftRestore(PETSC_NULL, &N, &tnames, &tvalues);
  return;
}

std::string getPETScOptionsUnused() {
  std::vector<std::string> names;
  std::vector<std::string> values;
  getPETScOptionsUnused(names, values);

  auto result = std::string();

  for (int i = 0; i < names.size(); ++i) {
    auto name = std::string(names[i]);
    auto value = std::string(values[i]);
    if (!name.empty()) {
      result += name;
      if (!value.empty()) {
        result += " ";
        result += value;
      }
      if (i < names.size() - 1) {
        result += " ";
      }
    } else {
      CUAS_WARN("option with empty name")
    }
  }

  return result;
}

void getPETScOptionsUsed(std::vector<std::string> &names, std::vector<std::string> &values) {
  // TODO getPETScOptionsAll and getPETScOptionsUnused --> used = all not in unused
  CUAS_WARN("{} not implemented yet", __PRETTY_FUNCTION__)
  return;
}

std::string getPETScOptionsUsed() {
  auto result = std::string();

  // TODO getPETScOptionsUsed, concatenate options
  CUAS_WARN("{} not implemented yet", __PRETTY_FUNCTION__)
  result = std::string(__PRETTY_FUNCTION__) + std::string(" not implemented yet");

  return result;
}
