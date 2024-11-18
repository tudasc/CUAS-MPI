/**
 * File: petscoptions.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_PETSCOPTIONSWRAPPER_H
#define CUAS_PETSCOPTIONSWRAPPER_H

#include "Logger.h"

#include "petsc.h"
#include "petscwrapperutils.h"

#include <string>

std::string getPETScOptionsAll() {
  char *copts;
  PetscOptionsGetAll(PETSC_NULLPTR, &copts);
  auto result = std::string(copts);
  PetscFree(copts);
  return result;
}

void getPETScOptionsAll(std::vector<std::string> &names, std::vector<std::string> &values) {
  // TODO getPETScOptionsAll --> split string
  CUAS_WARN("{} not implemented yet", __PRETTY_FUNCTION__)
  names.emplace_back(__PRETTY_FUNCTION__);
  values.emplace_back("not implemented yet");
}

void getPETScOptionsUnused(std::vector<std::string> &names, std::vector<std::string> &values) {
  PetscInt N;
  char **tnames;
  char **tvalues;
  PetscOptionsLeftGet(PETSC_NULLPTR, &N, &tnames, &tvalues);

  names.resize(N);
  values.resize(N);

  for (int i = 0; i < N; ++i) {
    auto name = std::string(tnames[i]);
    if (!name.empty()) {
      names[i] = "-" + name;
      if (tvalues[i] != nullptr) {
        auto value = std::string(tvalues[i]);
        values[i] = value;
      }
    } else {
      CUAS_WARN("option with empty name")
    }
  }

  PetscOptionsLeftRestore(PETSC_NULLPTR, &N, &tnames, &tvalues);
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
  names.emplace_back(__PRETTY_FUNCTION__);
  values.emplace_back("not implemented yet");
}

std::string getPETScOptionsUsed() {
  auto result = std::string();

  // TODO getPETScOptionsUsed, concatenate options
  CUAS_WARN("{} not implemented yet", __PRETTY_FUNCTION__)
  result = std::string(__PRETTY_FUNCTION__) + std::string(" not implemented yet");

  return result;
}

#endif  // CUAS_PETSCOPTIONSWRAPPER_H
