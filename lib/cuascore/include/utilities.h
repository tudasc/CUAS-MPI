/**
 * File: utilities.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_UTILITIES_H
#define CUAS_UTILITIES_H

#include <sstream>
#include <string>
#include <vector>

namespace CUAS {

/**
 * see https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
 */
std::string vformat(const char *format, ...);

/**
 * Provide string containing CUAS version number
 */
std::string version();

/**
 * Checks if a vector is strictly increasing.
 */
template <typename T>
bool isIncreasing(const std::vector<T> v) {
  for (typename std::vector<T>::size_type i = 0; i != v.size() - 1; ++i) {
    if (v[i] >= v[i + 1]) {
      return false;
    }
  }
  return true;
}

/**
 * string split
 * @param s input string to split
 * @param delim delimiter
 * @return a vector which contains the individual strings
 */
inline std::vector<std::string> split(std::string const &s, char delim) {
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;

  while (getline(ss, item, delim)) {
    result.push_back(item);
  }

  return result;
}

}  // namespace CUAS

#endif
