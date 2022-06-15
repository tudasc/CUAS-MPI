#ifndef CUAS_UTILITIES_H
#define CUAS_UTILITIES_H

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

}  // namespace CUAS

#endif
