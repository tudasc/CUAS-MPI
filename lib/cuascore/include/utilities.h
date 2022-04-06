#ifndef CUAS_UTILITIES_H
#define CUAS_UTILITIES_H

#include <string>

namespace CUAS {

/**
 * see https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
 */
std::string vformat(const char *format, ...);

/**
 * Provide string containing CUAS version number
 */
std::string version();

}  // namespace CUAS

#endif
