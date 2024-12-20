/**
 * File: timeparse.h
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#ifndef CUAS_TIMEPARSE_H
#define CUAS_TIMEPARSE_H

#include "Logger.h"

#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace CUAS {

typedef long timeSecs;

struct Time {
  std::vector<timeSecs> timeSteps;
  std::string units;     // --> "seconds since 01-01-01 00:00:00"
  std::string calendar;  // --> "365_day"
};

// this function parses the input string timeString and returns how many seconds are in the specified time-window
inline timeSecs parseTime(std::string const &timeString) {
  if (timeString.empty()) {
    return 0;
  }

  // map of time-units and the amount of seconds for a single time-unit
  std::map<std::string, timeSecs> timeUnits = {{"year", 60 * 60 * 24 * 365},
                                               {"week", 60 * 60 * 24 * 7},
                                               {"day", 60 * 60 * 24},
                                               {"hour", 60 * 60},
                                               {"minute", 60},
                                               {"second", 1}};
  // map for input times
  std::map<std::string, timeSecs> givenTimes;
  // vector for splitting the input timeString
  std::vector<std::string> splitString;

  // split the timeString by ' '
  std::string s;
  std::istringstream inputStringStream(timeString);
  timeSecs secs = 0;
  while (getline(inputStringStream, s, ' ')) {
    // remove the s at the end if the plural of a time-unit is used
    if (s[s.size() - 1] == 's') {
      s = s.substr(0, s.size() - 1);
    }
    splitString.push_back(s);
  }

  // check if the size is a multiple of 2. This has to be the case as there is always a quantifyer (integer) with the
  // time unit
  if (splitString.size() % 2 != 0) {
    CUAS_ERROR("timeparse.h: Wrong format of input string. Needs to be multiple of 2. Exiting.")
    exit(1);
  }

  // convert the timeString to secs using the timeUnits map and the split timeString
  std::string currentTimeUnit;
  std::string timeValue;
  // With a machine, where int is 32bit bit, the maximum number that can be converted is: 2147483647 (10 digits).
  const std::regex reg("^[0-9]{1,9}$");  // only integer 1 to 9 digits

  for (int i = 1; i < splitString.size(); i = i + 2) {
    // check if current element from splitString is in timeUnits
    if (timeUnits.count(splitString[i]) != 0) {
      // add the value to secs
      timeValue = splitString[i - 1];
      currentTimeUnit = splitString[i];
      if (std::regex_match(timeValue, reg)) {
        // Do what you want (process convert)
        secs += timeUnits[currentTimeUnit] * std::stoi(timeValue);
        // remove the used time-unit so that constructs like 2 years 1 year are impossible in the input string
        timeUnits.erase(currentTimeUnit);
      } else {
        CUAS_ERROR("timeparse.h: Wrong format for timeValue: '" + timeValue + "'! Exiting.")
        exit(1);
      }
    } else {
      CUAS_ERROR(
          "timeparse.h: Wrong format! You either used a non-existing time-unit or used a time-unit multiple times in "
          "the input string! Exiting.")
      exit(1);
    }
  }

  return secs;
}

inline std::string parseTime(timeSecs const secs) {
  // secs < 0 is invalid input
  if (secs < 0) {
    CUAS_ERROR("timeparse.h: Invalid input. secs cannot be less than 0. Exiting.")
    exit(1);
  }
  // second is the smallest timeUnit
  if (secs == 0) {
    return "0 second";  // set to zero instead of 1 to allow for runs with no time steps (input to output runs)
  }

  std::map<std::string, timeSecs> timeUnits = {{"year", 60 * 60 * 24 * 365},
                                               {"week", 60 * 60 * 24 * 7},
                                               {"day", 60 * 60 * 24},
                                               {"hour", 60 * 60},
                                               {"minute", 60},
                                               {"second", 1}};
  // this array is used to impose an order on the map.
  std::array<std::string, 6> timeUnitStrings = {"year", "week", "day", "hour", "minute", "second"};
  assert(timeUnits.size() == timeUnitStrings.size());

  std::string timeString;
  timeSecs secsToBeDistributed = secs;
  for (std::string const &timeUnit : timeUnitStrings) {
    auto timeValue = secsToBeDistributed / timeUnits[timeUnit];
    if (timeValue == 0) {
      continue;
    }

    // append timeValue to timeString
    // append " ", the last " " is removed at the end of the function
    timeString += std::to_string(timeValue) + " " + timeUnit + (timeValue > 1 ? "s" : "") + " ";
    // subtract the secs used for the current timeUnit from the seconds which are still left
    secsToBeDistributed -= (timeValue * timeUnits[timeUnit]);

    if (secsToBeDistributed <= 0) {
      break;
    }
  }
  // remove the last " " from the string
  timeString.pop_back();
  return timeString;
}

inline std::vector<timeSecs> getTimeStepArray(timeSecs startTime, timeSecs endTime, timeSecs dt) {
  // negative time is considered to be fatal
  if (startTime < 0) {
    CUAS_ERROR("{} called with invalid startTime = {} < 0. Exiting.",  //
               __PRETTY_FUNCTION__, startTime)
    exit(1);
  }
  if (endTime < 0) {
    CUAS_ERROR("{} called with invalid endTime = {} < 0. Exiting.",  //
               __PRETTY_FUNCTION__, endTime)
    exit(1);
  }
  if (dt < 0) {
    CUAS_ERROR("{} called with invalid dt = {} < 0. Exiting.",  //
               __PRETTY_FUNCTION__, dt)
    exit(1);
  }

  // end has to be larger than start
  if (startTime > endTime) {
    CUAS_ERROR("{} called with invalid startTime = {} > endTime = {}. Exiting.",  //
               __PRETTY_FUNCTION__, startTime, endTime)
    exit(1);
  }

  // create vector and insert default value startTime
  std::vector<timeSecs> timeSteps;
  auto currTime = startTime;
  timeSteps.push_back(currTime);

  // check if run is diagnostic: dt == 0 or start == end
  // diagnostic runs only use startTime
  if (dt == 0 && startTime == endTime) {
    return timeSteps;
  }
  // nothing to do, but this could be intentional
  if (dt == 0) {
    CUAS_WARN_RANK0("{} called with dt == 0. Ignoring startTime and endTime.", __PRETTY_FUNCTION__)
    return timeSteps;
  }
  // nothing to do, but this could be intentional
  if (startTime == endTime) {
    CUAS_WARN_RANK0("{} called with startTime equal to endTime. Ignoring time step length dt.", __PRETTY_FUNCTION__)
    return timeSteps;
  }

  // Ensures that endTime is included, even if endTime is not an integer multiple of dt ()
  if (endTime % dt != 0) {  // note: % 0 aka modulo(0) is invalid and thus dt == 0 is invalid here
    CUAS_WARN_RANK0("{} endTime is not an integer multiple of dt. Reducing the last time step to match timeEnd.",
                    __PRETTY_FUNCTION__)
  }
  while (currTime < endTime) {
    currTime = std::min(currTime + dt, endTime);  // don't overshoot the endTime
    timeSteps.push_back(currTime);
  }
  return timeSteps;
}

}  // namespace CUAS

#endif
