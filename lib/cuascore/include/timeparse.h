#ifndef CUAS_TIMEPARSE_H
#define CUAS_TIMEPARSE_H

#include "Logger.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace CUAS {

typedef long timeSecs;
// this function parses the input string timeString and returns how many seconds are in the specified time-window
inline timeSecs parseTime(std::string const &timeString) {
  // map of time-units and the amount of seconds for a single time-unit
  std::map<std::string, int> timeUnits = {{"year", 60 * 60 * 24 * 365},
                                          {"month", 60 * 60 * 24 * 30},
                                          {"week", 60 * 60 * 24 * 7},
                                          {"day", 60 * 60 * 24},
                                          {"hour", 60 * 60}};
  // map for input times
  std::map<std::string, int> givenTimes;
  // vector for splitted input timeString
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
    Logger::instance().error("timeparse.h: Wrong format of input string. Needs to be multiple of 2. Exiting.");
    exit(1);
  }

  // convert the timeString to secs using the timeUnits map and the split timeString
  std::string currentTimeUnit;
  std::string timeValue;
  for (int i = 1; i < splitString.size(); i = i + 2) {
    // check if current element from splitString is in timeUnits
    if (timeUnits.count(splitString[i]) != 0) {
      // add the value to secs
      timeValue = splitString[i - 1];
      currentTimeUnit = splitString[i];
      secs += timeUnits[currentTimeUnit] * stoi(timeValue);
      // remove the used time-unit so that constructs like 2 years 1 year are impossible in the input string
      timeUnits.erase(currentTimeUnit);
    } else {
      Logger::instance().error(
          "timeparse.h: Wrong format! You either used a non-existing time-unit or used a time-unit multiple times in "
          "the input string! Exiting.");
      exit(1);
    }
  }

  return secs;
};

inline std::string parseTime(timeSecs const secs) {
  // secs < 0 is invalid input
  if (secs < 0) {
    Logger::instance().error("timeparse.h: Invalid input. secs cannot be less than 0. Exiting.");
    exit(1);
  }
  // hour is the smallest timeUnit
  if (secs < 60 * 60) {
    return "0 hours";
  }
  std::map<std::string, int> timeUnits = {{"year", 60 * 60 * 24 * 365},
                                          {"month", 60 * 60 * 24 * 30},
                                          {"week", 60 * 60 * 24 * 7},
                                          {"day", 60 * 60 * 24},
                                          {"hour", 60 * 60}};
  // this vector is used to impose an order on the map.
  std::vector<std::string> timeUnitStrings = {"year", "month", "week", "day", "hour"};
  std::string timeString = "";
  timeSecs secsToBeDistributed = secs;
  std::string currentTimeUnit;
  for (std::string const &i : timeUnitStrings) {
    int timeValue = (int)(secsToBeDistributed / timeUnits[i]);
    if (timeValue == 0) {
      continue;
    } else {
      // append timeValue to timeString
      timeString += std::to_string(timeValue) + " " + i + (timeValue > 1 ? "s" : "") + " ";
    }
    // subtract the secs used for the current timeUnit from the seconds which are still left
    secsToBeDistributed -= (timeValue * timeUnits[i]);
    if (secsToBeDistributed <= 0) {
      break;
    }
  }
  // remove the last " " from the string
  timeString.pop_back();
  return timeString;
}

inline std::vector<timeSecs> getTimeStepArray(timeSecs startTime, timeSecs endTime, timeSecs dt) {
  std::vector<CUAS::timeSecs> timeSteps;
  auto currTime = startTime;
  timeSteps.push_back(currTime);
  while (currTime < endTime) {
    currTime += dt;
    timeSteps.push_back(currTime);
  }
  return timeSteps;
}

}  // namespace CUAS

#endif
