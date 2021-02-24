#include "timeparse.h"

#include "gtest/gtest.h"

TEST(timeparseTest, parse19y1m2w15d5h) {
  std::string inputTime = "10 years 1 month 2 weeks 15 days 5 hours";
  CUAS::timeSecs secs = CUAS::parseTime(inputTime);
  ASSERT_EQ(secs, 320475600);
}

TEST(timeparseTest, parse50y10h) {
  std::string inputTime = "50 years 10 hours";
  CUAS::timeSecs secs = CUAS::parseTime(inputTime);
  ASSERT_EQ(secs, 1576836000);
}

TEST(timeparseTest, parse3h) {
  std::string inputTime = "3 hours";
  CUAS::timeSecs secs = CUAS::parseTime(inputTime);
  ASSERT_EQ(secs, 10800);
}

TEST(timeparseTest, parse1y2y1m) {
  std::string inputTime = "1 year 2 years 1 month";
  CUAS::timeSecs secs = CUAS::parseTime(inputTime);
  ASSERT_EQ(secs, 0);
}
