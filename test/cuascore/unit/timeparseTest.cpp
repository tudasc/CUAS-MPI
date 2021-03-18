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

TEST(timeparseReverseTest, parse1y1mReverse) {
  std::string inputTime = "1 year 1 month";
  CUAS::timeSecs secs = 34128000;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, parse3hReverse) {
  std::string inputTime = "3 hours";
  CUAS::timeSecs secs = 10800;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, parse50y10hReverse) {
  std::string inputTime = "50 years 10 hours";
  CUAS::timeSecs secs = 1576836000;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, parse19y1m2w6d5hReverse) {
  std::string inputTime = "10 years 1 month 2 weeks 6 days 5 hours";
  CUAS::timeSecs secs = 319698000;
  std::string compareString = CUAS::parseTime(secs);
  ASSERT_EQ(inputTime, compareString);
}

TEST(timeparseReverseTest, wrongInput) {
  CUAS::timeSecs secs = -1;
  std::string timeString = CUAS::parseTime(secs);
  ASSERT_EQ(timeString, "");
}

TEST(timeparseReverseTest, lessThan1hour) {
  CUAS::timeSecs secs = 60;
  std::string timeString = CUAS::parseTime(secs);
  ASSERT_EQ(timeString, "0 hours");
}
