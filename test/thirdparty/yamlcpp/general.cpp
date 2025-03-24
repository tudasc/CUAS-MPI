/**
 * File: readyaml.cpp
 * License: Part of the CUAS-MPI project. Licensed under BSD 3 clause license. See LICENSE.txt file at
 * https://github.com/tudasc/CUAS-MPI/LICENSE.txt
 */

#include <yaml-cpp/yaml.h>

#include "gtest/gtest.h"

TEST(yamlcppTest, readSimpleYAML) {
  YAML::Node config = YAML::LoadFile("testconfig.yaml");

  // read values from file
  auto name = config["name"].as<std::string>();
  auto enabled = config["enabled"].as<bool>();
  auto id = config["attributes"]["id"].as<int>();
  auto value = config["attributes"]["value"].as<double>();
  auto features = config["attributes"]["features"].as<std::vector<std::string>>();
  auto input = config["input"].as<std::string>();
  auto output = config["output"].as<std::string>();
  auto random = config["random"].as<int>();

  // check values
  EXPECT_EQ(name, "Simple Config");
  EXPECT_TRUE(enabled);
  EXPECT_EQ(id, 123);
  EXPECT_DOUBLE_EQ(value, 3.14);
  EXPECT_EQ(features, std::vector<std::string>({"feature1", "feature3"}));
  EXPECT_EQ(input, "in.txt");
  EXPECT_EQ(output, "out.txt");
  EXPECT_EQ(random, 32523423);
}

TEST(yamlcppTest, readMissingValues) {
  YAML::Node config = YAML::LoadFile("testconfig.yaml");

  // try to read missing value from file
  int value = 0;
  if (config["missingValue"].IsDefined()) {
    value = config["missingValue"].as<int>();
  }
  // check value
  EXPECT_EQ(value, 0);

  // now read an existing value
  if (config["random"].IsDefined()) {
    value = config["random"].as<int>();
  }
  // check value
  EXPECT_EQ(value, 32523423);
}

// TODO check missing groups

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  return result;
}
