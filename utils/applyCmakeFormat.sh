#!/bin/bash

folder="lib tools test docs"

cmake-format -i "CMakeLists.txt"
find $folder -iname "CMakeLists.txt" -exec cmake-format -i {} \;
cmake-format -i cmake/*

