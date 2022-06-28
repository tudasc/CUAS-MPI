#!/bin/bash

exitCode=0

cd build

CHECKS="-*,\
readability-*,\
-readability-make-member-function-const,\
-readability-redundant-access-specifiers,\
-readability-qualified-auto,\
-readability-isolate-declaration,\
-readability-function-cognitive-complexity,\
-readability-isolate-declaration,\
"

clang-tidy --checks=$CHECKS --warnings-as-errors=$CHECKS ../lib/petscwrapper/include/*
if [ $? -ne 0 ]; then
  exitCode=$(($exitCode+1))
fi

clang-tidy --checks=$CHECKS --warnings-as-errors=$CHECKS ../lib/petscwrapper/src/*
if [ $? -ne 0 ]; then
  exitCode=$(($exitCode+1))
fi

exit $exitCode

#clang-tidy --fix -checks="-*,modernize-use-nodiscard" ../lib/petscwrapper/include/*

