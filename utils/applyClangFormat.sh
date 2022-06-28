#!/bin/bash

folder="lib tools test"
extensions='.*\.(cpp|h)'

find $folder -regextype posix-extended -regex $extensions -not -type d -exec clang-format -i {} \;

