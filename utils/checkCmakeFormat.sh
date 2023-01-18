#!/bin/bash

notFormated=0

folder="lib tools test"

cmake-format --check CMakeLists.txt
if [ $? -ne 0 ]; then
  notFormated=$(($notFormated+1))
fi

for file in $(find $folder -iname "CMakeLists.txt"); do
  cmake-format --check $file
  if [ $? -ne 0 ]; then
    notFormated=$(($notFormated+1))
  fi
done

cmake-format --check cmake/*
if [ $? -ne 0 ]; then
  notFormated=$(($notFormated+1))
fi

exit $notFormated

