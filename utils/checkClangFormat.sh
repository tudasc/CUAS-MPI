#!/bin/bash

notFormated=0

folder="lib tools test"
extensions='.*\.(cpp|h)'

for file in $(find $folder -regextype posix-extended -regex $extensions -not -type d); do
  clang-format --dry-run -Werror $file
  if [ $? -ne 0 ]; then
    notFormated=$(($notFormated+1))
  fi
done

exit $notFormated

