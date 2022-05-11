# tkleiner: found on https://www.mattkeeter.com/blog/2018-01-06-versioning/
# and modified for cuas
if(Git_FOUND)
execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:'%h' -n 1
        OUTPUT_VARIABLE GIT_REV
        ERROR_QUIET)
endif()

# Check whether we got any revision (which isn't
# always the case, e.g. when someone downloaded a zip
# file from Github instead of a checkout)
if ("${GIT_REV}" STREQUAL "")
    set(GIT_REV "N/A")
    set(GIT_DIFF "")
    set(GIT_TAG "N/A")
    set(GIT_BRANCH "N/A")
else()
    execute_process(
            COMMAND bash -c "${GIT_EXECUTABLE} diff --quiet --exit-code || echo +"
            OUTPUT_VARIABLE GIT_DIFF)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --exact-match --tags
            OUTPUT_VARIABLE GIT_TAG ERROR_QUIET)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            OUTPUT_VARIABLE GIT_BRANCH)

    string(STRIP "${GIT_REV}" GIT_REV)
    string(SUBSTRING "${GIT_REV}" 1 7 GIT_REV)
    string(STRIP "${GIT_DIFF}" GIT_DIFF)
    string(STRIP "${GIT_TAG}" GIT_TAG)
    string(STRIP "${GIT_BRANCH}" GIT_BRANCH)
endif()

set(VERSION "const char *GIT_REV = \"${GIT_REV}${GIT_DIFF}\";
const char *GIT_TAG = \"${GIT_TAG}\";
const char *GIT_BRANCH = \"${GIT_BRANCH}\";")

if(EXISTS ${PROJECT_SOURCE_DIR}/lib/cuascore/src/version.cpp)
    file(READ ${PROJECT_SOURCE_DIR}/lib/cuascore/src/version.cpp VERSION_)
else()
    set(VERSION_ "")
endif()

if (NOT "${VERSION}" STREQUAL "${VERSION_}")
    file(WRITE ${PROJECT_SOURCE_DIR}/lib/cuascore/src/version.cpp "${VERSION}")
endif()
