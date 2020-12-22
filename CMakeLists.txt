# General parameters
cmake_minimum_required(VERSION 3.0)
project(qpmad_benchmark VERSION 0.1.0)


set (CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")


# --------------
# Options
# --------------
set(CMAKE_VERBOSE_MAKEFILE  ON)
set(CMAKE_BUILD_TYPE        "Release")
# --------------


# --------------
# qpmad
# --------------
set(QPMAD_USE_HOUSEHOLDER   OFF CACHE BOOL "")
add_subdirectory("qpmad" "qpmad")
# --------------


# --------------
# qpOASES
# --------------
add_definitions("-D__SUPPRESSANYOUTPUT__")
add_subdirectory("qpOASES" "qpOASES")
# --------------


# --------------
# eiquadprog
# --------------
set(INSTALL_PKG_CONFIG_FILE OFF CACHE BOOL "")
set(LOGGING_FILENAME "/dev/null" CACHE STRING "")
add_subdirectory("eiquadprog" "eiquadprog")
include_directories(SYSTEM "eiquadprog/include")
# --------------


# --------------
# ariles
# --------------
set(ARILES_VISITORS_DEFAULT_MODE    OFF     CACHE BOOL "")
set(ARILES_CCACHE                   ON      CACHE BOOL "")
set(ARILES_BUILD_REGRESSION_TESTS   OFF     CACHE BOOL "")

set(ARILES_ENABLE_CORE              ON      CACHE BOOL "")
set(ARILES_VISITOR_rapidjson        ON      CACHE BOOL "")
set(ARILES_VISITOR_octave           ON      CACHE BOOL "")

set(BUILD_SHARED_LIBS   OFF)

add_subdirectory("ariles" "ariles")
# --------------


# --------------
find_package(Boost REQUIRED filesystem)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
# --------------


# --------------
enable_testing()
add_subdirectory("${PROJECT_SOURCE_DIR}/test")
# --------------
