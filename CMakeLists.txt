# General parameters
cmake_minimum_required(VERSION 3.0)
project(qpmad_benchmark VERSION 0.1.0)


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
add_subdirectory("qpOASES" "qpOASES")
# --------------


# --------------
# ariles
# --------------
set(ARILES_VISITORS_DEFAULT_MODE    OFF     CACHE BOOL "")
set(ARILES_CCACHE                   ON      CACHE BOOL "")
set(ARILES_BUILD_REGRESSION_TESTS   OFF     CACHE BOOL "")

set(ARILES_ENABLE_CORE              ON      CACHE BOOL "")
set(ARILES_VISITOR_rapidjson        ON      CACHE BOOL "")

set(BUILD_SHARED_LIBS   OFF)

add_subdirectory("ariles" "ariles")
# --------------

