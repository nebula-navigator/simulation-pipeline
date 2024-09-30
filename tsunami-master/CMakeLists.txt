cmake_minimum_required(VERSION 3.14)
project(TSUNAMI LANGUAGES CXX)

cmake_policy(VERSION 3.14...3.26)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)
# Set a default build type if none was specified
set(default_build_type "Release")
include(CMakeDependentOption)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_FLAGS_RELEASE " -W -Wall -O3 -march=native") #-flto=auto
set(CMAKE_CXX_FLAGS_DEBUG " -W -Wall -ggdb -O0 -fno-inline-small-functions -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_AMUSE " ")

set(CMAKE_VERBOSE_MAKEFILE ON)

# set cmake binary and library directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include_directories(lib)
set(CMAKE_SWIG_FLAGS)

#########################################################
######## OPTIONS ########################################
######################################################### Enable spins
option(spin "Enable spin for PN term and tides" OFF)
if(spin)
    message("-- Enabling spin")
    add_compile_definitions(SPINS)
    list(APPEND CMAKE_SWIG_FLAGS "-DSPINS")
endif(spin)

######################################################### Debug prints options
option(dbg_lf "Enable leapfrog debug" OFF)
if(dbg_lf)
    message("-- Enabling leapfrog debug")
    add_compile_definitions(DEBUGLF)
    list(APPEND CMAKE_SWIG_FLAGS "-DDEBUGLF")
endif(dbg_lf)

option(dbg_bs "Enable Bulirsh-Stoer debug" OFF)
if(dbg_bs)
    message("-- Enabling Bulirsh-Stoer debug")
    add_compile_definitions(DEBUGBS)
    list(APPEND CMAKE_SWIG_FLAGS "-DDEBUGBS")
endif(dbg_bs)

option(dbg_ch "Enable chain debug" OFF)
if(dbg_ch)
    message("-- Enabling chain debug")
    add_compile_definitions(DEBUGCH)
    list(APPEND CMAKE_SWIG_FLAGS "-DDEBUGCH")
endif(dbg_ch)

######################################################### Enable timing
option(timing "Enables timing" ON)
if(timing)
    message("-- Timing enabled")
    add_compile_definitions(TIMINGU)
else()
    message("-- Timing disabled")
endif(timing)

######################################################### Santai classification scheme
option(profile "Include time profiling of integration" OFF)
if(profile)
    message("-- Profiling enabled")
    add_compile_definitions(PROFILE)
    list(APPEND CMAKE_SWIG_FLAGS "-DPROFILE")
endif(profile)

######################################################### Chaos absorb
option(santai "Use santai classification scheme (only for N=3)" OFF)
CMAKE_DEPENDENT_OPTION(absorb "Stop triple when 'absorbed' has been measured" OFF "santai" OFF)
if(santai)
    message("-- Using santai classification scheme (DO NOT USE FOR N!=3)")
    target_compile_options(base PUBLIC -DSANTAI)
    if(absorb)
        message("-- Absorption measurement mode: triple will stop at first excursion")
        target_compile_options(base PUBLIC -DABSORB)
    endif(absorb)
endif(santai)

######################################################### TDETracker
option(tdetracker "Include tracking of TDEs" OFF)
if(tdetracker)
    message("-- TDETracker enabled")
    add_compile_definitions(TDETRACKER)
    list(APPEND CMAKE_SWIG_FLAGS "-DTDETRACKER")
endif(tdetracker)

####### END OPTIONS #####################################
#########################################################

option(cpp "Build C++ executable" ON)
if(cpp)
    add_library(tsunamilibstatic STATIC lib/chain.hpp lib/leapfrog_stepped.hpp lib/bulirsch.hpp src/IO.cpp src/keplerutils.cpp lib/tsunami.hpp src/classification.cpp)
    set_target_properties(tsunamilibstatic PROPERTIES OUTPUT_NAME tsunamilib)
    message("-- Compiling C++ executable")

    add_executable(tsunami.x ${CMAKE_SOURCE_DIR}/tests/tsunami/main.cpp)
    target_link_libraries(tsunami.x tsunamilibstatic)

    #copy input files to the bin folder
    add_custom_target(copy-input-dir ALL COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/input/ ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/input/)
    add_custom_target(copy-bench-dir ALL COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/bench/ ${CMAKE_BINARY_DIR}/bench/)
endif(cpp)


# Compile as python module using swig
option(python "Compile as Python module using SWIG" ON)
if(python)
    message("-- Compiling Python interface")
    include(python)
endif(python)
