cmake_minimum_required(VERSION 3.20)

#	project name
project (Analysis LANGUAGES C CXX)

# global settings
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libs vs static")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# use c++11 at a minimum
set(CMAKE_CXX_STANDARD 11)
if (APPLE)
    add_compile_options(-fvisibility=hidden)
endif()

# Path to correctionlib
set(CORRECTIONLIB_DIR "${CMAKE_SOURCE_DIR}/ext/correctionlib")

# Add built headers path (e.g., correctionlib/src or correctionlib/build/src)
include_directories(${CORRECTIONLIB_DIR}/build/src)

# Include headers
include_directories(${CORRECTIONLIB_DIR}/include)
include_directories(${CORRECTIONLIB_DIR}/build/include) 
include_directories(/afs/cern.ch/user/d/dblend/test-analysis/ext/correctionlib/build/include) 
link_directories(${CORRECTIONLIB_DIR}/build ${CORRECTIONLIB_DIR}/build/src)

# Add built library path
link_directories(${CORRECTIONLIB_DIR}/build/src)

# Optional: log for debugging
message(STATUS "Using correctionlib from ${CORRECTIONLIB_DIR}")

# NOTE: cmake must have a CTestTEstfile.cmkae generated
# at the root of the build folder... otherwise it does not find tests
include(${CMAKE_SOURCE_DIR}/cmake-tools/common.cmake)
include(${CMAKE_SOURCE_DIR}/cmake-tools/target.cmake)

set(USE_CONAN FALSE CACHE BOOL "include conan for external deps")

if (${USE_CONAN})
    # conan setup
    set(_conan_paths ${CMAKE_BINARY_DIR}/conan_paths.cmake)
    if (NOT EXISTS ${_conan_paths})
        tcpp_fail_with_error("file ${_conan_paths} does not exist - rerun conan!")
    endif()
    include(${_conan_paths})
endif()

tcpp_debug_var(CMAKE_PREFIX_PATH)

find_package(ROOT REQUIRED)
  
tcpp_debug_var(ROOT_LIBRARIES)
tcpp_debug_var(ROOT_DEFINITIONS)
tcpp_debug_var(ROOT_INCLUDE_DIRS)
tcpp_debug_var(ROOT_CXX_FLAGS)

message(STATUS "before")
tcpp_auto_addsubdirs(src/)
message(STATUS "after")

#	
#	find and add ROOT
#
#set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
#find_package(ROOT REQUIRED COMPONENTS RHTTP)
#include_directories("${ROOT_INCLUDE_DIR}")
#add_definitions(${ROOT_CXX_FLAGS})
