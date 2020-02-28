# Copyright 2019, Max Planck Society.
# Distributed under the "GNU GPL v3" licence.
# (See accompanying file LICENSE.md)

#
# This file contains the global compilation settings
#

# support for C++14
set(CMAKE_CXX_STANDARD 14)

# add version
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROJECT_VERSION='${CMAKE_PROJECT_VERSION}' ")
message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

# build type, by default to release (with optimisations)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# some compilation options
# show all warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra ")

# additional options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -fomit-frame-pointer ")
    set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE} -O3 -ffast-math -fomit-frame-pointer")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -fomit-frame-pointer")
    set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE} -O3 -ffast-math -fomit-frame-pointer")
elseif(MSVC)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(MSVC_Additional_flags "/GF /Oy /GT /Ox /Ob2 /Oi /Os") # /fp:fast
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MSVC_Additional_flags}")
endif()
