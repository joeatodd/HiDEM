cmake_minimum_required(VERSION 2.8)
project(METIS)

set(GKLIB_PATH "${PROJECT_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib")
set(METIS_SHARED FALSE CACHE BOOL "build a shared library")

set(METIS_INSTALL FALSE)

# Configure libmetis library.
if(METIS_SHARED)
  set(METIS_LIBRARY_TYPE SHARED)
else()
  set(METIS_LIBRARY_TYPE STATIC)
endif(METIS_SHARED)

include(${GKLIB_PATH}/GKlibSystem.cmake)
# Add include directories.
include_directories(${GKLIB_PATH})
include_directories(include)
# Recursively look for CMakeLists.txt in subdirs.
add_subdirectory("include")
add_subdirectory("libmetis")
#add_subdirectory("programs")
