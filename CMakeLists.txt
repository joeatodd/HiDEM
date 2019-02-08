cmake_minimum_required (VERSION 2.8)
project (HiDEM Fortran C)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_subdirectory (libs/metis-5.1.0)
add_subdirectory (src)

