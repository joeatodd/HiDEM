SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR x86_64)
SET(CMAKE_SYSTEM_VERSION 1)

# Specify the cross compilers (serial)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_Fortran_COMPILER gfortran)
SET(CMAKE_CXX_COMPILER g++)

# Specify the cross compilers (parallel) - don't think these are used
SET(MPI_C_COMPILER mpicc)
SET(MPI_CXX_COMPILER mpic++)
SET(MPI_Fortran_COMPILER mpif90)

# Compilation flags (i.e. with optimization)
SET(CMAKE_C_FLAGS_RELEASE "-O3 -fPIC" CACHE STRING "")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -fPIC" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fPIC" CACHE STRING "")

# Debug compilation flags
SET(CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=bounds,do" CACHE STRING "")
