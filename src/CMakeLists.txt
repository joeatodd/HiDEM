#Define the source code which makes up HiDEM
add_executable (HiDEM typedefs.f90 utils.f90 io.f90 octree.F90 melange.f90 glas.f90 effload.f90 amat.f ttmat.f90 kmat.f ranmar.f wave2.F90)

#Locate MPI and add relevant flags, links
find_package(MPI REQUIRED)

if(MPI_Fortran_COMPILE_FLAGS)
  set_target_properties(HiDEM PROPERTIES
    COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
endif()

if(MPI_Fortran_LINK_FLAGS)
  set_target_properties(HiDEM PROPERTIES
    LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
endif()

if(MPI_Fortran_INCLUDE_PATH)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
endif()

target_link_libraries(HiDEM ${MPI_Fortran_LIBRARIES})

#Link to previously compiled metis library
target_link_libraries(HiDEM metis)

INSTALL(TARGETS HiDEM RUNTIME DESTINATION ".")
