#!/bin/bash
module swap PrgEnv-gnu PrgEnv-cray
cd ./src/
ftn -o ../HiDEM -O3 typedefs.f90 utils.f90 io.f90 octree.F90 glas.f90 dist.f90 effload.f90 amat.f tmat.f ttmat.f kmat.f ranmar.f wave2.f90
cd ..
rm INOUT.mod
rm TYPEDEFS.mod
rm LATTICE.mod
rm UTILS.mod
rm OCTREE.mod
module swap PrgEnv-cray PrgEnv-gnu
