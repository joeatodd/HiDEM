#!/bin/bash
cd ./src/
mpif90 -O3 typedefs.f90 io.f90 dist.f90 circ.f90 effload.f90 amat.f tmat.f ttmat.f kmat.f glas.f90 ranmar.f dt.f90 wave2.f90 -I/usr/include/mpi/ -I/usr/include/ -L/usr/lib/x86_64-linux-gnu/ -lmetis -o ../HiDEM
rm inout.mod
rm typedefs.mod
rm lattice.mod
cd ..
