#!/bin/bash
mpif90 -O3  io.f90 dist.f90 circ.f90 effload.f90 amat.f tmat.f ttmat.f kmat.f glas.f90 ranmar.f dt.f90 wave2.f90 -I/usr/include/mpi/ -o HiDEM
