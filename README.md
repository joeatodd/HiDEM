Jan's HiDEM model.

##### Compilation ####

Compilation requires the cray ftn compiler:

module swap PrgEnv-gnu PrgEnv-cray

#### Getting data ready ####

Replace no data with zero

scale friction down 

adjust the waterline

Choose the number of cores

do something to the restitution coefficient

#### simulation blow up ####

if the simulation explodes (particles moving too far):

Change DMP and DMP2 - multiply by two? High values for damping would be 16E4 and 8E4

Add drag near the bed

Change the timestep

Change fric scale factor in input - changes bed friction

Avoid friction gradients

Avoid geometry gradients

Kill all motion every x timesteps



#### Files etc ####

These are part of the compiled code:

wave2.f90
dist.f
circ.f
effload.f
amat.f
tmat.f
ttmat.f
kmat.f
glas.f90
ranmar.f
dt.f90

ave*.f - these compute averages of something

*.job - job scripts for torque

comp.txt - the compilation command

rc2.f90 - compute the calved size distrib

glas.f90 - computing the FCC lattice, dense packing
dist.f - efficiently finding neighbouring particles which may interact
tmat.f, ttmat.f - rotation matrices
kmat.f - stiffness matrix computation
amat.f - called by effload, does integration
ranmar.f - random number generator
dt.f90 - called by glas.f90, finds and write connections to FSfiles

### input files ###

inp.dat
mass3.dat
param.dat - don't worry - something about opening files

### mass3.dat ####

mass3.dat is the input configuration, in format:

x, y, surface, base, bed. friction (Newton seconds per metre) 

x and y must start at zero, must have a (0,0) corner.

output transformation matrix which takes from Elmer domain to HiDEM domain.

make sure the bed is buffered beyond the edge of the ice, and define these regions by setting surf and base equal to bed.


"YN" is the number of partitions in the Y direction.
"NTOT" is the total number of partitions in the model
GRID - The grid size of the mass3.dat data, make sure it matches
SCL - diameter of each particle
RESTART - 1 = true



#### OUTPUT - jyr files and STR files ####

List position of all particles in x,y,z, every 2 seconds.
Read this in paraview quite easily.

STR file:

midpoint
strain

for each node connection in initial geometry (including broken bonds)


dtop* files - these show the total energy in the system for different parts

dtop00 - T,WENS,ENMS+ENMS0,KINS,MGHS-MGH0
dtop01 - T,DPES,DMPENS,PSUMS,GSUMS
dtopr  - T,system energy, damping energy?

Time
WENS - elastic energy of the spheres - imagine the particles overlap and are deformed against each other. Then they might bounce back apart.
ENMS+ENMS0 - elastic deformation energy - the energy held in a deformed/bent system
MGHS-MGH0 - potential energy
Kins - translational kinetic energy
Kins2 - rotational kinetic energy

DPES,DMPENS - guess energy lost to drag and damping?
PSUMS - something like pressure
GSUMS - energy of bed interaction


#### Processing ####

rc2.f90 computes the size distribution of calved blocks

size  count

1  452000
2  4560
3  985
...
largest_block 1

'maxi' is the same as jyr but for the largest intact block (e.g. the remaining glacier)

rh2.f90 - computes velocity
rh.f - computes strain


INFI1 - first JYR
INFI2 - second JYR
N - wc -l JYR0001.csv
max - 10?

#### PARAVIEW ####

First simply load csv, with blank space delimiter, merge delimiters, has no header
Then apply table to points filter
Then ensure enable ospray and shadows, and pick point size that makes them touch

For the sea: add a box, set size and centrepoints

For the bed: get bed.csv, read it in same as other CSVs
then apply Delaunay2D filter.
