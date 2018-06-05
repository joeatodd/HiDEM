# The Helsinki Discrete Element Model (HiDEM) #

Author: Jan Åström

## Compilation ##

HiDEM can be compiled using one of the compilation scripts in this directory.

`./compile_cray.sh`

## input files ##

On starting, HiDEM reads the name of an input file from HIDEM\_STARTINFO (e.g. 'inp.dat').

inp.dat - defines various parameters for the simulation

The user must also supply the initial geometry in gridded format: X,Y,SURF,BASE,BED,FRICTION

This initial geometry file is specified in the input file as 'Geometry File':

`Geometry File = "mass3.dat"`

head of a sample mass3.dat file:

```
15851
0.000000000000000000e+00	0.000000000000000000e+00	9.248271630519999462e+02	9.248271630519999462e+02	9.248271630519999462e+02	1.833796257881020755e+07  
5.000000000000000000e+01	0.000000000000000000e+00	9.248386672000000317e+02	9.248386672000000317e+02	9.248386672000000317e+02	2.784068307483682781e+07  
1.000000000000000000e+02	0.000000000000000000e+00	9.248503686919999609e+02	9.248503686919999609e+02	9.248503686919999609e+02	4.768270480512356758e+07
```

## Getting data ready ##

Replace no data with zero

scale friction down 

adjust the waterline

Choose the number of cores

do something to the restitution coefficient

Point to the input file (e.g. testinp.dat) using HIDEM_STARTINFO

## Running the model ##

e.g. `mpirun -n 70 HiDEM`

An example PBS job script is provided in example.job

## simulation blow up ##

if the simulation explodes (particles moving too far):

Change DMP and DMP2 - multiply by two? High values for damping would be 16E4 and 8E4

Add drag near the bed

Change the timestep

Change fric scale factor in input - changes bed friction

Avoid friction gradients

Avoid geometry gradients

Kill all motion every x timesteps

## Files etc ##

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

compile_*.sh - compilation scripts

rc2.f90, rc3.f90 - compute the calved size distrib  

wave2.f90 - the main program  
dist.f - efficiently finding neighbouring particles which may interact  
circ.f - confirms which particles are in contact and computes forces  
effload.f - computes elastic forces in particle beams  
amat.f - called by effload, does integration  
tmat.f, ttmat.f - rotation matrices  
kmat.f - stiffness matrix computation  
glas.f90 - computing the FCC lattice, dense packing  
ranmar.f - random number generator  
dt.f90 - called by glas.f90, finds and write connections to FSfiles  

### inp.dat parameters ###
| Parameter | VarName | Description |
| --------- | ------- | ----------- |
| Run Name  | RUNNAME | The name of the simulation (prepended to output filenames) |
| Work Directory | wrkdir | The subdirectory in which to store working files (must exist!) |
| Results Directory | resdir | The subdirection in which to store result files (must exist!) |
| Geometry File | geomfile | The name of the file which defines the geometry (see above) |
| Density | RHO | The density of the material (ice) |
| Water Density | RHOW | The density of the water in which the ice floats |
| Gravity | GRAV | Magnitude of gravity (positive!) |
| Backwall Pressure | PRESS | optional backwall pressure    |
| Submarine Melt    | MELT  | optional basal melt rate passed to fibg3 for altering domain shape   |
| UC                | UC    | optional frontal melt rate   |
| Timestep          | DT    | timestep size   |
| Width             | S     | beam width (relative to unit particle)   |
| Youngs Modulus    | EF0   | particle bond young's mod, describes interaction between butting particles (or bonded particles?)   |
| Size              | LS    | Something to do with the fcc lattice 'box size'   |
| Domain Inclination | SUB   | Angle of the domain vs gravity vector, not used    |
| Water Line   | WL    | Sea level for buoyancy calc   |
| Grounding Line  | GL    | 'grounding line' - not used   |
| Shear Line   | SLIN  | A distance below the surface where all bonds are broken?   |
| No Timesteps   | STEPS0| The number of steps   |
| Max Load     | MLOAD | Maximum load on bond - bonds break beyond this   |
| Friction Scale    | FRIC  | Scale factor for friction input   |
| Restart     | REST  | 1 = restart from prev, 0 = new run   |
| Scale         | SCL   | Scale factor for particle and beam size   |
| Grid        | GRID  | Resolution of mass3.dat input grid   |
| Porosity         | POR   | The proportion of initially broken bonds   |
| Random Seed       | SEEDI | Seed for random number generator   |
| Translational Damping       | DAMP1 | The damping coefficient for translation   |
| Rotational Damping       | DAMP2 | The damping coefficient for rotation   |
| Drag Coefficient        | DRAG  | The drag coefficient   |
| Output Interval      | OUTINT| The output interval (every OUTINT steps, write out CSV)   |
| Restart Output Interval   |RESOUTINT| The restart output interval (every RESOUTINT, write out restart files) <- Joe's addition   |
| Maximum Displacement       | MAXUT | The maximum velocity of particles (particles faster than this may be frozen if this is turned on)   |
| Fracture After Time | FRACTIME | Fracture is permitted after this time. |
| Bed Stiffness Constant | BedIntConst | The stiffness constant of the bed |
| Bed Z Only | BedZOnly | Whether to consider only the z component of bed interaction (rather than normal) |

### mass3.dat ####

mass3.dat is the input configuration, in format:

x, y, surface, base, bed. friction (Newton seconds per metre)  

x and y must start at zero, must have a (0,0) corner.  

output transformation matrix which takes from Elmer domain to HiDEM domain.  

make sure the bed is buffered beyond the edge of the ice, and define these regions by setting surf and base equal to bed.  


### Internal variables ###

| VarName | Description |
| --------- | ------- |
| YN | the number of partitions in the Y direction. |
| NTOT | the total number of partitions in the model |
| GRID | The grid size of the mass3.dat data, make sure it matches |
| SCL | diameter of each particle |
| RESTART | 1 = true |
| VDP | drag coeff |
| UT  | current displacement |
| UTM | previous displacement |
| UTP | next displacement |
| FRZ (FRY,FRX) | contact forces (particle-particle, particle-bed) |
| BOYZ, BOYY | buoyant forces |
| R | elastic forces |
| WSX,WSY | wall contact forces? |
| MN | mass of particles |
| MFIL | mass of particles - per particle |
| JS | moment of rotational inertia |
| NAN | list of connections between particles e.g. NAN(1:2,1) lists the two particle numbers which make up connection 1   |
| NRXF | initial position of this partition's particles   |
| NRXFL,... | initial position of particles in the partition to the left   |


## OUTPUT - jyr files and STR files ##

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


## Processing ##

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

## PARAVIEW ##

First simply load csv, with blank space delimiter, merge delimiters, has no header  
Then apply table to points filter  
Then ensure enable ospray and shadows, and pick point size that makes them touch  

For the sea: add a box, set size and centrepoints  

For the bed: get bed.csv, read it in same as other CSVs  
then apply Delaunay2D filter.  



## Changes made by Jan since commit ##

Got rid of: UTPP, VELO, D, MAXDT,XIND,YIND  
DAMP1,DAMP2,DRAG,BEDZONLY,MAXUT  <- but did I add these? yes  

The change means BEDZONLY=TRUE is implicit  

Using lbound,rbound fib00,tbed files, but NOT dtmax  

Changed ice density  

lots of things involving SCL seem to have changed power  

some change to setting porisity/predamage - no longer considers proximity to bed/base?  

EFLOAD1,2 -> EFFLOAD: Only EFFLOAD exists, and its passed more stuff  

Seems to get rid of UTP prediction, damping  


## TO DO ##

Allow user to specify restart prefix (i.e. look for files [prefix]_REST0*, etc)

