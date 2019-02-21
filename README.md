# The Helsinki Discrete Element Model (HiDEM) #

Developers: Jan Åström & Joe Todd

This is HiDEM, the Helsinki Discrete Element Model, a particle model for simulating elastic behaviour, fracture and calving at marine terminating glaciers. Due to its computational demands, this code is designed to be run on parallel HPC facilities.

<p align="center">
<img src="https://i.imgur.com/LqkPXJe.png" width="400" alt="Rink Glacier with Melange">
</p>

See the model in action [here](https://youtu.be/owUrbm_3zi0) and [here](https://youtu.be/wXMK0e2isM4).

The model physics is described in detail in [Åström et al. 2013](https://www.the-cryosphere.net/7/1591/2013/).

## Compilation ##

Configuration, compilation and installation is handled by CMake. Example installation scripts for Cray and Ubuntu systems are located in scripts/compilation. These scripts invoke cmake with toolchain files located in scripts/toolchains, which set default compilers etc.

By default, compilation produces a single binary 'HiDEM' in the top level of the 'build' directory.

If you generate your own toolchain/compilation scripts for different systems, please get in touch or make a pull request!

## Input Files ##

On starting, HiDEM reads the name of an input file from HIDEM\_STARTINFO (e.g. 'inp.dat').

inp.dat - defines various parameters for the simulation

The user must also supply the initial geometry in gridded format: X,Y,SURF,BASE,BED,FRICTION

This initial geometry file is specified in the input file as 'Geometry File':

`Geometry File = "mass3.dat"`

Head of a sample mass3.dat file:

```
15851
0.000000000000000000e+00	0.000000000000000000e+00	9.248271630519999462e+02	9.248271630519999462e+02	9.248271630519999462e+02	1.833796257881020755e+07  
5.000000000000000000e+01	0.000000000000000000e+00	9.248386672000000317e+02	9.248386672000000317e+02	9.248386672000000317e+02	2.784068307483682781e+07  
1.000000000000000000e+02	0.000000000000000000e+00	9.248503686919999609e+02	9.248503686919999609e+02	9.248503686919999609e+02	4.768270480512356758e+07
```

See more information in section 'Geometry File' 

Note: To avoid spurious single particles appearing at the edge of the domain, interpolation of 
particle locations is not permitted where not all 4 supporting points have valid surface and bed
values. To permit this, set `	    Strict Domain Interpolation = False` in inp.dat.

## Important Points ##

#### Typical simulation domains: ####

* 100m x 100m * 100m - typical particle size 2-10m
* 100km x 100km * 1km - typical particle size 30-100m

#### Boundary conditions: ####

* Domain in x-y plane
* Main flow in y-direction
* Forces acting on edge should be defined

#### Particle size: ####

* Should have at least 10 particles in thickness direction (preferably 20-30).
* Maximum ~10 million particles
* Max timestep size scales with particle size:
  * Too long timestep => instability
  * Too short timestep => long execution times
  

## Getting data ready ##

Replace no data with zero

Adjust the waterline (lowest point in domain must be non-negative)

Choose the number of cores

Do something to the restitution coefficient

Point to the input file (e.g. testinp.dat) using HIDEM_STARTINFO

## Test Case ##

A test case is provided in the `test` directory. This is a simple rectangular domain with a sloping upper surface. This simulation will produce fracturing behaviour within a 6 hour simulation on 130 CPUs. The user should edit `example.job` to provide a valid PBS budget account, and the correct number of nodes and cores (currently 6, 130) depending on system architecture (cores per node).

The behaviour of the test case can be changed by modifying `Max Load` and `Friction Scale` in inp.dat.

## Running the model ##

The model runs in parallel using MPI, so the simulation should be started using `mpirun` or `aprun`:

```
mpirun -n 70 HiDEM
```

An example PBS job script is provided in example.job

The user may specify a run name which is preprended on all output files:

`Run Name = Example_Simulation`

#### Restarting ####

It is possible to restart one run from another:

`Restart = 1`

If the run name has changed, specify the run from which to restart:

```
Run Name = Example_Simulation_Part2
Restart from Run Name = Example_Simulation
```

## Melange ##

It is possible to use HiDEM to investigate glacier/melange interaction. Melange is generated by running an initial 'melange simulation' until the domain calves to produce icebergs & brash. Then a subsequent simulation can read in this broken ice from the 'melange simulation'. This is controlled by the keyword:

```
Melange Run Name = "Your_previous_simulation_run_name"
```

The 'melange simulation' setup could be identical to the standard calving simulation, or one could choose to increase the porosity or lower the maximum load to encourage breakup. It is not necessary for the melange simulation and the subsequent simulation to run on the same number of cores, but it is expected that the particle size (SCL) would be the same.

Note: unlike a regular simulation restart, the melange is assumed to be 'at rest' when read into the new sim.

## Troubleshooting ##

If the simulation [explodes](https://www.youtube.com/watch?v=LZIixgvlF8U) (particles moving too far):

* Set a higher `Translational Damping` or `Rotational Damping` in the input file - multiply by two? Very high values for damping might be e.g. 16E4 and 8E4

* Set `Fracture After Time = 50.0` or so, to allow the simulation to settle before permitting fracture

* Change `friction scale` in input - this scales the bed friction

* Change the `timestep`

* Smooth gradients in basal friction or geometry (in `geometry file`)

* Kill all motion every x timesteps

## Files, Variables and Parameters ##

### Structure of the repository ###

* `src` - all the code required to compile the HiDEM binary.
* Within the `scripts` directory:
  - `analysis` - various fortran and python scripts for post-processing
  - `compilation` - example compilation scripts
  - `job_scripts` - example PBS scripts
  - `io` - fortran programs for modifying input
  - `paraview` - macros for working with HiDEM data in Paraview

### HiDEM Source Code ###

These are part of the compiled code:

| File | Description |
| --------- | ----------- |
|wave2.f90 | The main program |
|dist.f | Efficiently finding neighbouring particles which may interact |
|circ.f | Confirms which particles are in contact and computes forces |
|effload.f | Computes elastic forces in particle beams |
|amat.f | Called by effload, does integration |
|tmat.f, ttmat.f | Rotation matrices |
|kmat.f | Stiffness matrix computation |
|glas.f90 | Computing the FCC lattice, dense packing |
|ranmar.f | Random number generator |
|dt.f90 | Called by glas.f90, finds and write connections to FSfiles |

### Other Files ###

ave*.f - these compute averages of various outputs

rc2.f90, rc3.f90 - compute the calved size distrib  


### inp.dat parameters ###
| Parameter | VarName | Description | Default (SI units) |
| --------- | ------- | ----------- |---------|
| Run Name  | RUNNAME | The name of the simulation (prepended to output filenames) |  |
| Work Directory | wrkdir | The subdirectory in which to store working files (directory must exist!) | . |
| Results Directory | resdir | The subdirectory in which to store result files (directory must exist!) | . |
| Geometry File | geomfile | The name of the file which defines the geometry (see above) | False |
| Melange Run Name  | MelRunName | The name of the simulation from which to read melange. |  |
| Density | RHO | The density of the material (ice) | 900 |
| Water Density | RHOW | The density of the water in which the ice floats | 1030 |
| Gravity | GRAV | Magnitude of gravity (positive!) | 9.81 |
| Backwall Pressure | PRESS | optional backwall/water pressure (not yet properly implemented!)| 0 |
| Submarine Melt    | MELT  | optional basal melt rate passed to fibg3 for altering domain shape   | 0 |
| UC                | UC    | optional frontal melt rate   | 0 |
| Timestep          | DT    | timestep size   | 1.0e-4 |
| Width             | S     | beam width (relative to unit particle)   | 0.7 |
| Youngs Modulus    | EF0   | particle bond young's mod, describes interaction between connected particles   | 1.0e+9 |
| Size              | LS    | Something to do with the fcc lattice 'box size'   | 100 |
| Domain Inclination | SUB   | Angle of the domain vs gravity vector, not used    | 0 |
| Water Line   | WL    | Sea level for buoyancy calc   | |
| Grounding Line  | GL    | 'grounding line' - not used   | -100 |
| Shear Line   | SLIN  | A distance below the surface where all bonds are broken?   | 2000 |
| No Timesteps   | STEPS0| The number of steps   | |
| Max Load     | MLOAD | Maximum load on bond - bonds break beyond this   | 0.0002 |
| Friction Scale    | FRIC  | Scale factor for friction input   | 1 |
| Restart     | REST  | 1 = restart from prev, 0 = new run   | 0 |
| Restart From Run Name    | restname  | Specifies the name of the run from which to restart (defaults to same as 'Run Name')| |
| Scale         | SCL   | Scale factor for particle and beam size   | |
| Grid        | GRID  | Resolution of mass3.dat input grid   | |
| Porosity         | POR   | The proportion of initially broken bonds   | 0.1 |
| Random Seed       | SEEDI | Seed for random number generator   |  11695378 |
| Translational Damping       | DAMP1 | The damping coefficient for translation   | 1e4 |
| Rotational Damping       | DAMP2 | The damping coefficient for rotation   | 1e4 |
| Air Drag Coefficient        | DRAG_AIR  | The drag coefficient in air   | 1e1 |
| Water Drag Coefficient        | DRAG_WATER  | The drag coefficient in water   | 1e1 |
| Drag Coefficient        | DRAG_WATER, DRAG_AIR  | The drag coefficient in both air & water (alternative to previous 2)  | 1e1 |
| Viscous Distance        | ViscDist  | The SCLed particle proximity for viscous interaction   | 4e-2 |
| Viscous Force        | ViscForce  | The strength of viscous particle interaction | 1e4 |
| Output Interval      | OUTINT| The output interval (every OUTINT steps, write out CSV)   | 20000 |
| Restart Output Interval   |RESOUTINT| The restart output interval (every RESOUTINT, write out restart files) <- Joe's addition   | 20000 |
| Maximum Displacement       | MAXUT | The maximum displacement of particles - default 1.0e6 metres (particles further than this are frozen)   | 1e6 |
| Fracture After Time | FRACTIME | Fracture is permitted after this time (in s). | 40 |
| Bed Stiffness Constant | BedIntConst | The stiffness constant of the bed | 1e8 |
| Bed Damping Factor | BedDampFactor | Alters the damping of bed interaction (1.0 = critically damped) | 1.0 |
| Bed Z Only | BedZOnly | Whether to consider only the z component of bed interaction (rather than normal) | True |
| Strict Domain Interpolation | StrictDomain | Determines limit of interpolation w.r.t geometry input file. See note above | True |
| CSV Output | CSVOutput | If true, produce output in .csv format rather than binary (uses more disk space) | False |
| Double Precision Output | DoublePrec | If true, output data will be Float64 (as opposed to Float32, doubles output filesize) | False | 
| Geometry File Has Mask | GeomMasked | Specifies whether the geometry file includes a mask column (required for 'Fixed Lateral Margins' | False |
| Fixed Lateral Margins | FixLat | If true, particles near the lateral margins are not permitted to move in XY plane | False |
| Fixed Inflow Margin | FixBack | If true, particles near the inflow margin are not permitted to move in XY plane | True |


### Geometry file ####

The geometry file contains the input configuration, in format:

x, y, surface, base, bed, friction, geom_mask (optional)

x and y must be gridded and start at zero, must have a (0,0) corner.  

Friction has units of Newton seconds per metre

output transformation matrix which takes from Elmer domain to HiDEM domain.  

make sure the bed is buffered beyond the edge of the ice, and define these regions by setting surf and base equal to bed.  

User may optionally specify a 'geom_mask' column in geometry input file, which tells the model which regions are ice (=1), fjord (=2), bedrock(=0). This is required for imposing lateral boundary conditions (Fixed Lateral Margins).

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
| NRXF%M,%F | initial position of the particles  |


## Output - JYR and STR files ##

By default the model produces particle information in .vtu format (readable in Paraview (v5.5) - use HiDEM_load.py macro) and bond strain information in binary format (readable by the python script rh.py). For CSV output, use the 'CSV Output' option in the inp.dat (beware this produces much larger files which make visualisation time consuming).

JYR files list the position of all particles in x,y,z, every 2 seconds.  
Read this in paraview quite easily.  

STR file list the midpoint position and strain of each *bond* between two particles, for each node connection in initial geometry (including broken bonds)  

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

## Paraview ##

First simply load csv, with blank space delimiter, merge delimiters, has no header  
Then apply table to points filter  
Then ensure enable ospray and shadows, and pick point size that makes them touch  

For the sea: add a box, set size and centrepoints  

For the bed: get bed.csv, read it in same as other CSVs  
then apply Delaunay2D filter.  


## TO DO ##

 * BC strategy

 * Translate & Rotate input

 ## Notes ##

Domain needs to be orientated in XY because:

Boundary conditions are applied in WSY (, WSX?) components (need to compute normal? Ask 

