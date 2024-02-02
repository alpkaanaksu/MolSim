# PSE Molecular Dynamics WS23/24 - Group F

## Members

- Alp Kaan Aksu
- Berke Saylan
- Feryal Ezgi Aşkın

## Code

- Link:     https://github.com/alpkaanaksu/MolSim
- Branch:   master
- Revision: a67015b2fa1326f47820aab6ed309b2bfb4675a3

**Works with:**

- **Compiler:** clang 15.0.0
- CMake 3.27.7
- GNU Make 3.81

*Other compilers / versions would probably also work but we only tested with these.*

## Compiling and running the program

- You need `xerces-c` and `boost` (`program_options` and `filesystem`) to compile the program.
- If you need parallelism, you need install `OpenMP`.

```bash
mkdir build # if it does not exist
cd build
```

```bash
ccmake ..
```

*CMake will automatically fetch some files needed for additional libraries like `googletest`, `spdlog`
and `nlohmann_json`*

```bash
make
```

```bash
./MolSim ../input/rayleigh-taylor-instability/3d.json
```

---

Set the environment variable `OMP_NUM_THREADS` to the number of threads you want to use for parallelism.

### Doxygen

```bash
make doc_doxygen
```

An online version of the documentation can be found [here](https://alpkaanaksu.github.io/MolSim/).

## Simulation of a Membrane

### Generator class

In the Generator class we created a new method named `membrane`, which is responsible for generating particles in a
cuboid structure. The membrane method initializes particles with specific attributes, including position, velocity,
mass, and interaction parameters. Neighbors are determined based on the cuboid structure, distinguishing between direct
and diagonal neighbors. Furthermore in this method the particles are determined, which are pulled by an external force
for a time period.

In the `membrane` method a nested loop structure is employed to iterate over neighboring particles in a 2D grid. The
outer loop, for the x-axis iteration, and the inner loop, z-axis iteration, iterate through the possible relative
positions of neighboring particles. The conditional statement within the loops ensures that the current particle is
skipped (not considered as its own neighbor). The variables nx and nz are then used to calculate the coordinates of the
potential neighbor. Subsequent checks verify whether the calculated coordinates fall within the valid range of the grid.
If the neighbor is within the valid range, its index is determined, and further checks are performed to distinguish
between direct and diagonal neighbors. Depending on the conditions, the particles are set as the direct and diagonal
neighbors of the respective iterated particle, by making use of the id attribute of the particles, which is determined
already in the cuboid generation step uniquely.

Here we hardcoded the iteration for only z- and x-axis, because in the last meeting it was said that we were allowed to
hardcode the membrane simulation according to the given parameters, as this will be the last iteration of the practical
course, thus there will be no use case for a more generalized membrane simulation that works according to the given axis
etc. parameters.

### Simulation class and Input file

In this class we define the external pulling force and apply it to the particles, which have the attribute pulled set to
true, for a hardcoded time period. Then we check if the simulation is a membrane one or not, depending on that we apply
the `applyToAllPairsOnceMembrane` or `applyToAllPairsOnce` methods for the correct force calculation between unique
particle pairs in the particle container.

We have made small changes in the input file for a better simulation view. Initially we have created a simulation with
the given parameters and it functioned correctly, however the simulation view was not smooth, since the gravity is
applied along the y-axis, but pulling force is applied along the z-axis. The outcome was that after some time point the
membrane folded, which in our opinion didn't provide a perfect membrane simulation view. Therefore we set the pulling
force the same but now along the y-axis and we also set the cuboid accordingly (lying on a x-z plane flatly), plus
particles to which the pulling force apply are set correctly with hardcoding in the `membrane` method. These input
parameters result in a simulation of a "napkin" being pulled up and then falling down, which in our opinion represented
a fully functioning membrane simulation better than the other simulation.

### LinkedCellParticleContainer class

The `applyToAllPairsOnceMembrane` method is designed to handle particle interactions within a 3D grid structure. It
works analogously `applyToAllPairsOnce` to iterate through all the particle pairs in an efficient manner utilizing the
linked-cell data structure optimization.
This implementation however makes a distinction for particle pairs, it applies the Lennard-Jones potential truncated at
a distance of 2^(1/6) * σ (which is achieved by a newly implemented `membraneModel` in the `Model` header) for particle
pairs that are not neighbors of each other. If the condition of being a direct or diagonal neighbors is fulfilled within
particle pairs the respective harmonic potential calculation is carried out for the force calculation between
neighboring particles.

The harmonic potential is employed for both direct and diagonal neighbor interactions, adapting the force calculations
based on the specified stiffness constant and average bond length.

## Parallelism

During this practical course, we executed our program, often waiting for hours or even overnight for its completion. This hands-on experience, more than any other exercise, highlighted the significance of parallelism in enhancing efficiency and performance.

One of our parallelization ideas was to parallelize the cell and particle iteration 
in the ‘applyAllPairsOnce’ method for pairwise force calculation and the 
results showed that it enhanced the performance significantly. In the last worksheet, we saw that the force calculation is the most time-consuming part of the simulation. Therefore, we decided to parallelize the force calculation part of the simulation, force calculation also does not use any dynamic values that can be changed by other threads during the force calculation, this makes it really easy to avoid race conditions. We also use dynamic scheduling to ensure that the work is distributed evenly to threads, maximazing the performance gain.

We also tried other approaches for parallelization, e.g. parallelizing also 
‘applyToAll’ method etc. However all of these approaches caused one of the 
two following problems: either the overhead resulting from threads was so 
high (especially because of many critical regions) or the result wasn’t 
deterministic, as critical regioning wasn’t applied correctly. In all of the cases 
this dilemma forced us to actually not employ the respective parallelization 
approach. We also tried to parallelize the visit of neighbor cells instead of the processing of the cells, but this also didn't work out, because the overhead of creating threads was too high, it made our program even slower, but we still kept the code in the repository, because we had to implement 2 strategies.

You can specify the parallelization method in the input file by setting `parallelization_strategy` to `"cells"`, `"neighbors"` or `"none"`.

Our application performs as expected both with and without OpenMP, a validation confirmed by our suite of tests.

## Rayleigh-Taylor Instability in 3D

Essentially our implementation for the 2D case was designed to cover also the 3D case. However, we realized that within
halo particle creation we were handling only the case where only one axis is periodic. Therefore we had to rewrite the
logic for halo particle copying algorithm.
Now we are using 3 methods named `handleBoundariesOneAxis`, `handleBoundariesTwoAxes` and `handleBoundariesThreeAxes` to
handle the periodicity of the simulation. These methods are called in the `updateHaloCells` method, which is responsible
for updating the halo cells of the particle container. The `updateHaloCells` method is called in the `update` method of
the `LinkedCellParticleContainer` class, which is responsible for updating the particle container in each time step of
the simulation.
The implementation of `handleBoundariesOneAxis` is almost the same as the previous auxiliary
functions `lowerBoundaryToUpperHaloOneAxis` and `UpperBoundaryToLowerHaloOneAxis` from iteration 4. This method is
responsible for copying particles from a boundary cell to the corresponding halo cell along one specified axis (X, Y, or
Z). In order to handle corner border cell cases for multiple periodic axes, we have implemented
the `handleBoundariesTwoAxes` and `handleBoundariesThreeAxes` methods. These methods extend the functionality to handle
periodic conditions along two/three axes simultaneously, providing flexibility for 2D/3D scenarios with different
periodicity. Their implementation logic is the same as the auxiliary function handling only one axis, but they are
adjusted to function in multiple axes. `updateHaloCells` method calls these methods on boundary cells, if respective
periodicity conditions are fulfilled, to ensure the periodicity of the simulation, maintaining accurate boundary
conditions across all axes.

## Nano-Scale Flow Simulation

We opted for the nano-scale flow simulation in task 4. The addition of a new attribute,`fixed`, to the `Particle` class
helped define and clearly distinguish non-moving objects within our simulation. Since `fixed` operates at the
fundamental level
of the project (`Particle.cpp`), restructuring the existing codebase posed little challenge, given that every simulation
inherently involves particles.

Regarding the new thermostat implementation, the corresponding velocity scaling method is chosen based on whether the
input contains objects with fixed particles, cuboids, etc. If such objects exist, the method `scaleVelocitiesWithAvg` is
called; otherwise, `scaleVelocities` from the last worksheet is employed. This approach simplifies the code and
minimizes the need for modifications to accommodate future thermostats.

The striking part of this simulation was its velocity and density profiles. As the flow is downwards due to the gravity
factor and periodic boundaries, the y-component of particle velocities was of importance.

![Gravity and flow plot](nanoflow/gravity_plot.png)

For particles that are nearest to the boundaries, the downward velocity was generally lower compared to those located in
the middle of the flow. This didn't change for varying gravity values as one would intuitively expect. The interaction
between fixed wall particles and free flowing ones is more significant the closer they are to each-other. The movable
particles, however cannot exert force on the fixed particles. This inability of movable particles to exert
force on the fixed particles prevents a transfer of momentum that could otherwise counteract the localized suppression
of acceleration. As a result, the downward velocity of particles near the boundaries remains subdued, creating a
distinctive pattern in the particle dynamics within the simulation.

![Flow after some time plot](nanoflow/plot.png)

We notice similarities when comparing higher gravity values to the behavior observed after some time has elapsed. The
average fluid velocity is influenced by two primary factors: temperature, which remains relatively constant due to the
thermostat, and gravity. Similar to an object in free fall, the velocity increases over time.

The profiles reminded us of the [no-slip condition](https://en.wikipedia.org/wiki/No-slip_condition), which enforces
that at a solid boundary, a viscous fluid attains zero bulk velocity. This concept has been successful in solving
problems related to bigger flows for the past two centuries. However, when dealing with fluid transport at the micro-
and nanoscale, things get more complex according to
some [research](https://www.beilstein-journals.org/bjnano/articles/12/91). In these tiny systems, there's a possibility
that liquids don't adhere perfectly to solid surfaces, causing what's known as liquid slippage. It was fascinating to
observe how our nano-scale simulation behaves at boundaries.

Placing a fixed sphere helped us observe the turbulence that occurred and had this cool plot:

![Flow after some time plot](nanoflow/plot_obstacle.png)

We also attempted to investigate how different σ values affect the profiles using two different fluids within the same
simulation, similar to the approach in Rayleigh-Taylor instability. However, we didn't observe anything significant in
the profiles. Perhaps adhering to one variable factor at a time would have been a better approach.

## Crystallization of Argon 
We implemented the coolest part of the task with the smoothed LJ, to perform the experiments out of curiosity.

## Contest 2

Instructions are very similar to the first one, just create two new job scripts:

```bash
echo '#!/bin/bash
#SBATCH -J molsimf2
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=28
#SBATCH --mail-type=end
#SBATCH --mail-user=<your email here>
#SBATCH --export=NONE
#SBATCH --time=00:02:00

module load slurm_setup

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
 
./MolSim ../input/performance/contest1.json' > molsimf2.cmd
```

```bash
echo '#!/bin/bash
#SBATCH -J molsimf3
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=56
#SBATCH --mail-type=end
#SBATCH --mail-user=<your email here>
#SBATCH --export=NONE
#SBATCH --time=00:05:00

module load slurm_setup

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
 
./MolSim ../input/performance/contest2.json' > molsimf3.cmd
```

and submit them to cluster: ```sbatch molsimf2.cmd``` and ```sbatch molsimf3.cmd```.

We got the following results:

**2D:**
```
[2024-02-02 01:25:18.750] [info] Time: 4.996
[2024-02-02 01:25:18.750] [info] MUP/s: 2001416.202105
```

**3D:**
```
[2024-02-02 01:26:55.846] [info] Time: 71.911
[2024-02-02 01:26:55.846] [info] MUP/s: 1390590.604828
```
