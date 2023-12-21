# PSE Molecular Dynamics WS23/24 - Group F

## Members
- Alp Kaan Aksu
- Berke Saylan
- Feryal Ezgi Aşkın

## Code
- Link:     https://github.com/alpkaanaksu/MolSim
- Branch:   master
- Revision: TODO

**Works with:**
- **Compiler:** clang 15.0.0
- CMake 3.27.7
- GNU Make 3.81

*Other compilers / versions would probably also work but we only tested with these.*


## Compiling and running the program
- You need `xerces-c` and `boost` (`program_options` and `filesystem`) to compile the program.

```bash
mkdir build # if it does not exist
cd build
```

```bash
ccmake ..
```

*CMake will automatically fetch some files needed for additional libraries like `googletest`, `spdlog` and `nlohmann_json`*

```bash
make
```

```bash
./MolSim ../input/rayleigh-taylor-instability/rayleigh-taylor.json
```


The last line starts the program to run the simulation of the Rayleigh-Taylor instability. You can also run a smaller version of the simulation with `quick.json`

---

### Doxygen

```bash
make doc_doxygen
```

An online version of the documentation can be found [here](https://alpkaanaksu.github.io/MolSim/).

## Thermostat

TODO

## Periodic Boundaries

TODO

## Gravitational Force
We added a new static method called `verticalGravityForce` to our `Model` class. Given `m`, this method calculates the gravitational force for a particle. We add the resulting vector to the force vector of the particles in each step, together with the Lennard-Jones force. 

## Lorentz-Berthelot Mixing Rules
We added two new attributes to the Particle class: `epsilon` and `sigma`, they used to be global simulation parameters, but now they are specific particles. This makes it possible to simulate substances of different types.

We made a small change in the `force` function of the `lennardJonesModel` (see `Model.h`), before we do the usual calculation.

## Rayleigh-Taylor Instability

We implemented the Rayleigh-Taylor instability as described in the exercise sheet. Two cuboids represent the two fluids, one is heavier than the other. The heavier fluid is on top of the lighter one. The simulation starts with a small perturbation in the interface between the two fluids. The heavier fluid starts to fall down and the lighter fluid starts to rise up. The interface between the two fluids becomes more and more unstable and the fluids start to mix.

See ...

## Checkpointing

TODO

## Falling Drop

TODO

## Performance / Profiling

TODO