# Fermi

Fermi is a parallel Finite Element (FE) code to model the neutron diffusion in
nuclear reactors. The code was design to solve the steady-state, and also, the
time-dependent multi-group equations. Fermi uses PETSc library to assemble and
solve the system of equations in distributed architectures.

Some additional features are expected to be implemented:

- [x] Control rod simulation.
- [ ] Xenon evolution and poisoning.
- [ ] Coupling with other codes capabilities.

# Installation

## Install Dependencies

The only dependency of Fermi is the PETSc library for solving the systems of
equations. For this, download and install [PETSc](www.mcs.anl.gov/petsc) library
and declare the environmental variables:  

```bash
   export PETSC_DIR=<path>
   export PETSC_ARCH=<path>
```

## Build and compile

The values should be the same as those set during the compilation of the
libraries. In the main folder:

```bash
   mkdir build
   cd build
   cmake ..
   make
```
If you want a `Release` version with the optimized compilations flag `-O3` do

```bash
   cmake -DCMAKE_BUILD_TYPE=Release ..
```

in the other case the compilation will not be optimized. We suggest to define
different `build` directories if the user want different versions of the code
for debugging or compare results in case of code improvements.

# Input Definition

```bash
$Mesh
    mesh_file   fuel_1.msh
$EndMesh

$Mode
  timedep QSTATIC
  p0 2.5e6  
  t0 0.0
  tf 0.15
  dt 0.05
$EndMode

$Xs
  egn 1
  #          F D   XA  nXF  eXF CHI
  "FUEL_Z1"  0 1.5 0.2   0.5  0.4 1.0
  "FLUID_Z1" 0 1.5 0.005 0.0  0.0 1.0
$EndXs

$Boundary
  "TOP_SURF"    1 0
  "BOT_SURF"    2 0
  "LAT_SURF"    3 1
$EndBoundary

$Output
  # Power in physical entities
  kind 2
  file pow_phys.dat
  nphy 2
  "FUEL_Z1" "FLUID_Z1"
$EndOutput
```

## Docker

Build the docker image:

```bash
docker build -t fermi .
```

## Future Workflow

- [ ] Testing
- [ ] Parallelization and performance evaluation
- [ ] Benchmarking
- [ ] Documentation

## Contributing

Please read
[CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for
details on our code of conduct, and the process for submitting pull requests to
us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available,
see the [tags on this repository](https://github.com/your/project/tags).

## Contact

Guido Giuntoli - [gagiuntoli@gmail.com]
