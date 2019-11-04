# README

`Fermi` is a code made for studying the physical behavior of nuclear reactors and to perform a design of them.

It solves an eigenvalue problem corresponding to the steady state multi-group diffusion equation, which is an approximation of the Boltzmann Transport Equation. This approximation is only valid in materials were the relation between absorption and scattering collision of neutrons is small.

Macroscopic cross sections should be provided in order to solve the equations in the domain. They may become from a lattice cell homogenization, and possibly, a condensation to few groups in order to decrease the computational effort.

To perform calculations and solve the problem, `Fermi` take advantage of the finite element method to discretize the equations. With this technique, the code is capable of solving the problem over unstructured 3d meshes.

The design of `Fermi` was aim to be a very simple and easy understanding code. The input is a 1d, 2d or 3d mesh that should be generated with `gmsh` code and an ASCII file with a particular format which contains information about the cross sections of each material that exist in the domain. The output is a `VTK` file which contains the solution of the problem (the scalar flux) and the a file containing the eigenvalue of multiplication factor "Keff" of the problem.

## Instalation

## `PETSC` libraries

Download and install [PETSc](www.mcs.anl.gov/petsc) library and declare the
environmental variables:

```bash
   export PETSC_DIR=<path>
   export PETSC_ARCH=<path>
```
The values should be the same as those set during the compilation of the
libraries. In `fermi`'s main folder do:

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

in the other case the compilation will not be optimized. We suggest to define different `build` directories if the user want different versions of the code for debugging or compare results in case of code improvements.

## `Input Structure`:

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
## Future Workflow

- [ ] Testing
- [ ] Parallelization and performance evaluation
- [ ] Benchmarking
- [ ] Documentation

## Contact

Guido Giuntoli - [gagiuntoli@gmail.com]
