# Fermi

Fermi is a parallel Finite Element (FE) code to model the neutron diffusion in nuclear reactors. The code was design to
solve the steady-state, and also, the time-dependent multi-group equations. Fermi uses PETSc library to assemble and
solve the system of equations in distributed architectures.

# Build, Compile and Test

```shell
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

# Run with Docker

Build the docker image:

```bash
docker build -t fermi .
```

Run the container and connect with the local volume. This is useful for local development:

```bash
docker run -v $PWD:/fermi -it --rm fermi /bin/bash
```

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

## Push images to the GitHub Container Registry

Login, tag and push to ghcr:

```bash
echo $GITHUB_PAT | docker login ghcr.io -u USERNAME --password-stdin
docker tag IMAGE_ID ghcr.io/USERNAME/IMAGE_NAME:TAG
docker push ghcr.io/USERNAME/IMAGE_NAME:TAG
```

## Contact me

Guido Giuntoli - [gagiuntoli@gmail.com]