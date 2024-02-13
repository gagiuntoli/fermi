# Use an existing base image that includes MPI and PETSc dependencies
FROM ubuntu:latest

# Install additional dependencies if needed
RUN apt-get update && apt-get install -y petsc-dev cmake libopenmpi-dev

ENV PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real

# Copy project files
COPY . /fermi

# Set working directory
WORKDIR /fermi

# Build and test your project
RUN cmake . && make && make test