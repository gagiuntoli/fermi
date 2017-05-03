#!/bin/bash

gmsh -3 iaea.geo
gmsh2metis iaea.msh 3
mpmetis iaea_vol.metismesh $1
