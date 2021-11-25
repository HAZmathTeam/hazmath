#!/bin/bash
for i in ./stokes ./heat_equation ./amr_grids ./approximation ./eigen ./basic_elliptic ./reaction_diffusion ./elasticity ./solvers ./darcy ./mg_geometric ; do
make -C $i clean ; make -C $i
done
