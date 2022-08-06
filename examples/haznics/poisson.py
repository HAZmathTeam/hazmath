"""
\file examples/haznics/poisson.py
Created by Ana Budisa on 2021-08-30.
Example copyright by bitbucket.org/fenics-apps/cbc.block.

We solve

    u + div(grad(u)) = f

on a unit square domain. The preconditioner for this system is AMG.
Outer solver is Conjugate Gradients.
"""
from block.iterative import ConjGrad
from block.algebraic.hazmath import AMG
from dolfin import *

# Function spaces, elements
mesh = UnitCubeMesh(16, 16, 16)

V = FunctionSpace(mesh, "CG", 1)

f = Expression("sin(pi*x[0])", degree=2)
u, v = TrialFunction(V), TestFunction(V)

a = u*v*dx + dot(grad(u), grad(v))*dx
L = f*v*dx

A = assemble(a)
b = assemble(L)

# here we use hazmath AMG:
B = AMG(A,parameters={"max_levels":10,"print_level":10,"coarse_solver":32})

Ainv = ConjGrad(A, precond=B, tolerance=1e-10, show=2)

# solve
x = Ainv*b

u = Function(V)
u.vector()[:] = x[:]

# default solver in Dolfin 
u2 = Function(V)
solve(A, u2.vector(), b)

print("Max differences between the two solutions: ", (u.vector()-u2.vector()).max())


