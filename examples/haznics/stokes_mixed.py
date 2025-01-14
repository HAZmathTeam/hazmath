"""
\file examples/haznics/stokes_mixed.py
Created by Ana Budisa on 2021-08-30.
Example copyright by bitbucket.org/fenics-apps/cbc.block.

This demo solves the steady Stokes equations for a lid driven
cavity.  We solve

    - div(symmgrad(u) - pI) = f
                   - div(u) = 0

with Taylor-Hood elements on a unit square domain.
The algebraic system can be written as

  BB^ AA [u p]^T = BB^ [0 b]^T,

where AA is a 2x2 block system with zero in the (2,2) block

       | A   B |
  AA = |       |,
       | C   0 |

and BB^ approximates the inverse of the block operator

       | A   0 |
  BB = |       |,
       | 0   M |

where A is the Laplace operator and M the inner product in L2.
For A, an AMG preconditioner is used, i.e. A^ = AMG(A).
For M, we use Jacobi preconditioner.
"""
from dolfin import *
from block import *
from block.algebraic.petsc import Jacobi 
from block.algebraic.hazmath import AMG
from block.iterative import MinRes
import haznics

mesh = UnitSquareMesh(32, 32)

P2 = VectorElement("Lagrange", triangle, 2)
P1 = FiniteElement("Lagrange", triangle, 1)
TH = MixedElement([P2, P1])

W = FunctionSpace(mesh, TH)

f = Constant((0., 0.))

u, p = TrialFunctions(W)
v, q = TestFunctions(W)

a = inner(grad(u), grad(v)) * dx \
  - p * div(v) * dx \
  - q * div(u) * dx

b = inner(grad(u), grad(v)) * dx + inner(u, v) * dx \
  + p * q * dx

L = inner(f, v) * dx

bcs = [DirichletBC(W.sub(0), (0., 0.), "on_boundary&&(x[1]<1-DOLFIN_EPS)"),
       DirichletBC(W.sub(0), (1., 0.), "on_boundary&&(x[1]>1-DOLFIN_EPS)")]

# assemble as block matrices
A, _, rhs = block_assemble(a, L, bcs, symmetric=True)
B, _, _ = block_assemble(b, L, bcs, symmetric=True)

# build the preconditioner
params = {
   'AMG_type': haznics.SA_AMG,
   "aggregation_type": haznics.VMB,
   "max_levels":10,"print_level":10,"coarse_solver":32
}

from block.algebraic.hazmath import PETSc_to_dCSRmat
XXX = B[0,0].array()
import numpy as np
xnorm0=np.linalg.norm(XXX - XXX.T)
print("\nxnorm0=",xnorm0)
Ahaz=PETSc_to_dCSRmat(B[0,0])

print('Haz B[0, 0] symmetry', haznics.chk_symmetry(Ahaz))

P = block_mat([[AMG(B[0, 0], parameters=params), 0],
               [          0, Jacobi(B[1, 1])]])

Ainv = MinRes(A, precond=P, relativeconv=True, tolerance=1e-5, show=3)
x = Ainv * rhs

# plotting
V, Q = [sub_space.collapse() for sub_space in W.split()]
u, p = map(Function, [V, Q], x)

print('|u|_0', sqrt(abs(assemble(inner(u, u)*dx))))
print('|p|_0', sqrt(abs(assemble(inner(p, p)*dx))))

File('./results/stokes_mixed_u.pvd') << u
File('./results/stokes_mixed_p.pvd') << p
