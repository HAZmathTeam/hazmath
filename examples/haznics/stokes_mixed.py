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
For M, we use Jacobi.
"""
from dolfin import *
from block import *
from block.algebraic.petsc import Jacobi
from block.algebraic.hazmath import AMG
from block.iterative import MinRes
import haznics

mesh = UnitSquareMesh(32, 32)

V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
W = (V, Q)

f = Constant((0., 0.))

u, p = map(TrialFunction, W)
v, q = map(TestFunction, W)

a = [[inner(grad(u), grad(v)) * dx, -p * div(v) * dx],
     [-q * div(u) * dx,                            0]]

b = [[inner(grad(u), grad(v)) * dx + inner(u, v) * dx, 0],
     [0                                     , p * q * dx]]

L = [inner(f, v) * dx, 0]

Vbcs = [DirichletBC(V, (0., 0.), "on_boundary&&(x[1]<1-DOLFIN_EPS)"),
        DirichletBC(V, (1., 0.), "on_boundary&&(x[1]>1-DOLFIN_EPS)")]
bcs = block_bc([Vbcs, []], symmetric=True)

# assemble as block matrices
A, B, rhs = map(block_assemble, (a, b, L))
bcs.apply(A).apply(rhs)
bcs.apply(B)

# build the preconditioner
params = {
    'AMG_type': haznics.SA_AMG,
    "aggregation_type": haznics.VMB,
    "max_levels":10,"print_level":10,"coarse_solver":32
}
P = block_mat([[AMG(B[0, 0], parameters=params), 0],
               [          0, Jacobi(B[1, 1])]])

Ainv = MinRes(A, precond=P, relativeconv=True, tolerance=1e-5, show=3)
x = Ainv * rhs

from block.algebraic.hazmath import PETSc_to_dCSRmat
XXX = B[0,0].array()
import numpy as np
xnorm0=np.linalg.norm(XXX - XXX.T)
print("\nxnorm0=",xnorm0)
Ahaz=PETSc_to_dCSRmat(B[0,0])
# THIS MATRIX SEEMS TO BE NON_SYMMETRIC! haznics.dcsr_write_dcoo("AAA",Ahaz)
haznics.chk_symmetry(Ahaz)

# plotting
u, p = Function(V, x[0]), Function(Q, x[1])
File('./results/stokes_mixed_u.pvd') << u
File('./results/stokes_mixed_p.pvd') << p
