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
For M, we use SOR method, i.e. M^ = SOR(M).
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
A, rhs = block_assemble(a, L, bcs)
B, _ = block_assemble(b, L, bcs)

# build the preconditioner
params = {
    'AMG_type': haznics.SA_AMG,
    "aggregation_type": haznics.VMB,
}
P = block_mat([[AMG(B[0, 0], parameters=params), 0],
               [          0, Jacobi(B[1, 1])]])

# We don't want to solve too precisely since we have not accounted 
# for the constant pressure nullspace 
Ainv = MinRes(A, precond=P, relativeconv=True, tolerance=1e-5, show=3)
x = Ainv * rhs

# plotting
# V, Q = [sub_space.collapse() for sub_space in W.split()]
# u, p = list(map(Function, [V, Q], x))
# plot(u)
# plot(p)
