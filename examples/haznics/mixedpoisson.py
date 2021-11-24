"""
\file examples/haznics/mixedpoisson.py
Created by Ana Budisa on 2021-08-30.
Example copyright by bitbucket.org/fenics-apps/cbc.block.

We solve

    div(sigma) = f
         sigma = - K*grad(u)

on a unit square domain. The algebraic system to be solved can be written as

  BB^ AA [sigma u]^T = BB^ [0 b]^T,

where AA is a 2x2 block system with zero in the (2,2) block

       | A   B |
  AA = |       |,
       | C   0 |

and BB^ approximates the inverse of the block operator

       | A   0 |
  BB = |       |,
       | 0   L |

where L is the Laplace operator. Since the DG(0) approximation of L is zero, we
calculate it instead as L = C*B. For A, an ML multilevel preconditioner is used:

  A^ = ML(A)

For L, we create a composite operator so that

  x = L*v ==> w = B*v; x = C*w,

so we use an inner iterative solver:

  L^ = Richardson(L, precond=0.5, iter=40).
"""
from __future__ import division
from __future__ import print_function

from block import *
from block.iterative import MinRes, Richardson
from block.algebraic.hazmath import AMG
from dolfin import *

# Create mesh
mesh = UnitSquareMesh(32, 32)

# Define function spaces
BDM = FunctionSpace(mesh, "BDM", 1)
DG = FunctionSpace(mesh, "DG", 0)

# Define trial and test functions
tau, sigma = TestFunction(BDM), TrialFunction(BDM)
v, u = TestFunction(DG),  TrialFunction(DG)


# Define essential boundary
def boundary(x):
    return near(x[1], 0.0) or near(x[1], 1.0)


# Define function G such that G \cdot n = g
class BoundarySource(UserExpression):
    def __init__(self, mesh):
        super().__init__(self)
        self.mesh = mesh

    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        g = sin(5*x[0])
        values[0] = g*n[0]
        values[1] = g*n[1]

    def value_shape(self):
        return (2,)


G = BoundarySource(mesh)

# Define the blockwise boundary conditions -- a Dirichlet condition on the
# first block, and no conditions on the second block.
bcs_BDM = [DirichletBC(BDM, G, boundary)]
bcs = block_bc([bcs_BDM, None], True)

# Define source function
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=4)

# Define variational forms
a11 = dot(sigma, tau) * dx
a12 = div(tau) * u * dx
a21 = div(sigma) * v *dx
L2  = - f * v * dx

AA = block_assemble([[a11, a12],
                     [a21,  0 ]])

bb = block_assemble([0, L2])

bcs.apply(AA).apply(bb)

# Extract the individual submatrices
[[A, B],
 [C, _]] = AA

# Use multilevel preconditioner from hazmath (UA-AMG) for A
Ap = AMG(A)

# Create an approximate inverse of L=C*B using inner Richardson iterations
L = C*B
Lp = Richardson(L, precond=0.5, iter=40, name='L^')

# Define the block preconditioner
AAp = block_mat([[Ap, 0],
                 [0,  Lp]])

# Use MinRes as outer solver
AAinv = MinRes(AA, precond=AAp, show=2, name='AA^')

# =====================
# Solve the system
Sigma, U = AAinv * bb
# =====================

# Check that the norms of solutions are as expected
print(('norm Sigma:', Sigma.norm('l2')))
print(('norm U    :', U.norm('l2')))

if abs(1.213-Sigma.norm('l2')) > 1e-3 or abs(6.716-U.norm('l2')) > 1e-3:
    raise RuntimeError("Wrong value in norms -- please check!")

