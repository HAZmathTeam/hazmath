"""
\file examples/haznics/hodge_Hcurl_3d.py
Created by Ana Budisa on 2021-09-03.
Example copyright by bitbucket.org/fenics-apps/cbc.block.

This demo shows the use of a non-trivial block preconditioner for the Hodge
equations namely the curl-curl auxiliary space preconditioner is employed
for the Schur complement . It is adapted from the code described in the block
preconditioning chapter of the FENiCS book, by Kent-Andre Mardal <kent-and@simula.no>.

The block structure is as follows,

       | A   B |
  AA = |       |,
       | C  -D |

where C=B' and D is positive definite; hence, the system as a whole is
symmetric indefinite.

The block preconditioner is based on an approximation of the Schur complement
of the (0,0) block, L=D+B*A^*C:

        | L  0 |
  BB^ = |      |,
        | 0  S |

For L we use HXCurl auxiliary space preconditioner, i.e.

  L^ = HXCurl(L),

and for S we use AMG, i.e.

  S^ = AMG(S).
"""

from __future__ import division
from __future__ import print_function

from dolfin import *
from block import *
from block.iterative import MinRes
from block.algebraic.hazmath import AMG, HXCurl
import sys

set_log_level(30)

N = 4
dim = 3
# Parse command-line arguments like "N=6"
for s in sys.argv[1:]:
    exec(s)

mesh = UnitCubeMesh(N, N, N)

V = FunctionSpace(mesh, "N1curl", 1)
Q = FunctionSpace(mesh, "CG", 1)

v, u = TestFunction(V), TrialFunction(V)
q, p = TestFunction(Q), TrialFunction(Q)

A = assemble(dot(u,v)*dx + dot(curl(v), curl(u))*dx)
B = assemble(dot(grad(p),v)*dx)
C = assemble(dot(grad(q),u)*dx)
D = assemble(p*q*dx)
E = assemble(p*q*dx + dot(grad(p),grad(q))*dx)

AA = block_mat([[A,  B],
                [C, -D]])

gdim = mesh.geometry().dim()
b0 = assemble(inner(v, Constant((1, )*gdim))*dx)
              
b1 = assemble(inner(q, Constant(2))*dx)
bb = block_vec([b0, b1])

prec = block_mat([[HXCurl(A, V),  0  ],
                  [0,            AMG(E)]])

AAinv = MinRes(AA, precond=prec, tolerance=1e-9, maxiter=2000, show=2)

[Uh, Ph] = AAinv*bb

"""
if MPI.size(mesh.mpi_comm()) == 1 and gdim == 2:
    import matplotlib.pyplot as plt

    plt.subplot(121)
    plot(Function(V, Uh))

    plt.subplot(122)    
    plot(Function(Q,  Ph))

    plt.show()
"""
