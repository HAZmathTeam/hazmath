"""
\file examples/haznics/demo_sumfractional.py
Created by Miroslav Kuchta, Ana Budisa on 2021-07-26.

We solve a fractional operator problem

     (alpha * D^s + beta * D^t) u = f

on the boundary of a unit square or unit cube, where D is the Laplacian and
s,t in (-1,1) are the fractional powers. The system operator
(alpha * D^s + beta * D^t) is created solving the generalized eigenvalue problem

    A V = Lambda M V

where M is the mass matrix and A is the matrix representation of D. Then

    alpha * A^s + beta * A^t = M V (alpha * Lambda^s + beta * Lambda^t) (M V)^T

The preconditioner for this system is the rational approximation (RA) of the
function (alpha * x^s + beta * x^t)^(-1). We test the approximation and
efficiency power of RA for sum of fractionalities.
Outer solver is Conjugate Gradients.
"""
from scipy.sparse import csr_matrix
from petsc4py import PETSc
from dolfin import *
import numpy as np


# NOTE: This is only temporary.
# The eigenvalue solver is taken from hseig.py in HsMG.
def my_eigh(A, B):
    '''Au = lmbda Bu transforming to EVP'''
    # Transformation
    timer = Timer('Eigh Power')
    beta, U = np.linalg.eigh(B)
    Bnh = U.dot(np.diag(beta ** -0.5).dot(U.T))
    print('\tDone power in %s' % timer.stop())

    S = Bnh.dot(A.dot(Bnh))
    lmbda, V = np.linalg.eigh(S)
    # With transformed eigenvectors
    return lmbda, Bnh.dot(V)


def Hs_matrix(n, tdim=1, s=[(1.0, 0.5)]):
    """ [(a, s), (b, t)] -> Fractional Laplacian a*D^{s} + b*D^{t} """
    assert n > 0, " Make sure ndof > 0! "

    # Get boundary mesh (trace dim)
    if tdim == 1:
        mesh = UnitSquareMesh(n, n)
    else:
        mesh = UnitCubeMesh(n, n, n)
    mesh = BoundaryMesh(mesh, 'exterior')

    # Assemble inner product and div grad
    V = FunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(V), TestFunction(V)

    a = inner(grad(u), grad(v)) * dx + inner(u, v) * dx
    m = inner(u, v) * dx

    A, M = map(assemble, (a, m))

    # Solve generalized eigenvalue problem
    A_, M_ = A.array(), M.array()
    print('GEVP %d x %d' % (V.dim(), V.dim()))
    timer = Timer('Hs')
    lmbda, U = my_eigh(A_, M_)
    print('Min eig: ', lmbda[0], '\n Max eig: ', lmbda[-1])
    print('Done %s' % timer.stop())

    W = M_.dot(U)

    # Assemble fractional matrix
    if isinstance(s, (int, float)):
        s = [(1, s)]
    assert all(len(si) == 2 and -1 <= si[1] <= 1 for si in s), \
        " Make sure s in (-1, 1) and coefs > 0. "

    return (mesh, V,
            csr_matrix(W.dot(np.diag(sum(ai * lmbda ** si for ai, si in s)).dot(W.T))),
            DomainBoundary(), M, A)


def csr_to_dolfin(A):
    """ PETScMatrix from scipy.sparse.csr_matrix """
    Amat = PETSc.Mat().createAIJ(size=A.shape,
                                 csr=(A.indptr, A.indices, A.data))

    return PETScMatrix(Amat)


# --------------------------------------------------------------------


if __name__ == '__main__':
    import haznics
    from block.algebraic.hazmath import RA
    from block.iterative import MinRes, ConjGrad

    # Parameters
    s = -0.5
    t = 0.5
    alpha = 1E-6
    beta = 1E-2

    tdim = 2
    n = 2 ** 3

    # Get sum fractional matrix
    print("Alpha: ", alpha, "\ts: ", s)
    print("Beta: ", beta, "\tt: ", t)
    _, V, H, _, M, A = Hs_matrix(n, tdim=tdim, s=[(alpha, s), (beta, t)])  # Csr

    # Parameters for the preconditioner
    params = {'coefs': [alpha, beta], 'pwrs': [s, t],
              'print_level': 0,
              'AMG_type': haznics.SA_AMG,
              'cycle_type': haznics.V_CYCLE,
              "max_levels": 20,
              "tol": 1E-10,
              "smoother": haznics.SMOOTHER_GS,
              "relaxation": 1.2,  # Relaxation in the smoother
              "coarse_dof": 10,
              "aggregation_type": haznics.VMB,
              "strong_coupled": 0.0,  # threshold in SA
              "max_aggregation": 100
              }

    # Set up preconditioner
    B = RA(A, M, parameters=params)

    # Set up RHS
    b = Vector(MPI.comm_world, V.dim())
    b.set_local(np.random.rand(V.dim()))

    # Set up solver
    Hsmat = csr_to_dolfin(H)
    Ainv = ConjGrad(Hsmat, precond=B, tolerance=1E-10, show=2)

    # Solve
    x = Ainv * b

    # Results
    cvrg = Ainv.residuals

    print('*' * 32)
    print('CG iterations', len(cvrg), 'final norm', cvrg[-1], '#dofs', V.dim())
    print('*' * 32)
