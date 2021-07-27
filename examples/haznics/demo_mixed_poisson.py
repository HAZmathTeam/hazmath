"""
\file examples/haznics/demo_mixed_poisson.py
Created by Miroslav Kuchta, Ana Budisa on 2021-07-22.

We solve

    div(sigma) = f
         sigma = - K*grad(u)

on a unit square with sigma.n bcs enforced on left and top edge by Lagrange
multiplier and the pressure u is fixed set on the rest. The preconditioner for
this system is based on Riesz map wrt inner product
(1/K)*(I-grad div) x K*L^2 x K*H^{0.5}.
We use exact inverse for the first and second blocks and rational
approximation (RA) for the fractional block. Outer solver is MinRes.
"""

from block.algebraic.petsc import LU
from block.algebraic.hazmath import RA
from petsc4py import PETSc
from dolfin import *
import numpy as np
from xii import *
import haznics


def get_system(n, K, f, sigma0, u0):
    '''A, b, W, bcs'''
    K = Constant(K)
    
    mesh = UnitSquareMesh(n, n)
    #   4
    # 1   2
    #   3
    bdries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    CompiledSubDomain('near(x[0], 0)').mark(bdries, 1)
    CompiledSubDomain('near(x[0], 1)').mark(bdries, 2)
    CompiledSubDomain('near(x[1], 0)').mark(bdries, 3)
    CompiledSubDomain('near(x[1], 1)').mark(bdries, 4)

    # The LM is on 1 and 4
    bmesh = EmbeddedMesh(bdries, (1, 4))

    S = FunctionSpace(mesh, 'RT', 1)
    V = FunctionSpace(mesh, 'DG', 0)
    Q = FunctionSpace(bmesh, 'DG', 0)
    W = [S, V, Q]

    sigma, u, p = map(TrialFunction, W)
    tau, v, q = map(TestFunction, W)
    
    n_ = OuterNormal(bmesh, [0.5, 0.5])
    Ttau, Tsigma = Trace(tau, bmesh), Trace(sigma, bmesh)
    dx_ = Measure('dx', domain=bmesh)
    
    # We're building a 3x3 problem
    a = block_form(W, 2)
    a[0][0] = inner((1./K)*sigma, tau)*dx
    a[0][1] = -inner(u, div(tau))*dx
    a[0][2] = inner(p, dot(Ttau, n_))*dx_
    a[1][0] = -inner(v, div(sigma))*dx
    a[2][0] = inner(q, dot(Tsigma, n_))*dx_

    n = FacetNormal(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=bdries)
    
    L = block_form(W, 1)
    L[0] = -inner(u0, dot(tau, n))*ds(2) - inner(u0, dot(tau, n))*ds(3)
    L[1] = -inner(f, v)*dx
    L[2] = inner(dot(sigma0, n_), q)*dx_

    A, b = map(ii_assemble, (a, L))

    A = ii_convert(A, algorithm='')

    return A, b, W, mesh


def get_rational_preconditioner(W, K):
    '''Realize inv(H^{0.5}) by AAA'''
    S, V, Q = W
    
    KK = Constant(K)
    # (1/K)*(I-grad(div))
    sigma, tau = TrialFunction(S), TestFunction(S)
    a = (1/KK)*(inner(sigma, tau)*dx + inner(div(sigma), div(tau))*dx)
    B0 = LU(assemble(a))  # Exact inverse
    
    # K*I
    u, v = TrialFunction(V), TestFunction(V)
    a = KK*inner(u, v)*dx
    B1 = LU(assemble(a))

    # K*H^{0.5}
    p, q = TrialFunction(Q), TestFunction(Q)
    if Q.ufl_element().family() == 'Discontinuous Lagrange':
        assert Q.ufl_element().degree() == 0
        h = CellDiameter(Q.mesh())
        h_avg = avg(h)

        a = h_avg ** (-1) * dot(jump(p), jump(q)) * dS + \
            h ** (-1) * dot(p, q) * ds + \
            inner(p, q) * dx
    else:
        a = inner(grad(p), grad(q)) * dx + KK * inner(p, q) * dx

    m = inner(p, q) * dx
    params = {'coefs': [K, 0.0], 'pwrs': [0.5, 0.0],
              'print_level': 0,
              'AMG_type': haznics.SA_AMG,
              'cycle_type': haznics.V_CYCLE,
              "max_levels": 20,
              "tol": 1E-10,
              "smoother": haznics.SMOOTHER_GS,
              "relaxation": 1.2,            # Relaxation in the smoother
              "coarse_dof": 10,
              "aggregation_type": haznics.VMB,  # (VMB, MIS, MWM, HEC)
              "strong_coupled": 0.0,  # threshold
              "max_aggregation": 100
              }

    B2 = RA(assemble(a), assemble(m), parameters=params)

    return block_diag_mat([B0, B1, B2])


def solve_minres(AA, bb, BB, W, tol=1E-10):
    '''AA*x = bb with BB-preconditioned MinRes'''
    # opts = PETSc.Options()
    ksp = PETSc.KSP().create()

    reshist = []

    def monitor(ksp, its, rnorm):
        reshist.append(rnorm)
        print(its, rnorm, rnorm / reshist[0])

    ksp.setMonitor(monitor)

    # opts.setValue('-ksp_monitor_true_residual', True)
    
    ksp.setConvergenceHistory()
    ksp.setNormType(PETSc.KSP.NormType.NORM_PRECONDITIONED)
    ksp.setTolerances(rtol=tol, max_it=200)
    ksp.setType('minres')

    ksp.incrementTabLevel(1)

    ksp.setOperators(ii_PETScOperator(AA, None))  # Wrap block_mat
    ksp.setPC(ii_PETScPreconditioner(BB, ksp))  # Wrapped block_op
    ksp.setFromOptions()
    
    wh = ii_Function(W) 
    # Want the iterations to start from random
    wh.block_vec().randomize()
    ksp.solve(as_petsc_nest(bb), wh.petsc_vec())

    return wh, ksp.getConvergenceHistory()

# --------------------------------------------------------------------


if __name__ == '__main__':
    import sympy as sp
    from block.iterative import MinRes

    # Setup MMS
    def as_expression(thing):
        return Expression(sp.printing.ccode(thing), degree=4, K=1)
    
    x, y, K = sp.symbols('x[0] x[1] K')

    u = sp.sin(x**2 + y**2)
    sigma_x = -K*u.diff(x, 1)
    sigma_y = -K*u.diff(y, 1)
    f = sigma_x.diff(x, 1) + sigma_y.diff(y, 1)
    # On the left edge and top

    # Now wrap as expressions
    u, sigma_x, sigma_y, f = map(as_expression, (u, sigma_x, sigma_y, f))
    sigma = Expression(('sigma_x', 'sigma_y'), degree=5, sigma_x=sigma_x, sigma_y=sigma_y)
    # Just for quick updates (as loop) of K
    expressions = (u, sigma_x, sigma_y, sigma, f)

    # Parameters
    n = 64
    K = 1E0
    [setattr(e, 'K', K) for e in expressions]

    AA, bb, W, mesh = get_system(n, K=K, f=f, sigma0=sigma, u0=u)

    # DIRECT SOLVE
    # wh = ii_Function(W)
    # solve(ii_convert(AA), wh.vector(), ii_convert(bb))

    # Preconditioner from hazmath
    BB = get_rational_preconditioner(W, K)

    # Solve with MinRes
    AAinv = MinRes(AA, precond=BB, tolerance=1E-10, show=2)
    xx = AAinv * bb

    # Solution
    wh = ii_Function(W)
    ndofs = 0
    for i, xxi in enumerate(xx):
        wh[i].vector()[:] = xxi
        ndofs += xxi.size()

    # Results
    cvrg = AAinv.residuals

    print('*' * 32)
    print('MinRes iterations', len(cvrg), 'final norm', cvrg[-1],
          '#dofs', ndofs)
    print('Error sigma: ', errornorm(sigma, wh[0], 'Hdiv', degree_rise=2),
          '\nError u: ', errornorm(u, wh[1], 'L2', degree_rise=2))
    print('*' * 32)

    # Plot
    # File('solution_mixedp/solution_sigma.pvd') << wh[0]
    # File('solution_mixedp/solution_u.pvd') << wh[1]
