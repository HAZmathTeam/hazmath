"""
\file examples/haznics/demo_poisson.py
Created by Miroslav Kuchta, Ana Budisa on 2020-07-15.

We solve

    div(sigma) = f
         sigma = - K*grad(u)

on a unit square with u bcs enforced on left and top edge by
Lagrange multiplier and the sigma.n is fixed set on the rest.
The preconditioner for this system is based on Riesz map wrt inner product
K*H^1 x (1/K)*H^{-0.5}. Using exact inverse for the first block and
rational approximation (RA) for the fractional block. Outer solver is MinRes.
"""

from block.algebraic.petsc import LU
from block.algebraic.hazmath import RA
from dolfin import *
from xii import *
import haznics


def get_system(n, K, f, sigma0, u0):
    '''A, b, W, bcs'''
    K = Constant(K)

    mesh = UnitSquareMesh(n, n)
    #   4
    # 1   2
    #   3
    bdries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1, 0)
    CompiledSubDomain('near(x[0], 0)').mark(bdries, 1)
    CompiledSubDomain('near(x[0], 1)').mark(bdries, 2)
    CompiledSubDomain('near(x[1], 0)').mark(bdries, 3)
    CompiledSubDomain('near(x[1], 1)').mark(bdries, 4)

    # The LM is on 1 and 4
    bmesh = EmbeddedMesh(bdries, (1, 4))

    V = FunctionSpace(mesh, 'CG', 1)
    Q = FunctionSpace(bmesh, 'CG', 1)
    W = [V, Q]

    u, p = map(TrialFunction, W)
    v, q = map(TestFunction, W)

    n_ = OuterNormal(bmesh, [0.5, 0.5])
    Tu, Tv = Trace(u, bmesh), Trace(v, bmesh)
    dx_ = Measure('dx', domain=bmesh)

    # We're building a 4x4 problem
    a = block_form(W, 2)
    a[0][0] = inner(K * grad(u), grad(v)) * dx
    a[0][1] = inner(p, Tv) * dx_
    a[1][0] = inner(q, Tu) * dx_

    n = FacetNormal(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=bdries)

    L = block_form(W, 1)
    L[0] = inner(f, v) * dx - inner(dot(sigma0, n), v) * ds(2) - \
           inner(dot(sigma0, n), v) * ds(3)
    L[1] = inner(u0, q) * dx_

    A, b = map(ii_assemble, (a, L))

    return A, b, W


def get_rational_preconditioner(W, KK):
    '''Realize inv(H^{-0.5}) by AAA'''
    V, Q = W

    K = Constant(KK)
    u, v = TrialFunction(V), TestFunction(V)
    # K*H1
    a = inner(K * grad(u), grad(v)) * dx + inner(K * u, v) * dx
    B0 = LU(assemble(a))  # Exact inverse

    # K*H^{-0.5}
    p, q = TrialFunction(Q), TestFunction(Q)
    if Q.ufl_element().family() == 'Discontinuous Lagrange':
        assert Q.ufl_element().degree() == 0
        h = CellDiameter(Q.mesh())
        h_avg = avg(h)

        a = h_avg ** (-1) * dot(jump(p), jump(q)) * dS + \
            h ** (-1) * dot(p, q) * ds + \
            inner(p, q) * dx
    else:
        a = inner(grad(p), grad(q)) * dx + inner(p, q) * dx

    m = inner(p, q) * dx
    params = {'coefs': [1. / KK, 0.0], 'pwrs': [-0.5, 0.0],
              'print_level': 0,
              'AMG_type': haznics.SA_AMG,
              'cycle_type': haznics.V_CYCLE,
              "max_levels": 20,
              "tol": 1E-10,
              "smoother": haznics.SMOOTHER_GS,
              "relaxation": 1.2,  # Relaxation in the smoother
              "coarse_dof": 10,
              "aggregation_type": haznics.VMB,  # (VMB, MIS, MWM, HEC)
              "strong_coupled": 0.0,  # threshold
              "max_aggregation": 100
              }

    B1 = RA(assemble(a), assemble(m), parameters=params)

    return block_diag_mat([B0, B1])


# --------------------------------------------------------------------


if __name__ == '__main__':
    import sympy as sp
    from block.iterative import MinRes

    # Setup MMS
    def as_expression(thing):
        return Expression(sp.printing.ccode(thing), degree=4, K=1)

    x, y, K = sp.symbols('x[0] x[1] K')

    u = sp.sin(x ** 2 + y ** 2)
    sigma_x = -K * u.diff(x, 1)
    sigma_y = -K * u.diff(y, 1)
    f = sigma_x.diff(x, 1) + sigma_y.diff(y, 1)
    # On the left edge and top

    # Now wrap as expressions
    u, sigma_x, sigma_y, f = map(as_expression, (u, sigma_x, sigma_y, f))
    sigma = Expression(('sigma_x', 'sigma_y'), degree=5,
                       sigma_x=sigma_x, sigma_y=sigma_y)
    # Just for quick updates (as loop) of K
    expressions = (u, sigma_x, sigma_y, sigma, f)

    # Parameters
    n = 128
    K = 1E0
    [setattr(e, 'K', K) for e in expressions]

    AA, bb, W = get_system(n, K=K, f=f, sigma0=sigma, u0=u)

    # Direct solve
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
    print('Errors', errornorm(u, wh[0], 'H1', degree_rise=2))
    print('*' * 32)

    # Plot
    # File('solution_poisson/solution_u.pvd') << wh[0]
