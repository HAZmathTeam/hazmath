# We solve
#
#     div(sigma) = f
#          sigma = - K*grad(u)
#
# on a unit square with sigma.n bcs enforced on left and top edge by
# Lagrange multiplier and the pressure set on the rest. The preconditioner
# for this system is based on Riesz map wrt inner product
# (1/K)*(I-grad div) x K*H^{0.5}. Using exact inverse for the first block
# we want to see how different approximations of the fractional block
# translate to MinRes iterations
import sys
sys.path.append('../../../')
from block.algebraic.petsc import LU
from block.block_util import flatten
from petsc4py import PETSc
from hseig import HsNorm
from hsmg import HsNormAMG
from precond_haz import RA
import pyamg
from block import block_mat
from dolfin import *
import numpy as np
from xii import *
from dolfin import Matrix as dolfin_Matrix
import scipy.sparse as sps
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
    
    # We're building a 4x4 problem
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


def get_eigenvalue_preconditioner(AA, W, K):
    '''H^{0.5} by eigenvalues and take its exact inverse'''
    S, V, Q = W
    
    K = Constant(K)
    # (1/K)*(I-grad(div))
    sigma, tau = TrialFunction(S), TestFunction(S)
    a = (1/K)*(inner(sigma, tau)*dx + inner(div(sigma), div(tau))*dx)
    B0 = LU(assemble(a))  # Exact inverse
    
    # K*I
    u, v = TrialFunction(V), TestFunction(V)
    a = K*inner(u, v)*dx
    B1 = LU(assemble(a))

    # K*H^{0.5}
    B2 = HsNorm(Q, s=0.5, kappa=K, bcs=True)
    B2*B2.create_vec()

    B2 = LU(B2.matrix)

    return block_diag_mat([B0, B1, B2])


def get_hsmg_preconditioner(AA, W, K):
    '''Realize inv(H^{0.5}) by AMG'''
    S, V, Q = W
    
    K = Constant(K)
    # (1/K)*(I-grad(div))
    sigma, tau = TrialFunction(S), TestFunction(S)
    a = (1/K)*(inner(sigma, tau)*dx + inner(div(sigma), div(tau))*dx)
    B0 = LU(assemble(a))  # Exact inverse
    
    # K*I
    u, v = TrialFunction(V), TestFunction(V)
    a = K*inner(u, v)*dx
    B1 = LU(assemble(a))

    # K*H^{0.5}
    mg_params = {'pyamg_solver': lambda A: pyamg.ruge_stuben_solver(A),
                 'mass_inv': 'amg'}  # In BPL algo

    B2 = HsNormAMG(Q, s=0.5, bdry=DomainBoundary(),
                   mg_params=mg_params,
                   kappa=K,
                   neg_mg='bpl')

    return block_diag_mat([B0, B1, B2])


def get_rational_preconditioner(AA, W, K):
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
    # Setup a MMS
    import sympy as sp
    from block.iterative import MinRes

    as_expression = lambda thing: Expression(sp.printing.ccode(thing), degree=4, K=1)
    
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

    ### CONDUCTIVITY PARAMETER ###
    n = 64
    K = 1E0
    niters, ndofs_, cputime = [], [], []
    Kgrid = [1E-8, 1E-6, 1E-4, 1E-2, 1E0, 1E2, 1E4, 1E6, 1E8]
    ngrid = [4, 8, 16, 32, 64, 128, 256, 512]
    [setattr(e, 'K', K) for e in expressions]
    
    # for n in (4, 8, 16, 32, 64, 128, 256, 512):
    for K in Kgrid:
        AA, bb, W, mesh = get_system(n, K=K, f=f, sigma0=sigma, u0=u)

        # DIRECT SOLVE
        # wh = ii_Function(W)
        # solve(ii_convert(AA), wh.vector(), ii_convert(bb))

        # EIGEN PRECOND
        # BB = get_eigenvalue_preconditioner(AA, W, K)
        """
        # Miro's HSMG PRECOND
        BB_miro = get_hsmg_preconditioner(AA, W, K)
        
        AAinv_miro = MinRes(AA, precond=BB_miro, tolerance=1E-10, show=2)
        xx_miro = AAinv_miro*bb

        res_miro = bb - AA*xx_miro
        res_norm_miro = res_miro.norm()
        print("Residual norm after Minres - Miro's precond: ", res_norm_miro)

        cvrg_miro = AAinv_miro.residuals

        wh = ii_Function(W)
        ndofs = 0
        for i, xxi in enumerate(xx_miro):
            wh[i].vector()[:] = xxi
            ndofs += xxi.size()

        niters.append(len(cvrg_miro))
        ndofs_.append(K)
        cputime.append(AAinv_miro.cputime)

        # Check solution hsmg
        print('*' * 32)
        print('MinRes iterations eigen', len(cvrg_miro), 'final norm',
              cvrg_miro[-1], '#dofs',
              ndofs)
        print('Error sigma: ',
              errornorm(sigma, wh[0], 'Hdiv', degree_rise=2),
              '\nError u: ',
              errornorm(u, wh[1], 'L2', degree_rise=2))
        print('*' * 32)
        """


        # Hazmath RA PRECOND
        BB_haz = get_rational_preconditioner(AA, W, K)

        AAinv_haz = MinRes(AA, precond=BB_haz, tolerance=1E-10, show=2)
        xx_haz = AAinv_haz * bb

        res_haz = bb - AA*xx_haz
        res_norm_haz = res_haz.norm()
        print("Residual norm after Minres - Hazmath's precond: ", res_norm_haz)

        cvrg_haz = AAinv_haz.residuals

        wh = ii_Function(W)
        ndofs = 0
        for i, xxi in enumerate(xx_haz):
            wh[i].vector()[:] = xxi
            ndofs += xxi.size()

        niters.append(len(cvrg_haz))
        ndofs_.append(K)
        cputime.append(AAinv_haz.cputime)
        # Check solution hazmath
        print('*'*32)
        print('MinRes iterations haz', len(cvrg_haz), 'final norm', cvrg_haz[-1], '#dofs', ndofs)
        
        print('Error sigma: ',
              errornorm(sigma, wh[0], 'Hdiv', degree_rise=2),
              '\nError u: ',
              errornorm(u, wh[1], 'L2', degree_rise=2))
        print('*'*32)

        """
        V1, Q1, Lambda = W
        u1h, p1h, lmh = wh

        error_u1 = Function(V1)
        error_u1 = project(sigma - u1h, V1)
        error_u1_total = assemble(inner(error_u1, error_u1) * dx
                                  + inner(div(error_u1), div(error_u1)) * dx)
        sigma_proj = project(sigma, V1)
        norm_u1 = assemble(inner(sigma_proj, sigma_proj) * dx
                           + inner(div(sigma_proj), div(sigma_proj)) * dx)
        print("Rel error u1 in Hdiv: ", sqrt(error_u1_total) / sqrt(norm_u1))

        error_p1 = Function(Q1)
        error_p1 = project(u - p1h, Q1)

        # import pdb; pdb.set_trace()
        File('solution_poisson/error_u1.pvd') << error_u1
        File('solution_poisson/error_p1.pvd') << error_p1
        """

    from tabulate import tabulate

    table = [ndofs_, niters, cputime]
    table = np.array(table).T.tolist()
    # import pdb; pdb.set_trace()
    print('*' * 70)
    table_t = tabulate(table, headers=["K", "niters", "cputime"])
    print(table_t)
    print('*' * 70)
    open('table_hazmath.txt', 'w').write(table_t)
