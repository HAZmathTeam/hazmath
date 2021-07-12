# We solve on [0, 1]^2 split down middle Stokes|Darcy
#
# In Stokes
# -div(T(u1, p1)) = f1
# div(u1) = 0
#
# with T(u1, p1) = -p1*I + mu*sym(grad(u1))
#
#
# In Darcy
# K^-1 u2 + grad(p2) = f2
# div(u2)            = 0
#
# coupling u1*n1 + u2*n2 = gD
#         -(T.n1).n1 = p2 + gN
#          -(T.n1).tau = gT
#
# Bcs on the boundaries are chosen so that we need powers of the SAME
# operator in the multiplier preconditioner
import sys

sys.path.append('../../../')

from xii import *

from block.block_base import block_base
from block.object_pool import vec_pool

from block.algebraic.petsc import LU, AMG
from petsc4py import PETSc
from hseig import HsNorm
from hsmg import HsNormAMG
from precond_haz import RA
import pyamg
from block import block_mat
from dolfin import *
import numpy as np
import sympy as sp
import ulfy
import haznics


def setup_mms(mu_value=1, K_value=1):
    '''MMS problem'''
    mesh = UnitSquareMesh(2, 2)  # Dummy

    V = FunctionSpace(mesh, 'CG', 2)
    S = FunctionSpace(mesh, 'DG', 0)
    # Define constants as function to allow ufly substition
    mu = Function(S)
    K = Function(S)

    # Auxiliary function for defining Stokes velocity
    phi = Function(V)
    # Now the Stokes velocity is
    u1 = as_vector((phi.dx(1), -phi.dx(0)))  # To be divergence free
    p1 = Function(V)

    # Stokes stress
    stokes_stress = lambda u, p: -p * Identity(2) + mu * grad(u)
    stokes_traction = lambda n, u=u1, p=p1: dot(n,
                                                stokes_stress(u, p))  # Vector
    # Forcing for Stokes
    f1 = -div(stokes_stress(u1, p1))

    # Darcy pressure, velocity and force
    p2 = Function(V)
    u2 = -K * grad(p2)
    f2 = div(u2)

    R = Constant(((0, 1), (-1, 0)))
    # Normal and tangegent of the interface. NOTE: Stokes is master
    n = Constant((1, 0))
    tau = dot(R, n)

    # Coupling data, for piece of boundary
    gD = dot(u1, n) - dot(u2, n)
    gN = -p2 - dot(n, stokes_traction(n, u1, p1))
    gT = -dot(tau, stokes_traction(n, u1, p1))

    x, y, mu_, K_ = sp.symbols('x y mu K')

    phi_ = sp.sin(pi * (x + y))  # Aux expr

    p1_ = sp.sin(4 * pi * x) * sp.cos(pi * y)
    p1_ = p1_ + (x + y - 3. / 4)

    p2_ = sp.cos(2 * pi * x) * sp.cos(2 * pi * y)
    p2_ = p2_ + (x + y - 5. / 4)

    subs = {phi: phi_, p1: p1_, p2: p2_, mu: mu_, K: K_}

    as_expression = lambda f: ulfy.Expression(f, subs=subs, degree=4,
                                              mu=mu_value,
                                              K=K_value)

    data = {'solution': [as_expression(f) for f in (u1, p1, u2, p2)],
            'vol_f': [as_expression(f) for f in (f1, f2)],
            'iface_f': [as_expression(f) for f in (gD, gN, gT)],
            'dirichlet': [as_expression(f) for f in (u1, u2)],
            'neumann': [as_expression(f) for f in (stokes_stress(u1, p1), p2)]}

    return data


# ---

def get_system(n, K, mu, data):
    '''A, b, W, bcs'''
    K, mu = Constant(K), Constant(mu)

    mesh = UnitSquareMesh(n, n)
    #   4       4 
    # 1   2   1   2
    #   3       3
    subdomains = MeshFunction('size_t', mesh, 2, 2)
    CompiledSubDomain('x[0] < 0.5+DOLFIN_EPS').mark(subdomains, 1)

    # Stokes
    mesh1 = EmbeddedMesh(subdomains, 1)
    # Tag it
    bdries1 = MeshFunction('size_t', mesh1, mesh1.topology().dim() - 1, 0)
    CompiledSubDomain('near(x[0], 0)').mark(bdries1, 1)
    CompiledSubDomain('near(x[0], 0.5)').mark(bdries1, 2)
    CompiledSubDomain('near(x[1], 0)').mark(bdries1, 3)
    CompiledSubDomain('near(x[1], 1)').mark(bdries1, 4)

    # Darcy
    mesh2 = EmbeddedMesh(subdomains, 2)
    # Tag it
    bdries2 = MeshFunction('size_t', mesh2, mesh2.topology().dim() - 1, 0)
    CompiledSubDomain('near(x[0], 0.5)').mark(bdries2, 1)
    CompiledSubDomain('near(x[0], 1.0)').mark(bdries2, 2)
    CompiledSubDomain('near(x[1], 0)').mark(bdries2, 3)
    CompiledSubDomain('near(x[1], 1)').mark(bdries2, 4)

    # And interface
    bmesh = EmbeddedMesh(bdries1, 2)

    # Onto weak form
    # Stokes
    V1 = VectorFunctionSpace(mesh1, 'CG', 2)
    Q1 = FunctionSpace(mesh1, 'CG', 1)
    # Darcy
    V2 = FunctionSpace(mesh2, 'RT', 1)
    Q2 = FunctionSpace(mesh2, 'DG', 0)
    # The multiplier
    Q = FunctionSpace(bmesh, 'DG', 0)

    W = [V1, Q1, V2, Q2, Q]

    u1, p1, u2, p2, p = map(TrialFunction, W)
    v1, q1, v2, q2, q = map(TestFunction, W)
    # Stokes traces
    Tu1, Tv1 = Trace(u1, bmesh), Trace(v1, bmesh)
    # Darcy traces
    Tu2, Tv2 = Trace(u2, bmesh), Trace(v2, bmesh)

    R = Constant(((0, 1), (-1, 0)))
    # The line integral and orient iface
    dx_ = Measure('dx', domain=bmesh)
    n_ = Constant((1, 0))
    tau_ = dot(R, n_)

    mu, K = Constant(mu), Constant(K)

    a = block_form(W, 2)
    # Stokes
    a[0][0] = inner(mu * grad(u1), grad(v1)) * dx
    a[0][1] = -inner(p1, div(v1)) * dx
    a[0][4] = inner(p, dot(Tv1, n_)) * dx_
    # Darcy
    a[2][2] = K ** -1 * inner(u2, v2) * dx
    a[2][3] = -inner(p2, div(v2)) * dx
    a[2][4] = -inner(p, dot(Tv2, n_)) * dx_
    # Symmetrize
    a[1][0] = -inner(q1, div(u1)) * dx
    a[3][2] = -inner(q2, div(u2)) * dx

    a[4][0] = inner(q, dot(Tu1, n_)) * dx_
    a[4][2] = -inner(q, dot(Tu2, n_)) * dx_

    f1, f2 = data['vol_f']
    gD, gN, gT = data['iface_f']
    # Data for dirichlet bc is only
    u1_true = data['dirichlet'][0]  # Strongly on top and bottom
    # Data for neumann bcs; For stokes we set stress on left and for
    # darcy the pressure is imposed everywhere
    sigma1_true, p2_true = data['neumann']  # 

    # For Neuamnn bcs we'll need
    n1, n2 = FacetNormal(mesh1), FacetNormal(mesh2)

    ds1 = Measure('ds', domain=mesh1, subdomain_data=bdries1)
    ds2 = Measure('ds', domain=mesh2, subdomain_data=bdries2)

    L = block_form(W, 1)

    L[0] = (inner(f1, v1) * dx +
            inner(dot(sigma1_true, n1), v1) * ds1(1) +
            inner(gT, dot(Tv1, tau_)) * dx_)

    # Darcy forcing has the contrib due to coupling for normal tractions
    L[2] = (-inner(gN, dot(Tv2, n_)) * dx_ -
            inner(p2_true, dot(v2, n2)) * ds2(2) -
            inner(p2_true, dot(v2, n2)) * ds2(3) -
            inner(p2_true, dot(v2, n2)) * ds2(4))

    L[3] = -inner(f2, q2) * dx
    # Coupling of velocities
    L[4] = inner(gD, q) * dx_

    A, b = map(ii_assemble, (a, L))

    V1_bcs = [DirichletBC(V1, u1_true, bdries1, 3),
              DirichletBC(V1, u1_true, bdries1, 4)]
    bcs = [V1_bcs, [], [], [], []]

    A, b = apply_bc(A, b, bcs)

    return A, b, W


def get_eigenvalue_preconditioner(AA, W, mu, K, no_frac=False):
    '''H^{0.5} by eigenvalues and take its exact inverse'''
    mu, K = Constant(mu), Constant(K)

    V1, Q1, V2, Q2, Q = W

    # Velocity matrix extracted from the system
    B0 = AA[0][0]

    # Stokes pressure is
    p1, q1 = TrialFunction(Q1), TestFunction(Q1)
    B1 = assemble((1 / mu) * inner(p1, q1) * dx)

    # Darcy flux
    u2, v2 = TrialFunction(V2), TestFunction(V2)
    B2 = assemble((1 / K) * (inner(u2, v2) * dx + inner(div(u2), div(v2)) * dx))

    # Darcy pressure
    p2, q2 = TrialFunction(Q2), TestFunction(Q2)
    B3 = assemble(K * inner(p2, q2) * dx)

    # Multiplier
    # For other preconditioners we want to use this function but not
    # spent time on eigenvalue computations so skip
    if not no_frac:
        X = HsNorm(Q, s=-0.5, kappa=(1 / mu), bcs=True)
        X * X.create_vec()

        Y = HsNorm(Q, s=0.5, kappa=K, bcs=True)
        Y * Y.create_vec()

        B4 = ii_convert(X.matrix + Y.matrix)
    else:
        B4 = 0

    blocks = [LU(mat) for mat in (B0, B1, B2, B3)]
    if B4:
        blocks.append(LU(B4))
    else:
        blocks.append(0)

    return block_diag_mat(blocks)


class Hs0SumNormAMG(block_base):
    '''
    Preconditioner for 

      alpha*D^s + beta*D^(1+s)

    where s in (-1, 0) and D = -Delta + I on V. We will factorize out 
    negative fractional Laplacian
    '''

    def __init__(self, V, s, alpha, beta, mg_params, neg_mg='bpl'):
        assert between(s, (-1, 0))
        assert V.ufl_element().value_shape() == ()
        # FIXME:
        assert V.ufl_element().family() == 'Discontinuous Lagrange'
        assert V.ufl_element().degree() == 0

        self.size = V.dim()

        # We want to apply the decomposiotion
        # D^(s/2) [alpha + beta*D] D^(s/2)
        self.frac = HsNormAMG(V, s=0.5 * s, bdry=None, mg_params=mg_params,
                              neg_mg=neg_mg)
        # FIXME: this always puts on bcs!!!
        self.nlevels = self.frac.nlevels
        # The issue here is that H^s * H^s is not H^2s. What is true though
        # is H^s inv(M) H^s = H^2s we  therefore have
        #
        # H^{s/2}*[alpha*inv(M) + beta*inv(M)*A*inv(M)]*H^{s/2}
        #
        # However this is equal to
        #
        # H^{s/2}*[alpha*inv(M)*M*inv(M) + beta*inv(M)*A*inv(M)]*H^{s/2}
        #
        # so that
        #
        # H^{s/2}*inv(M)*[alpha*M + beta*A]*inv(M)*H^{s/2}
        # Let's build the inner operator.
        # FIXME: here we are again specific to DG0 and Dirichlet bcs!
        u, v = TrialFunction(V), TestFunction(V)
        h = CellDiameter(V.mesh())
        h_avg = avg(h)

        a = beta * h_avg ** (-1) * dot(jump(v), jump(u)) * dS + \
            beta * h ** (-1) * dot(u, v) * ds + \
            beta * inner(u, v) * dx

        m = alpha * inner(u, v) * dx

        N = assemble(a + m)
        # Exact?
        self.nofrac_inv = {'amg': AMG,
                           'lu': LU}[mg_params.get('nofrac_inv', 'amg')](N)

        self.M = assemble(inner(u, v) * dx)

    # Implementation of cbc.block API --------------------------------
    def matvec(self, b):
        x0 = self.frac * b
        x1 = self.M * x0
        x2 = self.nofrac_inv * x1
        x3 = self.M * x2
        x4 = self.frac * x3

        return x4

    @vec_pool
    def create_vec(self, dim=1):
        return Vector(mpi_comm(), self.size)


def get_hsmg_preconditioner(AA, W, K, mu):
    '''Realize inv(H^{0.5}+H^{-0.5}) by AMG'''
    # Get the block from eigenvalue
    BB = get_eigenvalue_preconditioner(AA, W, K=K, mu=mu, no_frac=True)
    # We only want to replace the last block
    mg_params = {'pyamg_solver': lambda A: pyamg.ruge_stuben_solver(A),
                 'mass_inv': 'splu',  # In BPL algo
                 'nofrac_inv': 'lu'}  # Inverting the middle part in the

    mu, K = Constant(mu), Constant(K)
    Q = W[-1]
    B4 = Hs0SumNormAMG(Q, s=-0.5, alpha=(1 / mu), beta=K, mg_params=mg_params)

    BB[4][4] = B4

    # This is a'la Layton, does not liek small mu
    # HsNormAMG(Q, s=0.5, kappa=K, bdry=DomainBoundary(), mg_params=mg_params)
    return BB


def get_rational_preconditioner(AA, W, K, mu):
    '''Realize inv(H^{0.5}+H^{-0.5}) by RA'''
    # Get the block from eigenvalue
    BB = get_eigenvalue_preconditioner(AA, W, K=K, mu=mu, no_frac=True)

    # Get L.m. matrices
    Q = W[-1]
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

    # parameters for RA and AMG for shifted laplacians
    params = {'coefs': [1. / mu, K], 'pwrs': [-0.5, 0.5],
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

    B4 = RA(assemble(a), assemble(m), parameters=params)

    BB[4][4] = B4

    return BB


def get_preconditioner_blocks(AA, W, mu, K):
    '''H^{0.5} by eigenvalues and take its exact inverse'''
    mu, K = Constant(mu), Constant(K)

    V1, Q1, V2, Q2, Q = W

    # Velocity matrix extracted from the system
    B0 = AA[0][0]

    # Stokes pressure is
    p1, q1 = TrialFunction(Q1), TestFunction(Q1)
    B1 = assemble((1 / mu) * inner(p1, q1) * dx)

    # Darcy flux
    u2, v2 = TrialFunction(V2), TestFunction(V2)
    B2 = assemble((1 / K) * (inner(u2, v2) * dx + inner(div(u2), div(v2)) * dx))

    # Darcy pressure
    p2, q2 = TrialFunction(Q2), TestFunction(Q2)
    B3 = assemble(K * inner(p2, q2) * dx)

    # Multiplier
    X = HsNorm(Q, s=-0.5, kappa=(1 / mu), bcs=True)
    X * X.create_vec()

    Y = HsNorm(Q, s=0.5, kappa=K, bcs=True)
    Y * Y.create_vec()

    B4 = ii_convert(X.matrix + Y.matrix)

    blocks = [B0, B1, B2, B3, B4]

    return block_diag_mat(blocks)


def solve_minres(AA, bb, BB, W, tol=1E-10):
    '''AA*x = bb with BB-preconditioned MinRes'''
    
    ksp = PETSc.KSP().create()

    reshist = []
    def monitor(ksp, its, rnorm):
        reshist.append(rnorm)
        print(its, rnorm, rnorm/reshist[0])
    ksp.setMonitor(monitor)
    
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
    from block.iterative import MinRes, LGMRES
    from scipy.linalg import eigvalsh

    ### VISCOSITY AND PERMEABILITY PARAMS ###
    mu, K = 1E0, 1E0
    data = setup_mms(mu_value=mu, K_value=K)

    # u1_true, p1_true, u2_true, p2_true = data['solution']
    # eu1, ep1, eu2, ep2 = [], [], [], []
    niters, ndofs_, cputime = [], [], []
    # ngrid = [16, 32, 64, 128, 256]
    # ngrid = [128]
    n = 64
    Kgrid = [1E0]
    mugrid = [1E-8]
    for mu in mugrid:
        for K in Kgrid:
            AA, bb, W = get_system(n, K=K, mu=mu, data=data)
            # This one is for reference (eigenvalue check)
            # BB = get_eigenvalue_preconditioner(AA, W, K=K, mu=mu)
            # iBB = block_diag_mat([BB[i][i].A for i in range(len(W))])
            # eigw = eigvalsh(ii_convert(AA).array(), ii_convert(iBB).array())
            # lmin, lmax = np.sort(np.abs(eigw))[[0, -1]]
            # print('>>>', lmin, lmax, lmax/lmin, '<<<')

            # hsmg preconditioner
            # BB = get_hsmg_preconditioner(AA, W, K=K, mu=mu)

            # Hazmath rational approximation preconditioner
            BB = get_rational_preconditioner(AA, W, K=K, mu=mu)
            print("mu = ", mu, "\t K = ", K)
            # AAinv = LGMRES(AA, precond=BB, tolerance=1E-10, show=4)
            AAinv = MinRes(AA, precond=BB, tolerance=1E-10, show=4)
            xx = AAinv * bb

            # wh, cvrg = solve_minres(AA, bb, BB, W, tol=1E-10)

            cvrg = AAinv.residuals
            cputime.append(AAinv.cputime)
            # cputime.append(0.)

            wh = ii_Function(W)
            ndofs = 0
            for i, xxi in enumerate(xx):
                wh[i].vector()[:] = xxi
                ndofs += xxi.size()
            # ndofs = 0.

            # Check solution
            # eu1.append(errornorm(u1_true, wh[0], 'H1', degree_rise=2))
            # ep1.append(errornorm(p1_true, wh[1], 'L2', degree_rise=2))
            # eu2.append(errornorm(u2_true, wh[2], 'Hdiv', degree_rise=2))
            # ep2.append(errornorm(p2_true, wh[3], 'L2', degree_rise=2))

            print('*' * 32, 'HAZMATH', '*' * 32)
            print('MinRes iterations haz', len(cvrg), 'final norm', cvrg[-1],
                  '#dofs', ndofs)
            niters.append(len(cvrg))
            # ndofs_.append(ndofs)
            # ndofs_.append(mu)
            ndofs_.append(K)

            print('*' * 70)

            # V1, Q1, V2, Q2, Lambda = W
            # u1h, p1h, u2h, p2h, lmh = wh

            # error_u1 = Function(V1)
            # error_u1 = project(u1_true - u1h, V1)
            # error_u1_total = assemble(inner(error_u1, error_u1) * dx
            #                           + inner(grad(error_u1), grad(error_u1)) * dx)

            # u1_true_proj = project(u1_true, V1)
            # norm_u1 = assemble(inner(u1_true_proj, u1_true_proj) * dx
            #                    + inner(grad(u1_true_proj), grad(u1_true_proj)) * dx)
            # print("Rel error u1 in H1: ", sqrt(error_u1_total) / sqrt(norm_u1))
            # eu1.append(sqrt(error_u1_total) / sqrt(norm_u1))
            """
            if n > 128:
                error_p1 = Function(Q1)
                error_p1 = project(p1_true - p1h, Q1)
    
                error_u2 = Function(V2)
                error_u2 = project(u2_true - u2h, V2)
    
                error_p2 = Function(Q2)
                error_p2 = project(p2_true - p2h, Q2)
    
                # import pdb; pdb.set_trace()
                File('solution_ds/error_u1.pvd') << error_u1
                File('solution_ds/error_p1.pvd') << error_p1
                File('solution_ds/error_u2.pvd') << error_u2
                File('solution_ds/error_p2.pvd') << error_p2
            """

    # print("mu = ", mu, "\t K = ", K)
    from tabulate import tabulate
    table = [ndofs_, niters, cputime]
    table = np.array(table).T.tolist()
    # import pdb; pdb.set_trace()
    print('*' * 70)
    table_t = tabulate(table, headers=["K", "niters", "cputime"])
    print(table_t)
    print('*' * 70)
    open('table_hazmath_%g.txt'%(1/mu), 'w').write(table_t)
