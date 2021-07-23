"""
\file examples/haznics/demo_darcy_stokes.py
Created by Miroslav Kuchta, Ana Budisa on 2021-07-23.

We solve on [0, 1]^2 split down middle Stokes|Darcy

(Stokes)
-div(T(u1, p1)) = f1
        div(u1) = 0

with T(u1, p1) = -p1*I + mu*sym(grad(u1))

(Darcy)
K^-1 u2 + grad(p2) = f2
           div(u2) = 0

(coupling)
u1*n1 + u2*n2 = gD
   -(T.n1).n1 = p2 + gN
  -(T.n1).tau = gT

The coupling is enforced by Lagrange multiplier. The preconditioner for this
system is based on Riesz map wrt inner product
mu*H^1 x (1/mu)*L^2 x (1/K)*(I-grad div) x K*L^2 x ((1/mu)*H^{-0.5} + K*H^{0.5}).
Bcs are chosen so that we need powers of the SAME operator in the multiplier
preconditioner. We use exact inverse for the Stokes and Darcy blocks and
rational approximation (RA) for the fractional block. Outer solver is MinRes.
"""
from xii import *
from block.algebraic.petsc import LU
from block.algebraic.hazmath import RA
from dolfin import *
import sympy as sp
import ulfy
import haznics


def setup_mms(mu_value=1, K_value=1):
    """MMS problem"""
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
    stokes_traction = lambda n, u=u1, p=p1: dot(n, stokes_stress(u, p))  # Vector
    # Forcing for Stokes
    f1 = -div(stokes_stress(u1, p1))

    # Darcy pressure, velocity and force
    p2 = Function(V)
    u2 = -K * grad(p2)
    f2 = div(u2)

    R = Constant(((0, 1), (-1, 0)))
    # Normal and tangent of the interface. NOTE: Stokes is master
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


def get_system(n, K, mu, data):
    """A, b, W, bcs"""
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


def get_rational_preconditioner(AA, W, K, mu):
    """ Realize inv((1/mu)*H^{-0.5} + K*H^{-0.5}) by RA """
    mu_, K_ = Constant(mu), Constant(K)

    V1, Q1, V2, Q2, Q = W

    # Velocity matrix extracted from the system
    B0 = AA[0][0]

    # Stokes pressure is
    p1, q1 = TrialFunction(Q1), TestFunction(Q1)
    B1 = assemble((1 / mu_) * inner(p1, q1) * dx)

    # Darcy flux
    u2, v2 = TrialFunction(V2), TestFunction(V2)
    B2 = assemble((1 / K_) * (inner(u2, v2) * dx + inner(div(u2), div(v2)) * dx))

    # Darcy pressure
    p2, q2 = TrialFunction(Q2), TestFunction(Q2)
    B3 = assemble(K_ * inner(p2, q2) * dx)

    # Multiplier
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

    # parameters for RA and AMG
    params = {'coefs': [1. / mu, K], 'pwrs': [-0.5, 0.5],  # for RA
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

    blocks = [LU(mat) for mat in (B0, B1, B2, B3)]
    blocks.append(B4)

    return block_diag_mat(blocks)

# --------------------------------------------------------------------


if __name__ == '__main__':
    from block.iterative import MinRes

    # Parameters
    n = 64
    mu, K = 1E-6, 1E-4
    print("Parameters: mu = ", mu, "\t K = ", K)

    # Setup a MMS
    data = setup_mms(mu_value=mu, K_value=K)
    u1_true, p1_true, u2_true, p2_true = data['solution']

    # Setup system
    AA, bb, W = get_system(n, K=K, mu=mu, data=data)

    # Hazmath rational approximation preconditioner
    BB = get_rational_preconditioner(AA, W, K=K, mu=mu)

    # Solve
    AAinv = MinRes(AA, precond=BB, tolerance=1E-10, show=4)
    xx = AAinv * bb

    cvrg = AAinv.residuals

    # Solution
    wh = ii_Function(W)
    ndofs = 0
    for i, xxi in enumerate(xx):
        wh[i].vector()[:] = xxi
        ndofs += xxi.size()

    # Results
    print('*' * 32)
    print('MinRes iterations ', len(cvrg), 'final norm', cvrg[-1],
          '#dofs', ndofs)
    print('Error u_stokes: ', errornorm(u1_true, wh[0], 'H1', degree_rise=2),
          '\nError p_stokes: ', errornorm(p1_true, wh[1], 'L2', degree_rise=2))
    print('Error u_darcy: ', errornorm(u2_true, wh[2], 'Hdiv', degree_rise=2),
          '\nError p_darcy: ', errornorm(p2_true, wh[3], 'L2', degree_rise=2))
    print('*' * 32)

    # Plot
    # File('solution_ds/solution_us.pvd') << wh[0]
    # File('solution_ds/solution_ps.pvd') << wh[1]
    # File('solution_ds/solution_ud.pvd') << wh[2]
    # File('solution_ds/solution_pd.pvd') << wh[3]
