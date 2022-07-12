from dolfin import *
from block.iterative import ConjGrad
from block.algebraic.hazmath import AMG as AMG_haz
from block.algebraic.petsc.precond import AMG as AMG_petsc
import haznics

# Function spaces, elements
# mesh = UnitSquareMesh(1024, 1024)
mesh = UnitCubeMesh(64, 64, 64)

V = FunctionSpace(mesh, "CG", 1)

f = Expression("sin(3.14*x[0])", degree=2)
u, v = TrialFunction(V), TestFunction(V)

a = u*v*dx + dot(grad(u), grad(v))*dx
L = f*v*dx

A = assemble(a)
b = assemble(L)

# choices for hazmath AMG parameters
params = {
            "print_level": 2,              # (0, 1, 2, 4, 8, 10), 0 - print nothing, 10 - print everything
            "AMG_type": haznics.UA_AMG,             # (UA, SA) + _AMG
            "cycle_type": haznics.NL_AMLI_CYCLE,            # (V, W, AMLI, NL_AMLI, ADD) + _CYCLE
            # "max_levels": 10,
            # "maxit": 1,
            # "tol": 1E-6,
            # "amli_degree": 2,             # Degree of amli polynomials
            # "nl_amli_krylov_type": 4,     # krylov method of nonlinear amli cycle (4 - FGMRES)
            "smoother": haznics.SMOOTHER_GS,         # SMOOTHER_ + (JACOBI, GS, SGS, CG, SOR, SSOR, GSOR, SGSOR, POLY, L1DIAG, FJACOBI, FGS, FSGS)
            # "relaxation": 1.0,            # Relaxation in the smoother
            # "presmooth_iter": 1,
            # "postsmooth_iter": 1,
            # "polynomial_degree": 2,
            # "coarse_dof": 10000,
            "coarse_solver": 32,    # (32 = DIRECT, 0 = ITERATIVE)
            # "coarse_scaling": haznics.ON,      # (OFF, ON)
            # "fpwr": 1.0,                  # fractional power (only in fractional AMG)
            "aggregation_type": haznics.VMB,    # (VMB, MIS, MWM, HEC)
            "strong_coupled": 0.0,       # threshold
            "max_aggregation": 100,
            # "tentative_smooth": 0.67,     # Set the relaxation weight for a particular level
            # "smooth_filter": 1,        # (OFF = 0, ON = 1)
            # "Schwarz_levels": 0,          # number for levels for Schwarz smoother
            # "Schwarz_mmsize": 200,
            # "Schwarz_maxlvl": 3,
            # "Schwarz_type": haznics.SCHWARZ_FORWARD,  # (SCHWARZ_FORWARD, SCHWARZ_BACKWARD, SCHWARZ_SYMMETRIC)
            # "Schwarz_blksolver": 0,       # type of Schwarz block solver, 0 - iterative
            # "damping_param": 0.0,
            # "BSR_alpha": -1000.0,         # weights for the BSR precond
            # "BSR_omega": -1000.0,         # weights for the BSR precond
            }

# choices for hypre parameters
# i don't know how to optimize this more
params_petsc = {
            "pc_hypre_boomeramg_cycle_type": "W", # (V,W)
            #"pc_hypre_boomeramg_max_levels": 25,
            #"pc_hypre_boomeramg_max_iter": 1,
            #"pc_hypre_boomeramg_tol": 0.0,
            "pc_hypre_boomeramg_truncfactor": 0.3,      # Truncation factor for interpolation (between 0 and 1)
            #"pc_hypre_boomeramg_P_max": 0,             # Max elements per row for interpolation
            "pc_hypre_boomeramg_agg_nl": 2,            # Number of levels of aggressive coarsening (0-1 for 2d, 2-5 for 3d)
            #"pc_hypre_boomeramg_agg_num_paths": 1,     # Number of paths for aggressive coarsening
            "pc_hypre_boomeramg_strong_threshold": 0.5,  #.25 for 2d and .5 for 3d,# Threshold for being strongly connected
            #"pc_hypre_boomeramg_max_row_sum": 0.9,
            #"pc_hypre_boomeramg_grid_sweeps_all": 1,   # Number of sweeps for the up and down grid levels
            #"pc_hypre_boomeramg_grid_sweeps_down": 1,
            #"pc_hypre_boomeramg_grid_sweeps_up":1,
            #"pc_hypre_boomeramg_grid_sweeps_coarse": 1,# Number of sweeps for the coarse level (None)
            #"pc_hypre_boomeramg_relax_type_all":  "symmetric-SOR/Jacobi", # (Jacobi, sequential-Gauss-Seidel, seqboundary-Gauss-Seidel,
                                                                          #  SOR/Jacobi, backward-SOR/Jacobi,  symmetric-SOR/Jacobi,
                                                                          #  l1scaled-SOR/Jacobi Gaussian-elimination, CG, Chebyshev,
                                                                          #  FCF-Jacobi, l1scaled-Jacobi)
            #"pc_hypre_boomeramg_relax_type_down": "symmetric-SOR/Jacobi",
            #"pc_hypre_boomeramg_relax_type_up": "symmetric-SOR/Jacobi",
            #"pc_hypre_boomeramg_relax_type_coarse": "Gaussian-elimination",
            #"pc_hypre_boomeramg_relax_weight_all": 1,   # Relaxation weight for all levels (0 = hypre estimates, -k = determined with k CG steps)
            #"pc_hypre_boomeramg_relax_weight_level": (1,1), # Set the relaxation weight for a particular level
            #"pc_hypre_boomeramg_outer_relax_weight_all": 1,
            #"pc_hypre_boomeramg_outer_relax_weight_level": (1,1),
            #"pc_hypre_boomeramg_no_CF": "",               # Do not use CF-relaxation
            #"pc_hypre_boomeramg_measure_type": "local",   # (local global)
            "pc_hypre_boomeramg_coarsen_type": "HMIS", # (Ruge-Stueben, modifiedRuge-Stueben, Falgout, PMIS, HMIS)
            "pc_hypre_boomeramg_interp_type": "ext+i-cc",# (classical, direct, multipass, multipass-wts, ext+i, ext+i-cc, standard, standard-wts, FF, FF1)
            "pc_hypre_boomeramg_print_statistics": 0,
            "pc_hypre_boomeramg_print_debug": 0,
            }

# from hazmath preconds
B = AMG_haz(A, parameters=params)
print("*" * 70)

# from block/petsc
B1 = AMG_petsc(A, parameters=params_petsc)
print("*" * 70)

Ainv = ConjGrad(A, precond=B, tolerance=1e-6, show=10)
Ainv_petsc = ConjGrad(A, precond=B1, tolerance=1e-6, show=10)

print("hazmath" + "*" * 70)
x1 = Ainv*b
print("hypre" + "*" * 70)
x2 = Ainv_petsc*b

# u1 = Function(V)
# u1.vector()[:] = x1[:]
# plot(u1, title="u, computed by hazmath [x=Ainv*b]")

# u2 = Function(V)
# u2.vector()[:] = x2[:]
# plot(u2, title="u, computed by petsc [x=Ainv*b]")