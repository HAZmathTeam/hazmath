"""
\file examples/haznics/poisson_dolfinx.py
dolfinx version of examples/haznics/poisson.py (no cbc.block needed).

We solve

    u - div(grad(u)) = f

on the unit cube with natural boundary conditions. The preconditioner
is hazmath AMG, applied through a PETSc python PC.
Outer solver is PETSc Conjugate Gradients.

Serial only: hazmath gets the whole matrix in CSR format.
Run with the python that matches the system petsc4py/dolfinx, e.g.
    PYTHONPATH=../../swig_files /usr/bin/python3 poisson_dolfinx.py
"""
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import ufl
from dolfinx import mesh, fem
from dolfinx.fem.petsc import assemble_matrix, assemble_vector
import haznics

# Function spaces, elements
msh = mesh.create_unit_cube(MPI.COMM_WORLD, 16, 16, 16)
V = fem.functionspace(msh, ("Lagrange", 1))

x = ufl.SpatialCoordinate(msh)
f = ufl.sin(ufl.pi*x[0])
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

a = fem.form(u*v*ufl.dx + ufl.dot(ufl.grad(u), ufl.grad(v))*ufl.dx)
L = fem.form(f*v*ufl.dx)

A = assemble_matrix(a)
A.assemble()
b = assemble_vector(L)
b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

# hazmath matrix from the PETSc CSR arrays (hazmath INT is 32 bit)
ia, ja, aa = A.getValuesCSR()
Ahaz = haznics.create_matrix(np.ascontiguousarray(aa, dtype=np.float64),
                             np.ascontiguousarray(ja, dtype=np.int32),
                             np.ascontiguousarray(ia, dtype=np.int32),
                             A.getSize()[1])

# here we use hazmath AMG:
amgparam = haznics.AMG_param()
haznics.param_amg_init(amgparam)
haznics.param_amg_set_dict({"max_levels": 10, "print_level": 10,
                            "coarse_solver": 32}, amgparam)
prec = haznics.create_precond_amg(Ahaz, amgparam)


class HazmathAMG:
    """PETSc python PC: y = B x with B the hazmath AMG V-cycle."""
    def __init__(self, prec):
        self.prec = prec

    def apply(self, pc, xin, yout):
        r = xin.getArray(readonly=True).copy()
        z = np.zeros_like(r)
        haznics.apply_precond(r, z, self.prec)
        yout.setArray(z)


ksp = PETSc.KSP().create(msh.comm)
ksp.setOperators(A)
ksp.setType("cg")
ksp.setTolerances(rtol=1e-10)
pc = ksp.getPC()
pc.setType("python")
pc.setPythonContext(HazmathAMG(prec))
ksp.setMonitor(lambda ksp, it, rnorm: print(f"  CG iter {it:3d}  |r| = {rnorm:.3e}"))

# solve
uh = fem.Function(V)
ksp.solve(b, uh.x.petsc_vec)
print("CG + hazmath AMG: iterations =", ksp.getIterationNumber(),
      " reason =", ksp.getConvergedReason())

# default sparse direct solver in PETSc, for comparison
lu = PETSc.KSP().create(msh.comm)
lu.setOperators(A)
lu.setType("preonly")
lu.getPC().setType("lu")
u2 = fem.Function(V)
lu.solve(b, u2.x.petsc_vec)

print("Max differences between the two solutions: ",
      np.max(np.abs(uh.x.array - u2.x.array)))
