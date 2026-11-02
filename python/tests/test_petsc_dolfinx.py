# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Tests of the PCs through dolfinx.fem.petsc."""

import numpy as np
import pytest

pytest.importorskip("dolfinx")
PETSc = pytest.importorskip("petsc4py.PETSc")

from mpi4py import MPI

import dolfinx.fem.petsc
import ufl
from dolfinx import fem, mesh

pytestmark = pytest.mark.skipif(
    np.dtype(PETSc.ScalarType) != np.float64, reason="HAZmath needs real double precision"
)


def poisson(msh: mesh.Mesh, degree: int, petsc_options: dict, prefix: str):
    """Lagrange Poisson problem -div grad u = 1, u = 0 on the boundary."""
    V = fem.functionspace(msh, ("Lagrange", degree))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    facets = mesh.exterior_facet_indices(msh.topology)
    bc = fem.dirichletbc(0.0, fem.locate_dofs_topological(V, tdim - 1, facets), V)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = ufl.inner(1.0, v) * ufl.dx
    return dolfinx.fem.petsc.LinearProblem(
        a, L, bcs=[bc], petsc_options=petsc_options, petsc_options_prefix=prefix
    )


def hazmath_options(comm: MPI.Comm) -> dict:
    """CG with HAZmath AMG, on each process's block in parallel."""
    options = {"ksp_type": "cg", "ksp_rtol": 1e-12, "ksp_max_it": 200}
    if comm.size == 1:
        options |= {"pc_type": "python", "pc_python_type": "haznicsx.HazmathAMG"}
    else:
        options |= {
            "pc_type": "bjacobi",
            "sub_pc_type": "python",
            "sub_pc_python_type": "haznicsx.HazmathAMG",
        }
    return options


def unit_mesh(comm: MPI.Comm, dim: int, n: int) -> mesh.Mesh:
    if dim == 2:
        return mesh.create_unit_square(comm, n, n)
    return mesh.create_unit_cube(comm, n, n, n)


@pytest.mark.parametrize("degree", [1, 2, 3])
@pytest.mark.parametrize(
    "dim, n, max_its",
    # serial iterations: 2D about 25, 27, 33 and 3D about 13, 19, 20 for
    # degrees 1, 2, 3; HEC aggregation is randomized
    [(2, 64, 50), (3, 12, 30)],
    ids=["2D", "3D"],
)
def test_poisson(dim, n, max_its, degree):
    comm = MPI.COMM_WORLD
    # about the same number of dofs for each degree
    msh = unit_mesh(comm, dim, n // degree)
    problem = poisson(msh, degree, hazmath_options(comm), "haznicsx_poisson_")
    uh = problem.solve()
    solver = problem.solver
    assert solver.getConvergedReason() > 0
    # block Jacobi weakens the preconditioner as the number of processes grows
    assert solver.getIterationNumber() < max_its * (1 if comm.size == 1 else comm.size + 1)

    lu = poisson(
        msh,
        degree,
        {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"},
        "haznicsx_poisson_lu_",
    )
    u_lu = lu.solve()
    error = np.abs(uh.x.array - u_lu.x.array).max()
    assert comm.allreduce(error, op=MPI.MAX) < 1e-8
