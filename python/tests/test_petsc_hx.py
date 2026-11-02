# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Tests of the HX PCs (HazmathHXCurl, HazmathHXDiv) with petsc4py and DOLFINx."""

import numpy as np
import pytest

import haznicsx

pytest.importorskip("dolfinx")
PETSc = pytest.importorskip("petsc4py.PETSc")

from mpi4py import MPI

import dolfinx.fem.petsc
import ufl
from dolfinx import fem, mesh

pytestmark = pytest.mark.skipif(
    np.dtype(PETSc.ScalarType) != np.float64, reason="HAZmath needs real double precision"
)


def unit_mesh(comm: MPI.Comm, dim: int, n: int) -> mesh.Mesh:
    if dim == 2:
        return mesh.create_unit_square(comm, n, n)
    return mesh.create_unit_cube(comm, n, n, n)


def operator(kind: str, msh: mesh.Mesh):
    """The H(curl) or H(div) space and PETSc matrix of the unit square or cube."""
    if kind == "curl":
        V = fem.functionspace(msh, ("N1curl", 1))
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        a = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx + ufl.inner(u, v) * ufl.dx
    else:
        V = fem.functionspace(msh, ("RT", 1))
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        a = ufl.inner(ufl.div(u), ufl.div(v)) * ufl.dx + ufl.inner(u, v) * ufl.dx
    A = dolfinx.fem.petsc.assemble_matrix(fem.form(a))
    A.assemble()
    return V, A


def context(kind: str, V: fem.FunctionSpace, **kwargs):
    if kind == "curl":
        return haznicsx.HazmathHXCurl(*haznicsx.hx_curl_operators(V), **kwargs)
    return haznicsx.HazmathHXDiv(*haznicsx.hx_div_operators(V), **kwargs)


def solve(A, ctx, ksp_type: str = "cg", prefix: str = ""):
    """Solve A x = 1 with the context as PC; returns the KSP and the relative residual."""
    ksp = PETSc.KSP().create(A.getComm())
    ksp.setOptionsPrefix(prefix)
    ksp.setOperators(A)
    ksp.setType(ksp_type)
    # with right preconditioning, rtol applies to the true residual
    ksp.setPCSide(PETSc.PC.Side.RIGHT if ksp_type == "gmres" else PETSc.PC.Side.LEFT)
    ksp.setTolerances(rtol=1e-10, max_it=300)
    pc = ksp.getPC()
    pc.setType("python")
    pc.setPythonContext(ctx)
    ksp.setFromOptions()
    b = A.createVecLeft()
    b.set(1.0)
    x = A.createVecRight()
    ksp.solve(b, x)
    r = b.duplicate()
    A.mult(x, r)
    r.axpy(-1.0, b)
    return ksp, r.norm() / b.norm()


@pytest.mark.parametrize("kind", ["curl", "div"])
@pytest.mark.parametrize("dim, n", [(2, 16), (3, 6)], ids=["2D", "3D"])
def test_hx_pc(kind, dim, n):
    V, A = operator(kind, unit_mesh(MPI.COMM_SELF, dim, n))
    ksp, residual = solve(A, context(kind, V))
    assert ksp.getConvergedReason() > 0
    assert ksp.getIterationNumber() < 40
    assert residual < 1e-8
    ksp.destroy()


@pytest.mark.parametrize("kind", ["curl", "div"])
def test_hx_pc_multiplicative(kind):
    V, A = operator(kind, unit_mesh(MPI.COMM_SELF, 3, 6))
    # not symmetric, so GMRES instead of CG
    ksp, residual = solve(A, context(kind, V, multiplicative=True), ksp_type="gmres")
    assert ksp.getConvergedReason() > 0
    assert ksp.getIterationNumber() < 60
    assert residual < 1e-8
    ksp.destroy()


def test_hx_pc_from_options():
    options = PETSc.Options()
    options["hx_pc_hazmath_cycle_type"] = "v"
    options["hx_pc_hazmath_smoother"] = "sgs"
    try:
        V, A = operator("curl", unit_mesh(MPI.COMM_SELF, 3, 4))
        ctx = context("curl", V)
        ksp, residual = solve(A, ctx, prefix="hx_")
        assert ctx.param.cycle_type == haznicsx.V_CYCLE
        assert ctx.param.smoother == haznicsx.SMOOTHER_SGS
        assert residual < 1e-8
        ksp.destroy()
    finally:
        options.delValue("hx_pc_hazmath_cycle_type")
        options.delValue("hx_pc_hazmath_smoother")


def test_hx_pc_petsc_auxiliary_matrices():
    """The auxiliary matrices may also be PETSc matrices."""
    msh = unit_mesh(MPI.COMM_SELF, 3, 4)
    V, A = operator("curl", msh)
    P1 = fem.functionspace(msh, ("Lagrange", 1))
    P1v = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    G = dolfinx.fem.petsc.discrete_gradient(P1, V)
    G.assemble()
    Pcurl = dolfinx.fem.petsc.interpolation_matrix(P1v, V)
    Pcurl.assemble()
    ksp, residual = solve(A, haznicsx.HazmathHXCurl(Pcurl, G))
    assert ksp.getConvergedReason() > 0
    assert residual < 1e-8
    ksp.destroy()
    G.destroy()
    Pcurl.destroy()


@pytest.mark.skipif(MPI.COMM_WORLD.size == 1, reason="needs more than one process")
def test_hx_pc_parallel_rejected():
    _, A = operator("curl", unit_mesh(MPI.COMM_WORLD, 3, 4))
    # no discarded PETSc objects here: garbage-collecting one while PETSc
    # handles the error below crashes
    V = fem.functionspace(unit_mesh(MPI.COMM_SELF, 3, 4), ("N1curl", 1))
    # destroy the distributed objects explicitly, on every process together
    ksp = PETSc.KSP().create(A.getComm())
    ksp.setOperators(A)
    pc = ksp.getPC()
    pc.setType("python")
    pc.setPythonContext(context("curl", V))
    with pytest.raises(PETSc.Error):
        ksp.setUp()
    ksp.destroy()
    A.destroy()
