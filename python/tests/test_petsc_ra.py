# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Tests of the RA preconditioner (HazmathRA, estimate_lambda_max) with petsc4py and slepc4py."""

import numpy as np
import pytest
import scipy.linalg
import scipy.sparse as sp

import haznicsx

PETSc = pytest.importorskip("petsc4py.PETSc")
pytest.importorskip("slepc4py.SLEPc")

pytestmark = [
    pytest.mark.skipif(
        np.dtype(PETSc.ScalarType) != np.float64, reason="HAZmath needs real double precision"
    ),
    pytest.mark.skipif(
        not haznicsx.HAS_LAPACK, reason="the RA preconditioner needs libhazmath with LAPACK"
    ),
]


s_values = [0.5, -0.5]


def fractional_operator(A: sp.csr_matrix, M: sp.csr_matrix, s: float) -> np.ndarray:
    """The matrix of A^s in the M inner product, M V diag(lambda^s) V^T M."""
    lam, V = scipy.linalg.eigh(A.toarray(), M.toarray())
    Md = M.toarray()
    return Md @ V @ np.diag(lam**s) @ V.T @ Md


def fgmres(operator, ctx, rtol: float = 1e-8):
    """FGMRES solve of operator x = 1 with ctx as PC; returns the KSP and the relative residual."""
    ksp = PETSc.KSP().create(PETSc.COMM_SELF)
    ksp.setOperators(*operator)
    ksp.setType("fgmres")
    ksp.setTolerances(rtol=rtol, max_it=100)
    pc = ksp.getPC()
    pc.setType("python")
    pc.setPythonContext(ctx)
    b = operator[0].createVecLeft()
    b.set(1.0)
    x = operator[0].createVecRight()
    ksp.solve(b, x)
    r = b.duplicate()
    operator[0].mult(x, r)
    r.axpy(-1.0, b)
    return ksp, r.norm() / b.norm()


def test_estimate_lambda_max(p1_1d, petsc_matrix):
    A, M = p1_1d(63)
    exact = scipy.linalg.eigh(A.toarray(), M.toarray(), eigvals_only=True)[-1]
    A_petsc, M_petsc = petsc_matrix(A), petsc_matrix(M)
    assert haznicsx.estimate_lambda_max(A_petsc, M_petsc) == pytest.approx(exact, rel=1e-4)
    A_petsc.destroy()
    M_petsc.destroy()


@pytest.mark.parametrize("s", s_values)
def test_ra_default_scaling(s, p1_1d, preconditioner_matrix):
    A, M = p1_1d(31)
    lam, V = scipy.linalg.eigh(A.toarray(), M.toarray())
    exact = V @ np.diag(lam**-s) @ V.T
    B = haznicsx.ra(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), s)
    Z = preconditioner_matrix(B)
    assert np.linalg.norm(Z - exact) / np.linalg.norm(exact) < 1e-4


def test_ra_pc_inverse(p1_1d, petsc_matrix):
    """With s = 1 the RA preconditioner approximates A^-1."""
    A, M = p1_1d(63)
    A_petsc, M_petsc = petsc_matrix(A), petsc_matrix(M)
    ksp, residual = fgmres((A_petsc,), haznicsx.HazmathRA(M_petsc, 1.0))
    assert ksp.getConvergedReason() > 0
    assert ksp.getIterationNumber() <= 5
    assert residual < 1e-8
    ksp.destroy()


@pytest.mark.parametrize("s", s_values)
def test_ra_pc_fractional_operator(s, p1_1d, petsc_matrix):
    """FGMRES for the dense fractional operator, preconditioned by RA with the sparse A and M."""
    A, M = p1_1d(63)
    Fs = PETSc.Mat().createDense(A.shape, array=fractional_operator(A, M, s), comm=PETSc.COMM_SELF)
    A_petsc, M_petsc = petsc_matrix(A), petsc_matrix(M)
    ksp, residual = fgmres((Fs, A_petsc), haznicsx.HazmathRA(M_petsc, s))
    assert ksp.getConvergedReason() > 0
    assert ksp.getIterationNumber() <= 5
    assert residual < 1e-8
    ksp.destroy()


def test_ra_pc_options(p1_1d, petsc_matrix):
    """AMGParam options apply to the AMG solves of the shifted problems."""
    options = PETSc.Options()
    options["ra_pc_hazmath_cycle_type"] = "v"
    try:
        A, M = p1_1d(31)
        ctx = haznicsx.HazmathRA(petsc_matrix(M), -0.5)
        ksp = PETSc.KSP().create(PETSc.COMM_SELF)
        ksp.setOptionsPrefix("ra_")
        ksp.setOperators(petsc_matrix(A))
        pc = ksp.getPC()
        pc.setType("python")
        pc.setPythonContext(ctx)
        ksp.setFromOptions()
        ksp.setUp()
        assert ctx.param.cycle_type == haznicsx.V_CYCLE
        ksp.destroy()
    finally:
        options.delValue("ra_pc_hazmath_cycle_type")


@pytest.mark.parametrize("s", s_values)
def test_ra_pc_bounded_under_refinement(s, petsc_matrix):
    """Fractional P1 problems in 2D: iteration counts do not grow under refinement."""
    pytest.importorskip("dolfinx")
    from mpi4py import MPI

    import ufl
    from dolfinx import fem, mesh

    its = []
    for n in (6, 12):
        msh = mesh.create_unit_square(MPI.COMM_SELF, n, n)
        V = fem.functionspace(msh, ("Lagrange", 1))
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        # stiffness plus mass, so that A is positive definite without boundary conditions
        A = fem.assemble_matrix(
            fem.form(ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx + ufl.inner(u, v) * ufl.dx)
        ).to_scipy()
        M = fem.assemble_matrix(fem.form(ufl.inner(u, v) * ufl.dx)).to_scipy()
        Fs = PETSc.Mat().createDense(
            A.shape, array=fractional_operator(A, M, s), comm=PETSc.COMM_SELF
        )
        A_petsc, M_petsc = petsc_matrix(A), petsc_matrix(M)
        ksp, residual = fgmres((Fs, A_petsc), haznicsx.HazmathRA(M_petsc, s))
        assert ksp.getConvergedReason() > 0
        assert residual < 1e-8
        its.append(ksp.getIterationNumber())
        ksp.destroy()
    assert max(its) <= 6


def test_petsc_copy(p1_1d):
    """The PETSc copy that the default scaling makes has PETSc's index type."""
    from haznicsx.precond import _petsc

    A, _ = p1_1d(7)
    A_petsc = _petsc(haznicsx.csr_from_scipy(A))
    indptr, indices, data = A_petsc.getValuesCSR()
    assert indptr.dtype == PETSc.IntType
    np.testing.assert_array_equal(indices, A.indices)
    np.testing.assert_allclose(data, A.data)
    A_petsc.destroy()


def test_default_scaling_needs_positive_eigenvalue(p1_1d):
    A, M = p1_1d(15)
    with pytest.raises(ValueError, match="needs a positive one"):
        haznicsx.default_scaling_a(haznicsx.csr_from_scipy(-A), haznicsx.csr_from_scipy(M))


def test_ra_pc_new_mass_matrix(tmp_path, p1_1d, petsc_matrix):
    """Replacing ctx.M is picked up at the next setup, by the scaling and by HAZmath."""
    A, M = p1_1d(31)
    A_petsc = petsc_matrix(A)
    ctx = haznicsx.HazmathRA(petsc_matrix(M), -0.5)
    ksp, _ = fgmres((A_petsc,), ctx)
    scaling = ctx._scaling_a_used
    # scaling M by 2 halves the eigenvalues of M^-1 A, so scaling_a doubles
    ctx.M = petsc_matrix(2.0 * M)
    ksp.setOperators(petsc_matrix(A))
    ksp.setUp()
    assert ctx._scaling_a_used == pytest.approx(2 * scaling, rel=1e-3)
    np.testing.assert_allclose(ctx._M.get(ctx.M).data, 2.0 * M.data)
    viewer = PETSc.Viewer().createASCII(str(tmp_path / "ksp.txt"), comm=PETSc.COMM_SELF)
    ksp.view(viewer)
    viewer.destroy()
    text = (tmp_path / "ksp.txt").read_text()
    assert "s = -0.5" in text
    assert "(computed)" in text
    ksp.destroy()


@pytest.mark.skipif(PETSc.COMM_WORLD.getSize() == 1, reason="needs more than one process")
def test_ra_pc_parallel_rejected(p1_1d, petsc_matrix):
    A, M = p1_1d(16)
    A_world = petsc_matrix(A, PETSc.COMM_WORLD)
    ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
    ksp.setOperators(A_world)
    pc = ksp.getPC()
    pc.setType("python")
    pc.setPythonContext(haznicsx.HazmathRA(haznicsx.csr_from_scipy(M), 0.5, scaling_a=1.0))
    with pytest.raises(PETSc.Error):
        ksp.setUp()
    ksp.destroy()
    A_world.destroy()
