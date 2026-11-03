# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Tests of the petsc4py preconditioners (haznicsx.pc) without DOLFINx."""

import collections.abc

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import haznicsx

PETSc = pytest.importorskip("petsc4py.PETSc")

pytestmark = pytest.mark.skipif(
    np.dtype(PETSc.ScalarType) != np.float64, reason="HAZmath needs real double precision"
)


def petsc_matrix(A: sp.csr_matrix, comm=None):
    """AIJ matrix with the rows of A split evenly over the processes of comm."""
    comm = PETSc.COMM_SELF if comm is None else comm
    rank, size = comm.getRank(), comm.getSize()
    n = A.shape[0]
    start, end = rank * n // size, (rank + 1) * n // size
    local = A[start:end]
    M = PETSc.Mat().createAIJ(
        size=((end - start, n), (end - start, n)),
        csr=(local.indptr, local.indices, local.data),
        comm=comm,
    )
    M.assemble()
    return M


def cg(A, configure: collections.abc.Callable, prefix: str = ""):
    """CG solve of A x = 1 with the PC set up by configure(ksp).

    Returns the KSP, the solution and the right-hand side; the caller destroys
    the KSP.
    """
    ksp = PETSc.KSP().create(A.getComm())
    ksp.setOptionsPrefix(prefix)
    ksp.setOperators(A)
    ksp.setType("cg")
    ksp.setTolerances(rtol=1e-10, max_it=200)
    configure(ksp)
    b = A.createVecLeft()
    b.set(1.0)
    x = A.createVecRight()
    ksp.solve(b, x)
    return ksp, x, b


def use_context(ctx):
    def configure(ksp):
        pc = ksp.getPC()
        pc.setType("python")
        pc.setPythonContext(ctx)
        ksp.setFromOptions()

    return configure


@pytest.fixture
def options():
    """The global options database, cleared of the options a test sets."""
    opts = PETSc.Options()
    before = set(opts.getAll())
    yield opts
    for name in set(opts.getAll()) - before:
        opts.delValue(name)


def test_amg_pc(tmp_path, poisson_2d):
    S = poisson_2d(32)
    A = petsc_matrix(S)
    ksp, x, _ = cg(A, use_context(haznicsx.HazmathAMG()))
    assert ksp.getConvergedReason() > 0
    assert ksp.getIterationNumber() < 25
    np.testing.assert_allclose(x.array, spla.spsolve(S.tocsc(), np.ones(S.shape[0])), rtol=1e-8)
    viewer = PETSc.Viewer().createASCII(str(tmp_path / "ksp.txt"), comm=PETSc.COMM_SELF)
    ksp.view(viewer)
    viewer.destroy()
    assert f"cycle_type = {haznicsx.W_CYCLE}" in (tmp_path / "ksp.txt").read_text()
    ksp.destroy()


def test_amg_pc_from_options(options, poisson_2d):
    options["hz_pc_type"] = "python"
    options["hz_pc_python_type"] = "haznicsx.HazmathAMG"
    options["hz_pc_hazmath_cycle_type"] = "v"
    options["hz_pc_hazmath_max_levels"] = 4
    options["hz_pc_hazmath_smoother"] = "SGS"
    options["hz_pc_hazmath_relaxation"] = 0.9
    A = petsc_matrix(poisson_2d(32))
    ksp, _, _ = cg(A, lambda ksp: ksp.setFromOptions(), prefix="hz_")
    param = ksp.getPC().getPythonContext().param
    assert param.cycle_type == haznicsx.V_CYCLE
    assert param.max_levels == 4
    assert param.smoother == haznicsx.SMOOTHER_SGS
    assert param.relaxation == 0.9
    assert ksp.getConvergedReason() > 0
    ksp.destroy()


def test_amg_param_from_options(options):
    param = haznicsx.AMGParam.from_dict({"max_levels": 5})
    options["x_pc_hazmath_cycle_type"] = "1"
    options["x_pc_hazmath_aggregation_type"] = "vmb"
    new = haznicsx.amg_param_from_options(param, PETSc.Options("x_"))
    assert new.cycle_type == 1
    assert new.aggregation_type == haznicsx.VMB
    assert new.max_levels == 5
    # param is not changed
    assert param.cycle_type == haznicsx.W_CYCLE
    options["x_pc_hazmath_cycle_type"] = "bogus"
    with pytest.raises(ValueError, match="cycle_type bogus"):
        haznicsx.amg_param_from_options(param, PETSc.Options("x_"))
    # cycle_type is a C short
    options["x_pc_hazmath_cycle_type"] = "70000"
    with pytest.raises(ValueError, match="out of range"):
        haznicsx.amg_param_from_options(param, PETSc.Options("x_"))


def test_amg_pc_view_non_ascii(tmp_path, poisson_2d):
    ksp, _, _ = cg(petsc_matrix(poisson_2d(8)), use_context(haznicsx.HazmathAMG()))
    viewer = PETSc.Viewer().createBinary(str(tmp_path / "ksp.dat"), "w", comm=PETSc.COMM_SELF)
    ksp.view(viewer)
    viewer.destroy()
    ksp.destroy()


def test_amg_pc_rebuilt_for_new_operator(poisson_2d):
    S = poisson_2d(16)
    ctx = haznicsx.HazmathAMG()
    ksp, x, b = cg(petsc_matrix(S), use_context(ctx))
    first = ctx._precond
    ksp.setOperators(petsc_matrix(2.0 * S))
    y = x.duplicate()
    ksp.solve(b, y)
    assert ctx._precond is not first
    np.testing.assert_allclose(y.array, x.array / 2.0, rtol=1e-8)
    ksp.destroy()


def test_famg_pc(poisson_2d):
    n = 32
    S = poisson_2d(n)
    M = petsc_matrix(sp.identity(S.shape[0], format="csr") / (n + 1) ** 2)
    param = haznicsx.AMGParam.from_dict({"smoother": haznicsx.SMOOTHER_FGS, "fpwr": 1.0})
    ctx = haznicsx.HazmathFAMG(M, param)
    ksp, x, b = cg(petsc_matrix(S), use_context(ctx))
    assert ksp.getConvergedReason() > 0
    assert ksp.getIterationNumber() < 30
    # M is converted once, not at every setup
    first = ctx._M
    ksp.setOperators(petsc_matrix(2.0 * S))
    ksp.solve(b, x)
    assert ctx._M is first
    ksp.destroy()


def test_famg_pc_needs_mass_matrix(poisson_2d):
    ksp = PETSc.KSP().create(PETSc.COMM_SELF)
    ksp.setOperators(petsc_matrix(poisson_2d(4)))
    pc = ksp.getPC()
    pc.setType("python")
    pc.setPythonContext(haznicsx.HazmathFAMG())
    # older petsc4py (e.g. 3.19) passes the ValueError on unwrapped
    with pytest.raises((PETSc.Error, ValueError)):
        ksp.setUp()
    ksp.destroy()


@pytest.mark.skipif(PETSc.COMM_WORLD.getSize() == 1, reason="needs more than one process")
def test_pc_parallel_rejected(poisson_2d):
    ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
    ksp.setOperators(petsc_matrix(poisson_2d(8), PETSc.COMM_WORLD))
    pc = ksp.getPC()
    pc.setType("python")
    pc.setPythonContext(haznicsx.HazmathAMG())
    with pytest.raises(PETSc.Error):
        ksp.setUp()
    ksp.destroy()


def test_block_jacobi(options, poisson_2d):
    """HAZmath AMG on each process's diagonal block, as the module docstring suggests."""
    options["bj_pc_type"] = "bjacobi"
    options["bj_sub_pc_type"] = "python"
    options["bj_sub_pc_python_type"] = "haznicsx.HazmathAMG"
    S = poisson_2d(32)
    A = petsc_matrix(S, PETSc.COMM_WORLD)
    ksp, x, b = cg(A, lambda ksp: ksp.setFromOptions(), prefix="bj_")
    assert ksp.getConvergedReason() > 0
    r = b.duplicate()
    A.mult(x, r)
    r.axpy(-1.0, b)
    assert r.norm() < 1e-8 * b.norm()
    ksp.destroy()
