# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Tests of csr_from_dolfinx, which needs DOLFINx but not PETSc."""

import sys

from mpi4py import MPI

import numpy as np
import pytest

import haznicsx

pytest.importorskip("dolfinx")

import ufl
from dolfinx import fem, mesh


def mass_matrix(comm: MPI.Comm, shape: tuple[int, ...] = (), dtype=np.float64):
    msh = mesh.create_unit_square(comm, 4, 4, dtype=np.real(np.zeros(1, dtype)).dtype)
    V = fem.functionspace(msh, ("Lagrange", 1, shape))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    return fem.assemble_matrix(fem.form(ufl.inner(u, v) * ufl.dx, dtype=dtype))


def test_csr_from_dolfinx():
    A = mass_matrix(MPI.COMM_SELF)
    B = haznicsx.csr_from_dolfinx(A)
    S = A.to_scipy()
    assert B.shape == S.shape
    np.testing.assert_array_equal(B.indptr, S.indptr)
    np.testing.assert_array_equal(B.indices, S.indices)
    np.testing.assert_array_equal(B.data, S.data)
    assert np.shares_memory(B.data, A.data)


def test_csr_from_dolfinx_blocked():
    A = mass_matrix(MPI.COMM_SELF, shape=(2,))
    assert list(A.block_size) == [2, 2]
    B = haznicsx.csr_from_dolfinx(A)
    S = A.to_scipy().tocsr()
    assert B.shape == S.shape == (2 * A.index_map(0).size_local,) * 2
    np.testing.assert_array_equal(B.indptr, S.indptr)
    np.testing.assert_array_equal(B.indices, S.indices)
    np.testing.assert_array_equal(B.data, S.data)


def test_csr_from_dolfinx_single_precision():
    with pytest.raises(TypeError, match="double precision"):
        haznicsx.csr_from_dolfinx(mass_matrix(MPI.COMM_SELF, dtype=np.float32))


@pytest.mark.skipif(MPI.COMM_WORLD.size == 1, reason="needs more than one process")
def test_csr_from_dolfinx_parallel_rejected():
    with pytest.raises(ValueError, match="one process"):
        haznicsx.csr_from_dolfinx(mass_matrix(MPI.COMM_WORLD))


# None in sys.modules makes the import fail, as if the package were not installed
def test_csr_from_dolfinx_blocked_without_scipy(monkeypatch: pytest.MonkeyPatch):
    A = mass_matrix(MPI.COMM_SELF, shape=(2,))
    monkeypatch.setitem(sys.modules, "scipy.sparse", None)
    with pytest.raises(ImportError, match="pip install scipy"):
        haznicsx.csr_from_dolfinx(A)
