# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Fixtures shared by the haznicsx tests."""

import collections.abc
import typing

import numpy as np
import numpy.typing as npt
import pytest
import scipy.sparse as sp


def _poisson_2d(n: int) -> sp.csr_matrix:
    T = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n))
    identity = sp.identity(n)
    return (sp.kron(T, identity) + sp.kron(identity, T)).tocsr()


@pytest.fixture
def poisson_2d() -> collections.abc.Callable[[int], sp.csr_matrix]:
    """Five-point Laplacian on an n x n grid, as a function of n."""
    return _poisson_2d


def _pcg(
    A: typing.Any,
    b: npt.NDArray[np.float64],
    precond: collections.abc.Callable,
    rtol: float = 1e-10,
    maxit: int = 1000,
) -> tuple[npt.NDArray[np.float64], int]:
    x = np.zeros_like(b)
    r = b.copy()
    z = np.zeros_like(b)
    precond(r, z)
    p = z.copy()
    rz = r @ z
    b_norm = np.linalg.norm(b)
    for k in range(1, maxit + 1):
        Ap = A @ p
        alpha = rz / (p @ Ap)
        x += alpha * p
        r -= alpha * Ap
        if np.linalg.norm(r) < rtol * b_norm:
            return x, k
        precond(r, z)
        rz_new = r @ z
        p = z + rz_new / rz * p
        rz = rz_new
    raise RuntimeError("PCG did not converge")


@pytest.fixture
def pcg() -> collections.abc.Callable:
    """Preconditioned CG, ``pcg(A, b, precond)``, returning the solution and iterations.

    ``precond(r, z)`` sets ``z`` to the preconditioned ``r``. PETSc is not
    needed.
    """
    return _pcg


def _p1_1d(n: int) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    h = 1.0 / (n + 1)
    A = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n)) / h
    M = sp.diags([1.0, 4.0, 1.0], [-1, 0, 1], shape=(n, n)) * (h / 6)
    return A.tocsr(), M.tocsr()


@pytest.fixture
def p1_1d() -> collections.abc.Callable[[int], tuple[sp.csr_matrix, sp.csr_matrix]]:
    """P1 stiffness and mass matrices on n interior nodes of (0, 1), as a function of n."""
    return _p1_1d


def _preconditioner_matrix(B: typing.Any) -> npt.NDArray[np.float64]:
    n = B.size
    # row j of Zt is B applied to e_j (a contiguous output array)
    Zt = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1.0
        B.apply(e, Zt[j])
    return Zt.T


@pytest.fixture
def preconditioner_matrix() -> collections.abc.Callable:
    """The dense matrix of a :class:`haznicsx.Precond`, as a function of it."""
    return _preconditioner_matrix


@pytest.fixture
def petsc_matrix() -> collections.abc.Iterator[collections.abc.Callable]:
    """``petsc_matrix(S, comm=None)``: AIJ matrix of a SciPy matrix.

    The rows are split evenly over the processes of ``comm`` (default
    ``PETSc.COMM_SELF``). The matrices are destroyed at the end of the test.
    """
    from petsc4py import PETSc

    created: list[PETSc.Mat] = []

    def create(S: sp.csr_matrix, comm: typing.Any = None) -> PETSc.Mat:
        comm = PETSc.COMM_SELF if comm is None else comm
        rank, size = comm.getRank(), comm.getSize()
        n = S.shape[0]
        start, end = rank * n // size, (rank + 1) * n // size
        local = S[start:end]
        M = PETSc.Mat().createAIJ(
            size=((end - start, n), (end - start, n)),
            csr=(local.indptr, local.indices, local.data),
            comm=comm,
        )
        M.assemble()
        created.append(M)
        return M

    yield create
    for M in created:
        M.destroy()
