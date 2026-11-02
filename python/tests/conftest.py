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
