# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Fixtures shared by the haznicsx tests."""

import collections.abc

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
