# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Tests of the HX preconditioners with DOLFINx operators; PETSc is not needed."""

import numpy as np
import pytest
import scipy.sparse.linalg as spla

import haznicsx

pytest.importorskip("dolfinx")

from mpi4py import MPI

import ufl
from dolfinx import fem, mesh


def unit_mesh(dim: int, n: int) -> mesh.Mesh:
    if dim == 2:
        return mesh.create_unit_square(MPI.COMM_SELF, n, n)
    return mesh.create_unit_cube(MPI.COMM_SELF, n, n, n)


def scipy_matrix(B: haznicsx.CSRMatrix):
    import scipy.sparse as sp

    return sp.csr_matrix((B.data, B.indices, B.indptr), shape=B.shape)


def hcurl_problem(msh: mesh.Mesh):
    V = fem.functionspace(msh, ("N1curl", 1))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx + ufl.inner(u, v) * ufl.dx
    return V, fem.assemble_matrix(fem.form(a))


def hdiv_problem(msh: mesh.Mesh):
    V = fem.functionspace(msh, ("RT", 1))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.div(u), ufl.div(v)) * ufl.dx + ufl.inner(u, v) * ufl.dx
    return V, fem.assemble_matrix(fem.form(a))


def interpolate(V: fem.FunctionSpace, expr) -> np.ndarray:
    f = fem.Function(V)
    f.interpolate(fem.Expression(expr, V.element.interpolation_points))
    return f.x.array


@pytest.mark.parametrize("dim", [2, 3])
def test_hx_curl_operators(dim):
    msh = unit_mesh(dim, 3)
    V, _ = hcurl_problem(msh)
    Pcurl, Grad = haznicsx.hx_curl_operators(V)
    x = ufl.SpatialCoordinate(msh)
    P1 = fem.functionspace(msh, ("Lagrange", 1))
    phi = fem.Function(P1)
    phi.interpolate(lambda x: x[0] + 2 * x[1] - x[dim - 1])
    np.testing.assert_allclose(
        scipy_matrix(Grad) @ phi.x.array, interpolate(V, ufl.grad(phi)), atol=1e-12
    )
    # a linear vector field is in vector P1
    P1v = fem.functionspace(msh, ("Lagrange", 1, (dim,)))
    w = fem.Function(P1v)
    w.interpolate(
        fem.Expression(ufl.as_vector([x[1]] + [x[0]] * (dim - 1)), P1v.element.interpolation_points)
    )
    np.testing.assert_allclose(scipy_matrix(Pcurl) @ w.x.array, interpolate(V, w), atol=1e-12)


@pytest.mark.parametrize("dim", [2, 3])
def test_hx_div_operators(dim):
    msh = unit_mesh(dim, 3)
    V, _ = hdiv_problem(msh)
    Pdiv, Curl, Pcurl = haznicsx.hx_div_operators(V)
    x = ufl.SpatialCoordinate(msh)
    P1v = fem.functionspace(msh, ("Lagrange", 1, (dim,)))
    w = fem.Function(P1v)
    w.interpolate(
        fem.Expression(ufl.as_vector([x[1]] + [x[0]] * (dim - 1)), P1v.element.interpolation_points)
    )
    np.testing.assert_allclose(scipy_matrix(Pdiv) @ w.x.array, interpolate(V, w), atol=1e-12)
    if dim == 2:
        assert Pcurl is None
        P1 = fem.functionspace(msh, ("Lagrange", 1))
        phi = fem.Function(P1)
        phi.interpolate(lambda x: x[0] ** 2 + np.sin(3 * x[1]))
        rot = ufl.as_vector((phi.dx(1), -phi.dx(0)))
        np.testing.assert_allclose(
            scipy_matrix(Curl) @ phi.x.array, interpolate(V, rot), atol=1e-12
        )
    else:
        Ned = fem.functionspace(msh, ("N1curl", 1))
        e = fem.Function(Ned)
        e.interpolate(
            fem.Expression(ufl.as_vector((x[1], x[2], x[0])), Ned.element.interpolation_points)
        )
        np.testing.assert_allclose(
            scipy_matrix(Curl) @ e.x.array, interpolate(V, ufl.curl(e)), atol=1e-12
        )
        np.testing.assert_allclose(scipy_matrix(Pcurl) @ w.x.array, interpolate(Ned, w), atol=1e-12)


@pytest.mark.parametrize(
    "family, degree, operators",
    [
        ("Lagrange", 1, haznicsx.hx_curl_operators),
        ("N1curl", 2, haznicsx.hx_curl_operators),
        ("N1curl", 1, haznicsx.hx_div_operators),
    ],
)
def test_hx_operators_wrong_space(family, degree, operators):
    V = fem.functionspace(unit_mesh(2, 2), (family, degree))
    with pytest.raises(ValueError, match="lowest order"):
        operators(V)


def hx_preconditioner(kind: str, dim: int, n: int, multiplicative: bool = False):
    """The matrix, as SciPy, and its HX preconditioner on the unit square or cube."""
    msh = unit_mesh(dim, n)
    if kind == "curl":
        V, A = hcurl_problem(msh)
        B = haznicsx.hx_curl(
            haznicsx.csr_from_dolfinx(A),
            *haznicsx.hx_curl_operators(V),
            multiplicative=multiplicative,
        )
    else:
        V, A = hdiv_problem(msh)
        B = haznicsx.hx_div(
            haznicsx.csr_from_dolfinx(A),
            *haznicsx.hx_div_operators(V),
            multiplicative=multiplicative,
        )
    return A.to_scipy(), B


@pytest.mark.parametrize("kind", ["curl", "div"])
@pytest.mark.parametrize("dim, sizes", [(2, (8, 16, 32)), (3, (4, 8))], ids=["2D", "3D"])
def test_hx_additive_bounded_under_refinement(kind, dim, sizes, pcg):
    its = []
    for n in sizes:
        A, B = hx_preconditioner(kind, dim, n)
        b = np.random.default_rng(1).random(A.shape[0])
        _, k = pcg(A, b, B.apply, rtol=1e-8)
        its.append(k)
    # about 14 to 27 iterations; HEC aggregation is randomized
    assert max(its) < 40
    assert its[-1] < its[0] + 20


@pytest.mark.parametrize("kind", ["curl", "div"])
@pytest.mark.parametrize("dim, n", [(2, 16), (3, 6)], ids=["2D", "3D"])
def test_hx_multiplicative(kind, dim, n):
    A, B = hx_preconditioner(kind, dim, n, multiplicative=True)
    b = np.random.default_rng(2).random(A.shape[0])
    # the multiplicative variant changes r in HAZmath; haznicsx passes a copy
    b0 = b.copy()
    z = np.zeros_like(b)
    B.apply(b, z)
    np.testing.assert_array_equal(b, b0)

    def apply(r):
        # LinearOperator may pass column vectors
        r = np.ascontiguousarray(r, dtype=np.float64).ravel()
        out = np.zeros_like(r)
        B.apply(r, out)
        return out

    # not symmetric, so GMRES instead of CG
    M = spla.LinearOperator(A.shape, matvec=apply)
    its = []
    x, info = spla.gmres(
        A, b, M=M, rtol=1e-8, restart=100, callback=its.append, callback_type="pr_norm"
    )
    assert info == 0
    assert len(its) < 60
    np.testing.assert_allclose(A @ x, b, rtol=0, atol=1e-6 * np.linalg.norm(b))
