# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Tests of the HAZmath bindings that need only NumPy and SciPy."""

import gc
import sys
import threading

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import haznicsx

# HAZmath's index type: int32, or int64 when configured with int64=yes
INDEX = haznicsx.index_dtype
OTHER_INDEX = np.dtype(np.int64) if INDEX == np.int32 else np.dtype(np.int32)


def with_index_dtype(A: sp.csr_matrix) -> sp.csr_matrix:
    """A with index arrays of HAZmath's index type, which CSRMatrix shares."""
    A.indptr = A.indptr.astype(INDEX)
    A.indices = A.indices.astype(INDEX)
    return A


def test_csr_round_trip():
    A = sp.random(30, 20, density=0.2, format="csr", random_state=1)
    B = haznicsx.csr_from_scipy(A)
    assert B.shape == (30, 20)
    assert B.nnz == A.nnz
    np.testing.assert_array_equal(B.indptr, A.indptr)
    np.testing.assert_array_equal(B.indices, A.indices)
    np.testing.assert_array_equal(B.data, A.data)


def test_csr_arrays_are_read_only_views(poisson_2d):
    B = haznicsx.csr_from_scipy(poisson_2d(4))
    indptr = B.indptr
    assert not indptr.flags.writeable
    # the view keeps the matrix alive
    del B
    gc.collect()
    assert indptr[-1] == poisson_2d(4).nnz


def test_csr_shares_arrays():
    indptr = np.array([0, 1, 2], dtype=INDEX)
    indices = np.array([0, 1], dtype=INDEX)
    data = np.array([1.0, 2.0])
    B = haznicsx.CSRMatrix(indptr, indices, data, (2, 2))
    assert np.shares_memory(B.indptr, indptr)
    assert np.shares_memory(B.indices, indices)
    assert np.shares_memory(B.data, data)
    data[0] = 3.0
    np.testing.assert_array_equal(B.data, [3.0, 2.0])
    # the matrix keeps the arrays alive
    del indptr, indices, data
    gc.collect()
    np.testing.assert_array_equal(B.data, [3.0, 2.0])


def test_csr_constructor_does_not_convert():
    indptr = np.array([0, 1, 2], dtype=INDEX)
    indices = np.array([0, 1], dtype=INDEX)
    data = np.array([1.0, 2.0])
    # a cast from int64 to int32 would wrap indices around without a check
    with pytest.raises(TypeError):
        haznicsx.CSRMatrix(indptr, indices.astype(OTHER_INDEX), data, (2, 2))
    with pytest.raises(TypeError):
        haznicsx.CSRMatrix(indptr, indices, data.astype(np.float32), (2, 2))
    with pytest.raises(TypeError):
        haznicsx.CSRMatrix(indptr, np.array([[0, 9], [1, 9]], dtype=INDEX)[:, 0], data, (2, 2))


def test_csr_converts_other_arrays():
    indptr = np.array([0, 1, 2], dtype=OTHER_INDEX)
    indices = np.array([[0, 9], [1, 9]], dtype=INDEX)[:, 0]  # not contiguous
    data = np.array([1.0, 2.0], dtype=np.float32)
    B = haznicsx.csr_from_arrays(indptr, indices, data, (2, 2))
    assert not np.shares_memory(B.data, data)
    # the matrix keeps the converted arrays alive
    del indptr, indices, data
    gc.collect()
    np.testing.assert_array_equal(B.indptr, [0, 1, 2])
    np.testing.assert_array_equal(B.indices, [0, 1])
    np.testing.assert_array_equal(B.data, [1.0, 2.0])


def test_csr_from_scipy_shares_arrays(poisson_2d):
    A = with_index_dtype(poisson_2d(4))
    B = haznicsx.csr_from_scipy(A)
    assert np.shares_memory(B.indptr, A.indptr)
    assert np.shares_memory(B.indices, A.indices)
    assert np.shares_memory(B.data, A.data)


def test_csr_from_scipy_canonical_format():
    # row 0 has a duplicate diagonal entry, row 1 unsorted columns
    A = sp.csr_matrix(([1.0, 1.0, 3.0, 4.0], [0, 0, 1, 0], [0, 2, 4]), shape=(2, 2))
    assert not A.has_canonical_format
    B = haznicsx.csr_from_scipy(A)
    np.testing.assert_array_equal(B.indptr, [0, 1, 3])
    np.testing.assert_array_equal(B.indices, [0, 0, 1])
    np.testing.assert_array_equal(B.data, [2.0, 4.0, 3.0])
    # A is not changed
    assert A.nnz == 4
    np.testing.assert_array_equal(A.indices, [0, 0, 1, 0])


def test_csr_checked_again_before_setup(poisson_2d):
    A = with_index_dtype(poisson_2d(4))
    B = haznicsx.csr_from_scipy(A)
    A.indices[0] = 1000
    with pytest.raises(ValueError, match="Column index out of range"):
        haznicsx.amg(B, haznicsx.AMGParam())
    # eliminate_zeros changes the shared indptr in place
    A = with_index_dtype(poisson_2d(4))
    B = haznicsx.csr_from_scipy(A)
    A.data[1] = 0.0
    A.eliminate_zeros()
    with pytest.raises(ValueError, match="indptr"):
        haznicsx.amg(B, haznicsx.AMGParam())


@pytest.mark.parametrize(
    "indptr, indices, shape, message",
    [
        ([0, 1], [0], (2, 2), "indptr must have length"),
        ([1, 1, 2], [0, 1], (2, 2), "indptr\\[0\\] must be 0"),
        ([0, 2, 1], [0, 1], (2, 2), "non-decreasing"),
        ([0, 1, 3], [0, 1], (2, 2), "indices and data must have length"),
        ([0, 1, 2], [0, 2], (2, 2), "Column index out of range"),
        ([0, 1, 2], [0, -1], (2, 2), "Column index out of range"),
        ([0, 2, 2], [1, 0], (2, 2), "sorted"),
        ([0, 2, 2], [1, 1], (2, 2), "duplicates"),
        ([0], [], (0, 2), "dimensions must be positive"),
    ],
)
def test_csr_invalid(indptr, indices, shape, message):
    data = np.ones(len(indices))
    with pytest.raises(ValueError, match=message):
        haznicsx.csr_from_arrays(indptr, indices, data, shape)


def test_csr_complex_rejected():
    with pytest.raises(TypeError, match="real"):
        haznicsx.csr_from_arrays([0, 1], [0], np.array([1j]), (1, 1))


def test_csr_index_overflow():
    too_big = int(np.iinfo(INDEX).max) + 1
    with pytest.raises(OverflowError, match="indptr"):
        haznicsx.csr_from_arrays(np.array([0, too_big], dtype=np.uint64), [0], [1.0], (1, 1))
    with pytest.raises(OverflowError, match="dimensions"):
        haznicsx.csr_from_arrays([0, 1], [0], [1.0], (1, too_big))
    if INDEX == np.int32:
        # would wrap around to 0
        with pytest.raises(OverflowError, match="indices"):
            haznicsx.csr_from_arrays([0, 1], np.array([-(2**32)]), [1.0], (1, 1))


def test_csr_int64_indices(poisson_2d):
    A = poisson_2d(5)
    B = haznicsx.csr_from_arrays(
        A.indptr.astype(np.int64), A.indices.astype(np.int64), A.data, A.shape
    )
    np.testing.assert_array_equal(B.indices, A.indices)


def test_csr_from_scipy_not_sparse():
    with pytest.raises(TypeError, match="SciPy sparse"):
        haznicsx.csr_from_scipy(np.eye(3))  # type: ignore[arg-type]


# None in sys.modules makes the import fail, as if the package were not installed
def test_csr_from_scipy_without_scipy(monkeypatch: pytest.MonkeyPatch, poisson_2d):
    A = poisson_2d(3)
    monkeypatch.setitem(sys.modules, "scipy.sparse", None)
    with pytest.raises(ImportError, match="pip install scipy"):
        haznicsx.csr_from_scipy(A)


def test_csr_from_petsc_without_petsc4py(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setitem(sys.modules, "petsc4py", None)
    with pytest.raises(ImportError, match=r"haznicsx\[petsc\]"):
        haznicsx.csr_from_petsc(None)  # type: ignore[arg-type]


def test_amg_param_defaults():
    param = haznicsx.AMGParam()
    assert param.AMG_type == haznicsx.UA_AMG
    assert param.cycle_type == haznicsx.W_CYCLE
    assert param.max_levels == 10
    assert param.maxit == 1
    assert param.print_level == haznicsx.PRINT_NONE


def test_amg_param_fields():
    fields = haznicsx.AMGParam.fields
    assert "max_levels" in fields
    assert "fpwr" in fields
    assert "from_dict" not in fields
    param = haznicsx.AMGParam()
    for name in fields:
        getattr(param, name)


def test_default_param(poisson_2d):
    A = haznicsx.csr_from_scipy(poisson_2d(8))
    for B in (haznicsx.amg(A), haznicsx.famg(A, A)):
        z = np.zeros(A.shape[0])
        B.apply(np.ones(A.shape[0]), z)
        assert np.isfinite(z).all()


def test_amg_param_from_dict():
    param = haznicsx.AMGParam.from_dict({"max_levels": 5, "strong_coupled": 0.1})
    assert param.max_levels == 5
    assert param.strong_coupled == pytest.approx(0.1)
    with pytest.raises(AttributeError):
        haznicsx.AMGParam.from_dict({"no_such_parameter": 1})


@pytest.mark.parametrize("cycle_type, max_its", [("V_CYCLE", 45), ("W_CYCLE", 30)])
def test_amg_pcg(cycle_type, max_its, poisson_2d, pcg):
    A = poisson_2d(64)
    b = np.ones(A.shape[0])
    param = haznicsx.AMGParam.from_dict({"cycle_type": getattr(haznicsx, cycle_type)})
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), param)
    assert B.size == A.shape[0]
    x, its = pcg(A, b, B.apply)
    _, its_cg = pcg(A, b, lambda r, z: np.copyto(z, r))
    assert its <= max_its < its_cg
    np.testing.assert_allclose(x, spla.spsolve(A.tocsc(), b), rtol=1e-7)


def test_amg_iterations_bounded_under_refinement(poisson_2d, pcg):
    its = []
    for n in (32, 64, 128):
        A = poisson_2d(n)
        B = haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam())
        its.append(pcg(A, np.ones(A.shape[0]), B.apply)[1])
    # the W-cycle iteration count grows slowly; plain CG doubles
    assert its[-1] < 1.6 * its[0]


def test_apply_does_not_modify_r(poisson_2d):
    A = poisson_2d(16)
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam())
    r = np.random.default_rng(2).random(A.shape[0])
    r0 = r.copy()
    z = np.zeros_like(r)
    B.apply(r, z)
    np.testing.assert_array_equal(r, r0)
    assert r @ z > 0
    # a read-only r is accepted
    r.flags.writeable = False
    z2 = np.zeros_like(r)
    B.apply(r, z2)
    np.testing.assert_array_equal(z2, z)
    # r and z may be the same array
    x = r0.copy()
    B.apply(x, x)
    np.testing.assert_array_equal(x, z)


def test_apply_arguments(poisson_2d):
    A = poisson_2d(8)
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam())
    n = A.shape[0]
    with pytest.raises(ValueError, match="length"):
        B.apply(np.ones(n), np.zeros(n + 1))
    # z is written in place, so it is never converted
    with pytest.raises(TypeError):
        B.apply(np.ones(n), np.zeros(n, dtype=np.float32))
    with pytest.raises(TypeError):
        B.apply(np.ones(n), np.zeros(2 * n)[::2])
    # r may be converted
    z = np.zeros(n)
    B.apply(np.ones(n, dtype=np.int32), z)
    assert np.isfinite(z).all()


def test_amg_non_square():
    A = haznicsx.csr_from_scipy(sp.random(4, 3, density=0.5, format="csr", random_state=0))
    with pytest.raises(ValueError, match="square"):
        haznicsx.amg(A, haznicsx.AMGParam())


def test_amg_unsupported_smoother(poisson_2d):
    A = haznicsx.csr_from_scipy(poisson_2d(4))
    param = haznicsx.AMGParam.from_dict({"smoother": haznicsx.SMOOTHER_FGS})
    with pytest.raises(ValueError, match="smoother"):
        haznicsx.amg(A, param)


@pytest.mark.parametrize(
    "options, message",
    [
        ({"max_levels": 0}, "max_levels"),
        ({"max_levels": -1}, "max_levels"),
        ({"aggregation_type": 2}, "aggregation_type"),
        ({"AMG_type": haznicsx.SA_AMG, "aggregation_type": 2}, "aggregation_type"),
        ({"cycle_type": haznicsx.AMLI_CYCLE, "amli_degree": -1}, "amli_degree"),
    ],
)
def test_amg_invalid_param(options, message, poisson_2d):
    A = haznicsx.csr_from_scipy(poisson_2d(4))
    with pytest.raises(ValueError, match=message):
        haznicsx.amg(A, haznicsx.AMGParam.from_dict(options))
    options = {"smoother": haznicsx.SMOOTHER_FGS, **options}
    with pytest.raises(ValueError, match=message):
        haznicsx.famg(A, A, haznicsx.AMGParam.from_dict(options))


@pytest.mark.parametrize("aggregation_type", ["VMB", "HEC", "HEM"])
def test_amg_aggregation_types(aggregation_type, poisson_2d, pcg):
    A = poisson_2d(16)
    param = haznicsx.AMGParam.from_dict({"aggregation_type": getattr(haznicsx, aggregation_type)})
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), param)
    _, its = pcg(A, np.ones(A.shape[0]), B.apply)
    assert its < 30


def test_amg_param_not_modified(poisson_2d):
    param = haznicsx.AMGParam.from_dict({"cycle_type": haznicsx.AMLI_CYCLE})
    haznicsx.amg(haznicsx.csr_from_scipy(poisson_2d(16)), param)
    assert param.cycle_type == haznicsx.AMLI_CYCLE
    assert param.max_levels == 10


@pytest.mark.parametrize("cycle_type", ["V_CYCLE", "W_CYCLE", "AMLI_CYCLE", "NL_AMLI_CYCLE"])
def test_amg_cycles(cycle_type, poisson_2d):
    A = poisson_2d(32)
    param = haznicsx.AMGParam.from_dict({"cycle_type": getattr(haznicsx, cycle_type)})
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), param)
    z = np.zeros(A.shape[0])
    B.apply(np.ones(A.shape[0]), z)
    assert np.isfinite(z).all()


def test_famg(poisson_2d, pcg):
    n = 32
    A = poisson_2d(n)
    M = sp.identity(A.shape[0], format="csr") / (n + 1) ** 2
    param = haznicsx.AMGParam.from_dict({"smoother": haznicsx.SMOOTHER_FGS, "fpwr": 1.0})
    B = haznicsx.famg(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), param)
    _, its = pcg(A, np.ones(A.shape[0]), B.apply)
    assert its < 30


def test_famg_unsupported_smoother(poisson_2d):
    A = haznicsx.csr_from_scipy(poisson_2d(4))
    with pytest.raises(ValueError, match="smoother"):
        haznicsx.famg(A, A, haznicsx.AMGParam())


def test_precond_outlives_matrix(poisson_2d):
    A = poisson_2d(16)
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam())
    gc.collect()
    z = np.zeros(A.shape[0])
    B.apply(np.ones(A.shape[0]), z)
    assert np.isfinite(z).all()
    # the setup copied the matrix, so changing A's arrays does not change B
    A.data[:] = 0.0
    del A
    gc.collect()
    z2 = np.zeros_like(z)
    B.apply(np.ones(z.size), z2)
    np.testing.assert_array_equal(z2, z)


def test_create_and_free_repeatedly(poisson_2d):
    A = poisson_2d(16)
    for cycle_type in (haznicsx.V_CYCLE, haznicsx.AMLI_CYCLE) * 10:
        param = haznicsx.AMGParam.from_dict({"cycle_type": cycle_type})
        B = haznicsx.amg(haznicsx.csr_from_scipy(A), param)
        z = np.zeros(A.shape[0])
        B.apply(np.ones(A.shape[0]), z)
        del B
    gc.collect()


def test_random_seed(poisson_2d):
    """HEC aggregation with the same seed gives the same preconditioner."""
    A = poisson_2d(16)
    r = np.random.default_rng(5).random(A.shape[0])
    param = haznicsx.AMGParam.from_dict({"aggregation_type": haznicsx.HEC, "random_seed": 7})
    assert param.random_seed == 7
    z = [np.zeros_like(r), np.zeros_like(r)]
    for zi in z:
        haznicsx.amg(haznicsx.csr_from_scipy(A), param).apply(r, zi)
    np.testing.assert_array_equal(z[0], z[1])


def test_threads(poisson_2d):
    """Applies from several threads, while other setups run, give the serial results."""
    A = poisson_2d(16)
    r = np.random.default_rng(3).random(A.shape[0])
    preconds = [haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam()) for _ in range(3)]
    expected = []
    for B in preconds:
        z = np.zeros_like(r)
        B.apply(r, z)
        expected.append(z)
    mismatches: list[int] = []

    def work(i: int) -> None:
        haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam())
        for _ in range(5):
            for B, z_expected in zip(preconds, expected, strict=True):
                z = np.zeros_like(r)
                B.apply(r, z)
                if not np.array_equal(z, z_expected):
                    mismatches.append(i)

    threads = [threading.Thread(target=work, args=(i,)) for i in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not mismatches


def random_csr(rows: int, cols: int) -> haznicsx.CSRMatrix:
    return haznicsx.csr_from_scipy(sp.random(rows, cols, density=0.5, format="csr", random_state=0))


def test_hx_curl_arguments(poisson_2d):
    A = haznicsx.csr_from_scipy(poisson_2d(4))
    n = A.shape[0]
    with pytest.raises(ValueError, match="Pcurl must have 16 rows"):
        haznicsx.hx_curl(A, random_csr(n - 1, 8), random_csr(n, 4))
    with pytest.raises(ValueError, match="Grad must have 16 rows"):
        haznicsx.hx_curl(A, random_csr(n, 8), random_csr(n + 1, 4))
    param = haznicsx.AMGParam.from_dict({"cycle_type": haznicsx.AMLI_CYCLE})
    with pytest.raises(ValueError, match="cycle_type"):
        haznicsx.hx_curl(A, random_csr(n, 8), random_csr(n, 4), param)


def test_hx_div_arguments(poisson_2d):
    A = haznicsx.csr_from_scipy(poisson_2d(4))
    n = A.shape[0]
    with pytest.raises(ValueError, match="Pdiv must have 16 rows"):
        haznicsx.hx_div(A, random_csr(n - 1, 8), random_csr(n, 6))
    with pytest.raises(ValueError, match="Curl must have 16 rows"):
        haznicsx.hx_div(A, random_csr(n, 8), random_csr(n - 1, 6))
    # 3D: Pcurl maps vector P1 into the space that Curl maps from
    with pytest.raises(ValueError, match="Pcurl must have 6 rows"):
        haznicsx.hx_div(A, random_csr(n, 8), random_csr(n, 6), random_csr(5, 8))
    param = haznicsx.AMGParam.from_dict({"cycle_type": haznicsx.NL_AMLI_CYCLE})
    with pytest.raises(ValueError, match="cycle_type"):
        haznicsx.hx_div(A, random_csr(n, 8), random_csr(n, 6), param=param)


needs_lapack = pytest.mark.skipif(
    not haznicsx.HAS_LAPACK, reason="the RA preconditioner needs libhazmath with LAPACK"
)


@pytest.mark.skipif(haznicsx.HAS_LAPACK, reason="libhazmath has LAPACK")
def test_ra_without_lapack(p1_1d):
    A, M = p1_1d(7)
    with pytest.raises(RuntimeError, match="LAPACK"):
        haznicsx.ra(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), 0.5, scaling_a=1.0)


@needs_lapack
@pytest.mark.parametrize(
    "s, t, alpha, beta, rtol",
    [
        (0.5, 0.0, 1.0, 0.0, 1e-6),
        # s < 0: the inner PCG solves to 1e-6 limit the accuracy
        (-0.5, 0.0, 1.0, 0.0, 1e-4),
        (0.5, 1.0, 1.0, 0.1, 1e-6),
    ],
)
def test_ra_inverse_of_fractional_operator(s, t, alpha, beta, rtol, p1_1d, preconditioner_matrix):
    """B is the inverse of alpha A^s + beta A^t in the M inner product."""
    import scipy.linalg

    A, M = p1_1d(31)
    lam, V = scipy.linalg.eigh(A.toarray(), M.toarray())
    exact = V @ np.diag(1 / (alpha * lam**s + beta * lam**t)) @ V.T
    B = haznicsx.ra(
        haznicsx.csr_from_scipy(A),
        haznicsx.csr_from_scipy(M),
        s,
        t,
        alpha,
        beta,
        scaling_a=1 / (1.1 * lam[-1]),
    )
    Z = preconditioner_matrix(B)
    assert np.linalg.norm(Z - exact) / np.linalg.norm(exact) < rtol


@needs_lapack
def test_ra_apply(p1_1d):
    A, M = p1_1d(15)
    B = haznicsx.ra(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), -0.5, scaling_a=1e-3)
    r = np.random.default_rng(4).random(B.size)
    r0 = r.copy()
    z = np.zeros_like(r)
    # HAZmath scales r in place; haznicsx passes a copy
    B.apply(r, z)
    np.testing.assert_array_equal(r, r0)
    # r and z may be the same array
    B.apply(r, r)
    np.testing.assert_allclose(r, z)


@needs_lapack
@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"alpha": -1.0}, "alpha and beta"),
        ({"alpha": 0.0, "beta": 0.0}, "alpha and beta"),
        ({"scaling_a": 0.0}, "scaling_a"),
        ({"scaling_m": -1.0}, "scaling_m"),
        ({"ra_tol": 0.0}, "ra_tol"),
        ({"s": 1.5}, "s and t"),
        ({"t": -2.0}, "s and t"),
        ({"alpha": float("nan")}, "finite"),
        ({"s": float("nan")}, "finite"),
        ({"scaling_a": float("inf")}, "scaling_a"),
        ({"scaling_m": float("nan")}, "finite"),
        ({"param": haznicsx.AMGParam.from_dict({"cycle_type": haznicsx.AMLI_CYCLE})}, "cycle_type"),
    ],
)
def test_ra_arguments(kwargs, message, p1_1d):
    A, M = p1_1d(7)
    kwargs = {"s": 0.5, "scaling_a": 1.0, **kwargs}
    with pytest.raises(ValueError, match=message):
        haznicsx.ra(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), **kwargs)


@needs_lapack
def test_ra_sizes(p1_1d):
    A, _ = p1_1d(7)
    with pytest.raises(ValueError, match="M must have 7 rows"):
        haznicsx.ra(
            haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(p1_1d(6)[1]), 0.5, scaling_a=1.0
        )


@needs_lapack
def test_ra_checked_before_scaling(p1_1d, monkeypatch: pytest.MonkeyPatch):
    """Invalid arguments are rejected before the eigenvalue solve for scaling_a."""
    A, M = p1_1d(7)
    # computing the scaling would raise ImportError
    monkeypatch.setitem(sys.modules, "slepc4py", None)
    with pytest.raises(ValueError, match="M must have 7 rows"):
        haznicsx.ra(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(p1_1d(6)[1]), 0.5)
    with pytest.raises(ValueError, match="s and t"):
        haznicsx.ra(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), 1.5)


@needs_lapack
@pytest.mark.parametrize("missing", ["petsc4py", "slepc4py"])
def test_ra_default_scaling_needs_slepc(missing, monkeypatch: pytest.MonkeyPatch, p1_1d):
    A, M = p1_1d(7)
    monkeypatch.setitem(sys.modules, missing, None)
    with pytest.raises(ImportError, match="pass scaling_a"):
        haznicsx.ra(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), 0.5)
    # the PC context checks when it is created
    with pytest.raises(ImportError, match="pass scaling_a"):
        haznicsx.HazmathRA(M, 0.5)
