"""Tests of the HAZmath bindings that need only NumPy and SciPy."""

import gc
import sys

import numpy as np
import numpy.typing as npt
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import haznicsx


def poisson_2d(n: int) -> sp.csr_matrix:
    """Five-point Laplacian on an n x n grid."""
    T = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n))
    identity = sp.identity(n)
    return (sp.kron(T, identity) + sp.kron(identity, T)).tocsr()


def pcg(A, b: npt.NDArray[np.float64], precond, rtol: float = 1e-10, maxit: int = 1000):
    """Preconditioned conjugate gradients. Returns the solution and iteration count."""
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


def test_csr_round_trip():
    A = sp.random(30, 20, density=0.2, format="csr", random_state=1)
    B = haznicsx.csr_from_scipy(A)
    assert B.shape == (30, 20)
    assert B.nnz == A.nnz
    np.testing.assert_array_equal(B.indptr, A.indptr)
    np.testing.assert_array_equal(B.indices, A.indices)
    np.testing.assert_array_equal(B.data, A.data)


def test_csr_arrays_are_read_only_views():
    B = haznicsx.csr_from_scipy(poisson_2d(4))
    indptr = B.indptr
    assert not indptr.flags.writeable
    # the view keeps the matrix alive
    del B
    gc.collect()
    assert indptr[-1] == poisson_2d(4).nnz


def test_csr_shares_arrays():
    indptr = np.array([0, 1, 2], dtype=np.int32)
    indices = np.array([0, 1], dtype=np.int32)
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


def test_csr_converts_other_arrays():
    indptr = np.array([0, 1, 2], dtype=np.int64)
    indices = np.array([[0, 9], [1, 9]], dtype=np.int32)[:, 0]  # not contiguous
    data = np.array([1.0, 2.0], dtype=np.float32)
    B = haznicsx.CSRMatrix(indptr, indices, data, (2, 2))
    assert not np.shares_memory(B.data, data)
    # the matrix keeps the converted arrays alive
    del indptr, indices, data
    gc.collect()
    np.testing.assert_array_equal(B.indptr, [0, 1, 2])
    np.testing.assert_array_equal(B.indices, [0, 1])
    np.testing.assert_array_equal(B.data, [1.0, 2.0])


def test_csr_from_scipy_shares_arrays():
    A = poisson_2d(4)
    B = haznicsx.csr_from_scipy(A)
    assert np.shares_memory(B.indices, A.indices)
    assert np.shares_memory(B.data, A.data)


@pytest.mark.parametrize(
    "indptr, indices, shape, message",
    [
        ([0, 1], [0], (2, 2), "indptr must have length"),
        ([1, 1, 2], [0, 1], (2, 2), "indptr\\[0\\] must be 0"),
        ([0, 2, 1], [0, 1], (2, 2), "non-decreasing"),
        ([0, 1, 3], [0, 1], (2, 2), "indices and data must have length"),
        ([0, 1, 2], [0, 2], (2, 2), "Column index out of range"),
        ([0, 1, 2], [0, -1], (2, 2), "Column index out of range"),
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
    big = np.array([0, 2**31], dtype=np.int64)
    with pytest.raises(OverflowError, match="indptr"):
        haznicsx.csr_from_arrays(big, [0], [1.0], (1, 1))
    # would wrap around to 0
    with pytest.raises(OverflowError, match="indices"):
        haznicsx.csr_from_arrays([0, 1], np.array([-(2**32)]), [1.0], (1, 1))
    with pytest.raises(OverflowError, match="dimensions"):
        haznicsx.csr_from_arrays([0, 1], [0], [1.0], (1, 2**31))


def test_csr_int64_indices():
    A = poisson_2d(5)
    B = haznicsx.csr_from_arrays(
        A.indptr.astype(np.int64), A.indices.astype(np.int64), A.data, A.shape
    )
    np.testing.assert_array_equal(B.indices, A.indices)


def test_csr_from_scipy_not_sparse():
    with pytest.raises(TypeError, match="SciPy sparse"):
        haznicsx.csr_from_scipy(np.eye(3))  # type: ignore[arg-type]


# None in sys.modules makes the import fail, as if the package were not installed
def test_csr_from_scipy_without_scipy(monkeypatch: pytest.MonkeyPatch):
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


def test_amg_param_from_dict():
    param = haznicsx.AMGParam.from_dict({"max_levels": 5, "strong_coupled": 0.1})
    assert param.max_levels == 5
    assert param.strong_coupled == pytest.approx(0.1)
    with pytest.raises(AttributeError):
        haznicsx.AMGParam.from_dict({"no_such_parameter": 1})


@pytest.mark.parametrize("cycle_type, max_its", [("V_CYCLE", 45), ("W_CYCLE", 30)])
def test_amg_pcg(cycle_type, max_its):
    A = poisson_2d(64)
    b = np.ones(A.shape[0])
    param = haznicsx.AMGParam.from_dict({"cycle_type": getattr(haznicsx, cycle_type)})
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), param)
    assert B.size == A.shape[0]
    x, its = pcg(A, b, B.apply)
    _, its_cg = pcg(A, b, lambda r, z: np.copyto(z, r))
    assert its <= max_its < its_cg
    np.testing.assert_allclose(x, spla.spsolve(A.tocsc(), b), rtol=1e-7)


def test_amg_iterations_bounded_under_refinement():
    its = []
    for n in (32, 64, 128):
        A = poisson_2d(n)
        B = haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam())
        its.append(pcg(A, np.ones(A.shape[0]), B.apply)[1])
    # the W-cycle iteration count grows slowly; plain CG doubles
    assert its[-1] < 1.6 * its[0]


def test_apply_does_not_modify_r():
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


def test_apply_arguments():
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


def test_amg_unsupported_smoother():
    A = haznicsx.csr_from_scipy(poisson_2d(4))
    param = haznicsx.AMGParam.from_dict({"smoother": haznicsx.SMOOTHER_FGS})
    with pytest.raises(ValueError, match="smoother"):
        haznicsx.amg(A, param)


def test_amg_param_not_modified():
    param = haznicsx.AMGParam.from_dict({"cycle_type": haznicsx.AMLI_CYCLE})
    haznicsx.amg(haznicsx.csr_from_scipy(poisson_2d(16)), param)
    assert param.cycle_type == haznicsx.AMLI_CYCLE
    assert param.max_levels == 10


@pytest.mark.parametrize("cycle_type", ["V_CYCLE", "W_CYCLE", "AMLI_CYCLE", "NL_AMLI_CYCLE"])
def test_amg_cycles(cycle_type):
    A = poisson_2d(32)
    param = haznicsx.AMGParam.from_dict({"cycle_type": getattr(haznicsx, cycle_type)})
    B = haznicsx.amg(haznicsx.csr_from_scipy(A), param)
    z = np.zeros(A.shape[0])
    B.apply(np.ones(A.shape[0]), z)
    assert np.isfinite(z).all()


def test_famg():
    n = 32
    A = poisson_2d(n)
    M = sp.identity(A.shape[0], format="csr") / (n + 1) ** 2
    param = haznicsx.AMGParam.from_dict({"smoother": haznicsx.SMOOTHER_FGS, "fpwr": 1.0})
    B = haznicsx.famg(haznicsx.csr_from_scipy(A), haznicsx.csr_from_scipy(M), param)
    _, its = pcg(A, np.ones(A.shape[0]), B.apply)
    assert its < 30


def test_famg_unsupported_smoother():
    A = haznicsx.csr_from_scipy(poisson_2d(4))
    with pytest.raises(ValueError, match="smoother"):
        haznicsx.famg(A, A, haznicsx.AMGParam())


def test_precond_outlives_matrix():
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


def test_create_and_free_repeatedly():
    A = poisson_2d(16)
    for cycle_type in (haznicsx.V_CYCLE, haznicsx.AMLI_CYCLE) * 10:
        param = haznicsx.AMGParam.from_dict({"cycle_type": cycle_type})
        B = haznicsx.amg(haznicsx.csr_from_scipy(A), param)
        z = np.zeros(A.shape[0])
        B.apply(np.ones(A.shape[0]), z)
        del B
    gc.collect()
