"""Conversion of PETSc matrices. Runs with every PETSc build (PETSC_ARCH)."""

import numpy as np
import pytest
import scipy.sparse as sp

import haznicsx

PETSc = pytest.importorskip("petsc4py.PETSc")

real64 = np.dtype(PETSc.ScalarType) == np.float64


def petsc_matrix(A: sp.csr_matrix, mat_type: str = "aij", block_size: int = 1):
    """Assembled sequential PETSc matrix with the entries of A."""
    M = PETSc.Mat().create(comm=PETSc.COMM_SELF)
    M.setSizes(A.shape, bsize=block_size)
    M.setType(mat_type)
    M.setUp()
    C = A.tocoo()
    for i, j, value in zip(C.row, C.col, C.data):
        M.setValue(int(i), int(j), value)
    M.assemble()
    return M


def laplace_1d(n: int) -> sp.csr_matrix:
    return sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n), format="csr")


@pytest.mark.skipif(not real64, reason="needs real double precision PETSc")
@pytest.mark.parametrize("mat_type, block_size", [("aij", 1), ("baij", 2), ("sbaij", 2)])
def test_petsc_round_trip(mat_type, block_size):
    A = laplace_1d(10)
    if mat_type == "sbaij":
        # SBAIJ stores the upper triangle only
        M = petsc_matrix(sp.triu(A).tocsr(), mat_type, block_size)
    else:
        M = petsc_matrix(A, mat_type, block_size)
    B = haznicsx.csr_from_petsc(M)
    assert B.shape == A.shape
    C = sp.csr_matrix((B.data, B.indices, B.indptr), shape=B.shape)
    np.testing.assert_allclose(C.toarray(), A.toarray())


@pytest.mark.skipif(not real64, reason="needs real double precision PETSc")
def test_petsc_amg():
    A = laplace_1d(50)
    B = haznicsx.amg(haznicsx.csr_from_petsc(petsc_matrix(A)), haznicsx.AMGParam())
    z = np.zeros(50)
    B.apply(np.ones(50), z)
    assert np.isfinite(z).all()


@pytest.mark.skipif(real64, reason="for complex and single precision PETSc builds")
def test_petsc_scalar_type_rejected():
    M = PETSc.Mat().createAIJ((2, 2), nnz=1, comm=PETSc.COMM_SELF)
    M.setValue(0, 0, 1.0)
    M.setValue(1, 1, 1.0)
    M.assemble()
    with pytest.raises(TypeError, match="real double precision"):
        haznicsx.csr_from_petsc(M)


@pytest.mark.skipif(PETSc.COMM_WORLD.getSize() == 1, reason="needs more than one process")
def test_petsc_parallel_rejected():
    M = PETSc.Mat().createAIJ((4 * PETSc.COMM_WORLD.getSize(),) * 2, nnz=1)
    M.assemble()
    with pytest.raises(ValueError, match="one process"):
        haznicsx.csr_from_petsc(M)
