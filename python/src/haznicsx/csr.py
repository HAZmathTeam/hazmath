"""Conversion of sparse matrices to HAZmath's CSR format.

HAZmath uses 32-bit indices and double precision values. The functions here
convert the arrays only when needed and check that the indices fit. The
:class:`haznicsx.CSRMatrix` shares arrays that are already int32 and float64
and C-contiguous, so changing them afterwards changes the matrix.
"""

from __future__ import annotations

import typing

import numpy as np
import numpy.typing as npt

from haznicsx._haznicsx import CSRMatrix

if typing.TYPE_CHECKING:
    from petsc4py import PETSc

__all__ = ["csr_from_arrays", "csr_from_petsc", "csr_from_scipy"]

_INT32 = np.iinfo(np.int32)


class _SupportsToCSR(typing.Protocol):
    """A SciPy sparse matrix or array."""

    def tocsr(self) -> typing.Any: ...


def _as_int32(a: npt.ArrayLike, name: str) -> npt.NDArray[np.int32]:
    a = np.asarray(a)
    if a.dtype != np.int32 and a.size > 0 and (a.max() > _INT32.max or a.min() < _INT32.min):
        raise OverflowError(f"{name} does not fit in the 32-bit integers that HAZmath uses")
    return np.ascontiguousarray(a, dtype=np.int32)


def csr_from_arrays(
    indptr: npt.ArrayLike,
    indices: npt.ArrayLike,
    data: npt.ArrayLike,
    shape: tuple[int, int],
) -> CSRMatrix:
    """HAZmath matrix from CSR arrays.

    The matrix shares the arrays that need no conversion. The column indices
    must be sorted in each row, without duplicates.

    Args:
        indptr: Row pointers, length ``shape[0] + 1``.
        indices: Column indices.
        data: Real values; they are converted to float64.
        shape: Number of rows and columns.

    Raises:
        TypeError: If the values are complex.
        OverflowError: If an index or dimension does not fit in 32 bits.
        ValueError: If the arrays are not a valid CSR matrix.
    """
    values = np.asarray(data)
    if np.iscomplexobj(values):
        raise TypeError("HAZmath supports real matrices only")
    num_rows, num_cols = (int(s) for s in shape)
    if max(num_rows, num_cols) > _INT32.max:
        raise OverflowError("The matrix dimensions do not fit in 32-bit integers")
    return CSRMatrix(
        _as_int32(indptr, "indptr"),
        _as_int32(indices, "indices"),
        np.ascontiguousarray(values, dtype=np.float64),
        (num_rows, num_cols),
    )


def csr_from_scipy(A: _SupportsToCSR) -> CSRMatrix:
    """HAZmath matrix from a SciPy sparse matrix or array (any format).

    Duplicate entries are summed and the column indices sorted, in a copy, if
    the matrix is not in canonical format.

    Raises:
        ImportError: If SciPy is not installed.
        TypeError: If ``A`` is not a SciPy sparse matrix or array.
    """
    try:
        import scipy.sparse
    except ImportError as err:
        raise ImportError("csr_from_scipy needs SciPy: pip install scipy") from err

    if not scipy.sparse.issparse(A):
        raise TypeError(f"Expected a SciPy sparse matrix or array, got {type(A).__name__}")
    B = A.tocsr()
    if not B.has_canonical_format:
        # tocsr can return A itself, which must not change
        B = B.copy()
        B.sum_duplicates()
    return csr_from_arrays(B.indptr, B.indices, B.data, B.shape)


def csr_from_petsc(A: PETSc.Mat) -> CSRMatrix:
    """HAZmath matrix from an assembled PETSc matrix on one process.

    Matrices that are not AIJ (for example BAIJ or SBAIJ) are converted to
    AIJ first.

    Raises:
        ImportError: If petsc4py is not installed.
        TypeError: If PETSc is not built with real double precision scalars.
        ValueError: If the matrix is distributed over more than one process.
    """
    try:
        from petsc4py import PETSc
    except ImportError as err:
        raise ImportError("csr_from_petsc needs petsc4py: pip install haznicsx[petsc]") from err

    if np.dtype(PETSc.ScalarType) != np.float64:
        raise TypeError(
            "HAZmath needs a PETSc build with real double precision scalars; "
            f"this petsc4py uses {np.dtype(PETSc.ScalarType)}"
        )
    num_processes = A.getComm().getSize()
    if num_processes > 1:
        raise ValueError(f"HAZmath needs a matrix on one process; this one is on {num_processes}")
    if A.getType() != PETSc.Mat.Type.SEQAIJ:
        A = A.convert(PETSc.Mat.Type.SEQAIJ)
    indptr, indices, data = A.getValuesCSR()
    return csr_from_arrays(indptr, indices, data, A.getSize())
