# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Conversion of sparse matrices to HAZmath's CSR format.

HAZmath uses double precision values and indices of type ``INT``: 32-bit, or
64-bit when HAZmath is configured with ``int64=yes`` (see
:data:`haznicsx.index_dtype`). The functions here convert the arrays only when
needed and check that the indices fit. The :class:`haznicsx.CSRMatrix` shares
arrays that already have the index dtype and float64 and are C-contiguous, so
changing them afterwards changes the matrix.
"""

from __future__ import annotations

import typing

import numpy as np
import numpy.typing as npt

from haznicsx._haznicsx import INT_SIZE, CSRMatrix

if typing.TYPE_CHECKING:
    from petsc4py import PETSc

    import dolfinx.la

__all__ = [
    "csr_from_arrays",
    "csr_from_dolfinx",
    "csr_from_petsc",
    "csr_from_scipy",
    "index_dtype",
]

#: NumPy dtype of HAZmath's index type INT
index_dtype = np.dtype(np.int64) if INT_SIZE == 8 else np.dtype(np.int32)
_INDEX = np.iinfo(index_dtype)


class _SupportsToCSR(typing.Protocol):
    """A SciPy sparse matrix or array."""

    def tocsr(self) -> typing.Any: ...


def _as_index(a: npt.ArrayLike, name: str) -> npt.NDArray[np.integer]:
    a = np.asarray(a)
    # compare as Python ints: NumPy < 2 compares uint64 with int64 in float64
    if a.dtype != index_dtype and a.size > 0:
        if int(a.max()) > _INDEX.max or int(a.min()) < _INDEX.min:
            raise OverflowError(
                f"{name} does not fit in the {_INDEX.bits}-bit integers that HAZmath uses"
            )
    return np.ascontiguousarray(a, dtype=index_dtype)


def _check_one_process(num_processes: int) -> None:
    if num_processes > 1:
        raise ValueError(f"HAZmath needs a matrix on one process; this one is on {num_processes}")


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

    Returns:
        The HAZmath matrix.

    Raises:
        TypeError: If the values are complex.
        OverflowError: If an index or dimension does not fit in HAZmath's ``INT``.
        ValueError: If the arrays are not a valid CSR matrix.
    """
    values = np.asarray(data)
    if np.iscomplexobj(values):
        raise TypeError("HAZmath supports real matrices only")
    num_rows, num_cols = (int(s) for s in shape)
    if max(num_rows, num_cols) > _INDEX.max:
        raise OverflowError(f"The matrix dimensions do not fit in {_INDEX.bits}-bit integers")
    return CSRMatrix(
        _as_index(indptr, "indptr"),
        _as_index(indices, "indices"),
        np.ascontiguousarray(values, dtype=np.float64),
        (num_rows, num_cols),
    )


def csr_from_scipy(A: _SupportsToCSR) -> CSRMatrix:
    """HAZmath matrix from a SciPy sparse matrix or array (any format).

    Duplicate entries are summed and the column indices sorted, in a copy, if
    the matrix is not in canonical format.

    Args:
        A: The matrix.

    Returns:
        The HAZmath matrix.

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

    Args:
        A: The matrix.

    Returns:
        The HAZmath matrix, with copies of the PETSc arrays.

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
    _check_one_process(A.getComm().getSize())
    if A.getType() != PETSc.Mat.Type.SEQAIJ:
        A = A.convert(PETSc.Mat.Type.SEQAIJ)
    indptr, indices, data = A.getValuesCSR()
    return csr_from_arrays(indptr, indices, data, A.getSize())


def csr_from_dolfinx(A: dolfinx.la.MatrixCSR) -> CSRMatrix:
    """HAZmath matrix from a DOLFINx matrix on one process.

    For a scalar space the matrix shares the values of ``A``, and the index
    arrays that have HAZmath's index dtype. Blocked matrices (vector spaces)
    are expanded to scalar entries, which needs SciPy.

    Args:
        A: The matrix.

    Returns:
        The HAZmath matrix.

    Raises:
        ImportError: If ``A`` is blocked and SciPy is not installed.
        TypeError: If the values are not real double precision.
        ValueError: If the matrix is distributed over more than one process.
    """
    if A.data.dtype != np.float64:
        raise TypeError(
            f"HAZmath needs real double precision values; this matrix has {A.data.dtype}"
        )
    # with more than one process the columns are process-local indices
    _check_one_process(A.index_map(0).comm.size)
    if list(A.block_size) == [1, 1]:
        shape = (A.index_map(0).size_local, A.index_map(1).size_local)
        return csr_from_arrays(A.indptr, A.indices, A.data, shape)
    # the arrays of a blocked matrix index blocks; to_scipy() gives a BSR matrix
    try:
        S = A.to_scipy()
    except ImportError as err:
        raise ImportError(
            "csr_from_dolfinx needs SciPy for blocked matrices: pip install scipy"
        ) from err
    return csr_from_scipy(S)
