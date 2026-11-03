# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""AMG, fractional AMG, auxiliary space (HX) and rational approximation preconditioners."""

from __future__ import annotations

import typing

import numpy as np

from haznicsx import _haznicsx as _hz
from haznicsx._haznicsx import AMGParam, CSRMatrix, Precond

if typing.TYPE_CHECKING:
    from petsc4py import PETSc

__all__ = [
    "AMGParam",
    "CSRMatrix",
    "Precond",
    "amg",
    "default_scaling_a",
    "estimate_lambda_max",
    "famg",
    "hx_curl",
    "hx_div",
    "ra",
]

#: Safety factor on the estimated largest eigenvalue when ra() chooses the
#: scaling, so that the scaled spectrum stays inside [0, 1]
RA_SCALING_MARGIN = 1.1


def amg(A: CSRMatrix, param: AMGParam | None = None) -> Precond:
    """AMG preconditioner for a square matrix.

    Args:
        A: The matrix. HAZmath copies it during setup.
        param: AMG parameters. ``None`` means the HAZmath defaults.

    Returns:
        The preconditioner.

    Raises:
        ValueError: If ``A`` is not square or not valid CSR, or ``param``
            has a value that HAZmath cannot handle.
        RuntimeError: If the HAZmath setup fails.
    """
    return _hz.amg(A, AMGParam() if param is None else param)


def famg(A: CSRMatrix, M: CSRMatrix, param: AMGParam | None = None) -> Precond:
    """Fractional AMG preconditioner, for the power ``param.fpwr`` of ``A``.

    Args:
        A: The matrix. HAZmath copies it during setup.
        M: Mass matrix, of the same size as ``A``.
        param: AMG parameters. ``param.smoother`` must be a fractional
            smoother (``SMOOTHER_FJACOBI``, ``SMOOTHER_FGS`` or
            ``SMOOTHER_FSGS``). ``None`` means the HAZmath defaults with
            ``SMOOTHER_FGS``.

    Returns:
        The preconditioner.

    Raises:
        ValueError: If ``A`` or ``M`` is not square or not valid CSR, their
            sizes differ, or ``param`` has a value that HAZmath cannot
            handle.
        RuntimeError: If the HAZmath setup fails.
    """
    if param is None:
        param = AMGParam.from_dict({"smoother": _hz.SMOOTHER_FGS})
    return _hz.famg(A, M, param)


def hx_curl(
    A: CSRMatrix,
    Pcurl: CSRMatrix,
    Grad: CSRMatrix,
    param: AMGParam | None = None,
    multiplicative: bool = False,
) -> Precond:
    """HX (auxiliary space) preconditioner for an H(curl) matrix.

    Combines smoothing in the space of ``A`` with AMG for the vector
    Laplacian ``Pcurl^T A Pcurl`` and the scalar Laplacian
    ``Grad^T A Grad``.

    Args:
        A: The H(curl) matrix, e.g. for ``curl u . curl v + u . v`` with
            Nedelec elements.
        Pcurl: Interpolation from vector P1 into the space of ``A``.
        Grad: Discrete gradient from P1 into the space of ``A``.
        param: AMG parameters for both Laplacians. The cycle must be V or W.
            ``None`` means the HAZmath defaults.
        multiplicative: Combine the corrections multiplicatively instead of
            additively.

    Returns:
        The preconditioner.

    Raises:
        ValueError: If the matrices do not fit together, or ``param`` has a
            value that HAZmath cannot handle.
        RuntimeError: If the HAZmath setup fails.
    """
    return _hz.hx_curl(A, Pcurl, Grad, AMGParam() if param is None else param, multiplicative)


def hx_div(
    A: CSRMatrix,
    Pdiv: CSRMatrix,
    Curl: CSRMatrix,
    Pcurl: CSRMatrix | None = None,
    param: AMGParam | None = None,
    multiplicative: bool = False,
) -> Precond:
    """HX (auxiliary space) preconditioner for an H(div) matrix.

    In 3D (``Pcurl`` given), the curl sub-problem is handled through the
    vector Laplacian on Nedelec's auxiliary space; in 2D (``Pcurl`` is
    ``None``), ``Curl`` maps P1 into the space of ``A`` and the scalar
    Laplacian is used.

    Args:
        A: The H(div) matrix, e.g. for ``div u div v + u . v`` with
            Raviart-Thomas elements.
        Pdiv: Interpolation from vector P1 into the space of ``A``.
        Curl: In 3D, the discrete curl from Nedelec into the space of ``A``;
            in 2D, the discrete rot from P1.
        Pcurl: In 3D, the interpolation from vector P1 into Nedelec.
        param: AMG parameters for the Laplacians. The cycle must be V or W.
            ``None`` means the HAZmath defaults.
        multiplicative: Combine the corrections multiplicatively instead of
            additively.

    Returns:
        The preconditioner.

    Raises:
        ValueError: If the matrices do not fit together, or ``param`` has a
            value that HAZmath cannot handle.
        RuntimeError: If the HAZmath setup fails.
    """
    param = AMGParam() if param is None else param
    return _hz.hx_div(A, Pdiv, Curl, Pcurl, param, multiplicative)


def estimate_lambda_max(A: PETSc.Mat, M: PETSc.Mat, rtol: float = 1e-4) -> float:
    """Largest eigenvalue of ``A x = lambda M x``, computed with SLEPc.

    Uses SLEPc's default eigensolver for the generalized Hermitian problem,
    with ``M`` solved by CG preconditioned with Jacobi.

    Args:
        A: Symmetric matrix.
        M: Symmetric positive definite matrix, of the same size.
        rtol: Relative tolerance of the eigenvalue.

    Returns:
        The largest eigenvalue.

    Raises:
        ImportError: If slepc4py is not installed.
        RuntimeError: If SLEPc does not converge.
    """
    try:
        from slepc4py import SLEPc
    except ImportError as err:
        raise ImportError(
            "estimate_lambda_max needs slepc4py: pip install haznicsx[slepc]"
        ) from err

    eps = SLEPc.EPS().create(A.getComm())
    try:
        eps.setOperators(A, M)
        eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
        eps.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
        eps.setDimensions(1)
        eps.setTolerances(rtol)
        ksp = eps.getST().getKSP()
        ksp.setType("cg")
        ksp.getPC().setType("jacobi")
        ksp.setTolerances(rtol=1e-10)
        eps.solve()
        if eps.getConverged() < 1:
            raise RuntimeError("SLEPc did not find the largest eigenvalue")
        return float(eps.getEigenvalue(0).real)
    finally:
        eps.destroy()


def _petsc(A: CSRMatrix) -> PETSc.Mat:
    from petsc4py import PETSc

    # HAZmath's INT and PetscInt may differ (make config int64=yes)
    if max(*A.shape, int(A.indptr[-1])) > np.iinfo(PETSc.IntType).max:
        raise ValueError("The matrix is too large for PETSc's 32-bit indices")
    return PETSc.Mat().createAIJ(
        size=A.shape,
        csr=(
            A.indptr.astype(PETSc.IntType, copy=False),
            A.indices.astype(PETSc.IntType, copy=False),
            A.data,
        ),
        comm=PETSc.COMM_SELF,
    )


def _check_slepc() -> None:
    """Raise ``ImportError`` unless the default RA scaling can be computed."""
    try:
        import petsc4py.PETSc  # noqa: F401

        import slepc4py.SLEPc  # noqa: F401
    except ImportError as err:
        raise ImportError(
            "Computing the RA scaling needs petsc4py and slepc4py "
            "(pip install haznicsx[slepc]); or pass scaling_a"
        ) from err


def _check_ra_available(scaling_a: float | None) -> None:
    """Check that :func:`ra` can run, before any matrix is at hand.

    Args:
        scaling_a: The ``scaling_a`` that :func:`ra` will get.

    Raises:
        RuntimeError: If libhazmath has no LAPACK.
        ImportError: If ``scaling_a`` is ``None`` and petsc4py or slepc4py
            cannot be imported.
    """
    if not _hz.HAS_LAPACK:
        raise RuntimeError("ra needs libhazmath built with LAPACK (USE_LAPACK=ON)")
    if scaling_a is None:
        _check_slepc()


def default_scaling_a(A: typing.Any, M: typing.Any, scaling_m: float = 1.0) -> float:
    """The scaling of ``A`` that :func:`ra` uses when none is given.

    ``scaling_m / (RA_SCALING_MARGIN * lambda_max)``, with ``lambda_max`` the
    largest eigenvalue of ``M^-1 A`` from :func:`estimate_lambda_max`.

    Args:
        A: The matrix; a PETSc matrix or :class:`haznicsx.CSRMatrix`.
        M: Mass matrix; a PETSc matrix or :class:`haznicsx.CSRMatrix`.
        scaling_m: Scaling of ``M``.

    Returns:
        The scaling.

    Raises:
        ImportError: If petsc4py or slepc4py is not installed.
        ValueError: If the largest eigenvalue is not positive.
    """
    _check_slepc()
    # destroy the PETSc copies made here, also when an error occurs
    created: list[PETSc.Mat] = []
    try:
        if isinstance(A, CSRMatrix):
            A = _petsc(A)
            created.append(A)
        if isinstance(M, CSRMatrix):
            M = _petsc(M)
            created.append(M)
        lambda_max = estimate_lambda_max(A, M)
    finally:
        for X in created:
            X.destroy()
    if not lambda_max > 0:
        raise ValueError(
            f"The largest eigenvalue of M^-1 A is {lambda_max}; the RA preconditioner "
            "needs a positive one"
        )
    return scaling_m / (RA_SCALING_MARGIN * lambda_max)


def ra(
    A: CSRMatrix,
    M: CSRMatrix,
    s: float,
    t: float = 0.0,
    alpha: float = 1.0,
    beta: float = 0.0,
    scaling_a: float | None = None,
    scaling_m: float = 1.0,
    ra_tol: float = 1e-10,
    param: AMGParam | None = None,
) -> Precond:
    """Rational approximation preconditioner for a fractional operator.

    Approximates the inverse of ``alpha A^s + beta A^t`` in the ``M`` inner
    product: ``B = (alpha X^s + beta X^t)^-1 M^-1`` with ``X = M^-1 A``. A
    rational approximation of ``1 / (alpha x^s + beta x^t)`` on ``[0, 1]``
    (AAA algorithm) turns this into one shifted problem ``A - p M`` per pole
    ``p``, each solved by PCG with AMG. The inner solves make ``B`` slightly
    nonlinear, so use a flexible Krylov method (FGMRES).

    Args:
        A: Symmetric positive (semi-)definite matrix, e.g. a stiffness matrix.
        M: Mass matrix, of the same size.
        s: First fractional power, in ``[-1, 1]``.
        t: Second fractional power, in ``[-1, 1]``.
        alpha: Coefficient of ``A^s``; non-negative.
        beta: Coefficient of ``A^t``; non-negative, not both zero.
        scaling_a: Scaling of ``A``, so that ``scaling_a / scaling_m`` times
            the largest eigenvalue of ``M^-1 A`` is at most 1. ``None``
            uses :func:`default_scaling_a`, which computes that eigenvalue
            with SLEPc and leaves a margin of :data:`RA_SCALING_MARGIN`.
        scaling_m: Scaling of ``M``.
        ra_tol: Tolerance of the AAA algorithm.
        param: AMG parameters for the shifted problems. The cycle must be V
            or W. ``None`` means the HAZmath defaults.

    Returns:
        The preconditioner.

    Raises:
        ValueError: If the matrices do not fit together, or a parameter is
            out of range.
        ImportError: If ``scaling_a`` is ``None`` and petsc4py or slepc4py
            is not installed.
        RuntimeError: If libhazmath has no LAPACK (needed for the AAA
            algorithm), or the HAZmath setup fails.
    """
    param = AMGParam() if param is None else param
    # the cheap checks first: the default scaling_a needs an eigenvalue solve
    _hz.check_ra(A, M, s, t, alpha, beta, scaling_m, ra_tol, param)
    if scaling_a is None:
        scaling_a = default_scaling_a(A, M, scaling_m)
    return _hz.ra(A, M, s, t, alpha, beta, scaling_a, scaling_m, ra_tol, param)
