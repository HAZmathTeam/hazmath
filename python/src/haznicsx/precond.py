# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""AMG, fractional AMG and auxiliary space (HX) preconditioners."""

from __future__ import annotations

from haznicsx import _haznicsx as _hz
from haznicsx._haznicsx import AMGParam, CSRMatrix, Precond

__all__ = ["AMGParam", "CSRMatrix", "Precond", "amg", "famg", "hx_curl", "hx_div"]


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
