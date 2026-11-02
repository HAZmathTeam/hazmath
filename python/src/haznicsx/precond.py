# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""AMG and fractional AMG preconditioners."""

from __future__ import annotations

from haznicsx import _haznicsx as _hz
from haznicsx._haznicsx import AMGParam, CSRMatrix, Precond

__all__ = ["AMGParam", "CSRMatrix", "Precond", "amg", "famg"]


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
