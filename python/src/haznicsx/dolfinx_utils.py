# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""Auxiliary space matrices for the HX preconditioners, from DOLFINx.

The matrices come from DOLFINx's discrete operators on the mesh of the given
space, and are returned as :class:`haznicsx.CSRMatrix` (on one process).
"""

from __future__ import annotations

import typing

from haznicsx.csr import csr_from_arrays, csr_from_dolfinx

if typing.TYPE_CHECKING:
    import dolfinx.fem

    from haznicsx.precond import CSRMatrix

__all__ = ["hx_curl_operators", "hx_div_operators"]


def _check_space(V: dolfinx.fem.FunctionSpace, family: str, name: str) -> None:
    import basix

    element = V.element.basix_element
    if element.family != getattr(basix.ElementFamily, family) or element.degree != 1:
        raise ValueError(f"{name} needs a lowest order {family} space")


def hx_curl_operators(V: dolfinx.fem.FunctionSpace) -> tuple[CSRMatrix, CSRMatrix]:
    """Auxiliary matrices of :func:`haznicsx.hx_curl` for a Nedelec space.

    Args:
        V: Lowest order Nedelec (first kind) space, in 2D or 3D.

    Returns:
        ``Pcurl``, the interpolation from vector P1 into ``V``, and ``Grad``,
        the discrete gradient from P1 into ``V``.

    Raises:
        ValueError: If ``V`` is not a lowest order Nedelec space.
    """
    from dolfinx import fem

    _check_space(V, "N1E", "hx_curl_operators")
    msh = V.mesh
    P1 = fem.functionspace(msh, ("Lagrange", 1))
    P1v = fem.functionspace(msh, ("Lagrange", 1, (msh.geometry.dim,)))
    Pcurl = csr_from_dolfinx(fem.interpolation_matrix(P1v, V))
    Grad = csr_from_dolfinx(fem.discrete_gradient(P1, V))
    return Pcurl, Grad


def hx_div_operators(
    V: dolfinx.fem.FunctionSpace,
) -> tuple[CSRMatrix, CSRMatrix, CSRMatrix | None]:
    """Auxiliary matrices of :func:`haznicsx.hx_div` for a Raviart-Thomas space.

    Args:
        V: Lowest order Raviart-Thomas space, in 2D or 3D.

    Returns:
        ``Pdiv``, the interpolation from vector P1 into ``V``; ``Curl``, the
        discrete curl from lowest order Nedelec into ``V`` in 3D, or the
        discrete rot ``(d/dy, -d/dx)`` from P1 into ``V`` in 2D; and in 3D
        ``Pcurl``, the interpolation from vector P1 into Nedelec (``None``
        in 2D).

    Raises:
        ValueError: If ``V`` is not a lowest order Raviart-Thomas space.
    """
    from dolfinx import fem

    _check_space(V, "RT", "hx_div_operators")
    msh = V.mesh
    P1v = fem.functionspace(msh, ("Lagrange", 1, (msh.geometry.dim,)))
    Pdiv = csr_from_dolfinx(fem.interpolation_matrix(P1v, V))
    Ned = fem.functionspace(msh, ("N1curl", 1))
    if msh.topology.dim == 3:
        Curl = csr_from_dolfinx(fem.discrete_curl(Ned, V))
        return Pdiv, Curl, csr_from_dolfinx(fem.interpolation_matrix(P1v, Ned))
    # In 2D, Nedelec and Raviart-Thomas share their degrees of freedom (one per
    # edge), and the matrix of rot into Raviart-Thomas is minus that of grad
    # into Nedelec
    G = csr_from_dolfinx(fem.discrete_gradient(fem.functionspace(msh, ("Lagrange", 1)), Ned))
    return Pdiv, csr_from_arrays(G.indptr, G.indices, -G.data, G.shape), None
