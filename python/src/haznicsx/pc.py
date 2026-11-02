# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of HAZmath (https://github.com/HAZmathTeam/hazmath)
#
# SPDX-License-Identifier:    MIT
"""HAZmath preconditioners for petsc4py (``PCPYTHON`` contexts).

A context is attached to a PETSc ``PC`` either directly::

    pc.setType("python")
    pc.setPythonContext(haznicsx.HazmathAMG())

or through the options database, for example in the ``petsc_options`` of
``dolfinx.fem.petsc.LinearProblem``::

    -pc_type python -pc_python_type haznicsx.HazmathAMG -pc_hazmath_cycle_type w

Every :class:`haznicsx.AMGParam` field can be set with
``-pc_hazmath_<field>`` (lower case, after the PC's options prefix). The
fields that select a method also accept names, see :data:`OPTION_CHOICES`.

HAZmath is serial: the preconditioning matrix must be on one process. In
parallel, use ``HazmathAMG`` as the sub-preconditioner of block Jacobi or
additive Schwarz (``-pc_type bjacobi -sub_pc_type python
-sub_pc_python_type haznicsx.HazmathAMG``). The HX contexts need their
auxiliary matrices, so they are created in code, e.g.
``HazmathHXCurl(*haznicsx.hx_curl_operators(V))``.
"""

from __future__ import annotations

import typing

from haznicsx import _haznicsx as _hz
from haznicsx.csr import csr_from_dolfinx, csr_from_petsc, csr_from_scipy
from haznicsx.precond import AMGParam, CSRMatrix, Precond, amg, famg, hx_curl, hx_div

if typing.TYPE_CHECKING:
    from petsc4py import PETSc

__all__ = [
    "OPTION_CHOICES",
    "HazmathAMG",
    "HazmathFAMG",
    "HazmathHXCurl",
    "HazmathHXDiv",
    "amg_param_from_options",
]

#: Names accepted in the options database, besides the integer values, for
#: the AMGParam fields that select a method
OPTION_CHOICES: dict[str, dict[str, int]] = {
    "amg_type": {"ua": _hz.UA_AMG, "sa": _hz.SA_AMG},
    "cycle_type": {
        "v": _hz.V_CYCLE,
        "w": _hz.W_CYCLE,
        "amli": _hz.AMLI_CYCLE,
        "nl_amli": _hz.NL_AMLI_CYCLE,
        "add": _hz.ADD_CYCLE,
    },
    "smoother": {
        "jacobi": _hz.SMOOTHER_JACOBI,
        "gs": _hz.SMOOTHER_GS,
        "sgs": _hz.SMOOTHER_SGS,
        "cg": _hz.SMOOTHER_CG,
        "sor": _hz.SMOOTHER_SOR,
        "ssor": _hz.SMOOTHER_SSOR,
        "gsor": _hz.SMOOTHER_GSOR,
        "sgsor": _hz.SMOOTHER_SGSOR,
        "l1diag": _hz.SMOOTHER_L1DIAG,
        "fjacobi": _hz.SMOOTHER_FJACOBI,
        "fgs": _hz.SMOOTHER_FGS,
        "fsgs": _hz.SMOOTHER_FSGS,
    },
    "aggregation_type": {"vmb": _hz.VMB, "hec": _hz.HEC, "hem": _hz.HEM},
    "coarse_solver": {"default": _hz.SOLVER_DEFAULT, "umfpack": _hz.SOLVER_UMFPACK},
    "print_level": {
        "none": _hz.PRINT_NONE,
        "min": _hz.PRINT_MIN,
        "some": _hz.PRINT_SOME,
        "more": _hz.PRINT_MORE,
        "most": _hz.PRINT_MOST,
        "all": _hz.PRINT_ALL,
    },
}

_OPTION_PREFIX = "pc_hazmath_"


def _copy(param: AMGParam) -> AMGParam:
    return AMGParam.from_dict({name: getattr(param, name) for name in AMGParam.fields})


def _option_value(options: PETSc.Options, key: str, field: str, current: typing.Any) -> typing.Any:
    choices = OPTION_CHOICES.get(field.lower())
    if choices is None:
        return options.getReal(key) if isinstance(current, float) else options.getInt(key)
    value = options.getString(key).strip().lower()
    if value in choices:
        return choices[value]
    try:
        return int(value)
    except ValueError:
        raise ValueError(
            f"-{options.prefix or ''}{key} {value}: use one of {', '.join(choices)} or an integer"
        ) from None


def amg_param_from_options(
    param: AMGParam, options: PETSc.Options, prefix: str = _OPTION_PREFIX
) -> AMGParam:
    """Update AMG parameters from the options database.

    Args:
        param: AMG parameters, not modified.
        options: Options database, with the PC's options prefix.
        prefix: Prefix of the options, followed by the lower-case field name.

    Returns:
        A copy of ``param`` with the ``<prefix><field>`` options applied.

    Raises:
        ValueError: If an option has a value the field does not accept.
    """
    result = _copy(param)
    for field in AMGParam.fields:
        key = prefix + field.lower()
        if not options.hasName(key):
            continue
        value = _option_value(options, key, field, getattr(param, field))
        try:
            setattr(result, field, value)
        except TypeError as err:
            raise ValueError(f"-{options.prefix or ''}{key} {value}: out of range") from err
    return result


def _to_hazmath(A: PETSc.Mat, hint: str = "") -> CSRMatrix:
    """HAZmath matrix of A, with ``hint`` added to the error for a distributed A."""
    try:
        return csr_from_petsc(A)
    except ValueError as err:
        if hint and A.getComm().getSize() > 1:
            raise ValueError(f"{err}. {hint}") from err
        raise


def _auxiliary(M: typing.Any) -> CSRMatrix:
    """HAZmath matrix of an auxiliary matrix: HAZmath, PETSc, DOLFINx or SciPy."""
    if isinstance(M, CSRMatrix):
        return M
    if hasattr(M, "getValuesCSR"):
        return csr_from_petsc(M)
    if hasattr(M, "to_scipy") and hasattr(M, "block_size"):
        return csr_from_dolfinx(M)
    return csr_from_scipy(M)


class HazmathAMG:
    """HAZmath AMG preconditioner (:func:`haznicsx.amg`) as a ``PCPYTHON`` context.

    The preconditioner is built in ``setUp`` from the PC's preconditioning
    matrix, and rebuilt whenever PETSc sets the PC up again.
    """

    def __init__(self, param: AMGParam | None = None) -> None:
        """Create the context.

        Args:
            param: AMG parameters. ``None`` means the HAZmath defaults.
                Options from the options database are applied on top.
        """
        self.param = AMGParam() if param is None else _copy(param)
        self._precond: Precond | None = None

    def setFromOptions(self, pc: PETSc.PC) -> None:
        """Apply the ``-pc_hazmath_<field>`` options (PETSc callback)."""
        from petsc4py import PETSc

        self.param = amg_param_from_options(self.param, PETSc.Options(pc.getOptionsPrefix()))

    def setUp(self, pc: PETSc.PC) -> None:
        """Build the preconditioner from the preconditioning matrix (PETSc callback)."""
        _, P = pc.getOperators()
        # free the previous preconditioner before building the next
        self._precond = None
        self._precond = self._build(P)

    def _build(self, P: PETSc.Mat) -> Precond:
        hint = (
            "Use HazmathAMG as a sub-preconditioner, e.g. -pc_type bjacobi -sub_pc_type python "
            "-sub_pc_python_type haznicsx.HazmathAMG"
        )
        return amg(_to_hazmath(P, hint), self.param)

    def apply(self, pc: PETSc.PC, x: PETSc.Vec, y: PETSc.Vec) -> None:
        """Compute ``y = B x`` (PETSc callback)."""
        if self._precond is None:
            raise RuntimeError(f"{type(self).__name__} is not set up")
        # no copies: HAZmath reads x's array and writes y's in place
        with x.getBuffer(readonly=True) as r, y.getBuffer() as z:
            self._precond.apply(r, z)

    def reset(self, pc: PETSc.PC) -> None:
        """Free the preconditioner (PETSc callback)."""
        self._precond = None

    def view(self, pc: PETSc.PC, viewer: PETSc.Viewer) -> None:
        """Print the parameters to an ASCII viewer (PETSc callback)."""
        from petsc4py import PETSc

        if viewer.getType() != PETSc.Viewer.Type.ASCII:
            return
        viewer.printfASCII(f"  {type(self).__name__}:\n")
        for field in AMGParam.fields:
            viewer.printfASCII(f"    {field} = {getattr(self.param, field)}\n")


class HazmathFAMG(HazmathAMG):
    """HAZmath fractional AMG (:func:`haznicsx.famg`) as a ``PCPYTHON`` context.

    Preconditions the power ``param.fpwr`` of the operator, using a mass
    matrix. HAZmath is serial and the mass matrix cannot be given through
    the options database, so there is no parallel (sub-preconditioner) use.
    """

    def __init__(self, M: PETSc.Mat | None = None, param: AMGParam | None = None) -> None:
        """Create the context.

        Args:
            M: Mass matrix, on one process and of the same size as the
                preconditioning matrix. It is converted at the first setup.
            param: AMG parameters. ``None`` means the HAZmath defaults with
                the fractional Gauss-Seidel smoother; see :func:`haznicsx.famg`.
        """
        if param is None:
            param = AMGParam.from_dict({"smoother": _hz.SMOOTHER_FGS})
        super().__init__(param)
        self.M = M
        self._M: CSRMatrix | None = None
        self._M_source: PETSc.Mat | None = None

    def _build(self, P: PETSc.Mat) -> Precond:
        if self.M is None:
            raise ValueError("HazmathFAMG needs the mass matrix: HazmathFAMG(M)")
        # convert M once, not at every setup
        if self._M is None or self._M_source is not self.M:
            self._M = _to_hazmath(self.M)
            self._M_source = self.M
        return famg(_to_hazmath(P), self._M, self.param)


class HazmathHXCurl(HazmathAMG):
    """HAZmath HX preconditioner for H(curl) (:func:`haznicsx.hx_curl`).

    A ``PCPYTHON`` context. The auxiliary matrices are converted once, when
    the context is created; :func:`haznicsx.hx_curl_operators` builds them
    for a DOLFINx Nedelec space.
    """

    def __init__(
        self,
        Pcurl: typing.Any,
        Grad: typing.Any,
        param: AMGParam | None = None,
        multiplicative: bool = False,
    ) -> None:
        """Create the context.

        Args:
            Pcurl: Interpolation from vector P1 into the H(curl) space; a
                :class:`haznicsx.CSRMatrix`, PETSc, DOLFINx or SciPy matrix.
            Grad: Discrete gradient from P1 into the H(curl) space.
            param: AMG parameters for the auxiliary Laplacians. ``None``
                means the HAZmath defaults.
            multiplicative: Combine the corrections multiplicatively; the
                result is not symmetric, so use it with e.g. GMRES, not CG.
        """
        super().__init__(param)
        self._Pcurl = _auxiliary(Pcurl)
        self._Grad = _auxiliary(Grad)
        self.multiplicative = multiplicative

    def _build(self, P: PETSc.Mat) -> Precond:
        return hx_curl(_to_hazmath(P), self._Pcurl, self._Grad, self.param, self.multiplicative)


class HazmathHXDiv(HazmathAMG):
    """HAZmath HX preconditioner for H(div) (:func:`haznicsx.hx_div`).

    A ``PCPYTHON`` context. The auxiliary matrices are converted once, when
    the context is created; :func:`haznicsx.hx_div_operators` builds them
    for a DOLFINx Raviart-Thomas space.
    """

    def __init__(
        self,
        Pdiv: typing.Any,
        Curl: typing.Any,
        Pcurl: typing.Any | None = None,
        param: AMGParam | None = None,
        multiplicative: bool = False,
    ) -> None:
        """Create the context.

        Args:
            Pdiv: Interpolation from vector P1 into the H(div) space; a
                :class:`haznicsx.CSRMatrix`, PETSc, DOLFINx or SciPy matrix.
            Curl: In 3D, the discrete curl from Nedelec into the H(div)
                space; in 2D, the discrete rot from P1.
            Pcurl: In 3D, the interpolation from vector P1 into Nedelec;
                ``None`` in 2D.
            param: AMG parameters for the auxiliary Laplacians. ``None``
                means the HAZmath defaults.
            multiplicative: Combine the corrections multiplicatively; the
                result is not symmetric, so use it with e.g. GMRES, not CG.
        """
        super().__init__(param)
        self._Pdiv = _auxiliary(Pdiv)
        self._Curl = _auxiliary(Curl)
        self._Pcurl = None if Pcurl is None else _auxiliary(Pcurl)
        self.multiplicative = multiplicative

    def _build(self, P: PETSc.Mat) -> Precond:
        return hx_div(
            _to_hazmath(P), self._Pdiv, self._Curl, self._Pcurl, self.param, self.multiplicative
        )
