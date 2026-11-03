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
-sub_pc_python_type haznicsx.HazmathAMG``).

``HazmathAMG`` and ``HazmathFAMG`` (with ``M`` set afterwards) can be created
from the options database (``-pc_python_type``). The HX and RA contexts need
matrices or parameters that the options database cannot provide, so they are
created in code, e.g. ``HazmathHXCurl(*haznicsx.hx_curl_operators(V))`` or
``HazmathRA(M, s=0.5)``; their AMG parameters can still be set with
``-pc_hazmath_<field>``.
"""

from __future__ import annotations

import typing

from haznicsx import _haznicsx as _hz
from haznicsx.csr import csr_from_dolfinx, csr_from_petsc, csr_from_scipy
from haznicsx.precond import (
    AMGParam,
    CSRMatrix,
    Precond,
    _check_ra_available,
    amg,
    default_scaling_a,
    famg,
    hx_curl,
    hx_div,
    ra,
)

if typing.TYPE_CHECKING:
    from petsc4py import PETSc

__all__ = [
    "OPTION_CHOICES",
    "HazmathAMG",
    "HazmathFAMG",
    "HazmathHXCurl",
    "HazmathHXDiv",
    "HazmathRA",
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


class _Converted:
    """HAZmath matrix of a context's auxiliary matrix, kept between setups.

    Converted when the context is created, and again only after the context's
    matrix is replaced.
    """

    def __init__(self, M: typing.Any) -> None:
        self._source = M
        self._matrix = None if M is None else _auxiliary(M)

    def get(self, M: typing.Any) -> CSRMatrix:
        """The HAZmath matrix of ``M``, converted if ``M`` is new."""
        if self._matrix is None or M is not self._source:
            self._matrix = _auxiliary(M)
            self._source = M
        return self._matrix


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

    def __init__(self, M: typing.Any = None, param: AMGParam | None = None) -> None:
        """Create the context.

        Args:
            M: Mass matrix, of the same size as the preconditioning matrix; a
                PETSc matrix on one process, :class:`haznicsx.CSRMatrix`,
                DOLFINx or SciPy matrix. It is converted when the context is
                created, and again at a setup after ``ctx.M`` is replaced.
            param: AMG parameters. ``None`` means the HAZmath defaults with
                the fractional Gauss-Seidel smoother; see :func:`haznicsx.famg`.
        """
        if param is None:
            param = AMGParam.from_dict({"smoother": _hz.SMOOTHER_FGS})
        super().__init__(param)
        self.M = M
        self._M = _Converted(M)

    def _build(self, P: PETSc.Mat) -> Precond:
        if self.M is None:
            raise ValueError("HazmathFAMG needs the mass matrix: HazmathFAMG(M)")
        return famg(_to_hazmath(P), self._M.get(self.M), self.param)


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


class HazmathRA(HazmathAMG):
    """HAZmath rational approximation preconditioner (:func:`haznicsx.ra`).

    A ``PCPYTHON`` context for the fractional operator ``alpha A^s + beta A^t``
    of the PC's preconditioning matrix ``A``, in the inner product of the mass
    matrix ``M``. Its inner solves make it slightly nonlinear, so use it with
    a flexible Krylov method (``-ksp_type fgmres``).
    """

    def __init__(
        self,
        M: typing.Any,
        s: float,
        t: float = 0.0,
        alpha: float = 1.0,
        beta: float = 0.0,
        scaling_a: float | None = None,
        scaling_m: float = 1.0,
        ra_tol: float = 1e-10,
        param: AMGParam | None = None,
    ) -> None:
        """Create the context.

        Args:
            M: Mass matrix; a PETSc, :class:`haznicsx.CSRMatrix`, DOLFINx or
                SciPy matrix. It is converted when the context is created,
                and again at a setup after ``ctx.M`` is replaced.
            s: First fractional power, in ``[-1, 1]``.
            t: Second fractional power, in ``[-1, 1]``.
            alpha: Coefficient of ``A^s``.
            beta: Coefficient of ``A^t``.
            scaling_a: See :func:`haznicsx.ra`. ``None`` computes it at every
                setup with :func:`haznicsx.default_scaling_a` (SLEPc).
            scaling_m: See :func:`haznicsx.ra`.
            ra_tol: Tolerance of the AAA algorithm.
            param: AMG parameters for the shifted problems. ``None`` means the
                HAZmath defaults.

        Raises:
            RuntimeError: If libhazmath has no LAPACK.
            ImportError: If ``scaling_a`` is ``None`` and petsc4py or
                slepc4py is not installed.
        """
        _check_ra_available(scaling_a)
        super().__init__(param)
        self.M = M
        self._M = _Converted(M)
        self.s, self.t, self.alpha, self.beta = s, t, alpha, beta
        self.scaling_a, self.scaling_m, self.ra_tol = scaling_a, scaling_m, ra_tol
        self._scaling_a_used: float | None = None

    def _build(self, P: PETSc.Mat) -> Precond:
        A = _to_hazmath(P)
        M = self._M.get(self.M)
        scaling_a = self.scaling_a
        if scaling_a is None:
            # SLEPc takes a PETSc M as it is; the others go through M
            M_petsc = self.M if hasattr(self.M, "getValuesCSR") else M
            scaling_a = default_scaling_a(P, M_petsc, self.scaling_m)
        self._scaling_a_used = scaling_a
        return ra(
            A,
            M,
            self.s,
            self.t,
            self.alpha,
            self.beta,
            scaling_a,
            self.scaling_m,
            self.ra_tol,
            self.param,
        )

    def view(self, pc: PETSc.PC, viewer: PETSc.Viewer) -> None:
        """Print the RA and AMG parameters to an ASCII viewer (PETSc callback)."""
        from petsc4py import PETSc

        if viewer.getType() != PETSc.Viewer.Type.ASCII:
            return
        if self.scaling_a is not None:
            scaling_a = f"{self.scaling_a}"
        elif self._scaling_a_used is not None:
            scaling_a = f"{self._scaling_a_used} (computed)"
        else:
            scaling_a = "computed at setup"
        super().view(pc, viewer)
        viewer.printfASCII(
            f"    s = {self.s}, t = {self.t}, alpha = {self.alpha}, beta = {self.beta}\n"
            f"    scaling_a = {scaling_a}, scaling_m = {self.scaling_m}, "
            f"ra_tol = {self.ra_tol}\n"
        )
