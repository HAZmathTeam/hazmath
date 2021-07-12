from block.block_base import block_base
from builtins import str
import os
import ctypes as ct
import numpy as np
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import norm as sparse_norm
from dolfin import *

# ------------------ link swig files ---------------- #
# run "source setup.rc" in ../
# --------------------------------------------------- #
import haznics


class Precond(block_base):
    """
    Class of general preconditioners from HAZmath using SWIG

    """

    def __init__(self, A, prectype, parameters={}, precond=None):
        # haznics.dCSRmat* type (assert?)
        self.A = A
        # python dictionary of parameters
        self.parameters = parameters

        # init and set preconditioner (precond *)
        if precond:
            self.precond = precond
        else:
            import warnings
            warnings.warn("!! Preconditioner not specified !! Creating default UA-AMG precond... ", RuntimeWarning)
            # change data type for the matrix (to dCSRmat pointer)
            petsc_mat = as_backend_type(A).mat()

            # NB! store copies for now
            csr0 = petsc_mat.getValuesCSR()[0]
            csr1 = petsc_mat.getValuesCSR()[1]
            csr2 = petsc_mat.getValuesCSR()[2]

            A_ptr = haznics.create_matrix(csr2, csr1, csr0)

            # initialize amg parameters (AMG_param pointer)
            amgparam = haznics.amg_param_alloc(1)

            # print (relevant) amg parameters
            haznics.param_amg_print(amgparam)

            self.precond = haznics.create_precond(A_ptr, amgparam)

            # if fail, setup returns null
            if not precond:
                raise RuntimeError(
                    "AMG levels failed to set up (null pointer returned) ")

        # preconditioner type (string)
        self.prectype = prectype

    def matvec(self, b):
        from dolfin import GenericVector
        if not isinstance(b, GenericVector):
            return NotImplemented

        x = self.A.create_vec(dim=1)
        if len(x) != len(b):
            raise RuntimeError('incompatible dimensions for matvec, %d != %d'
                               % (len(x), len(b)))

        # convert rhs and dx to numpy arrays
        b_np = b[:]
        x_np = x[:]

        # apply the preconditioner (solution dx saved in x_np)
        haznics.apply_precond(b_np, x_np, self.precond)

        # convert dx to GenericVector
        x.set_local(x_np)

        return x

    # noinspection PyMethodMayBeStatic
    def down_cast(self):
        return NotImplemented

    def __str__(self):
        return '<%s prec of %s>' % (self.__class__.__name__, str(self.A))


class AMG(Precond):
    """
    AMG preconditioner from the HAZmath Library with SWIG

    """

    def __init__(self, A, parameters={}):
        # change data type for the matrix (to dCSRmat pointer)
        petsc_mat = as_backend_type(A).mat()

        # NB! store copies for now
        csr0 = petsc_mat.getValuesCSR()[0]
        csr1 = petsc_mat.getValuesCSR()[1]
        csr2 = petsc_mat.getValuesCSR()[2]

        A_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # initialize amg parameters (AMG_param pointer)
        amgparam = haznics.amg_param_alloc(1)

        # set extra amg parameters
        if parameters:
            for key in parameters:
                if isinstance(parameters[key], str):
                    exec("amgparam.%s = \"%s\"" % (key, parameters[key]))
                # elif isinstance(parameters[key], function):
                #     haznics.py_callback_setup(parameters[key], amgparam)
                else:
                    exec("amgparam.%s = %s" % (key, parameters[key]))

        # print (relevant) amg parameters
        haznics.param_amg_print(amgparam)

        # set AMG preconditioner
        precond = haznics.create_precond_amg(A_ptr, amgparam)

        # if fail, setup returns null
        if not precond:
            raise RuntimeError(
                "AMG levels failed to set up (null pointer returned) ")

        Precond.__init__(self, A, "amg", parameters, precond)


class FAMG(Precond):
    """
    AMG preconditioner from the HAZmath Library

    """

    def __init__(self, A, M, parameters={'fpwr': 0.5, 'smoother': 'fjacobi'}):
        # change data type for the A matrix (to dCSRmat pointer)
        petsc_mat_A = as_backend_type(A).mat()

        # NB! store copies for now
        csr0 = petsc_mat_A.getValuesCSR()[0]
        csr1 = petsc_mat_A.getValuesCSR()[1]
        csr2 = petsc_mat_A.getValuesCSR()[2]

        A_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # change data type for the M matrix (to dCSRmat pointer)
        petsc_mat_M = as_backend_type(M).mat()

        # NB! store copies for now
        csr0 = petsc_mat_M.getValuesCSR()[0]
        csr1 = petsc_mat_M.getValuesCSR()[1]
        csr2 = petsc_mat_M.getValuesCSR()[2]

        M_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # initialize amg parameters (AMG_param pointer)
        amgparam = haznics.amg_param_alloc(1)

        # set extra amg parameters
        if parameters:
            for key in parameters:
                if isinstance(parameters[key], str):
                    exec("amgparam.%s = \"%s\"" % (key, parameters[key]))
                else:
                    exec("amgparam.%s = %s" % (key, parameters[key]))

        # print (relevant) amg parameters
        haznics.param_amg_print(amgparam)

        # set AMG preconditioner
        precond = haznics.create_precond_famg(A_ptr, M_ptr, amgparam)

        # if fail, setup returns null
        if not precond:
            raise RuntimeError(
                "FAMG levels failed to set up (null pointer returned) ")

        Precond.__init__(self, A, "famg", parameters, precond)


class SumFAMG(Precond):
    """
    Fractional AMG preconditioner from the HAZmath Library
    for a sum of fractional operators alpha * A^s + beta * A^(1+s)
    """

    def __init__(self, A, M, parameters={'coefs': [1.0, 1.0], 'pwrs': [-0.5, 0.5],
                                         'smoother': 'fjacobi'}):
        # initialize amg parameters
        lib.amg_param_alloc.restype = ct.POINTER(hazmath.AMG_param)
        amgparam_p = lib.amg_param_alloc(ct.c_short(1))

        # set extra amg parameters
        amgparam = hazmath.set_amg_param(parameters)
        lib.param_amg_cp(ct.byref(amgparam), amgparam_p)

        # print (relevant) amg parameters
        lib.param_amg_print(amgparam_p)

        # change data type for the matrix (to dCSRmat pointer)
        petsc_mat_A = as_backend_type(A).mat()

        # NB! store copies for now
        csr0 = petsc_mat_A.getValuesCSR()[0]
        csr1 = petsc_mat_A.getValuesCSR()[1]
        csr2 = petsc_mat_A.getValuesCSR()[2]

        A_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # change data type for the M matrix (to dCSRmat pointer)
        petsc_mat_M = as_backend_type(M).mat()

        # NB! store copies for now
        csr0 = petsc_mat_M.getValuesCSR()[0]
        csr1 = petsc_mat_M.getValuesCSR()[1]
        csr2 = petsc_mat_M.getValuesCSR()[2]

        M_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # initialize amg parameters (AMG_param pointer)
        amgparam = haznics.amg_param_alloc(1)

        # set extra amg parameters
        if parameters:
            for key in parameters:
                if isinstance(parameters[key], str):
                    exec("amgparam.%s = \"%s\"" % (key, parameters[key]))
                else:
                    exec("amgparam.%s = %s" % (key, parameters[key]))

        # print (relevant) amg parameters
        haznics.param_amg_print(amgparam)

        # other parameters
        alpha, beta = parameters['coefs']
        s, t = parameters['pwrs']

        # set AMG preconditioner #  fixme: create this function
        precond = haznics.create_precond_sum_famg(A_ptr, M_ptr, alpha, beta, s, t, amgparam)

        # if fail, setup returns null
        if not precond:
            raise RuntimeError(
                "SUM FAMG levels failed to set up (null pointer returned) ")

        # set up preconditioner function pointer  # fixme: set this up in precond create function
        try:
            cycle_type = parameters['cycle_type']
        except KeyError:
            cycle_type = ''

        if cycle_type in ['add', 'ADD']:
            precond.fct = haznics.precond_sum_famg_add2
        else:
            # default should be V or W cycle, but I don't think it's
            # implemented in hazmath yet
            precond.fct = haznics.precond_sum_famg_add2

        Precond.__init__(self, A, "sum_famg", parameters, precond)


class RA(Precond):
    """
    Rational approximation preconditioner from the HAZmath library

    """

    def __init__(self, A, M, dim=2,
                 parameters={'coefs': [1.0, 0.0], 'pwrs': [0.5, 0.0]}):

        # change data type for the matrix (to dCSRmat pointer)
        petsc_mat_A = as_backend_type(A).mat()

        # NB! store copies for now
        csr0 = petsc_mat_A.getValuesCSR()[0]
        csr1 = petsc_mat_A.getValuesCSR()[1]
        csr2 = petsc_mat_A.getValuesCSR()[2]

        A_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # change data type for the M matrix (to dCSRmat pointer)
        petsc_mat_M = as_backend_type(M).mat()

        # NB! store copies for now
        csr0 = petsc_mat_M.getValuesCSR()[0]
        csr1 = petsc_mat_M.getValuesCSR()[1]
        csr2 = petsc_mat_M.getValuesCSR()[2]

        M_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # initialize amg parameters (AMG_param pointer)
        amgparam = haznics.amg_param_alloc(1)

        # set extra amg parameters
        if parameters:
            for key in parameters:
                if isinstance(parameters[key], str):
                    exec("amgparam.%s = \"%s\"" % (key, parameters[key]))
                else:
                    exec("amgparam.%s = %s" % (key, parameters[key]))

        # print (relevant) amg parameters
        haznics.param_amg_print(amgparam)

        # get scalings
        scaling_a = 1. / A.norm("linf")  # / sparse_norm(Acsr, np.inf)
        scaling_m = 1. / petsc_mat_M.getDiagonal().min()[1]  # / np.min(Mcsr.diagonal())  # 32.0

        # get coefs and powers
        alpha, beta = parameters['coefs']
        s_power, t_power = parameters['pwrs']
        """
        print("----------------------------------------------------------")
        print("Before create precond RA: ")
        print("A col, row, nnz: ", A_ptr.col, A_ptr.row, A_ptr.nnz)
        print("M col, row, nnz: ", M_ptr.col, M_ptr.row, M_ptr.nnz)
        print("alpha, beta: ", alpha, beta)
        print("s, t: ", s_power, t_power)
        print("scalings a, m: ", scaling_a, scaling_m)
        print("----------------------------------------------------------")
        """
        # set AMG preconditioner #
        precond = haznics.create_precond_ra(A_ptr, M_ptr, s_power, t_power, alpha, beta, scaling_a, scaling_m, amgparam)

        # if fail, setup returns null
        if not precond:
            raise RuntimeError(
                "Rational Approximation data failed to set up (null pointer returned) ")
        # import pdb; pdb.set_trace()
        # self.npoles = haznics.get_poles_no(precond)

        Precond.__init__(self, A, "RA", parameters, precond)


class HXCurl(Precond):
    """
    HX preconditioner from the HAZmath library for the curl-curl inner product
    NB! only for 3D problems
    
    """

    def __init__(self, Acurl, Pcurl, Grad, parameters={}):
        # change data type for the Acurl matrix (to dCSRmat pointer)
        petsc_mat_A = as_backend_type(Acurl).mat()

        # NB! store copies for now
        csr0 = petsc_mat_A.getValuesCSR()[0]
        csr1 = petsc_mat_A.getValuesCSR()[1]
        csr2 = petsc_mat_A.getValuesCSR()[2]

        Acurl_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # change data type for the Pcurl matrix (to dCSRmat pointer)
        petsc_mat_P = as_backend_type(Pcurl).mat()

        # NB! store copies for now
        csr0 = petsc_mat_P.getValuesCSR()[0]
        csr1 = petsc_mat_P.getValuesCSR()[1]
        csr2 = petsc_mat_P.getValuesCSR()[2]

        Pcurl_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # change data type for the Grad matrix (to dCSRmat pointer)
        petsc_mat_G = as_backend_type(Grad).mat()

        # NB! store copies for now
        csr0 = petsc_mat_G.getValuesCSR()[0]
        csr1 = petsc_mat_G.getValuesCSR()[1]
        csr2 = petsc_mat_G.getValuesCSR()[2]

        Grad_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # initialize amg parameters (AMG_param pointer)
        amgparam = haznics.amg_param_alloc(1)

        # set extra amg parameters
        if parameters:
            for key in parameters:
                if isinstance(parameters[key], str):
                    exec("amgparam.%s = \"%s\"" % (key, parameters[key]))
                else:
                    exec("amgparam.%s = %s" % (key, parameters[key]))

        # print (relevant) amg parameters
        haznics.param_amg_print(amgparam)

        # add or multi
        try:
            prectype = parameters['prectype']
        except KeyError:
            prectype = haznics.PREC_HX_CURL_A

        # set AMG preconditioner
        precond = haznics.create_precond_hxcurl(Acurl_ptr, Pcurl_ptr, Grad_ptr, prectype, amgparam)

        # if fail, setup returns null
        if not precond:
            raise RuntimeError(
                "HXcurl data failed to set up (null pointer returned) ")

        try:
            prectype = parameters['prectype']
        except KeyError:
            prectype = ''

        if prectype in ["add", "Add", "ADD", "additive", "ADDITIVE"]:
            precond.fct = haznics.precond_hx_curl_additive
        elif prectype in ["multi", "MULTI", "Multi", "multiplicative", "MULTIPLICATIVE"]:
            precond.fct = haznics.precond_hx_curl_multiplicative
        else:  # default is additive
            precond.fct = haznics.precond_hx_curl_additive

        Precond.__init__(self, Acurl, "HXCurl_add", parameters, precond)


class HXDiv(Precond):
    """
    HX preconditioner from the HAZmath library for the div-div inner product

    """

    def __init__(self, Adiv, Pdiv, Curl, Pcurl=None, parameters={'dimension': 2}):
        # change data type for the Adiv matrix (to dCSRmat pointer)
        petsc_mat_A = as_backend_type(Adiv).mat()

        # NB! store copies for now
        csr0 = petsc_mat_A.getValuesCSR()[0]
        csr1 = petsc_mat_A.getValuesCSR()[1]
        csr2 = petsc_mat_A.getValuesCSR()[2]

        Adiv_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # change data type for the Pdiv matrix (to dCSRmat pointer)
        petsc_mat_P = as_backend_type(Pdiv).mat()

        # NB! store copies for now
        csr0 = petsc_mat_P.getValuesCSR()[0]
        csr1 = petsc_mat_P.getValuesCSR()[1]
        csr2 = petsc_mat_P.getValuesCSR()[2]

        Pdiv_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # change data type for the Curl matrix (to dCSRmat pointer)
        petsc_mat_C = as_backend_type(Curl).mat()

        # NB! store copies for now
        csr0 = petsc_mat_C.getValuesCSR()[0]
        csr1 = petsc_mat_C.getValuesCSR()[1]
        csr2 = petsc_mat_C.getValuesCSR()[2]

        Curl_ptr = haznics.create_matrix(csr2, csr1, csr0)

        # initialize amg parameters (AMG_param pointer)
        amgparam = haznics.amg_param_alloc(1)

        # set extra amg parameters
        if parameters:
            for key in parameters:
                if isinstance(parameters[key], str):
                    exec("amgparam.%s = \"%s\"" % (key, parameters[key]))
                else:
                    exec("amgparam.%s = %s" % (key, parameters[key]))

        # print (relevant) amg parameters
        haznics.param_amg_print(amgparam)

        # get dimension and type of HX precond application
        try:
            dim = parameters['dimension']
        except KeyError:
            dim = 2

        # add or multi
        try:
            prectype = parameters['prectype']
        except KeyError:
            prectype = haznics.PREC_HX_DIV_A

        if dim == 3:
            # check Pcurl
            assert Pcurl, "For 3D case, Pcurl operator is needed!"

            # change data type for the Pcurl matrix (to dCSRmat pointer)
            petsc_mat_P = as_backend_type(Pcurl).mat()

            # NB! store copies for now
            csr0 = petsc_mat_P.getValuesCSR()[0]
            csr1 = petsc_mat_P.getValuesCSR()[1]
            csr2 = petsc_mat_P.getValuesCSR()[2]

            Pcurl_ptr = haznics.create_matrix(csr2, csr1, csr0)

            # set AMG preconditioner
            precond = haznics.create_precond_hxdiv_3D(Adiv_ptr, Pdiv_ptr, Curl_ptr, Pcurl_ptr, prectype, amgparam)

            # if fail, setup returns null
            if not precond:
                raise RuntimeError(
                    "HXdiv data failed to set up (null pointer returned) ")
            """
            if prectype in ["add", "Add", "ADD", "additive", "ADDITIVE"]:
                precond.fct = haznics.precond_hx_div_additive

            elif prectype in ["multi", "MULTI", "Multi", "multiplicative", "MULTIPLICATIVE"]:
                precond.fct = haznics.precond_hx_div_multiplicative

            else:
                # default is additive
                precond.fct = haznics.precond_hx_div_additive
            """
            Precond.__init__(self, Adiv, "HXDiv_add", parameters, precond)

        else:
            # set AMG preconditioner
            precond = haznics.create_precond_hxdiv_2D(Adiv_ptr, Pdiv_ptr, Curl_ptr, prectype, amgparam)

            # if fail, setup returns null
            if not precond:
                raise RuntimeError(
                    "HXdiv data failed to set up (null pointer returned) ")
            """
            if prectype in ["add", "Add", "ADD", "additive", "ADDITIVE"]:
                precond.fct = haznics.precond_hx_div_additive_2D

            elif prectype in ["multi", "MULTI", "Multi", "multiplicative", "MULTIPLICATIVE"]:
                precond.fct = haznics.precond_hx_div_multiplicative_2D

            else:
                # default is additive
                precond.fct = haznics.precond_hx_div_additive_2D
            """
            Precond.__init__(self, Adiv, "HXDiv_add", parameters, precond)


# ----------------------------------- EOF ----------------------------------- #
