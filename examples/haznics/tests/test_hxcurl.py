from dolfin import *
from block.algebraic.hazmath import Pcurl, discrete_gradient, PETSc_to_dCSRmat
from block.iterative import ConjGrad
from block.algebraic.petsc import HypreAMS
import haznics
import numpy as np

N = 5

# set up 3d curl-curl problem
mesh = UnitCubeMesh(N, N, N)
V = FunctionSpace(mesh, "N1curl", 1)
W = VectorFunctionSpace(mesh, "CG", 1)

p, q = TrialFunction(W), TestFunction(W)

aa = inner(grad(p), grad(q))

v, u = TestFunction(V), TrialFunction(V)

a = inner(curl(v), curl(u)) * dx
m = inner(u, v) * dx
# A = assemble(a)
AM = assemble(a + m)

gdim = mesh.geometry().dim()
b = assemble(inner(v, Constant((1,) * gdim)) * dx)
x = np.zeros(AM.size(1))

# solve with petsc
# B = HypreAMS(AM, V)
# Ainv = ConjGrad(AM, precond=B, tolerance=1E-6, show=2)
# xx = Ainv * b

# parameters for hazmath
inparam = haznics.input_param("input.dat")
itparam = haznics.linear_itsolver_param()
haznics.param_linear_solver_set(itparam, inparam)
itparam.linear_precond_type = haznics.PREC_HX_CURL_A
haznics.param_linear_solver_print(itparam)

amgparam = haznics.AMG_param()
haznics.param_amg_set(amgparam, inparam)
haznics.param_amg_print(amgparam)

# get auxiliary operators for HX precond
Pc = Pcurl(mesh)
Grad = discrete_gradient(mesh)

# convert data to hazmath types
Acurl_ptr = PETSc_to_dCSRmat(AM)
Pcurl_ptr = PETSc_to_dCSRmat(Pc)
Grad_ptr = PETSc_to_dCSRmat(Grad)
bb = b.get_local()
b_ptr = haznics.create_dvector(bb)
x_ptr = haznics.create_dvector(x)

print("A: ", Acurl_ptr.row, Acurl_ptr.col, Acurl_ptr.nnz)
print("Pcurl: ", Pcurl_ptr.row, Pcurl_ptr.col, Pcurl_ptr.nnz)
print("Grad: ", Grad_ptr.row, Grad_ptr.col, Grad_ptr.nnz)


# import pdb; pdb.set_trace()
status = haznics.linear_solver_dcsr_krylov_hx_curl(Acurl_ptr,
                                                   b_ptr,
                                                   x_ptr,
                                                   itparam,
                                                   amgparam,
                                                   Pcurl_ptr,
                                                   Grad_ptr)

# print("The end.")
