from numpy import random
import haznics
from dolfin import *


mesh = UnitIntervalMesh(10) 
V = FunctionSpace(mesh, "Lagrange", 1)  
u, v = TrialFunction(V), TestFunction(V)
A = assemble(u*v*dx + inner(grad(u), grad(v))*dx) 
b = assemble(Constant(1)*v*dx())
x = assemble(Constant(1)*v*dx())
x[:] = random.random(x.size()) 

petsc_mat = as_backend_type(A).mat()
print ("nnz ", len(petsc_mat.getValuesCSR()[1]))
print ("col ", len(petsc_mat.getValuesCSR()[0]))

# store copies for now
csr0  = petsc_mat.getValuesCSR()[0]
csr1  = petsc_mat.getValuesCSR()[1]
csr2 =  petsc_mat.getValuesCSR()[2]
print (csr0) 
print (csr1) 
print (csr2)

print(" ------------- CREATE HAZMATH MATRIX --------------")
AA = haznics.create_matrix(csr2, csr1, csr0, A.size(1))
print(" ---------------------- END -----------------------")

print(" ------------- CREATE HAZMATH VECTOR --------------")
bbb = b[:]
xxx = x[:]

bb = haznics.create_dvector(bbb)
xx = haznics.create_dvector(xxx)
print(" ---------------------- END -----------------------")

"""
print(" ----------- CREATE HAZMATH AMG PARAM -------------")
amgparam = haznics.AMG_param()
haznics.param_amg_init(amgparam)
haznics.param_amg_print(amgparam)
print(" ---------------------- END -----------------------")

print(" ------------- SOLVE HAZMATH DIRECT ---------------")
# haznics.directsolve_UMF(AA, bb, xx, 12)
print(" ---------------------- END -----------------------")

print(" --------------- SOLVE HAZMATH AMG ----------------")
ret = haznics.linear_solver_amg(AA, bb, xx, amgparam)
print("Returns: \n", ret)
print(" ---------------------- END -----------------------")

print(" ------------ CREATE HAZMATH IT PARAM -------------")
it_param = haznics.linear_itsolver_param()
haznics.param_linear_solver_init(it_param)
haznics.param_linear_solver_print(it_param)
print(" ---------------------- END -----------------------")

print(" ------------ SOLVE HAZMATH KRYLOV AMG ------------")
ret = haznics.linear_solver_dcsr_krylov_amg(AA, bb, xx, it_param, amgparam)
print("Returns: \n", ret)
print(" ---------------------- END -----------------------")
"""
"""
prec = haznics.create_precond_amg(AA, amgparam)
print (type(prec)) 
print (dir(prec)) 
print ("before apply ") 
print (prec.precond_data())
print (dir(prec.precond_data()))
print ("the numbers below seem strange, wrongly initialized parameters?????")
print ("precond_data smoother ", prec.precond_data().smoother)
print ("precond_data maxit ", prec.precond_data().maxit)
print ("precond_data max_levels ", prec.precond_data().max_levels)
print ("precond_data tol ", prec.precond_data().tol)
print ("the numbers change from run to run ") 

# print ("This function produces a crash") 
prec.apply(bb, xx)
"""
