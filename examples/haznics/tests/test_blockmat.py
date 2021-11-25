import haznics
from dolfin import *


mesh = UnitIntervalMesh(10) 
V = FunctionSpace(mesh, "Lagrange", 1)  
u, v = TrialFunction(V), TestFunction(V)
A = assemble(u*v*dx + inner(grad(u), grad(v))*dx) 

b = assemble(Constant(1.0)*v*dx) 
x = assemble(Constant(1.0)*v*dx) 

# store copies for now
petsc_mat = as_backend_type(A).mat()
csr0  = petsc_mat.getValuesCSR()[0]
csr1  = petsc_mat.getValuesCSR()[1]
csr2 =  petsc_mat.getValuesCSR()[2]
print (csr0) 
print (csr1) 
print (csr2) 
AA = haznics.create_matrix(csr2, csr1, csr0, A.size(1))


AAA = haznics.block_dCSRmat() 
AAA.init(2,2)
AAA.debugPrint()


type(AAA)
# import pdb; pdb.set_trace()

AAA.set(i=0,j=0,mat=AA)
AAA.set(i=0,j=1,mat=AA)
AAA.set(i=1,j=0,mat=AA)
AAA.set(i=1,j=1,mat=AA)

xxx = x[:]
bbb = b[:]

bb = haznics.create_dvector(bbb)
xx = haznics.create_dvector(xxx)

# do  block mat vec 




