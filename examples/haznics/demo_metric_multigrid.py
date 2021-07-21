from fenics import * 
import haznics
import block 
from block.algebraic.hazmath import AMG as hazAMG
from block.algebraic.petsc import collapse 
from block.iterative import ConjGrad 


def run(eps, N):
    print ("running N=%d eps=%e" % (N, eps))
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, "CG", 1) 
    u = TrialFunction(V) 
    v = TestFunction(V)
    M = assemble(u*v*dx) 
    A = assemble(inner(grad(u), grad(v))*dx) 

    AA = block.block_mat([[M+eps*A, M], [M, M+eps*A]])
    BB = block.block_mat([[hazAMG(collapse(M+eps*A)), 0], [0, hazAMG(collapse(M+eps*A))]])

    xx = AA.create_vec()
    xx.randomize()
    
    bb = AA.create_vec() 
    bb.zero() 

    AAinv = ConjGrad(AA, precond=BB, initial_guess=xx, tolerance=1.0e-12, show=2) 

    xx = AAinv*bb
    


if __name__ == "__main__": 
    run(16, 1)
    run(32, 1)
    run(64, 1)

    run(16, 0.001)
    run(32, 0.001)
    run(64, 0.001)





    
