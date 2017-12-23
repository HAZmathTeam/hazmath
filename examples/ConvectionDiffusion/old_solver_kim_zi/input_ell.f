C=====================================================================
      subroutine input
C=====================================================================
      implicit real*8(a-h,o-z)
      character*100 meshf(0:20),mtrxf(20)
      common /fname/ meshf,mtrxf
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  This subroutine initializes various parameters.
C...
C...  Parameter:
C...  
C...  Mesh file names:
C...    meshf(k) - level k-th mesh for MG. 
C...    mtrxf(k) - level k-th matrix A and load vector B. 
C...
C...  General preference:
C...    ich(1)  - The value is 1 if the stiffness matrix A and the load 
C...              vector B have to be computed from the given mesh; it
C...              is 2 if the stiffness matrix A and the load vector B 
C...              are given as inputs in formatted form; it is 3 if the
C...              stiffness matrix A and the load vector B are given
C...              as inputs in unformatted form. If ich(1) > 1, then
C...              skip assigning the values of ich(2), ich(3), ich(5), 
C...              ich(11), ich(12), and ich(13). 
C...    ich(2)  - choice of the boundary condition. If Robin or Neumann
C...              boundary condition is identically zero or there is no
C...              such condition on the boundary at all, then set 
C...              ich(2) = 0, otherwise, set ich(2) = 1.
C...    ich(3)  - choice of finite elements, i.e., degrees of freedom:
C...              3 = linear element (default value)
C...    ich(5)  - mesh format: 1 = formatted, 2 = unformatted
C...
C...  Quadrature for integrations:
C...    ich(11) - number of the quadrature points in the triangle: 
C...              1, 3, or 7. The default value is 3.
C...    ich(12) - number of the Gauss quadrature points on each Neumann 
C...              boundary edge: 1, 2, 3, 4, 5, 6, or 7.
C...              The default value is 2.
C...    ich(13) - number of the Gauss quadrature points on the general 
C...              line segment: 1, 2, 3, 4, 5, 6, or 7.
C...              The default value is 2.
C...
C...  Functions:
C...    ich(22) - choice of the coefficients of the leading term, 
C...              ALPHA(i,j,x,y) (It has a matrix form.)
C...    ich(26) - choice of the coefficient of the zero order term, 
C...              GAMMA(x,y).
C...    ich(27) - choice of the right hand side, RHS(x,y).
C...    ich(28) - choice of the coefficient for a Robin or Neumann
C...              boundary data, SIGMA(x,y). If SIGMA(x,y) = 0 for all 
C...              x and y, then it becomes a Neumann boundary condition.
C...    ich(29) - choice of the right hand side for a Robin or Neumann 
C...              boundary condition, ZETA(x,y).
C...    ich(30) - choice for the Dirichlet boundary data, G(x,y).
C...
C...  Multigrid:
C...    ich(31) - choice of MG: 1 = V-cycle, 2 = \-cycle
C...              3 = Variable V-cycle
C...    ich(32) - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    ich(33) - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    ich(34) - number of pre-smoothings: 1,2,3, etc
C...    ich(35) - number of post-smoothings: 1,2,3, etc
C...    ich(36) - level of the finest mesh: 2,3,4, etc
C...    ich(37) - level of the coarsest mesh: 1,2,3, etc
C...    ich(38) - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...    ich(39) - maximum number of MG iterations
C...
C---------------------------------------------------------------------
      meshf(0) = 'matrices/regu0'
      meshf(1) = 'matrices/regu1'
      meshf(2) = 'matrices/regu2'
      meshf(3) = 'matrices/regu3'
      meshf(4) = 'matrices/regu4'
      meshf(5) = 'matrices/regu5'
      meshf(6) = 'matrices/regu6'
      meshf(7) = 'matrices/regu7'
      meshf(8) = 'matrices/regu8'
      meshf(9) = 'matrices/regu9'
      meshf(10) = 'matrices/regu10'
      meshf(11) = 'matrices/regu11'
C
      mtrxf(1) = 'matrices/mtrx1'
      mtrxf(2) = 'matrices/mtrx2'
      mtrxf(3) = 'matrices/mtrx3'
      mtrxf(4) = 'matrices/mtrx4'
      mtrxf(5) = 'matrices/mtrx5'
      mtrxf(6) = 'matrices/mtrx6'
      mtrxf(7) = 'matrices/mtrx7'
      mtrxf(8) = 'matrices/mtrx8'
      mtrxf(9) = 'matrices/mtrx9'
      mtrxf(10) = 'matrices/mtrx10'
C
C...  General preference:
C
C...  If iread = 1, then read the input data from a separate data file;
C...  if it is 2, then assign the appropriate values.
C
      iread = 2
      go to (100, 200) iread
C
 100  continue
      read(*,*) ich(1)
      read(*,*) ich(31)
      read(*,*) ich(32)
      read(*,*) ich(33)
      read(*,*) ich(34)
      read(*,*) ich(35)
      read(*,*) ich(36)
      read(*,*) ich(37)
      read(*,*) ich(38)
      read(*,*) ich(39)
      go to 300
C
 200  continue
      ich(1)  = 1     !To compute A and B or to read them from inputs:
      ich(31) = 1     !choice of MG:
      ich(32) = 1     !choice of pre-smoothing: 
      ich(33) = 1     !choice of post-smoothing: 
      ich(34) = 1     !number of pre-smoothings: 1,2,3, etc
      ich(35) = 1     !number of post-smoothings: 1,2,3, etc
      ich(36) = 8     !level of the finest mesh: 2,3,4, etc
      ich(37) = 1     !level of the coarsest mesh: 1,2,3, etc
      ich(38) = 3     !choice of the coarsest grid solver: 
      ich(39) = 50    !maximum number of MG iterations
C
 300  continue
      if (ich(1) .gt. 1) go to 400
C
      ich(2)  = 0     !Boundary condition for Robin or Neumann boundary 
      ich(3)  = 3     !Elements (degrees of freedom)
      ich(5)  = 2     !Mesh format: 1 = formatted, 2 = unformatted
C
      ich(11) = 3     !Quadrature in each triangle
      ich(12) = 2     !Quadrature in each Neumann edge
      ich(13) = 2     !Quadrature in general line segment
C
C...  Functions:
C
      ich(22) = 1     !Leading term: a matrix function, ALPHA(i,j,x,y) 
      ich(26) = 1     !Zero order term: GAMMA(x,y)
      ich(27) = 1     !Right hand side: RHS(x,y)
      ich(28) = 1     !Coefficient of Robin or Neumann boundary, SIGMA(x,y)
      ich(29) = 1     !Right side of Robin or Neumann baundary , ZETA(x,y)
      ich(30) = 2     !Dirichlet boundary data, G(x,y)
C
 400  return
      end
C=====================================================================
