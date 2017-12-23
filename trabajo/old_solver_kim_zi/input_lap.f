C=====================================================================
      subroutine input_data
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  This subroutine initializes various parameters.
C...
C...  Parameter:
C...  
C...  General preference:
C...    ich(2)  - choice of the boundary condition. If Robin or Neumann
C...              boundary condition is identically zero or there is no
C...              such condition on the boundary at all, then set 
C...              ich(2) = 0, otherwise, set ich(2) = 1.
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
C...    ich(27) - choice of the right hand side, RHS(x,y).
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
C...              2 = bwd G-S, 3 = sym G-S
C...    ich(39) - maximum number of MG iterations
C...
C---------------------------------------------------------------------
C
C...  If iread = 1, then read the input data from a separate data file;
C...  if it is 2, then assign the appropriate values.
C
      iread = 1
      go to (100, 200) iread
C
 100  continue
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
      ich(31) = 1     !choice of MG:
      ich(32) = 1     !choice of pre-smoothing: 
      ich(33) = 1     !choice of post-smoothing: 
      ich(34) = 1     !number of pre-smoothings: 1,2,3, etc
      ich(35) = 1     !number of post-smoothings: 1,2,3, etc
      ich(36) = 8     !level of the finest mesh: 2,3,4, etc
      ich(37) = 1     !level of the coarsest mesh: 1,2,3, etc
      ich(38) = 1     !choice of the coarsest grid solver: 
      ich(39) = 50    !maximum number of MG iterations
C
 300  continue
C
      ich(2)  = 0     !Boundary condition for Robin or Neumann boundary 
C
      ich(11) = 3     !Quadrature in each triangle
      ich(12) = 2     !Quadrature in each Neumann edge
      ich(13) = 2     !Quadrature in general line segment
C
C...  Functions:
C
      ich(27) = 1     !Right hand side: RHS(x,y)
      ich(30) = 2     !Dirichlet boundary data, G(x,y)
C
 400  return
      end
C=====================================================================
