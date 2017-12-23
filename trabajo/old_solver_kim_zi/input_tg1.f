C=====================================================================
      subroutine input_data
C=====================================================================
      implicit real*8(a-h,o-z)
      character*100 meshf(0:20)
      common /fname/ meshf
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  This subroutine initializes various parameters.
C...
C...  Parameter:
C...  
C...  Mesh file names:
C...    meshf(0) - original mesh.
C...    meshf(k) - level k-th mesh for MG. 
C...
C...  General preference:
C...    ich(1)  - choice of the PDE.
C...    ich(2)  - choice of the boundary condition. If Robin or Neumann
C...              boundary condition is identically zero or there is no
C...              such condition on the boundary at all, then set 
C...              ich(2) = 0, otherwise, set ich(2) = 1.
C...    ich(3)  - choice of finite elements, i.e., degrees of freedom
C...              3 = linear element, 6 = quadratic element,
C...              4 = bilinear element,
C...    ich(4)  - choice of the finite element methods.
C...              1 = Standard Galerkin,    2 = EAFE, 
C...              3 = Streamline diffusion, 4 = Upwind, 5 = SUPG,
C...              6 = Mixed with 3 and 4,   7 = 
C...    ich(5)  - mesh format: 1 = formatted, 2 = unformatted
C...    ich(6)  - 
C...    ich(7)  - choice of the solver of the linear system.
C...              1 = Gaussian elimination, 2 = Gauss-Seidel, 3 = MG,
C...              4 = Linear iterative method with some preconditioners.
C...    ich(8)  - the number of subdomains, 1, 2, or 4.
C...    ich(9)  - the coarse level of the two-grid algorithm, 3,4,5,6,7.
C...
C...  Quadrature for integrations:
C...    ich(11) - number of the quadrature points in the triangle: 
C...              1, 3, or 7.
C...    ich(12) - number of the Gauss quadrature points on each Neumann 
C...              boundary edge: 1, 2, 3, 4, 5, 6, or 7.
C...    ich(13) - number of the Gauss quadrature points on the general 
C...              line segment: 1, 2, 3, 4, 5, 6, or 7.
C...
C...  Functions:
C...    ich(20) - Exact solution if it is known.
C...    ich(21) - choice of the coefficient of the leading term, 
C...              EPS(x,y).
C...    ich(22) - choice of the coefficients of the leading term, 
C...              ALPHA(i,j,x,y) (This has a matrix form.)
C...    ich(23) - format of the convection term, 
C...              1 = usual vector format, i.e., 
C...                 BETA(x,y) = <BETA(1,x,y), BETA(2,x,y)>,
C...              2 = given by the gradient of potential function, i.e.,
C...                 BETA(x,y) = (grad PSI)(x,y)
C...    ich(24) - choice of the convection term, BETA(x,y), when 
C...              ich(23) = 1.
C...    ich(25) - choice of the convection term, BETA(x,y), when
C...              ich(23) = 2, i.e., the gradient of potential function.
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
C...  Multigrid: Set the following values if ich(7) = 3.
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
C...    ich(38) - choice of the coarsest grid solver: 
C...              1 = fwd G-S 
C...              2 = bwd G-S 
C...              3 = sym G-S 
C...              4 = Gauss elim 
C...    ich(39) - maximum number of MG iterations
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
C...  General preference:
C
      ich(1) = 6      !PDE
      ich(2) = 0      !Boundary condition for Robin or Neumann boundary 
      ich(3) = 3      !Elements (degrees of freedom)
      ich(4) = 1      !Finite element methods 
      ich(5) = 2      !Mesh format: 1 = formatted, 2 = unformatted
      ich(6) = 1      !Transposed EAFE = -1
      ich(7) = 3      !Solver of the linear system 
      ich(8) = 1      !the number of subdomains, 1, 2, or 4.
      ich(9) = 4      !the coarse level of the two-grid algorithm, 3,4,5,6,7.
C
C...  Quadrature for integrations:
C
      ich(11) = 3     !Quadrature in each triangle
      ich(12) = 2     !Quadrature in each Neumann edge
      ich(13) = 2     !Quadrature in general line segment
C
C...  Functions:
C
      ich(20) = 7     !Exact solution if it is known
      ich(21) = 7     !Leading term: a scalar function, EPS(x,y)
      ich(22) = 1     !Leading term: a matrix function, ALPHA(i,j,x,y) 
      ich(23) = 1     !Format of the convection term
      ich(24) = 1     !Convection term: BETA, when ich(23) = 1
      ich(25) = 1     !Convection term: BETA = grad PSI, when ich(23) = 2
      ich(26) = 1     !Zero order term: GAMMA(x,y)
      ich(27) = 6     !Right hand side: RHS(x,y)
      ich(28) = 2     !Coefficient of Robin or Neumann boundary, SIGMA(x,y)
      ich(29) = 1     !Right side of Robin or Neumann baundary , ZETA(x,y)
      ich(30) = 2     !Dirichlet boundary data, G(x,y)
C     
C...  Multigrid:  If ich(7) = 3, then set the following values.
C
      ich(31) = 1     !choice of MG:
      ich(32) = 1     !choice of pre-smoothing: 
      ich(33) = 1     !choice of post-smoothing: 
      ich(34) = 1     !number of pre-smoothings: 1,2,3, etc
      ich(35) = 1     !number of post-smoothings: 1,2,3, etc
      ich(36) = 8     !level of the finest mesh: 2,3,4, etc
      ich(37) = 1     !level of the coarsest mesh: 1,2,3, etc
      ich(38) = 1     !choice of the coarsest grid solver: 
      ich(39) = 20    !maximum number of MG iterations
C
      iread = 1
      if(iread .eq. 1) then
         read(*,*) ich(9)
         read(*,*) ich(32)
         read(*,*) ich(33)
         read(*,*) ich(34)
         read(*,*) ich(35)
         read(*,*) ich(36)
         read(*,*) ich(37)
         read(*,*) ich(38)
         read(*,*) ich(39)
      else if(iread.eq.2) then
         open(123,file='inpt1',status='unknown')
         read(123,*) ich(9)
         read(123,*) ich(32)
         read(123,*) ich(33)
         read(123,*) ich(34)
         read(123,*) ich(35)
         read(123,*) ich(36)
         read(123,*) ich(37)
         read(123,*) ich(38)
         read(123,*) ich(39)
         close (123)
      end if
C
      return
      end
C=====================================================================

