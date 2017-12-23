C=====================================================================
      real*8 function alpha(i,j,x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The 2 X 2 symmetric matrix function ALPHA(i,j,x,y) is the 
C...  coefficients of the leading term of the second order linear PDEs 
C...  of the form
C...
C...     - div (alpha(x,y) grad(u)) + gamma0(x,y) u = f(x,y).
C...
C---------------------------------------------------------------------
      go to (10,20,30), ich(22)
      stop
C
 10   go to (11,14), i
      stop
 11   go to (12,13), j
      stop
 12   alpha = 1.d00
      return
 13   alpha = 0.0d00
      return
 14   go to (15,16), j
      stop
 15   alpha = 0.d00
      return
 16   alpha = 1.0d00
      return
C
 20   go to (21,24), i
      stop
 21   go to (22,23), j
      stop
 22   alpha = 1.d00
      return
 23   alpha = 0.0d00
      return
 24   go to (25,26), j
      stop
 25   alpha = 0.d00
      return
 26   alpha = 1.0d00
      return
C
 30   go to (31,34), i
      stop
 31   go to (32,33), j
      stop
 32   alpha = 1.d00
      return
 33   alpha = 0.0d00
      return
 34   go to (35,36), j
      stop
 35   alpha = 0.d00
      return
 36   alpha = 1.0d00
      return
      end
C=====================================================================
      real*8 function rhs(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function RHS(x,y) is the right hand side of the given
C...  PDE:
C...     - div (alpha(x,y) grad(u)) + gamma0(x,y) u = rhs(x,y).
C---------------------------------------------------------------------
      go to (10,20,30,40,50), ich(27)
      stop
C
 10   rhs = 0.d0
      return
C
 20   rhs = x*x - 5*dsin(4*pi*y)
      return
C
 30   rhs = 1.d0
      return
C
 40   continue
      rhs = dexp(x-y)
      return
C
 50   continue
      rhs = 0.0d00
C
      return
      end
C=====================================================================
      real*8 function gamma0(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function GAMMA0(x,y) is the coefficient of the zero 
C...  order term of the second order linear PDEs of the form
C...
C...     - div (alpha(x,y) grad(u)) + gamma0(x,y) u = f(x,y).
C...
C---------------------------------------------------------------------
      go to (10,20,30,40,50), ich(26)
      stop
C
 10   gamma0 = 0.d0
      return
C
 20   gamma0 = 1.d0
      return
C
 30   gamma0 = 0.d0
      return
C
 40   gamma0 = 0.d0
      return
C
 50   gamma0 = 2.d0
      return
      end
C=====================================================================
      real*8 function sigma(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function SIGMA(x,y) is the coefficient of a Robin or
C...  Neumann boundary condition of the form
C...
C...     {alpha(x,y) grad(u)} \cdot normal + sigma(x,y)u = zeta(x,y),
C...
C...  where ALPHA(x,y) is a matrix or scalar function, NORMAL is the
C...  outward normal vector, and SIGMA(x,y) and ZETA(x,y) are scalar
C...  functions. If SIGMA(x,y) = 0 for all x and y, then it becomes 
C...  Neumann boundary.
C---------------------------------------------------------------------
      go to (10,20,30,40,50), ich(28)
      stop
C
 10   sigma = 0.d0
      return
C
 20   continue
      sigma = sin(pi/18.d00)
      return
C
 30   sigma = 0.d0
      return
C
 40   sigma = 0.d0
      return
C
 50   sigma = 0.d0
      return
      end
C=====================================================================
      real*8 function zeta(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function ZETA(x,y) is the right hand side for a Robin
C...  or Neumann boundary condition: 
C...
C...     {alpha(x,y) grad(u)} \cdot normal + sigma(x,y)u = zeta(x,y),
C...
C---------------------------------------------------------------------
      go to (10,20,30,40,50), ich(29)
      stop
C
 10   zeta = 0.d0
      return
C
 20   zeta = 0.d0
      return
C
 30   zeta = 0.d0
      return
C
 40   zeta = 0.d0
      return
C
 50   zeta = 0.d0
      return
      end
C=====================================================================
      real*8 function g(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function G(x,y) is the Dirichlet boundary data. 
C---------------------------------------------------------------------
      go to (10,20,30,40,50,60), ich(30)
      stop
C
 10   g = 1.d0
      return
C
 20   g = (5*x*x-x)*(2*y+y*y)
      return
C
 30   g = 0.d0
      return
C
 40   continue
      if(dabs(x) .lt. zero .and. y .gt. 0.5d00-zero) then
         g = 0.d00
      else
         g = 1.0d00
      end if
      return
C
 50   continue
      g0 = 0.d00
      xx = x
      yy = y
      if(1-x-y .lt. 0.d00)then
         xx = 1.d00-x
         yy = 1.d00-y
         g0 = 2.1
      end if
      rho = xx*xx + yy*yy
      if((dabs(xx) .lt. zero .and. yy .le. 0.25+zero) .or. 
     >   (dabs(yy) .lt. zero .and. xx .le. 0.25+zero)) g = g0
      return
C
 60   continue
      g0 = 0.d00
      xx = x
      yy = y
      if(1-x-y .lt. 0.d00)then
         xx = 1.d00-x
         yy = 1.d00-y
         g0 = 2.1
      end if
      rho = xx*xx + yy*yy
      if((dabs(xx) .lt. zero .and. yy .le. 0.25+zero) .or. 
     >         (yy .lt. zero .and. xx .le. 0.25+zero)) g = 2-x-y
      return
      end
C=====================================================================


