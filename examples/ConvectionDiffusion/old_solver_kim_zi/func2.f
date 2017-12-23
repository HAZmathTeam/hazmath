C=====================================================================
      real*8 function eps(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function EPS(x,y) is the coefficient of the leading 
C...  term of the second order linear PDEs of the form
C...
C...        - div(eps(x,y) grad(u)) + lower_order_terms  = f(x,y).
C...
C---------------------------------------------------------------------
      go to (10,20,30,40,50,60,70,80,90), ich(21)
      stop
C
 10   eps = 1.d-1
      return
C
 20   eps = 1.d-2
      return
C
 30   eps = 1.d-3
      return
C
 40   eps = 1.d-4
      return
C
 50   eps = 1.d-5
      return
C
 60   eps = 1.d-6
      return
C
 70   eps = 1.d0
      return
C
 80   continue
      circle = (x-0.5d0)**2 + (y-0.5d0)**2
      if (circle .le. 0.1111d0) then
         eps = 1.d-3
cc         eps = 1.d0
      else 
cc         eps = 1.d-6
         eps = 1.d-6
      end if
      return
C
 90   continue
      circle = x
      if (circle .le. 0.5d0) then
         eps = 1.0d0
      else 
         eps = 1.0d-1
      end if
      return
C
      end
C=====================================================================
      real*8 function alpha(i,j,x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The 2 X 2 matrix function ALPHA(i,j,x,y) is the coefficients of
C...  the leading term of the second order linear PDEs of the form
C...
C...     - div (alpha(x,y) grad(u)) + lower_order_terms  = f(x,y).
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
C...  PDEs.
C---------------------------------------------------------------------
      go to (10,20,30,40,50,60,70,80,90,100), ich(27)
      stop
C
 10   rhs = 1.d0
      return
C
 20   rhs = 0.d0
      return
C
 30   rhs = 0.d0
      xmx2  = x - x*x
      pi3   = 3*pi
      pi29  = 9*pi*pi
      pi3y  = pi3*y
      s3piy = dsin(pi3y)
      c3piy = dcos(pi3y)
      rhs   = eps(x,y)*s3piy*(2.d0 + pi29*xmx2)
      rhs   = rhs - beta(1,x,y)*(1.d0 - 2*x)*s3piy
      rhs   = rhs - beta(2,x,y)*pi3*xmx2*c3piy
      rhs   = 10*rhs
      return
C
 40   rhs = 0.d0
      xmx2  = x - x*x
      pi3   = 3*pi
      pi29  = 9*pi*pi
      pi3y  = pi3*y
      s3piy = dsin(pi3y)
      c3piy = dcos(pi3y)
      rhs   = eps(x,y)*s3piy*(2.d0 + pi29*xmx2)
      rhs   = rhs + beta(1,x,y)*(1.d0 - 2*x)*s3piy
      rhs   = rhs + beta(2,x,y)*pi3*xmx2*c3piy
      rhs   = 10*rhs
      return
C
 50   continue
      rhs = 0.0d00
      rhs = 2*eps(x,y)*(x - x*x + y - y*y) 
      rhs = rhs - beta(1,x,y)*(1.d0 - 2*x)*(y - y*y)
      rhs = rhs - beta(2,x,y)*(1.d0 - 2*y)*(x - x*x)
      rhs = 100*rhs
      return
C
 60   continue
      rhs = x - dcos(2*pi*y)
      rhs = 10*rhs
      return
C
 70   rhs = 0.d0
      p2 = pi*pi
      sx = dsin(4*pi*x)
      ey = dexp(y)
      eys = dexp(y*y)
      ey1 = ey - eys
      ey3 = ey - (2 + 4*y*y)*eys
      rhs = (16*p2*ey1-ey3)*sx 
      rhs = 10*rhs
      return
C
 80   rhs = (x - y*y)*dexp(y)*50d0
      return
C
 90   rhs = 70*dlog((x + 0.1d0)*(dsin(pi*y) + 1))
      return
C
 100  rhs = 50*(x - y*y)*dexp(y)*dsin(5*pi*x)*dcos(6*pi*x*y)
      return
C
      end
C=====================================================================
      real*8 function beta(i,x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The vector function BETA(i,x,y) is the coefficients of the
C...  convection term of the form
C...
C...       - div(beta(x,y)u)     or     beta(x,y) \cdot grad(u).
C...
C...  If the convection term is given by the gradient of potential
C...  function, i.e., BETA(x,y) = (grad PSI)(x,y), then go to the 
C...  function PSI to define it.
C---------------------------------------------------------------------
      go to (10,20,30,40,50,60,70,80,90,100), ich(24)
      stop
C
 10   go to (11,12), i
      stop
 11   beta =  3*x*x*y - 2*y*y
      return
 12   beta =  5*dsin(pi*x) - 3*dexp(y)
      return
C
 20   go to (21,22), i
      stop
 21   beta = 4*(x*x-y)-y*y
      return
 22   beta = 10*(2*x-y)*dsin(8*(x-y))
      return
C
 30   go to (31,32), i
      stop
 31   beta = 0.0d0
      return
 32   beta = 0.0d0
      return
C
 40   go to (41,42), i
      stop
 41   beta = 4.0d0
      return
 42   beta = -3.0d0
      return
C
 50   continue
      rho = sqrt(x*x + y*y)
      if (rho .gt. -zero .and. rho .lt. 0.8) then
         beta = 0.0d00
         return
      end if
      go to (51,52), i
 51   if (rho .lt. 0.9) then
         beta = 4*x/rho
      else
         beta = 0.0d00
      end if
      return
 52   if(rho .lt. 0.9) then
         beta = 4*y/rho
      else
         beta = 0.0d00
      end if
      return
C
 60   go to (61,62), i
      stop
 61   beta = -1.d0
      return
 62   beta = 0.d0
      return
C
 70   continue
      rho = dsqrt((x-0.5d0)**2+(y-0.5d0)**2)
      if(rho .lt. 1.d-10) rho = rho + 1.d-10
      go to (71,72), i
      stop
 71   beta = -(y-0.5d0)/rho
      return
 72   beta =  (x-0.5d0)/rho
      return
C
 80   go to (81,82), i
      stop
 81   beta = -dexp(y) + 2*x
C      beta =  10*beta
      return
 82   beta =  3*y*dsin(pi*x)
C      beta =  -15*beta
      return
C
 90   go to (91,92), i
      stop
 91   beta = - 2*x*y + 2.d0
cc      beta =  1.d3*beta
      return
 92   beta =  5*dcos(3*pi*x)*dexp(y)
cc      beta =  1.d3*beta
      return
C
 100  go to (101,102), i
      stop
 101  circle = (x-0.5d0)**2 + (y-0.5d0)**2
      if (circle .le. 0.1111d0) then
         beta = (x*x - 2.d0*y)*y*y
cc         beta = - 10.d0*(x*x - y)
      else 
         beta = 2000.d0*dcos(2*pi*x*y)
cc         beta = - 400.d0*y
      end if
      return
 102  circle = (x-0.5d0)**2 + (y-0.5d0)**2
      if (circle .le. 0.1111d0) then
         beta = 1000.d0*(dexp(x) + 2.d0* y)
cc         beta = 50.d0*(dcos(5*pi*x*y))
      else 
         beta = 2.d0*x*dsin(8.d0*(x-y))
cc         beta = - 7.d0*x
      end if
      return
C
      end
C=====================================================================
      real*8 function psi(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function PSI(x,y) is a potential function, i.e.,
C...  the convection term BETA(x,y) of the form
C...
C...       - div(beta(x,y)u)     or     beta(x,y) \cdot grad(u),
C...
C...  is given by the gradient of a potential function, i.e., 
C...  BETA = grad PSI. See the description of the function BETA above.
C---------------------------------------------------------------------
      go to (10,20,30,40,50), ich(25)
      stop
C
 10   psi = -5*x+3*y
      return
C
 20   psi = 1.d0
      return
C
 30   psi = 0.d0
      return
C
 40   psi = 0.d0
      return
C
 50   continue
      xx = x
      yy = y
      if (1-x-y .lt. zero) then
cc         xx = 1.d00-x
cc         yy = 1.d00-y
      end if
      rho = xx*xx + yy*yy
      rho = dsqrt(rho) + x
      if(rho .gt. -zero .and. rho .lt. 0.65) then
         psi = 0.0d00
      else
         if(rho .lt. 0.55) then
             psi = 2.0d00*(rho-0.55)
         else
            psi = 0.2
         end if
      end if
C
      return
      end
C=====================================================================
      real*8 function gamma(x,y)
C=====================================================================
      implicit real*8(a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function GAMMA(x,y) is the coefficient of the zero 
C...  order term of the second order linear PDEs of the form
C...
C...        higher_order_terms + gamma(x,y)u = f(x,y).
C...
C---------------------------------------------------------------------
      go to (10,20,30,40,50), ich(26)
      stop
C
 10   gamma = x*y*y
      return
C
 20   gamma = 0.d0
      return
C
 30   gamma = 0.d0
      return
C
 40   gamma = 0.d0
      return
C
 50   gamma = 2.d0
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
C...  or
C...     {alpha(x,y) grad(u) + beta(x,y)u} \cdot normal + sigma(x,y)u 
C...                                                     = zeta(x,y),
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
C...  or Neumann boundary condition. 
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
 20   g = 0.d0
      return
C
 30   g = (5*x*x-x)*(2*y+y*y)
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
      real*8 function exact_sol(x,y)
C=====================================================================
      implicit real*8 (a-h,o-z)
      common /menu/ ich(200)
      data pi/3.141592653589793d00/, zero/1.d-13/
C---------------------------------------------------------------------
C...  The scalar function EXACT_SOL(x,y) is the exact solution of the
C...  PDE.
C---------------------------------------------------------------------
      go to (10,20,30,40,50,60,70), ich(20)
      stop
C
 10   exact_sol = 100*(x-x*x)*(y-y*y) 
      return
C
 20   exact_sol = 100*(x-x*x)*(y-y*y)*dexp(x+y)
      return
C
 30   exact_sol = 0.d0
      return
C
 40   exact_sol = 0.d0
      return
C
 50   exact_sol = 0.d0
      return
C
 60   exact_sol = 0.d0
      return
C
 70   exact_sol = 10*dsin(4*pi*x)*(dexp(y) - dexp(y*y))

      return
      end
C=====================================================================


