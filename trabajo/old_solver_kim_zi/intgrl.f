C=====================================================================
      subroutine intgrl1d(xe,ye,nd,result,fxy)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),xx(7),yy(7)
      common /quad_gs/ w(7),z(7),lquad_gs
C---------------------------------------------------------------------
C...  This subroutine computes the line integral of the function
C...  fxy(x,y) on the line segment joining two points (XE(1),YE(1)) 
C...  and (XE(2),YE(2)) using Gauss quadrature. If the line integral
C...  is performed on the interval J = [a, b], then we set XE(1) = a, 
C...  XE(2) = b, YE(1) = YE(2) = 0.
C...
C...  Input:
C...    XE, YE - x and y coordinates of the end points of the line 
C...             segment. If the line integral is performed on the
C...             interval J = [a, b], then XE(1) = a, XE(2) = b, 
C...             YE(1) = YE(2) = 0.
C...    ND     - ND = 2.
C...
C...  Output
C...    result - the rusulting value of the line integral.
C...
C...  Parameter:
C...    TOL    - to decide whether function is zero or not.
C...
C...  Working space:
C...    XX, YY - x and y coordinates of the points on the line segment
C...             which correspond to the Gauss points.
C...
C...  Call module: 
C...    QUAD_GS  - This is a common data set which is defined in the 
C...               subroutine QUAD_DATA. 
C...    FXY(x,y) - an external function: the integrand. If the line 
C...               integral is performed on the interval J = [a, b], 
C...               then we still set FXY(x,y) = G(x) for the integrand 
C...               G(x).
C---------------------------------------------------------------------
C
      tol    = 1.d-20
      result = 0.d0
C
      xd = xe(2) - xe(1)
      yd = ye(2) - ye(1)
C
      if (abs(yd) .lt. tol) then
         d = abs(xd)
         do 10 k = 1, lquad_gs
            yy(k) = ye(1) 
 10      continue
      else
         d = sqrt(xd*xd + yd*yd)
         do 20 k = 1, lquad_gs
            yy(k) = z(k)*yd + ye(1) 
 20      continue
      end if  
      d2 = d*0.5d0     
C
      do 30 k = 1, lquad_gs
         xx(k) = z(k)*xd + xe(1)
 30   continue
C
      do 100 k = 1,lquad_gs
         x = xx(k)
         y = yy(k)
         f = fxy(x,y)
         if (abs(f) .lt. tol) go to 100
	 result = result  + f*w(k)
 100  continue
      result = result*d2
C
      return
      end
C=====================================================================
      subroutine intgrl2d(xe,ye,nd,result,fxy)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),xx(7),yy(7)
      common /quad_gs/ w(7),z(7),lquad_gs
C---------------------------------------------------------------------
C...  This subroutine computes the double integral of the function
C...  fxy(x,y) on the rectangle with vertices (XE(1),YE(1)), 
C...  (XE(2),YE(1)), (XE(2),YE(2)), and (XE(1),YE(2)) using Gauss
C...  quadrature. 
C...
C...  Input:
C...    XE, YE - x and y coordinates of the given rectangle.
C...             Note that XE(1) .le. XE(2) and YE(1) .le. YE(2). 
C...    ND     - ND = 2.
C...
C...  Output
C...    result - the rusulting value of the double integral.
C...
C...  Parameter:
C...    TOL    - to decide whether function is zero or not.
C...
C...  Working space:
C...    XX, YY - x and y coordinates of the points in the rectangle
C...             which correspond to the Gauss points.
C...
C...  Call module: 
C...    QUAD_GS  - This is a common data set which is defined in the 
C...               subroutine QUAD_DATA. 
C...    FXY(x,y) - an external function: the integrand. 
C---------------------------------------------------------------------
C
      tol    = 1.d-20
      result = 0.d0
C
      xd  = xe(2) - xe(1)
      yd  = ye(2) - ye(1)
      xd2 = xd*0.5d0
      yd2 = yd*0.5d0
C
      do 10 k = 1, lquad_gs
         xx(k) = z(k)*xd + xe(1)
         yy(k) = z(k)*yd + ye(1) 
 10   continue
C
      do 100 k = 1,lquad_gs
         x    = xx(k)
         temp = 0.d0
         do 50 j = 1,lquad_gs
            y    = yy(j)
            f    = fxy(x,y)
            if (abs(f) .lt. tol) go to 50
            temp = temp + f*w(j)
 50      continue
         temp   = temp*yd2
         result = result + temp*w(k)
 100  continue
      result = result*xd2
C
      return
      end
C=====================================================================










