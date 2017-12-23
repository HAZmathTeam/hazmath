*
* functions.f: the functions used by harmonic average method
* every function takes (x,y) as input and returns a real value.
* all other parameters should be passed through common block.
*

*     alpha1: currently it is eps
      function alpha(x,y)
      implicit real*8(a-h,o-z)
      common /params/lbeta,lalpha,lmethod
      common /coeff/eps, beta(2), gamma, delta, frhs

      goto (10,20,30) lalpha
C
c      print *,'Error: no such alpha'
      alpha = 1.0d00
      return
C
 10   alpha = eps
      return
 20   alpha = 10.d00
      return
 30   print *,'Error: no such alpha'
      alpha = 1.d00
      return
      end

*     beta_x: the x-component of beta vector function
      function beta_x(x,y)
      implicit real*8(a-h,o-z)
      common /params/lbeta,lalpha,lmethod
      common /coeff/eps, beta(2), gamma, delta, frhs

      goto (10,20,30,40,50,60) lbeta

c      print *,'Error: no such beta'
      beta_x = 0.0d00
      return

 10   beta_x = 1.0
      return

*     C. Johnson's
 20   beta_x = 0.98480775301221
      return

*     Quarter circle
 30   beta_x = -y
      return

*     Stagnation point
 40   continue 
      r = sqrt((x-.5)**2+(y-.5)**2)
      if (r .gt. 1e-12) then
         beta_x = (0.5-y)/r
      else
         beta_x = 0d0
      endif
      return

*     Annular domain
 50   beta_x = - y / (x*x+y*y)
      return

 60   print *,'Error: no such beta'
      beta_x = 0.0d00
      return
      end


*     beta_y: the y-component of beta vector function
      function beta_y(x,y)
      implicit real*8(a-h,o-z)
      common /params/lbeta,lalpha,lmethod
      common /coeff/eps, beta(2), gamma, delta, frhs

      goto (10, 20, 30, 40, 50,60) lbeta
c      print *,'Error: no such beta'
      beta_y = 0.0d00
      return
*
 10   beta_y = 0.0
      return
 20   beta_y = 0.17364817766693
      return
 30   beta_y = x
      return

 40   continue 
      r = sqrt((x-.5)**2+(y-.5)**2)
      if (r .gt. 1e-12) then
         beta_y = (x-0.5)/r
      else
         beta_y = 0d0
      endif
      return

 50   beta_y = y / (x*x+y*y)
      return

 60   print *,'Error: no such beta'
      beta_y = 0.0d00
      return
      end
      real*8 function bcond(x,y)
      implicit real*8(a-h,o-z)
      if(abs(x-1.d00) .lt. 1.d-10) then
         bcond  = 1.0d00
      else
         bcond  = 0.0d00
      end if
      return
      end
C===================================================================
      real*8 function bernoulli(x)
      implicit real*8(a-h,o-z)
* THESE NUMBERS CORRESPOND TO ALPHA DEC STATION. 
      data 
     *      x1 /-37.42994775023704d+00/,
     *      x2 /-0.7629386345735235d-05/,
     *      x3 / 0.7620686773800261d-06/,
     *      x4 / 37.42994775023704d+00/,
     *      x5 / 37.42994775025371d+00/
      if( x .le. x1) then
         bernoulli = -x
      else
         if(x .lt. x2) then
            bernoulli = x/(dexp(x)-1.d+00)
         else
            if(x .le. x3) then
               bernoulli = 1.d00 - 0.5d+00*x
            else
               if(x .lt. x4) then
                  bernoulli = (x*dexp(-x))/(1.d+00-dexp(-x))
               else
                  if(x .lt. x5) then
                     bernoulli = (x*dexp(-x))
                  else
                     bernoulli = 0.0d+00
                  end if
               end if
            end if
         end if
      end if
      return
      end




