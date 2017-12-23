C======================================================================
      subroutine elt_stf_eafe(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd),aed(3,3)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs/ wg(7),z(7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the EAFE (Edge Average
C...  Finite Element) method, which is applicable to the following PDE:
C...
C...      - div (eps(x,y) grad(u) + beta(x,y)u) = f(x,y),
C...
C...  where eps(x,y) is a scalar function, beta(x,y) is a vector
C...  valued function of the convection term which can be given as the 
C...  gradient of a potentail function PSI(x,y), and f(x,y) is the 
C...  right hand side.
C...
C...  Input:
C...    XE     - x coordinates of the given element.
C...    YE     - y coordinates of the given element.
C...    ND     - the degree of the freedom in each element.
C...
C...  Output
C...    AE, BE - the elementary stiffness matrix and the load vector.
C...
C...  Working space:
C...    B(2,2) - This matrix characterizes the transform from the 
C...             reference element to the current element.
C...    C(2,2) - C = B^(-t), i.e., the transpose of the inverse of B.
C...    PHX(3) - the inner product between the 1st row of C and the
C...             gradient vector of a reference basis function PHI,
C...             that is, <C(1,:), (PHIX,PHIY)>.
C...    PHY(3) - the inner product between the 2nd row of C and the
C...             gradient vector of a reference basis function PHI,
C...             that is, <C(2,:), (PHIX,PHIY)>.
C...
C...  Parameter:
C...
C...  Call module: 
C...    QUAD    - This is a common data set which is defined in the 
C...              subroutine QUAD_ELT. 
C...    QUAD_GS - This is a common data set which is defined in the 
C...              subroutine QUAD_DATA1. 
C...    MENU    - This is a common data set which is defined in the 
C...              subroutine INPUT_DATA. 
C...    EPS     - a function: the coefficient of the leading term of 
C...              the given PDE. This function should be given by 
C...              the user beforehand.
C...    BETA    - a function: the set of two functions which are the
C...              coefficients of the convection term of the given 
C...              PDE. Usual vector format, i.e., 
C...                 BETA(x,y) = <BETA(1,x,y), BETA(2,x,y)>.
C...    PSI     - a function: a potential function whose gradient is
C...              the convection term, i.e.,
C...                 BETA(x,y) = (grad PSI)(x,y)
C...    GAMMA   - a function: the coefficient of the zero order term 
C...              of the given PDE.
C...    RHS     - a function: the right hand side of the given PDE.
C----------------------------------------------------------------------
      lquad    = ich(11)
      lquad_gs = ich(13)
      lpsi     = ich(23)
      zero     = 1.d-12
C
C...  Initialize for ae and be to be zero.
C
      call nullv(ae,nd*nd)
cc      call nullv(aed,nd*nd)
      call nullv(be,nd)
C
C...  Computation of Laplacian term.
C
      b(1,1) = xe(2) - xe(1)
      b(1,2) = xe(3) - xe(1)
      b(2,1) = ye(2) - ye(1)
      b(2,2) = ye(3) - ye(1)
      det    = b(1,1)*b(2,2) - b(1,2)*b(2,1)
      adet   = dabs(det)
      det    = 1.d00/det
      c(1,1) =   b(2,2)*det
      c(1,2) = - b(2,1)*det
      c(2,1) = - b(1,2)*det
      c(2,2) =   b(1,1)*det
C
      do 200 k = 1,lquad
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
         rh  = rhs(x,y)
         wei = w(k)
C     
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         do 40 j = 1,nd
            do 30 i = 1,j-1
               ae(i,j) = ae(i,j) + (phx(j)*phx(i)+phy(j)*phy(i))*wei
cc               aed(i,j) = aed(i,j) + 
cc     >              eps(x,y)*(phx(j)*phx(i)+phy(j)*phy(i))*wei
 30         continue
            be(j) = be(j) + rh*phi(j,k)*wei
 40      continue
 200  continue
C
      do 220 j = 1,nd
         do 210 i = 1,j-1
            ae(i,j) = ae(i,j)*adet
cc            aed(i,j) = aed(i,j)*adet
 210     continue
         be(j) = be(j)*adet
 220  continue
C
      do 240 j = 1,nd-1
         ii = j + 1
         do 230 i = ii,nd
            ae(i,j) = ae(j,i)
cc            aed(i,j) = aed(j,i)
 230     continue
 240  continue
C
cc      aed(1,1) = - aed(2,1) - aed(3,1)
cc      aed(2,2) = - aed(1,2) - aed(3,2)
cc      aed(3,3) = - aed(1,3) - aed(2,3)
C
C...  End of computation of Laplacian term.
C
      do 330 i = 1,nd
         j = i + 1
         if (j .eq. nd+1) j = 1
         taux  = xe(i) - xe(j)
         tauy  = ye(i) - ye(j)
C
cc         beta_norm = 1.d00
         if (lpsi .eq. 1) then
            xmid  = (xe(i) + xe(j))*0.5d0
            ymid  = (ye(i) + ye(j))*0.5d0
            bt1 = beta(1,xmid,ymid)
            bt2 = beta(2,xmid,ymid)
cc            beta_norm = dsqrt(bt1*bt1+bt2*bt2)
            t_e = bt1*taux + bt2*tauy
         else if (lpsi .eq. 2) then
            t_e = psi(xe(i),ye(i)) - psi(xe(j),ye(j))
         end if	 
C
         alp_e = 0.d0
cc         emax = 0.0d00
         do 310 k = 1,lquad_gs
            x     = z(k)*taux + xe(j)
            y     = z(k)*tauy + ye(j)
            alp   = eps(x,y)
cc            if (beta_norm .gt. zero) then
cc               emax = dmax1(alp/beta_norm,emax)
cc            end if
            f     = 1.d0/alp
            alp_e = alp_e + f*wg(k)
 310     continue
C
         alp_e   = alp_e*0.5d0
         bx      = t_e*alp_e
C
         bern1   = bernoulli(bx)
         bern2   = bernoulli(-bx)
         ae(i,j) = ae(i,j)*bern1/alp_e
         ae(j,i) = ae(j,i)*bern2/alp_e
 330  continue
C
      ae(1,1) = - ae(2,1) - ae(3,1)
      ae(2,2) = - ae(1,2) - ae(3,2)
      ae(3,3) = - ae(1,3) - ae(2,3)
C
      call nullv(be,nd)
      call source_element_vector(xe,ye,be)
C
cc      if(emax .lt. 0.5d-02) then
cc         do k = 1 , nd
cc            do j = 1 , nd
cc               ae(j,k) = ae(j,k) + aed(j,k)
cc            end do
cc         end do
cc      end if
C
      return
      end
C======================================================================
      subroutine elt_eafe_bdry(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) corresponding to the edges
C...  of Robin or Neumann boundary. Boundary condition is of the 
C...  following form
C...
C...     {alpha(x,y) grad(u) + beta(x,y)u} \cdot normal + sigma(x,y)u 
C...                                                     = zeta(x,y),
C...
C...  where ALPHA(x,y) is a matrix or scalar function, NORMAL is the
C...  outward normal vector, and SIGMA(x,y) and ZETA(x,y) are scalar
C...  functions. If SIGMA(x,y) = 0 for all x and y, then it becomes 
C...  Neumann boundary.
C...
C...  Input:
C...    XE     - x coordinates of the end points of the given edge.
C...    YE     - y coordinates of the end points of the given edge.
C...    ND     - the degree of the freedom in each edge.
C...
C...  Output
C...    AE, BE - the elementary stiffness matrix and the load vector
C...             corresponding to the edge of Neumann boundary.
C...
C...  Call module: 
C...    SIGMA  - a function: the coefficient of the Robin boundary
C...             condition. If SIGMA(x,y) = 0 for all x and y, then it
C...             becomes Neumann boundary.
C...    ZETA   - a function: the right hand side of a Robin or Neumann 
C...             boundary condition.
C----------------------------------------------------------------------
C
      xd = xe(2) - xe(1)
      yd = ye(2) - ye(1)
      d  = dsqrt(xd*xd + yd*yd)
      d2 = d*0.5d0
C
      ae(1,1) = sigma(xe(1),ye(1))*d2
      ae(2,2) = sigma(xe(2),ye(2))*d2
C     
      be(1)   = zeta(xe(1),ye(1))*d2
      be(2)   = zeta(xe(2),ye(2))*d2
C
      return
      end
C======================================================================
      real*8 function bernoulli(x)
C======================================================================
      implicit real*8(a-h,o-z)
C----------------------------------------------------------------------
C...  These numbers correspond to alpha dec station. 
C----------------------------------------------------------------------
      data 
     >      x1 /-37.42994775023704d+00/,
     >      x2 /-0.7629386345735235d-05/,
     >      x3 / 0.7620686773800261d-06/,
     >      x4 / 37.42994775023704d+00/,
     >      x5 / 37.42994775025371d+00/
C
      if (x .le. x1) then
         bernoulli = -x
      else
         if (x .lt. x2) then
            bernoulli = x/(dexp(x)-1.d+00)
         else
            if (x .le. x3) then
               bernoulli = 1.d00 - 0.5d+00*x
            else
               if (x .lt. x4) then
                  bernoulli = (x*dexp(-x))/(1.d+00-dexp(-x))
               else
                  if (x .lt. x5) then
                     bernoulli = (x*dexp(-x))
                  else
                     bernoulli = 1.d-12
                  end if
               end if
            end if
         end if
      end if
      return
      end
C======================================================================
C      real*8 function bernoulli1(x)
C======================================================================
c      implicit real*8(a-h,o-z)
ccccc      real*8 n(21),id(21)
c      dimension coef(19)
cc      data n /1.d00,-1.d00,1.d00,-1.d00,1.d00,-1.d00,
cc     >     5.d00,-691.d00,7.d00,-3617.d00,43867.d00,-174611.d00,
cc     >     854513.d00,-236364091.d00,8553103.d00,
cc     >     -23749461029.d00,8615841276005.d00,
cc     >     -7709321041217.d00,2577687858367.d00,
cc     >     -26315271553053477373.d00,2929993913841559.d00/
cc      data id /1.d00,2.d00,6.d00,30.d00,
cc     >     42.d00,30.d00,66.d00,2730.d00,6.d00,
cc     >     510.d00,798.d00,330.d00,138.d00,
cc     >     2730.d00,6.d00,870.d00,14322.d00,510.d00,6.d00,
cc     >     1919190.d00,6.d00/
c      data coef/
c    >     0.166666666666666657d00,    
c     >     -0.333333333333333329d-01,
c     >     0.238095238095238082d-01,
c     >     -0.333333333333333329d-01,
c     >     0.757575757575757597d-01,
c     >     -0.253113553113553102d00,    
c     >     1.16666666666666674d00,    
c     >     -7.09215686274509771d00,    
c     >     54.9711779448621556d00,    
c     >     -529.124242424242425d00,    
c     >     6192.12318840579701d00,    
c     >     -86580.2531135531171d00,    
c     >     1425517.16666666674d00,    
c     >     -27298231.0678160936d00,    
c     >     601580873.900642395d00,    
c     >     -15116315767.0921574d00,    
c     >     429614643061.166687d00,    
c     >     -13711655205088.3340d00,    
c     >     488332318973593.188d00    /
C
c      k2 = 2
ccc      x = 2.d00
cccccccccc      zz = x*x / dble(k2)
c      ber = 0.0d00
c      if(dabs(x) .lt. 1.d-10) then
c        kk = 2
c         zz = 0.0d00
c         bernoulli1 = 1.d00
c         return
c      else
c         kk = ifix(sngl(log10(dabs(x))))
c         zz = 1.d00/ dble(k2)/x**(kk-2)
c      end if
c      do k = 1 , 19
c         ber = ber + zz*coef(k)
c         zz = (zz*(x*x)/(dble(k2+1)*dble(k2+2)))
c         k2 = k2 + 2
c      end do
c      ber = ber*x**kk + 1.d00-x*0.5d00
c      bernoulli1 = ber
c      return
c      end
C======================================================================










