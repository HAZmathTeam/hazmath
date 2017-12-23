C======================================================================
      subroutine elt_stf_sld(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the streamline diffusion
C...  method which is applicable to the following PDE:
C...
C...     - div(eps(x,y) grad(u)) + beta(x,y) \cdot grad(u) 
C...                             + gamma0(x,y)u = f(x,y),
C...
C...  where eps(x,y) is a scalar functiion, beta(x,y) is a vector
C...  valued function of the convection term, gamma0(x,y) is the 
C...  coefficient of the zero order term, and f(x,y) is the right hand 
C...  side. If the convection term is given by the gradient of a
C...  potential function, i.e., beta(x,y) = (grad psi)(x,y), then use 
C...  the function PSI to define it.
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
C...    TOL    - to decide whether coefficients are zero or not.
C...
C...  Call module: 
C...    QUAD   - This is a common data set which is defined in the 
C...             subroutine QUAD_ELT. 
C...    MENU   - This is a common data set which is defined in the 
C...             subroutine INPUT_DATA. 
C...    EPS    - a function: the coefficient of the leading term of 
C...             the given PDE.
C...    BETA   - a function: the set of two functions which are the
C...             coefficients of the convection term of the given PDE.
C...             If the convection term is given by the gradient of 
C...             a potential function, i.e., BETA = grad PSI, then
C...             use the function PSI to define it.
C...    PSI    - BETA = grad PSI (see BETA).
C...    GAMMA0  - a function: the coefficient of the zero order term of
C...             the given PDE.
C...    RHS    - a function: the right hand side of the given PDE.
C----------------------------------------------------------------------
      lquad = ich(11)
      lpsi  = ich(23)
C
      tol = 1.d-20
C
C...  Initialize for AE and BE to be zero.
C
      call nullv(ae,nd*nd)
      call nullv(be,nd)
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
C ... Compute delta...
C
      btmax = -1.d15
      epsmn =  1.d15
      hh    = dsqrt(adet)
C     
      do 10 k = 1 , nd
         bet1  = beta(1,xe(k),ye(k))
         bet2  = beta(2,xe(k),ye(k))
         epp   = eps(xe(k),ye(k))
         btmax = dmax1(dsqrt(bet1*bet1+bet2*bet2),btmax)
         epsmn = dmin1(epsmn,epp)
 10   continue
C     
      if(epsmn .gt. btmax*hh) then
         delta = 0.0d00
      else
         delta = hh/btmax
      end if

      delta = delta*.5d0

C      write(*,*) 'epsmn, btmax, hh, btmax*hh, delta',
C     >            epsmn, btmax, hh, btmax*hh, delta
C     
      do 200 k = 1,lquad
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
         ep  = eps(x,y)
         gam = gamma0(x,y)
         rh  = rhs(x,y)
         wei = w(k)
C 
         if (lpsi .eq. 1) then
            bt1   = beta(1,x,y)
            bt2   = beta(2,x,y)
         else if (lpsi .eq. 2) then
            bt1 = 0.0d00
            bt2 = 0.0d00
            do 30 kb = 1,nd
               bt1 = bt1 + psi(xe(kb),ye(kb))*phx(kb)
               bt2 = bt2 + psi(xe(kb),ye(kb))*phy(kb)
 30         continue
         end if
C
         do 50 j = 1,nd
            do 40 i = 1,nd
               ae(i,j) = ae(i,j) + ep*(phx(j)*phx(i) +
     >                   phy(j)*phy(i))*wei
 40         continue
 50      continue
C
         do 120 j = 1,nd
            do 110 i = 1,nd
               ae(i,j) = ae(i,j) + (bt1*phx(j) + bt2*phy(j))*
     >                   (phi(i,k)+delta*(bt1*phx(i)+bt2*phy(i)))*wei
 110        continue
 120     continue
C     
         do 160 j = 1,nd
            do 150 i = 1,nd
               ae(i,j) = ae(i,j) + gam*phi(j,k)*
     >                   (phi(i,k)+delta*(bt1*phx(i)+bt2*phy(i)))*wei
 150        continue
 160     continue
C   
         if (dabs(rh) .lt. tol) go to 200
         do 190 j = 1,nd
            be(j) = be(j) + rh*
     >              (phi(j,k)+delta*(bt1*phx(j)+bt2*phy(j)))*wei
 190     continue
 200  continue
C
      do 220 j = 1,nd
         do 210 i = 1,nd
            ae(i,j) = ae(i,j)*adet
 210     continue
         be(j) = be(j)*adet
 220  continue
C
      return
      end
C=====================================================================
