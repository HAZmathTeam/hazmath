C======================================================================
      subroutine elt_stf_1(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the following Poisson Eq.:
C...
C...      - Laplacian(u) = f(x,y),
C...
C...  where f(x,y) is the right hand side.
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
C...  Call module: 
C...    QUAD   - This is a common data set which is defined in the 
C...             subroutine QUAD_ELT. 
C...    MENU   - This is a common data set which is defined in the 
C...             subroutine INPUT_DATA. 
C...    RHS    - a function: the right hand side of the given PDE.
C----------------------------------------------------------------------
      lquad = ich(11)
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
      do 200 k = 1,lquad
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
C
         rh  = rhs(x,y)
         wei = w(k)
C
         do 40 j = 1,nd
            do 30 i = 1,j
               ae(i,j) = ae(i,j) + (phx(j)*phx(i)+phy(j)*phy(i))*wei
 30         continue
            be(j) = be(j) + rh*phi(j,k)*wei
 40      continue
 200  continue
C
      do 220 j = 1,nd
         do 210 i = 1,j
            ae(i,j) = ae(i,j)*adet
 210     continue
         be(j) = be(j)*adet
 220  continue
C
      do 240 j = 1,nd-1
         ii = j + 1
         do 230 i = ii,nd
            ae(i,j) = ae(j,i)
 230     continue
 240  continue
C
      return
      end
C======================================================================
      subroutine elt_stf_2(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the following symmetric PDE:
C...
C...      - Laplacian(u) + gamma0(x,y)u = f(x,y),
C...
C...  where gamma0(x,y) is the coefficient of the lower order term and
C...  f(x,y) is the right hand side.
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
C...  Call module: 
C...    QUAD   - This is a common data set which is defined in the 
C...             subroutine QUAD_ELT. 
C...    MENU   - This is a common data set which is defined in the 
C...             subroutine INPUT_DATA. 
C...    GAMMA0  - a function: the coefficient of the zero order term of
C...             the given PDE.
C...    RHS    - a function: the right hand side of the given PDE.
C----------------------------------------------------------------------
      lquad = ich(11)
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
      do 200 k = 1,lquad
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
C     
         gam = gamma0(x,y)
         rh  = rhs(x,y)
         wei = w(k)
C
         do 40 j = 1,nd
            do 30 i = 1,j
               ae(i,j) = ae(i,j) + (phx(j)*phx(i)+phy(j)*phy(i) 
     >                           + gam*phi(j,k)*phi(i,k))*wei
 30         continue
            be(j) = be(j) + rh*phi(j,k)*wei
 40      continue
 200  continue
C
      do 220 j = 1,nd
         do 210 i = 1,j
            ae(i,j) = ae(i,j)*adet
 210     continue
         be(j) = be(j)*adet
 220  continue
C
      do 240 j = 1,nd-1
         ii = j + 1
         do 230 i = ii,nd
            ae(i,j) = ae(j,i)
 230     continue
 240  continue
C
      return
      end
C======================================================================
      subroutine elt_stf_3(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the following PDE:
C...
C...      - div (eps(x,y) grad(u) + beta(x,y)u) = f(x,y),
C...
C...  where eps(x,y) is a scalar function, beta(x,y) is a vector
C...  valued function of the convection term, and f(x,y) is the right 
C...  hand side. If the convection term is given by the gradient of 
C...  a potential function, i.e., beta(x,y) = (grad psi)(x,y), then use 
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
C...    RHS    - a function: the right hand side of the given PDE.
C----------------------------------------------------------------------
      lquad = ich(11)
      lpsi  = ich(23)
C
      tol = 1.d-14
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
      do 200 k = 1,lquad    
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
         ep  = eps(x,y)
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
               ae(i,j) = ae(i,j) + bt1*phi(j,k)*phx(i)*wei
 110        continue
 120     continue
C     
         do 140 j = 1,nd
            do 130 i = 1,nd
               ae(i,j) = ae(i,j) + bt2*phi(j,k)*phy(i)*wei
 130        continue
 140     continue
C     
         if (dabs(rh) .lt. tol) go to 200
         do 190 j = 1,nd
            be(j) = be(j) + rh*phi(j,k)*wei
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
C======================================================================
      subroutine elt_stf_4(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the following PDE:
C...
C...     - div(eps(x,y) grad(u) + beta(x,y)u) + gamma0(x,y)u = f(x,y),
C...
C...  where eps(x,y) is a scalar function, beta(x,y) is a vector
C...  valued function of the convection term, gamma0(x,y) is the 
C...  coefficient of the zero order term, and f(x,y) is the right 
C...  hand side.  If the convection term is given by the gradient of 
C...  a potential function, i.e., beta(x,y) = (grad psi)(x,y), then use 
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
      tol = 1.d-14
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
               ae(i,j) = ae(i,j) + bt1*phi(j,k)*phx(i)*wei
 110        continue
 120     continue
C     
         do 140 j = 1,nd
            do 130 i = 1,nd
               ae(i,j) = ae(i,j) + bt2*phi(j,k)*phy(i)*wei
 130        continue
 140     continue
C  
         do 160 j = 1,nd
            do 150 i = 1,nd
               ae(i,j) = ae(i,j) + gam*phi(j,k)*phi(i,k)*wei
 150        continue
 160     continue
C   
         if (dabs(rh) .lt. tol) go to 200
         do 190 j = 1,nd
            be(j) = be(j) + rh*phi(j,k)*wei
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
C======================================================================
      subroutine elt_stf_5(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the following PDE:
C...
C...     - div(eps(x,y) grad(u)) + beta(x,y) \cdot grad(u) = f(x,y),
C...
C...  where eps(x,y) is a scalar functiion, beta(x,y) is a vector
C...  valued function of the convection term, and f(x,y) is the right 
C...  hand side. If the convection term is given by the gradient of 
C...  a potential function, i.e., beta(x,y) = (grad psi)(x,y), then use 
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
C...    RHS    - a function: the right hand side of the given PDE.
C----------------------------------------------------------------------
      lquad = ich(11)
      lpsi  = ich(23)
C
      tol = 1.d-14
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
      do 200 k = 1,lquad
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
         ep  = eps(x,y)
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
               ae(i,j) = ae(i,j) + bt1*phx(j)*phi(i,k)*wei
 110        continue
 120     continue
C     
         do 140 j = 1,nd
            do 130 i = 1,nd
               ae(i,j) = ae(i,j) + bt2*phy(j)*phi(i,k)*wei
 130        continue
 140     continue
C  
         if (dabs(rh) .lt. tol) go to 200
         do 190 j = 1,nd
            be(j) = be(j) + rh*phi(j,k)*wei
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
C======================================================================
      subroutine elt_stf_6(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the following PDE:
C...
C...     - div(eps(x,y) grad(u)) + beta(x,y) \cdot grad(u) 
C...                             + gamma0(x,y)u = f(x,y),
C...
C...  where eps(x,y) is a scalar functiion, beta(x,y) is a vector
C...  valued function of the convection term, gamma0(x,y) is the 
C...  coefficient of the zero order term, and f(x,y) is the right hand 
C...  side. If the convection term is given by the gradient of 
C...  a potential function, i.e., beta(x,y) = (grad psi)(x,y), then use 
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
      tol = 1.d-14
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
               ae(i,j) = ae(i,j) + bt1*phx(j)*phi(i,k)*wei
 110        continue
 120     continue
C     
         do 140 j = 1,nd
            do 130 i = 1,nd
               ae(i,j) = ae(i,j) + bt2*phy(j)*phi(i,k)*wei
 130        continue
 140     continue
C  
         do 160 j = 1,nd
            do 150 i = 1,nd
               ae(i,j) = ae(i,j) + gam*phi(j,k)*phi(i,k)*wei
 150        continue
 160     continue
C   
         if (dabs(rh) .lt. tol) go to 200
         do 190 j = 1,nd
            be(j) = be(j) + rh*phi(j,k)*wei
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
C======================================================================
      subroutine elt_stf_7(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the PDE of the following most 
C...  general form:
C...
C...     - div (alpha(x,y) grad(u)) + beta(x,y) \cdot grad(u) 
C...                                + gamma0(x,y)u  = f(x,y),
C...
C...  where alpha(x,y) is a 2 X 2 matrix functiion, beta(x,y) is a 
C...  vector valued function of the convection term, gamma0(x,y) is the 
C...  coefficient of the lower order term, and f(x,y) is the right hand 
C...  side. If the convection term is given by the gradient of 
C...  a potential function, i.e., beta(x,y) = (grad psi)(x,y), then use 
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
C...    ALPHA  - a function: the set of four functions which are the
C...             coefficients of the leading term of the given PDE.
C...             This function should be given by the user beforehand.
C...    BETA   - a function: the set of two functions which are the
C...             coefficients of the convection term of the given PDE.
C...             This function should be given by the user beforehand.
C...             If the convection term is given by the gradient of 
C...             a potential function, i.e., BETA = grad PSI, then
C...             use the function PSI to define it.
C...    PSI    - BETA = grad PSI (see BETA).
C...    GAMMA0  - a function: the coefficient of the zero order term of
C...             the given PDE.
C...             This function should be given by the user beforehand.
C...    RHS    - a function: the right hand side of the given PDE.
C...             This function should be given by the user beforehand.
C----------------------------------------------------------------------
      lquad = ich(11)
      tol   = 1.d-14
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
      do 200 k = 1,lquad
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
         alp11 = alpha(1,1,x,y)
         alp12 = alpha(1,2,x,y)
         alp21 = alpha(2,1,x,y)
         alp22 = alpha(2,2,x,y)
         bt1   = beta(1,x,y)
         bt2   = beta(2,x,y)
         gam   = gamma0(x,y)
         rh    = rhs(x,y)
         wei   = w(k)
C     
         if (dabs(alp11) .lt. tol) go to 45
         do 40 j = 1,nd
            do 30 i = 1,nd
               ae(i,j) = ae(i,j) + alp11*phx(j)*phx(i)*wei
 30         continue
 40      continue
C     
 45      if (dabs(alp12) .lt. tol) go to 65
         do 60 j = 1,nd
            do 50 i = 1,nd
               ae(i,j) = ae(i,j) + alp12*phy(j)*phx(i)*wei
 50         continue
 60      continue
C     
 65      if (dabs(alp21) .lt. tol) go to 85
         do 80 j = 1,nd
            do 70 i = 1,nd
               ae(i,j) = ae(i,j) + alp21*phx(j)*phy(i)*wei
 70         continue
 80      continue
C
 85      if (dabs(alp22) .lt. tol) go to 105
         do 100 j = 1,nd
            do 90 i = 1,nd
               ae(i,j) = ae(i,j) + alp22*phy(j)*phy(i)*wei
 90         continue
 100     continue
C     
 105     if (dabs(bt1) .lt. tol) go to 125
         do 120 j = 1,nd
            do 110 i = 1,nd
               ae(i,j) = ae(i,j) + bt1*phx(j)*phi(i,k)*wei
 110        continue
 120     continue
C     
 125     if (dabs(bt2) .lt. tol) go to 145
         do 140 j = 1,nd
            do 130 i = 1,nd
               ae(i,j) = ae(i,j) + bt2*phy(j)*phi(i,k)*wei
 130        continue
 140     continue
C     
 145     if (dabs(gam) .lt. tol) go to 165
         do 160 j = 1,nd
            do 150 i = 1,nd
               ae(i,j) = ae(i,j) + gam*phi(j,k)*phi(i,k)*wei
 150        continue
 160     continue
C     
 165     continue
C
         if (dabs(rh) .lt. tol) go to 200
         do 190 j = 1,nd
            be(j) = be(j) + rh*phi(j,k)*wei
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
C======================================================================
      subroutine elt_stf_8(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) for the standard Galerkin
C...  method which is applicable to the PDE of the following most 
C...  general symmetric form:
C...
C...     - div (alpha(x,y) grad(u)) + gamma0(x,y)u  = f(x,y),
C...
C...  where alpha(x,y) is a 2 X 2 symmetric matrix functiion, gamma0(x,y) 
C...  is the coefficient of the lower order term, and f(x,y) is the 
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
C...  Call module: 
C...    QUAD   - This is a common data set which is defined in the 
C...             subroutine QUAD_ELT. 
C...    MENU   - This is a common data set which is defined in the 
C...             subroutine INPUT_DATA. 
C...    ALPHA  - a function: the set of four functions which are the
C...             coefficients of the leading term of the given PDE.
C...             This function should be given by the user beforehand.
C...    GAMMA0  - a function: the coefficient of the zero order term of
C...             the given PDE.
C...             This function should be given by the user beforehand.
C...    RHS    - a function: the right hand side of the given PDE.
C...             This function should be given by the user beforehand.
C----------------------------------------------------------------------
      lquad = ich(11)
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
      do 200 k = 1,lquad
         do 20 i = 1,nd
            phx(i) = c(1,1)*phix(i,k) + c(1,2)*phiy(i,k)
            phy(i) = c(2,1)*phix(i,k) + c(2,2)*phiy(i,k)
 20      continue
C
         x = b(1,1)*quadpt(1,k) + b(1,2)*quadpt(2,k) + xe(1)
         y = b(2,1)*quadpt(1,k) + b(2,2)*quadpt(2,k) + ye(1)
         alp11 = alpha(1,1,x,y)
         alp12 = alpha(1,2,x,y)
         alp22 = alpha(2,2,x,y)
         gam   = gamma0(x,y)
         rh    = rhs(x,y)
         wei   = w(k)
C     
         do 40 j = 1,nd
            do 30 i = 1,j
               ae(i,j) = ae(i,j) + (alp11*phx(j)*phx(i)
     >                           + alp22*phy(j)*phy(i)
     >                           + alp12*(phy(j)*phx(i)+phx(j)*phy(i)) 
     >                           + gam*phi(j,k)*phi(i,k))*wei
 30         continue
            be(j) = be(j) + rh*phi(j,k)*wei
 40      continue
 200  continue
C
      do 220 j = 1,nd
         do 210 i = 1,j
            ae(i,j) = ae(i,j)*adet
 210     continue
         be(j) = be(j)*adet
 220  continue
C
      do 240 j = 1,nd-1
         ii = j + 1
         do 230 i = ii,nd
            ae(i,j) = ae(j,i)
 230     continue
 240  continue
C
      return
      end
C======================================================================
      subroutine elt_ntrl_bdry(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      common /quad_gs_ntrl/ z(7),zw(5,7)
      common /menu/ ich(200)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) and the load vector BE(nd) corresponding to the edges
C...  of Robin or Neumann boundary. Boundary condition is of the 
C...  following form
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
C...  Parameter:
C...    TOL    - to decide whether functions are zero or not.
C...
C...  Call module: 
C...    QUAD_GS_NTRL - This is a common data set which is defined in  
C...                   the subroutine QUAD_DATA. 
C...    MENU   - This is a common data set which is defined in the 
C...             subroutine INPUT_DATA. 
C...    SIGMA  - a function: the coefficient of the Robin boundary
C...             condition. If SIGMA(x,y) = 0 for all x and y, then it
C...             becomes Neumann boundary.
C...    ZETA   - a function: the right hand side of a Robin or Neumann 
C...             boundary condition.
C----------------------------------------------------------------------
      lquad_gs_ntrl = ich(12)
      tol           = 1.d-14
C
C...  Initialize for AE and BE to be zero.
C
      call nullv(ae,nd*nd)
      call nullv(be,nd)
C
      xd = xe(2) - xe(1)
      yd = ye(2) - ye(1)
      d  = dsqrt(xd*xd + yd*yd)
      d2 = d*0.5d0
C
      do 100 k = 1,lquad_gs_ntrl
         x = z(k)*xd + xe(1)
         y = z(k)*yd + ye(1)
         sig   = sigma(x,y)
         zet   = zeta(x,y)
C
         if (dabs(sig) .lt. tol) go to 10
	 ae(1,1) = ae(1,1) + sig*zw(3,k)
	 ae(1,2) = ae(1,2) + sig*zw(4,k)
	 ae(2,2) = ae(2,2) + sig*zw(5,k)
C
 10      if (dabs(zet) .lt. tol) go to 100
         be(1)   = be(1)   + zet*zw(1,k)
	 be(2)   = be(2)   + zet*zw(2,k)
 100  continue
      ae(2,1) = ae(1,2)
C
      do 120 j = 1,2
	 do 110 i = 1,2
	    ae(i,j) = ae(i,j)*d2
 110     continue
	 be(j) = be(j)*d2
 120  continue
C
      return
      end
C======================================================================
