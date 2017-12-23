C=====================================================================
      subroutine mtrx_sym(
     I     nel,n,nd,ie,je,x,y,
     O     ia,ja,a,nnz,b,
     W     ned,nedge,idir,ip)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1),ia(1),ja(1),idir(1),nedge(1),ip(1)
      dimension a(1),b(1),x(1),y(1)
      dimension ae(3,3),be(3),xe(3),ye(3),jedg(2),ab(2,2),bb(2)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  Compute (or assemble) stiffness matrix A and load vector b. 
C...  
C...  Input:  
C...    NEL    - number of elements.
C...    N      - number of nodes or dimension of the problem.
C...    ND     - order of AE and BE.
C...    IA, JA - structure of matrix A in RRCU. The order of A is N.
C...    IDIR   - Array which identifies Dirichlet nodes; it contains
C...             n+1 for a Dirichlet node, 0 otherwise.
C...    IE, JE - reference numbers of the nodes associated with the
C...             rows and columns of AE, i.e., global reference numbers 
C...             of nodes of a triangle.
C...    X, Y   - x and y coordinates of all the nodes.
C...    NED    - number of boundary edges.
C...    NEDGE  - global reference numbering of end points of each edge.
C...    NNZ    - the dimension of JA.
C...
C...  Output: 
C...    A      - numerical values of nonzeros of A in RRCU.
C...    B      - right-hand vector of system of linear equations.
C...
C...  Local memory:
C...    AE, BE - element matrix and vector, respectively, to be
C...             assembled in A, and B.
C...    XE,YE  - x and y coordinates of an element.
C...    JEDG   - global reference numbering of end points of an edge.
C...    AB,BB  - element matrix and vector corresponding to Robin or 
C...             Neumann boundary to be assembled in A and B.
C...
C...  Working space:
C...    IP     - of demension N, initialized to 0. IP is used and 
C...             then reset to 0. IP is the expanded integer array of
C...             pointers.
C...
C...  Note:
C...    The prescribed values of the Dirichlet unknowns should be 
C...    stored in the corresponding positions of B. In other words if
C...    i is a Dirichlet node with prescribed value C_i, then, before
C...    using this algorithm, set IDIR(i) = n+1, B(i) = C_i.
C--------------------------------------------------------------------
C
C...  Initialize to be zero.
C
      call inullv(ip,n)
      call nullv(b,n)
      call nullv(a,nnz)
C
C...  Element by element loop:
C
      do 100 i=1,nel
         j = 0
         do 20 it = ie(i), ie(i+1)-1
            j     = j + 1
            i1    = je(it)
            xe(j) = x(i1)
            ye(j) = y(i1)
 20      continue
C
         call elt_stf_3(xe,ye,nd,ae,be)
C     
         do 30 j = 1,nd
            k     = ie(i) + j - 1
            knode = je(k)
            if(idir(knode) .gt. n) then
               b(knode) = g(xe(j),ye(j))
            end if
 30      continue
         call assmbg(ia,ja,idir,ae,be,je(ie(i)),nd,a,b,ip,n)
 100  end do 
      write(*,*) ' Element loop ended.'   !End of element by element loop. 
C
C...  Check if Robin or Neumann boundary is identically zero or no such
C...  boundary exist at all.
C   
      if (ich(2) .eq. 0) go to 300
C
C...  Edge by edge loop.
C
      nedel = 2
      do 200 i = 1,ned
         jt1 = nedge(2*i-1)
         jt2 = nedge(2*i)
C...     Check if the edge is D-D. If it is, do nothing. 
         if (idir(jt1) .le. n .or. idir(jt2) .le. n) then
            jedg(1) = jt1
            jedg(2) = jt2
            xe(1) = x(jt1)
            ye(1) = y(jt1)
            xe(2) = x(jt2)
            ye(2) = y(jt2)
C     
            call elt_ntrl_bdry(xe,ye,nedel,ab,bb)     
            call assmbg(ia,ja,idir,ab,bb,jedg,nedel,a,b,ip,n)
         end if
 200  continue
C     
C...  End of edge by edge loop. 
C     
 300  do 330 i = 1, n
         if (idir(i) .gt. n) then
            do 310 k = ia(i),ia(i+1)-1
               j = ja(k)
               if(i .eq. j) then
                  a(k) = 1.d00
                  go to 320
               end if
 310        end do
 320        continue
         end if
 330  end do
C     
      return
      end
C=====================================================================
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
      tol           = 1.d-13
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
