C=====================================================================
      subroutine lapl_mass(lm,nel,n,nd,ie,je,x,y,ia,ja,a,nnz,idir,ip)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1),ia(1),ja(1),idir(1),ip(1)
      dimension a(1),x(1),y(1)
      dimension ae(3,3),xe(3),ye(3)
C---------------------------------------------------------------------
C...  Compute (or assemble) the stiffness matrix A or the mass matrix M
C...  for Laplace equation with Neumann or/and Dirichlet conditions.
C...  
C...  Input:
C...    LM     - 1 for the stiffness matrix, 2 for the mass matrix M
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
C...    NNZ    - the dimension of JA.
C...
C...  Output: 
C...    A      - numerical values of nonzeros of A in RRCU.
C...
C...  Local memory:
C...    AE     - element matrix to be assembled in A.
C...    XE,YE  - x and y coordinates of an element.
C...
C...  Working space:
C...    IP     - of demension N, initialized to 0. IP is used and 
C...             then reset to 0. IP is the expanded integer array of
C...             pointers.
C--------------------------------------------------------------------
C
C...  Initialize to be zero.
C
      call inullv(ip,n)
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
         if (lm .eq. 1) then
            call lapl_local(xe,ye,nd,ae)
         else if (lm .eq. 2) then
            call mass_local(xe,ye,nd,ae)
         else 
            write(*,*) ' Error: Check for stiffness or mass matrix.'
         end if
C
         call assmbg_only_matrix(ia,ja,idir,ae,je(ie(i)),nd,a,ip,n)
 100  end do 
cc      write(*,*) ' Element loop ended.'   !End of element by element loop. 
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
C====================================================================
      subroutine assmbg_only_matrix(ia,ja,idir,ae,jep,nn,an,ip,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),jep(1),ae(1),an(1)
      dimension ip(1),idir(1)
C--------------------------------------------------------------------
C...  Numerical assembly of an element matrix AE into the nodal  
C...  assembly matrix A: General case.
C...
C...  Input:  
C...    IA, JA - structure of matrix A in RRCU. The order of A is N.
C...    IDIR   - Array which identifies Dirichlet nodes; it contains
C...             n+1 for a Dirichlet node, 0 otherwise.
C...    AE     - element matrix and vector, respectively, to be
C...             assembled in AN, and B.
C...    JEP    - reference numbers of the nodes associated with the
C...             rows and columns of AE, i.e., global reference numbers 
C...             of nodes of a triangle.
C...    NN     - order of AE.
C...
C...  Output: 
C...    AN     - numerical values of nonzeros of A in RRCU.
C...
C...  Working space:
C...    IP     - of demension N, initialized to 0. IP is used and 
C...             then reset to 0. IP is the expanded integer array of
C...             pointers.
C--------------------------------------------------------------------
      do 40 l = 1, nn
         i = jep(l)
         if (idir(i) .gt. n) go to 40
         k = l - nn
C
         do 20 ll = 1, nn
            k = k + nn
            j = jep(ll)
            if(idir(j) .gt. n) go to 10
            ip(j) = k
            go to 20
 10         continue
 20      continue
C
         iaa = ia(i)
         iab = ia(i+1) - 1
         kkk = 0
         do 30 j = iaa, iab
            k = ip(ja(j))
            if(k .eq. 0) go to 30
            an(j) = an(j) + ae(k)
            ip(ja(j)) = 0
            kkk = kkk + 1
            if(kkk .eq. nn) go to 40
 30      continue
 40   continue
C     
      return
      end
C======================================================================
      subroutine lapl_local(xe,ye,nd,ae)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary stiffness matrix 
C...  AE(nd,nd) only for the standard Galerkin method which is 
C...  applicable to the following Poisson Eq.:
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
C...    AE     - the elementary stiffness matrix.
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
C----------------------------------------------------------------------
      lquad = 3
C
      call nullv(ae,nd*nd)
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
         wei = w(k)
C
         do 40 j = 1,nd
            do 30 i = 1,j
               ae(i,j) = ae(i,j) + (phx(j)*phx(i)+phy(j)*phy(i))*wei
 30         continue
 40      continue
 200  continue
C
      do 220 j = 1,nd
         do 210 i = 1,j
            ae(i,j) = ae(i,j)*adet
 210     continue
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
      subroutine mass_local(xe,ye,nd,ae)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd)
      dimension phx(3),phy(3),b(2,2),c(2,2)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
C----------------------------------------------------------------------
C...  This subroutine computes the elementary mass matrix 
C...  AE(nd,nd) only for the standard Galerkin method.
C...
C...  Input:
C...    XE     - x coordinates of the given element.
C...    YE     - y coordinates of the given element.
C...    ND     - the degree of the freedom in each element.
C...
C...  Output
C...    AE     - the elementary stiffness matrix.
C...
C...  Working space:
C...    B(2,2) - This matrix characterizes the transform from the 
C...             reference element to the current element.
C...
C...  Call module: 
C...    QUAD   - This is a common data set which is defined in the 
C...             subroutine QUAD_ELT. 
C...    MENU   - This is a common data set which is defined in the 
C...             subroutine INPUT_DATA. 
C----------------------------------------------------------------------
      lquad = 3
C
      call nullv(ae,nd*nd)
C
      b(1,1) = xe(2) - xe(1)
      b(1,2) = xe(3) - xe(1)
      b(2,1) = ye(2) - ye(1)
      b(2,2) = ye(3) - ye(1)
      det    = b(1,1)*b(2,2) - b(1,2)*b(2,1)
      adet   = dabs(det)
C
      do 200 k = 1,lquad
         wei = w(k)
C
         do 40 j = 1,nd
            do 30 i = 1,j
               ae(i,j) = ae(i,j) + phi(j,k)*phi(i,k)*wei
 30         continue
 40      continue
 200  continue
C
      do 220 j = 1,nd
         do 210 i = 1,j
            ae(i,j) = ae(i,j)*adet
 210     continue
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
C=====================================================================
      subroutine l2_h1s_norm(lh,u,n,sl2_h1s,nel,nd,ie,je,x,y)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1),x(1),y(1),u(1)
C---------------------------------------------------------------------
C...  Compute l^2 or H^1 seminorm of the vector u using the L^2 or H^1
C...  scalar product, i.e., using the stiffness matrix A or the mass 
C...  matrix M.
C...  Input:  
C---------------------------------------------------------------------
      call l2_h1s_product(lh,u,u,n,sl2_h1s,nel,nd,ie,je,x,y)
      sl2_h1s = dsqrt(sl2_h1s)
C
      return
      end
C=====================================================================
      subroutine l2_h1s_product(lh,u,v,n,sl2_h1s,nel,nd,ie,je,x,y)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1), x(1),y(1),u(1),v(1)
      dimension ae(3,3),xe(3),ye(3)
C---------------------------------------------------------------------
C...  Compute L^2 or H^1 scalar product of U and V, i.e., (Au,v) or 
C...  (u,v), by computing the stiffness matrix A or the mass matrix M
C...  locally. It doesn't use the global stiffness matrix A and the 
C...  mass matrix M.
C...
C...  Input:
C...    LH     - 0 for L^2 norm and 1 for H^1 seminorm.
C...    NEL    - number of elements.
C...    N      - number of nodes or dimension of the problem.
C...    ND     - order of AE.
C...    IE, JE - reference numbers of the nodes associated with the
C...             rows and columns of AE, i.e., global reference numbers 
C...             of nodes of a triangle.
C...    X, Y   - x and y coordinates of all the nodes.
C...
C...  Output: 
C...    SL2_H1S - L^2 or H^1 scalar product of U and V.
C...
C...  Local memory:
C...    AE     - element matrix.
C...    XE,YE  - x and y coordinates of an element.
C---------------------------------------------------------------------
      sl2_h1s = 0.0d00
C
      do 100 nli=1,nel
         j = 0
         do 20 it = ie(nli), ie(nli+1)-1
            j     = j + 1
            i1    = je(it)
            xe(j) = x(i1)
            ye(j) = y(i1)
 20      continue
C     
         if (lh .eq. 0) then
            call mass_local(xe,ye,nd,ae)
         else if (lh .eq. 1) then
            call lapl_local(xe,ye,nd,ae)
         end if
C     
         do i = ie(nli),ie(nli+1)-1
            ik = i - ie(nli) + 1
            numbi = je(i)
            do j = ie(nli),ie(nli+1)-1
               jk = j - ie(nli) + 1
               numbj = je(j)
               sl2_h1s = sl2_h1s + ae(ik,jk)*u(numbi)*v(numbj)
            end do
         end do
 100  end do 
C     
      return
      end
C=====================================================================
      subroutine subd_norm(lh,
     >     ns,sol_loc,n,sol,work,u_exact,sl2_h1s,
     >     ksub,ielsub,nel,ie,je,nd,x,y,jsub_all)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension work(1),u_exact(1),sol_loc(1),sol(1),
     >     ielsub(1),ie(1),je(1),x(1),y(1),jsub_all(1),
     >     xe(3),ye(3),ae(3,3)
C---------------------------------------------------------------------
C...  Compute the L^2 norm or H^1 seminorm of errors on each subdomain.
C...
C...  Parameter:
C...    LH - 0 for L^2 norm; 1 for H^1 seminorm
C---------------------------------------------------------------------      
      call nullv(work,n)
      do jk = 1,ns
         k = jsub_all(jk)
         work(k) = u_exact(k) - (sol(k) + sol_loc(jk))
      end do
C
      sl2_h1s = 0.0d00
C
      do nli = 1, nel
         if(ielsub(nli) .ne. ksub) go to 110
         j = 0
         do 20 it = ie(nli), ie(nli+1)-1
            j     = j + 1
            i1    = je(it)
            xe(j) = x(i1)
            ye(j) = y(i1)
 20      continue
C     
         if (lh .eq. 0) then
            call mass_local(xe,ye,nd,ae)
         else if (lh .eq. 1) then
            call lapl_local(xe,ye,nd,ae)
         end if
C     
         do i = ie(nli),ie(nli+1)-1
            ik = i - ie(nli) + 1
            numbi = je(i)
            do j = ie(nli),ie(nli+1)-1
               jk = j - ie(nli) + 1
               numbj = je(j)
               sl2_h1s = sl2_h1s + ae(ik,jk)*work(numbi)*work(numbj)
            end do
         end do
 110     continue
      end do 
C     
      return
      end
C=====================================================================


