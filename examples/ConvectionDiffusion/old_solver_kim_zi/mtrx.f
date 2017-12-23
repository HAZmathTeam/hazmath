C=====================================================================
      subroutine pde_choice(istiff,
     I     nel,n,nd,ie,je,x,y,
     O     ia,ja,a,nnz,b,
     W     ned,nedge,idir,ip)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1),ia(1),ja(1),idir(1),nedge(1),ip(1)
      dimension a(1),b(1),x(1),y(1)
      common /menu/ ich(200)
      external elt_ntrl_bdry, elt_eafe_bdry,
     >         elt_stf_1,     elt_stf_2,
     >         elt_stf_3,     elt_stf_4,
     >         elt_stf_5,     elt_stf_6,
     >         elt_stf_7,     elt_stf_8,
     >         elt_stf_eafe,  elt_stf_upwind,
     >         elt_stf_sld,   elt_stf_supg
C---------------------------------------------------------------------
C...  Choose a PDE and compute the stiffness matrix and the load vector.
C...
C...  Parameter:
C...    istiff - choice of stiffness matrices.
C...             1 = Standard Galerkin    
C...             2 = EAFE
C...             3 = Streamline diffusion
C...             4 = Upwind 
C...             5 = SUPG
C---------------------------------------------------------------------
      go to (10,20,30,40,50,60,70,80), ich(1)
      write(*,*) ' **NO SUCH METHOD** '
      stop 5
C
 10   continue
      call mtrx(elt_stf_1,elt_ntrl_bdry,
     I          nel,n,nd,ie,je,x,y,
     O          ia,ja,a,nnz,b,
     W          ned,nedge,idir,ip)
      return
C
 20   continue
      call mtrx(elt_stf_2,elt_ntrl_bdry,
     I          nel,n,nd,ie,je,x,y,
     O          ia,ja,a,nnz,b,
     W          ned,nedge,idir,ip)
      return
C
 30   continue
      if (istiff .eq. 1) then
         call mtrx(elt_stf_3,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else if (istiff .eq. 2) then
         call mtrx(elt_stf_eafe,elt_eafe_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else
         write(*,*) ' Check the value of ISTIFF.'
         stop
      end if
C
 40   continue
      call mtrx(elt_stf_4,elt_ntrl_bdry,
     I          nel,n,nd,ie,je,x,y,
     O          ia,ja,a,nnz,b,
     W          ned,nedge,idir,ip)
      return
C
 50   continue
      if (istiff .eq. 1) then
         call mtrx(elt_stf_5,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else if (istiff .eq. 3) then
         call mtrx(elt_stf_sld,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else if (istiff .eq. 4) then
         call mtrx(elt_stf_upwind,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else if (istiff .eq. 5) then
         call mtrx(elt_stf_supg,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else
         write(*,*) ' Check the value of ISTIFF.'
         stop
      end if
C
 60   continue
      if (istiff .eq. 1) then
         call mtrx(elt_stf_6,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else if (istiff .eq. 3) then
         call mtrx(elt_stf_sld,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else if (istiff .eq. 4) then
         call mtrx(elt_stf_upwind,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else if (istiff .eq. 5) then
         call mtrx(elt_stf_supg,elt_ntrl_bdry,
     I        nel,n,nd,ie,je,x,y,
     O        ia,ja,a,nnz,b,
     W        ned,nedge,idir,ip)
         return
      else
         write(*,*) ' Check the value of ISTIFF.'
         stop
      end if
C
 70   continue
      call mtrx(elt_stf_7,elt_ntrl_bdry,
     I     nel,n,nd,ie,je,x,y,
     O     ia,ja,a,nnz,b,
     W     ned,nedge,idir,ip)
C
 80   continue
      call mtrx(elt_stf_8,elt_ntrl_bdry,
     I     nel,n,nd,ie,je,x,y,
     O     ia,ja,a,nnz,b,
     W     ned,nedge,idir,ip)
C
      return
      end
C=====================================================================
      subroutine mtrx(local,local_bdry,
     I     nel,n,nd,ie,je,x,y,
     O     ia,ja,a,nnz,b,
     W     ned,nedge,idir,ip)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1),ia(1),ja(1),idir(1),nedge(1),ip(1)
      dimension a(1),b(1),x(1),y(1)
      dimension ae(3,3),be(3),xe(3),ye(3),jedg(2),ab(2,2),bb(2)
      external local,local_bdry
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
C...  External:
C...    LOCAL      - choice of PDE and FE.
C...    LOCAL_BDRY - Robin or Neumann boundary data.
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
         call local(xe,ye,nd,ae,be)
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
cc      write(*,*) ' Element loop ended.'   !End of element by element loop. 
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
            call local_bdry(xe,ye,nedel,ab,bb)     
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



