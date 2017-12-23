C====================================================================
      subroutine rdmesh(lread,mt,ie,je,ned,nedge,idir,
     >     x,y,nel,n,nd)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1),x(1),y(1),idir(1),nedge(1)
C--------------------------------------------------------------------
C...  Read the given mesh: Elements and corresponding nodes and 
C...  coordinates of each node are given in the file reference mt. 
C...
C...  Input:
C...    LREAD  - format of mesh: 1-formatted, 2-unformatted.
C...    MT     - reference for file name containing mesh.
C...    NEL    - number of elements.
C...    N      - number of nodes.
C...    NED    - number of boundary edges.
C...
C...  Parameter:
C...    ND     - degrees of freedom, for example, 3: linear elements;
C...             4: bilinear elements; 6: quadratic elements.
C...
C...  Output:
C...    NEDGE  - global reference numbering of end points of each edge.
C...    IE, JE - mesh connectivity matrix in RRCU.
C...    X, Y   - x and y coordinates of each node.
C...    IDIR   - array of dimension N which contains the value  
C...             N+1 in the positions corresponding to Dirichlet 
C...             nodes and 0 elsewhere.
C-------------------------------------------------------------------- 
      call inullv(idir,n)
C
C...  Read the mesh...
C
      if (lread .eq. 1) then
         read(mt,*) (je(i),i=1,3*nel)
         read(mt,*) (x(i),y(i),i=1,n)
         read(mt,*) (nedge(2*k-1),nedge(2*k),ie(k), k=1,ned)
      else
         read(mt) (je(i),i=1,3*nel)
         read(mt) (x(i),y(i),i=1,n)
         read(mt) (nedge(2*k-1),nedge(2*k),ie(k), k=1,ned)
      end if
C
      do 30 k = 1, ned
         if (ie(k) .gt. 0) then
            idir(nedge(2*k-1)) = n + 1
            idir(nedge(2*k)) = n + 1
         end if
 30   continue
C     
      ie(1) = 1
      do 40 k = 1, nel
         ie(k+1) = ie(k) + nd
 40   continue
C
      return
      end
C====================================================================
      subroutine rdmesh_perm(lread,mt,ie,je,ned,nedge,idir,iperm,
     >     x,y,nel,n,nd)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension ie(1),je(1),x(1),y(1),idir(1),nedge(1),iperm(1)
C--------------------------------------------------------------------
C...  Read the given mesh: Elements and corresponding nodes and 
C...  coordinates of each node are given in the file reference mt. 
C...  The node numbering is randomly chosen using the subroutine
C...  RANDOM and the ramdom numbering is stored the permutation array
C...  IPERM.
C...
C...  Input:
C...    LREAD  - format of mesh: 1-formatted, 2-unformatted.
C...    MT     - reference for file name containing mesh.
C...    NEL    - number of elements.
C...    N      - number of nodes.
C...    NED    - number of boundary edges.
C...
C...  Parameter:
C...    ND     - degrees of freedom, for example, 3: linear elements;
C...             4: bilinear elements; 6: quadratic elements.
C...    NPERM  - 1: if considering the random node numbering;
C...             0: if considering the original node numbering.
C...
C...  Output:
C...    NEDGE  - global reference numbering of end points of each edge.
C...    IE, JE - mesh connectivity matrix in RRCU.
C...    X, Y   - x and y coordinates of each node.
C...    IDIR   - array of dimension N+1 which contains the value  
C...             N+1 in the positions corresponding to Dirichlet 
C...             nodes and 0 elsewhere.
C...    IPERM  - a permutation of N integers {1,2,...,N}.
C-------------------------------------------------------------------- 
      call inullv(idir,n)
C
      nperm = 0 !1 !0
      if (nperm .eq. 1) then
         call random(iperm,n)
      else if (nperm .eq. 0) then
         call iseqv(iperm,n)
      else
         print*, 'Check the value of NPERM in the subroutine RDMESH.'
      end if
C
C...  Read the mesh...
C
      if (lread .eq. 1) then
         read(mt,*) (je(i),i=1,3*nel)
         read(mt,*) (x(iperm(i)),y(iperm(i)),i=1,n)
         read(mt,*) (nedge(2*k-1),nedge(2*k),ie(k), k=1,ned)
      else
         read(mt) (je(i),i=1,3*nel)
         read(mt) (x(iperm(i)),y(iperm(i)),i=1,n)
         read(mt) (nedge(2*k-1),nedge(2*k),ie(k), k=1,ned)
      end if
C
      do 10 k = 1, 3*nel
         je(k) = iperm(je(k))
 10   continue
C
      do 20 k = 1, 2*ned
         nedge(k) = iperm(nedge(k))
 20   continue
C
      do 30 k = 1, ned
         if (ie(k) .gt. 0) then
            idir(nedge(2*k-1)) = n + 1
            idir(nedge(2*k)) = n + 1
         end if
 30   continue
C     
      ie(1) = 1
      do 40 k = 1, nel
         ie(k+1) = ie(k) + nd
 40   continue
C
      return
      end
C====================================================================
      subroutine read_u(mt,ja,au,rhs,n,nnz)
C====================================================================
      implicit none
      integer ja(1),i,n,nnz,mt
      double precision au(1),rhs(1)
C--------------------------------------------------------------------
C...  Read the matrix and the right hand side vector in unformatted
C...  form.
C--------------------------------------------------------------------
C...  The first two records have been read already ... 
C
      read(mt)(ja(i),i=1,nnz)
      read(mt)(au(i),i=1,nnz)
      read(mt)(rhs(i),i=1,n)
C
      close(mt)
      return
      end
C====================================================================
      subroutine read_f(mt,ja,au,rhs,n,nnz)
C====================================================================
      implicit none
      integer ja(1),i,n,nnz,mt
      double precision au(1),rhs(1)
C--------------------------------------------------------------------
C...  Read the matrix and the right hand side vector in formatted
C...  form.
C--------------------------------------------------------------------
C...  The first two records have been read already ... 
C
      read(mt,*)(ja(i),i=1,nnz)
      read(mt,*)(au(i),i=1,nnz)
      read(mt,*)(rhs(i),i=1,n)
C
      close(mt)
      return
      end
C====================================================================
      subroutine lpri_more(ip,jp,n)
C====================================================================
      integer*4 ip(1),jp(1),n,k,kk
C--------------------------------------------------------------------
C     Print IP and JP.
C--------------------------------------------------------------------
      do k = 1, n
         write(*,*) ' row ', k, ' entries ',ip(k+1)-ip(k)
         write(*,*)' elements: ', (jp(kk),kk=ip(k),ip(k+1)-1)
cc         if(mod(k,10) .eq. 1) read(*,*)
      end do
      return
      end
C====================================================================
      subroutine lpri(ip,jp,n)
C====================================================================
      integer*4 ip(1),jp(1),n,k,kk
C--------------------------------------------------------------------
C     Print IP and JP.
C--------------------------------------------------------------------
      do k = 1, n
         write(*,*) ' k: ', k
         write(*,*)(jp(kk),kk=ip(k),ip(k+1)-1)
cc        if(mod(k,10) .eq. 1) read(*,*)
      end do
      return
      end
C====================================================================
      subroutine lprr(ip,jp,p,n)
C====================================================================
      real*8 p(1)
      integer*4 ip(1),jp(1),n,k,kk
C--------------------------------------------------------------------
C     Print IP, JP, and P.
C--------------------------------------------------------------------
      do k = 1, n
         write(*,*) ' row : ', k
         write(*,*)(jp(kk),p(kk),kk=ip(k),ip(k+1)-1)
cc         if(mod(k,10) .eq. 1) read(*,*)
      end do
      return
      end
C====================================================================
      subroutine output_sol(ie,je,x,y,u,n,nel)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension x(1),y(1),u(1),ie(1),je(1)
C--------------------------------------------------------------------
C...  Print out the x- and y-coordinates.
C...  Print out the solution vector.
C...
C...  Parameter:
C...    N      - number of nodes.
C...    NEL    - number of elements.
C...
C...  Input:
C...    IE, JE - mesh connectivity matrix in RRCU.
C...    X, Y   - x and y coordinates of each node.
C...    U      - solution vector.
C-------------------------------------------------------------------- 
      do 20 k = 1, nel
         ie1 = je(ie(k))
         ie2 = je(ie(k)+1)
         ie3 = je(ie(k)+2)
         write(800,'(3e20.12$)')x(ie1),x(ie2),x(ie3)
 20   continue
      write(800,*)
C
      do 21 k = 1, nel
         ie1 = je(ie(k))
         ie2 = je(ie(k)+1)
         ie3 = je(ie(k)+2)
         write(800,'(3e20.12$)')y(ie1),y(ie2),y(ie3)
 21   continue
      write(800,*)
C
      do 22 k = 1, nel
         ie1 = je(ie(k))
         ie2 = je(ie(k)+1)
         ie3 = je(ie(k)+2)
         write(800,'(3e20.12$)')u(ie1),u(ie2),u(ie3)
 22   continue
      write(800,*)
C
C...  Output for the graph plot...
C
      do 88 k = 1, n
         write(199,'(e20.12$)')x(k)
 88   continue
      write(199,*)
      do 99 k = 1, n
         write(199,'(e20.12$)')y(k)
 99   continue
      write(199,*)
      return
      end
C====================================================================
      subroutine wradj_r(xadj,adjncy,a,n,lpp)
C====================================================================
      integer xadj(1),adjncy(1),n,lpp,i,j,k
      real*8 a(1)
C--------------------------------------------------------------------
C...  Write the structure of A and A itself.
C--------------------------------------------------------------------
      if(n .gt. 3000) return
C      
      rewind lpp
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(i7$)') i
         end do
      end do
      write(lpp,*)
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(i7$)') j
         end do
      end do
      write(lpp,*)
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            write(lpp,'(e25.17$)') a(k)
         end do
      end do
      write(lpp,*)
      close(lpp)
C
      return
      end   
C====================================================================
      subroutine outmat(xadj,adjncy,n,lpp,a,rh)
C====================================================================
      integer xadj(1),adjncy(1),n,lpp,i,j,k,l
      real*8 a(1),rh(1)
C---------------------------------------------------------------------
C...  Write the structure of A, A itself, and right hand side RH. 
C---------------------------------------------------------------------
      rewind lpp
      l = 1
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(i7$)') i
         end do
      end do
      write(lpp,*)
C
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(i7$)') j
         end do
      end do
      write(lpp,*)
C
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(e20.12$)') a(k)
         end do
      end do
      write(lpp,*)
C
      do i = 1, n
         write(lpp+1,'(e20.12$)') rh(i)
      end do
      write(lpp,*)
      close(lpp)
C
      return
      end   
C====================================================================
      subroutine out_graph(xadj,adjncy,n,lpp)
C====================================================================
      integer xadj(1),adjncy(1),n,lpp,i,j,k,l
C---------------------------------------------------------------------
C...  Write the structure of A.
C---------------------------------------------------------------------
      rewind lpp
      l = 1
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(i7$)') i
         end do
      end do
      write(lpp,*)
C
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(i7$)') j
         end do
      end do
      write(lpp,*)
C
      do i = 1, n
         do k = xadj(i),xadj(i+1)-1
            j = adjncy(k)
            write(lpp,'(i7$)') 1
         end do
      end do
C
      return
      end   
C====================================================================
      subroutine outmat1(ia,ja,a,n,nnz)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),a(1)
C--------------------------------------------------------------------
C...  Print out the number of nonzero entries of the matrix A.
C...  Print out the nonzero column indices of each row of A.
C...
C...  Parameter:
C...    N      - number of nodes.
C...    NNZ    - number of nonzero entries of the matrix A.
C...
C...  Input:
C...    IA, JA - structure of the nodal assembly matrix of 
C...             dimension NxN, in RRCU form.
C-------------------------------------------------------------------- 
      write(*,*) ' Number of nonzeros :', nnz
C
      do 10 k = 1, n
         write(*,*) ' Row :', k
         write(*,*) (ja(kk),kk=ia(k),ia(k+1)-1)
         write(*,*) ' Entries of A :'
         write(*,*) (a(kk),kk=ia(k),ia(k+1)-1)
 10   continue
C
      return
      end
C====================================================================
      subroutine output(lpp,x,y,u,n)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension x(1),y(1),u(1)
C--------------------------------------------------------------------
C...  Print out the solution.
C-------------------------------------------------------------------- 
      write(*,'(A)') ' Solution:'
      do 21 k = 1, n
         if(dabs(u(k)) .lt. 1d-15) u(k) = 0d0
         write(lpp,'(e27.17$)') u(k)
 21   continue
      write(lpp,*)
      do 22 k = 1, n
         if(dabs(x(k)) .lt. 1d-15) x(k) = 0d0
         write(lpp,'(e27.17$)') x(k)
 22   continue
      write(lpp,*)
      do 23 k = 1, n
         if(dabs(y(k)) .lt. 1d-15) y(k) = 0d0
         write(lpp,'(e27.17$)') y(k)
 23   continue
      write(lpp,*)
C
      write(*,*) ' Solution has been written.'
      return
      end
C====================================================================
      subroutine wrxy(x,y,n,lpp)
C====================================================================
      integer i,n
      real*8 x(1),y(1)
C--------------------------------------------------------------------
C...  Write out vectors X and Y.
C--------------------------------------------------------------------
      if(n .gt. 3000) return
C
      rewind lpp
      do i = 1, n
         write(lpp,'(e15.7$)') x(i)
      end do
      write(lpp,*)
      do i = 1, n
         write(lpp,'(e15.7$)') y(i)
      end do
      write(lpp,*)
      close(lpp)
      return
      end   
C====================================================================
      subroutine wrsol01(sol,n,lpp)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension  sol(1)
C--------------------------------------------------------------------
C...  Write sol in the file reference lpp.
C--------------------------------------------------------------------
      write(*,*) ' Writing the solution.'
      rewind(lpp)
      do k = 1 , n
         write(lpp,'(e17.8$)')sol(k)
      end do
      write(lpp,*)
      write(*,*) ' Solution has been written.'
      return
      end
C====================================================================
      subroutine wruv_plus(u,v,n,lpp)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension  u(1),v(1)
C--------------------------------------------------------------------
C...  Write u + v in the file reference lpp.
C--------------------------------------------------------------------
      rewind(lpp)
      do k = 1 , n
         write(lpp,'(e17.8$)')u(k)+v(k)
      end do
      write(lpp,*)
      return
      end
C====================================================================
      subroutine write_soln(lvl,soln,n,iform)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension soln(1)
      character*100 fname
C--------------------------------------------------------------------
C...  Write the solution vector in the file 'soln_lvl'.
C...
C...  Parameter:
C...    IFORM  - formatted or unformatted,  1 or 2.
C--------------------------------------------------------------------
      if (lvl .eq. 2) then
         fname  = 'soln2'
      else if (lvl .eq. 3) then
         fname  = 'soln3'
      else if (lvl .eq. 4) then
         fname  = 'soln4'
      else if (lvl .eq. 5) then
         fname  = 'soln5'
      else if (lvl .eq. 6) then
         fname  = 'soln6'
      else if (lvl .eq. 7) then
         fname  = 'soln7'
      else if (lvl .eq. 8) then
         fname  = 'soln8'
      else if (lvl .eq. 9) then
         fname  = 'soln9'
      else if (lvl .eq. 10) then
         fname  = 'soln10'
      else if (lvl .eq. 11) then
         fname  = 'soln11'
      else
         fname  = 'soln'
      end if
C
      if (iform .eq. 1) then
         open(20,file=fname,form='formatted',status='unknown')
      else
         open(20,file=fname,form='unformatted',status='unknown')
      end if
C
      if (iform .eq. 1) then
          write(20,*) (soln(l), l = 1, n)
      else
         write(20) (soln(l), l = 1, n)
      end if
C
      write(*,*) ' Solution has been written in the file solt_i. '
C
      endfile(20)
      close(20)
C
      return
      end
C====================================================================
      subroutine write_matrix_rhs(lvl,n,ia,ja,a,rhs,nnz,iform)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),a(1),rhs(1)
      character*100 matname
C--------------------------------------------------------------------
C...  Write the matrix and right_hand_side in the file 'matname_i'.
C...
C...  Parameter:
C...    IFORM  - formatted or unformatted,  1 or 2.
C--------------------------------------------------------------------
      if (lvl .eq. 2) then
         matname  = 'matrices/matrix_rhs.unfcd2'
      else if (lvl .eq. 3) then
         matname  = 'matrices/matrix_rhs.unfcd3'
      else if (lvl .eq. 4) then
         matname  = 'matrices/matrix_rhs.unfcd4'
      else if (lvl .eq. 5) then
         matname  = 'matrices/matrix_rhs.unfcd5'
      else if (lvl .eq. 6) then
         matname  = 'matrices/matrix_rhs.unfcd6'
      else if (lvl .eq. 7) then
         matname  = 'matrices/matrix_rhs.unfcd7'
      else if (lvl .eq. 8) then
         matname  = 'matrices/matrix_rhs.unfcd8'
      else if (lvl .eq. 9) then
         matname  = 'matrices/matrix_rhs.unfcd9'
      else if (lvl .eq. 10) then
         matname  = 'matrices/matrix_rhs.unfcd10'
      else if (lvl .eq. 11) then
         matname  = 'matrices/matrix_rhs.unfcd11'
      else
         matname  = 'matrices/matrix_rhs.unf'
      end if
C
      if(iform .eq. 1) then
         open(19,file=matname,status='unknown',form='formatted')
      else
         open(19,file=matname,status='unknown',form='unformatted')
      end if
C
      if (iform .eq. 1) then
         write(19,*) n
         write(19,*) (ia(ii),ii=1,n+1)
         write(19,*) (ja(ii),ii=1,nnz)
         write(19,*) ( a(ii),ii=1,nnz)
         write(19,*) (rhs(ii),ii=1,n)
      else if (iform .eq. 2) then
         write(19) n
         write(19) (ia(ii),ii=1,n+1)
         write(19) (ja(ii),ii=1,nnz)
         write(19) ( a(ii),ii=1,nnz)
         write(19) (rhs(ii),ii=1,n)
      end if
C
      write(*,*) 'Printing the matrix and RHS is done.'
C
      endfile(19)
      close(19)
C
      return
      end
C====================================================================
      subroutine read_soln(lvl,soln,n,iform)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension soln(1)
      character*100 fname
C--------------------------------------------------------------------
C...  read the solution vector in the file 'soln_lvl'.
C...
C...  Parameter:
C...    IFORM  - formatted or unformatted,  1 or 2.
C--------------------------------------------------------------------
      if (lvl .eq. 2) then
         fname  = 'matrices/soln_cd2'
      else if (lvl .eq. 3) then
         fname  = 'matrices/soln_cd3'
      else if (lvl .eq. 4) then
         fname  = 'matrices/soln_cd4'
      else if (lvl .eq. 5) then
         fname  = 'matrices/soln_cd5'
      else if (lvl .eq. 6) then
         fname  = 'matrices/soln_cd6'
      else if (lvl .eq. 7) then
         fname  = 'matrices/soln_cd7'
      else if (lvl .eq. 8) then
         fname  = 'matrices/soln_cd8'
      else if (lvl .eq. 9) then
         fname  = 'matrices/soln_cd9'
      else if (lvl .eq. 10) then
         fname  = 'matrices/soln_cd10'
      else if (lvl .eq. 11) then
         fname  = 'matrices/soln_cd11'
      else
         fname  = 'matrices/soln'
      end if
C
      if (iform .eq. 1) then
         open(20,file=fname,form='formatted',status='unknown')
      else
         open(20,file=fname,form='unformatted',status='unknown')
      end if
C
      if (iform .eq. 1) then
          read(20,*) (soln(l), l = 1, n)
      else
         read(20) (soln(l), l = 1, n)
      end if
C
      write(*,*) ' Solution has been read from the file soln_cd_i.'
C
      endfile(20)
      close(20)
C
      return
      end
C====================================================================
      subroutine read_param(lsolve,tol,maxitr,iformatted)
C====================================================================
      implicit none
      include 'global.h'      
      include 'paramsmg.h'
      integer k,lsolve,maxitr,iformatted
      double precision tol
C--------------------------------------------------------------------
C     Read parameters.
C--------------------------------------------------------------------
      do k = 1, 100
         matrix_file(k:k) = ' '
         soln_file(k:k) = ' ' 
      end do
C
      open(11,file='solver.in',status='unknown')
      read(11,'(a)') matrix_file
      call input_name(matrix_file)
      read(11,'(a)') soln_file
      call input_name(soln_file)
      write(*,*) 'Matrix and RHS input file: ', matrix_file
      write(*,*) 'Computed exact solution file: ', soln_file
      read(11,*) iformatted
      read(11,*) lsolve
      read(11,*) tol
      read(11,*) maxitr
      read(11,*) iprsm
      read(11,*) ipssm
      read(11,*) nprsm
      read(11,*) npssm
      read(11,*) max_pr_smooths
      read(11,*) icycle
      read(11,*) isolv_coarse
      read(11,*) ich_ord
C
      close(11)
      write(*,*) 
      return
      end
C====================================================================
      subroutine read_param_multi(lsolve,tol,maxitr,iformatted,kk,jj)
C====================================================================
      implicit none
      include 'global.h'      
      include 'paramsmg.h'
      integer k,lsolve,maxitr,iformatted,kk,jj
      double precision tol
C--------------------------------------------------------------------
C     Read parameters.
C--------------------------------------------------------------------
      do k = 1, 100
         matrix_file(k:k) = ' '
         soln_file(k:k) = ' '
      end do
C
      if (kk .eq. 4) then
         matrix_file = 'matrix_rhs.unfcd4'
         soln_file = 'soln_cd4'
      else if (kk .eq. 5) then
         matrix_file = 'matrix_rhs.unfcd5'
         soln_file = 'soln_cd5'
      else if (kk .eq. 6) then
         matrix_file = 'matrix_rhs.unfcd6'
         soln_file = 'soln_cd6'
      else if (kk .eq. 7) then
         matrix_file = 'matrix_rhs.unfcd7'
         soln_file = 'soln_cd7'
      else if (kk .eq. 8) then
         matrix_file = 'matrix_rhs.unfcd8'
         soln_file = 'soln_cd8'
      else if (kk .eq. 9) then
         matrix_file = 'matrix_rhs.unfcd9'
         soln_file = 'soln_cd9'
      else if (kk .eq. 10) then
         matrix_file = 'matrix_rhs.unfcd10'
         soln_file = 'soln_cd10'
      else if (kk .eq. 11) then
         matrix_file = 'matrix_rhs.unfcd11'
         soln_file = 'soln_cd11'
      end if 
C
      if (jj .eq. 1) then
         open(11,file='solv_multi.in1',status='unknown')
      else if (jj .eq. 2) then
         open(11,file='solv_multi.in2',status='unknown')
      else if (jj .eq. 3) then
         open(11,file='solv_multi.in3',status='unknown')
      else if (jj .eq. 4) then
         open(11,file='solv_multi.in4',status='unknown')
      else if (jj .eq. 5) then
         open(11,file='solv_multi.in5',status='unknown')
      else if (jj .eq. 6) then
         open(11,file='solv_multi.in6',status='unknown')
      else if (jj .eq. 7) then
         open(11,file='solv_multi.in7',status='unknown')
      else if (jj .eq. 8) then
         open(11,file='solv_multi.in8',status='unknown')
      else if (jj .eq. 9) then
         open(11,file='solv_multi.in9',status='unknown')
      else if (jj .eq. 10) then
         open(11,file='solv_multi.in10',status='unknown')
      end if 
C
      call input_name(matrix_file)
      call input_name(soln_file)
      write(*,*) 'Matrix and RHS input file: ', matrix_file
      write(*,*) 'Computed exact solution file: ', soln_file
      read(11,*) iformatted
      read(11,*) lsolve
      read(11,*) tol
      read(11,*) maxitr
      read(11,*) iprsm
      read(11,*) ipssm
      read(11,*) nprsm
      read(11,*) npssm
      read(11,*) max_pr_smooths
      read(11,*) icycle
      read(11,*) isolv_coarse
      read(11,*) ich_ord
C
      close(11)
      write(*,*) 
      return
      end
C====================================================================
      subroutine soln_file_name(soln_file,kk)
C====================================================================
      implicit none
      character*100 soln_file
      integer k,kk
C--------------------------------------------------------------------
C     Find the file name of corresponding exact solution at level KK.
C--------------------------------------------------------------------
      do k = 1, 100
         soln_file(k:k) = ' '
      end do
C
      if (kk .eq. 4) then
         soln_file = 'soln_cd4'
      else if (kk .eq. 5) then
         soln_file = 'soln_cd5'
      else if (kk .eq. 6) then
         soln_file = 'soln_cd6'
      else if (kk .eq. 7) then
         soln_file = 'soln_cd7'
      else if (kk .eq. 8) then
         soln_file = 'soln_cd8'
      else if (kk .eq. 9) then
         soln_file = 'soln_cd9'
      else if (kk .eq. 10) then
         soln_file = 'soln_cd10'
      else if (kk .eq. 11) then
         soln_file = 'soln_cd11'
      end if 
C
      call input_name(soln_file)
      write(*,*) 'Computed exact solution file: ', soln_file
      write(*,*) 
C
      return
      end
C====================================================================
      subroutine input_name(name1)
C====================================================================
      character*100 name1, name2
C--------------------------------------------------------------------
C     Find the file name of an input matrix. 
C--------------------------------------------------------------------
      name2 = name1
      l1 = len('matrices/')
      l2 = index(name1,' ')
      name1(l1+1:l1+l2)=name2(1:l2-1)
      name1(1:l1)='matrices/'
      return
      end
C====================================================================
