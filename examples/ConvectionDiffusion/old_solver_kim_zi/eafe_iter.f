C=====================================================================
      subroutine eafe_iter(iao,jao,idiro,m,sol,rhs,ra,r,n,nnz) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension iao(1),jao(1),idiro(1),sol(1),rhs(1),ra(1)
      dimension m(1),r(1)
      character*100 meshf(0:20)
      common /menu/ ich(200)
      common /fname/ meshf
C---------------------------------------------------------------------
C...  Iterative method using ordered Gauss-Seidel method on EAFE (Edge 
C...  Average Finite Element) as a preconditioner.
C...
C...  Parameter:
C...    MAXIT  - maximum number of MG iterations
C...    TOL    - stopping iteration procedure.
C...    IFLAG  - useful for mixed MG; it is 0 in main program, but it
C...             becomes 1 in the subroutine MG.
C---------------------------------------------------------------------
      mt         = 16
      lread      = ich(5)
      iprecn     = ich(46)
      maxit      = ich(47)
      max_sweeps = ich(48)
      norder     = ich(50)
C
      tol    = 0.5d-6
      omega  = 0.1d00
      zero   = 0.d-13
C
      go to (3,4,6), ich(3)
      write(*,*) ' Error in the value of ICH(3).'
      stop
 3    jdf = 3           !Linear finite elements: degrees of freedom = 3.
      go to 10
 4    jdf = 4
      go to 10
 6    jdf = 6
 10   continue
C
C...  Generating the stiffness matrix A and the load vector B by EAFE.
C
      if(lread .eq. 1) then
         open(mt,file = meshf(0),status='old',form='formatted')
         rewind(mt)
         read(mt,*) n1,n2
         read(mt,*) nel,n,ned
      else
         open(mt,file = meshf(0),status='old',form='unformatted')
         rewind(mt)
         read(mt) n1,n2
         read(mt) nel,n,ned
      end if
      write(*,*)
     >     ' No. of elements, nodes, & bdry edges:',nel,n,ned
C
C...  Compute the addresses of each array in the arrays m and r.
C
      ia      = 1
      ja      = ia + n + 1
      idir    = ja + 2*(nel + n + 5) + n
      lasti   = idir + n + 5
      ie      = lasti
      je      = ie + nel + 1
      inedg   = je + nel*3
      iet     = inedg + ned * 2
      jet     = iet + n + 1
      jat     = jet + nel*3
      iend    = jat + 2*(nel + n + 5) + n
C
      ku      = 1
      kb      = ku + n
      ka      = kb + n
      lastr   = ka + 2*(nel + n + 5) + n
      kx      = lastr
      ky      = kx + n
      krat    = ky + n
      kend    = krat + 2*(nel + n + 5) + n
C     
      call rdmesh(lread,mt,m(ie),m(je),ned,m(inedg),
     >     m(idir),r(kx),r(ky),nel,n,jdf)
C
      call iit(m(ie),m(je),nel,n,m(iet),m(jet))
C
      nnz = 0
      call smbasg(m(ie),m(je),m(iet),m(jet),m(ia),m(ja),
     >     n,nnz,m(idir))
C
      iflag = 2
      iip   = iet
      call pde_choice(iflag,
     I     nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     O     m(ia),m(ja),r(ka),nnz,r(kb),
     W     ned,m(inedg),m(idir),m(iip))
      close(mt)
C
C...  To determine whether Tarjan's ordering is considered or not.
C...  If ordering is considered, then set NORDER = 1.
C
      if (norder .eq. 1) then
         iord   = idir
         iaw    = iord + n + 1
         jaw    = iaw + n + 1 + 1
         isub   = jaw + nnz 
         iwork1 = isub + n + 1
         iwork2 = iwork1 + n + 1
         iwork3 = iwork2 + n + 1
         lasti  = iwork3 + n + 1 + 1
C
         call cut_off(
     >        n,m(ia),m(ja),r(ka),m(iaw),m(jaw),m(iord),
     >        m(isub),nblk,m(iwork1),m(iwork2),m(iwork3))
C
      else
         iord   = idir
         call iseqv(m(iord),n)
      end if
C
      call init_guess(sol,rhs,idiro,n)
      call copyv(rhs,r(kb),n)
      call zero_dir(r(kb),idiro,n)
C
      call scpro(r(kb),r(kb),err_res0,n)
      err_res0 = dsqrt(err_res0)
      if (err_res0 .lt. zero)  err_res0 = 1.d0
C
      write(*,*)
      write(*,*) '=============================', 
     >           '=================================================='
      write(*,*) '     METHOD :   EAFE preconditioner '
C
      write(*,*) '     Residual :   r := b - Au; '
      write(*,'(a,a,e12.4)')
     >     '      Convergence :     ||r||_0/||r_0||_0  &',
     >     '  ||Br||_0/||r_0||_0 <', tol
      write(*,*) '-----------------------------',
     >           '--------------------------------------------------'
      write(*,*) 'Iteration  ||r||/||r_0||    ||r||  ',
     >                  '   ||Br||/||r_0|| ||Br||/||u||     ||Br||'
      write(*,*) '-----------------------------',
     >           '--------------------------------------------------'
C
      kk = 0
      write(*,'(2X,i5,2X,5(3X,e11.4))') 
     >     kk,1.d0,err_res0,1.d0,1.d0,err_res0
C
C...  Iterate until convergence:  sol = sol + B*b = sol + B(rhs - A*sol).
C
      do 200 kk = 1, maxit
         call nullv(r(ku),n)
C
C...     Iterator B in the iteration u = B*b = B(rhs - A*sol).
C     
C...     Doing Jacobi iteration for the original system.
C
         ijacob = 0 
         if (ijacob .eq. 1) then
            call jacobi_wr(r(ku),iao,jao,ra,r(kb),n,2,0,r(kx),omega)
         end if
C
         iw = 1                  !If it is 1, write a convergence history.
         if (iprecn .eq. 1) then
            call fwd_gs_ord(r(ku),m(ia),m(ja),
     >           r(ka),r(kb),m(iord),n,max_sweeps,iw)
         else if (iprecn .eq. 2) then
            call bwd_gs_ord(r(ku),m(ia),m(ja),
     >           r(ka),r(kb),m(iord),n,max_sweeps,iw)
         else if (iprecn .eq. 3) then
            call sym_gs_ord(r(ku),m(ia),m(ja),
     >           r(ka),r(kb),m(iord),n,max_sweeps,iw)
         else if (iprecn .eq. 4) then
            call jacobi(r(ku),m(ia),m(ja),
     >           r(ka),r(kb),n,max_sweeps,iw,r(kx),omega)
         end if
C
C...     New solution:  sol = sol + u = sol + B(rhs - A*sol).
         call uupluv(sol,r(ku),n)
C
C...     To compute the relative error.
C
	 call scpro(r(ku),r(ku),xdiffn,n)
         call scpro(sol,sol,xnewn,n)
         xdiffn  = dsqrt(xdiffn)
         xnewn   = dsqrt(xnewn)
         if (xnewn .gt. zero) err_rel = xdiffn/xnewn
         err_brr = xdiffn/err_res0
C
C...     To compute the residual: b = rhs - A*sol.
C
         call abyvam(rhs,iao,jao,ra,sol,n,r(kb))
C
C...     To compute the residual error.
C
         call scpro(r(kb),r(kb),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
         write(*,'(2X,i5,2X,5(3X,e11.4))')
     >        kk,err_relres,err_res,err_brr,err_rel,xdiffn
         if (err_relres .lt. tol .and. err_brr .lt. tol) go to 333
 200  continue
C
      write(*,*), ' Iteration limit exceeded:', maxit
 333  continue
C
      write(*,*) '=============================', 
     >           '=================================================='
C
      arfac  = (err_relres)**(1./dble(kk))
      arfac1 = (xdiffn/err_res0)**(1./dble(kk))
      write(*,*) 
      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
C
      niter = min(kk, maxit)
      write(*,*)
C
 500  return
      end
C=====================================================================
