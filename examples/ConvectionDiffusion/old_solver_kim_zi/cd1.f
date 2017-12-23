C=====================================================================
      program CDFEM
C=====================================================================
cc      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
cc      parameter (nlast = 16 000 000, nrlast = 8 000 000) !level = 9
      parameter (nlast = 70 000 000, nrlast = 100 000 000) !level = 10
C
      implicit real*8(a-h,o-z)
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      character*100 meshf(0:20),soln_file
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
      common /fname/ meshf
C
      integer lmax
      parameter (lmax = 20)
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
C
      external msolve, matvec
      common /prenormal/ ipre_normal
C
      parameter (length = 1025*1025)
      double precision soln
      COMMON /DSLBLK/ SOLN(length)
C
      character*7 mbell
      common /ring/mbell
C
      real dtime,timet(2),t_setup,t_iter
      external dtime
      common /time0/t_setup,t_iter
C---------------------------------------------------------------------
C...  This program solves the 2nd order linear PDEs which are possibly 
C...  convection dominated. The following finite element method can be 
C...  applied: Standard Galerkin, EAFE (Edge Average Finite Element),
C...  Streamline diffusion. For the solver of linear systems the 
C...  following method can be used: Gaussian Elimination, Gauss-Seidel,
C...  Multigrid.
C...
C...  Parameter:
C...    istiff - choice of stiffness matrices.
C...             1 = Standard Galerkin, 2 = EAFE, 
C...             3 = Streamline diffusion, 4 = Upwind, 5 = SUPG
C...    IEXACT - to determine whether reading or writing the computed
C...             exact solution is considered or not;
C...             0 - no reading and writing the computed exact soln;
C...             1 - writing the computed soln;
C...             2 - reading the computed exact soln;
C...             3 - writing the computed exact soln;
C...    IWR_A_RHS - to determine whether print out the matrix A and 
C...                the load vector RHS.
C...             0 - do not print them
C...             1 - print them
C---------------------------------------------------------------------
      iwr_a_rhs = 0 !1
      iexact    = 0 !3 !2
C
      tol = 1.d-6
      if (iexact .eq. 3) tol = 0.2d-11
cc      if (iexact .eq. 3) tol = 0.5d-12
C
      call input_data
      call quad_elt
      call quad_data
      call quad_data1
C 
      mt      = 16
      ms      = 17
C
      mbell(1:1) = char(7)
      mbell(2:7) = 'DONE.'
C
      jdf     = ich(3)
      lmethod = ich(4)
      lformat = ich(5)
      lsolve  = ich(7)
      lfmg    = ich(36)
      lcmg    = ich(37)
      norder  = ich(41)
C
cc      write(*,1000) lformat,ich(1),lmethod,ich(23),ich(11),lsolve
 1000 format(2x,' lformat = ', i7 / 2x,
     >          ' lpde    = ', i7 / 2x,
     >          ' lmethod = ', i7 / 2x,
     >          ' lconv_f = ', i7 / 2x,
     >          ' lquad   = ', i7 / 2x,
     >          ' lsolve  = ', i7)
C
 10   continue
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      la = lcmg
      do 980 ll = 1, 1, 2       !1,5,2
         lc = ll
         if(iexact .eq. 3) ll0 = 10
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         t_setup = dtime(timet)
C
         lf   = lfmg
         if (lc .ge. lf) lf = lc + 1
 50      continue 
C
      write(*,*)
      write(*,*), ' LF, LC :  ', lf, lc
C
      if(lformat .eq. 1) then
         open(mt,file=meshf(lf),status='unknown',form='formatted')
cc         rewind(mt)
         read(mt,*) n1(lf),n2(lf)
         read(mt,*) nel,n,ned
      else
         open(mt,file=meshf(lf),status='unknown',form='unformatted')
cc         rewind(mt)
         read(mt) n1(lf),n2(lf)
         read(mt) nel,n,ned
      end if
      nd(lf) = n 
cc      write(*,*) ' No. of elements, nodes, & bdry edges:', nel,n,ned
C
C...  Compute the addresses of each array in the arrays m and r.
C
      ie     = 1
      je     = ie + nel + 1
      iao    = je + nel*3
      jao    = iao + n + 1
      idiro  = jao + 2*(nel + n + 5) + n
      inedg  = idiro + n + 5
      iperm(lf) = inedg + ned * 2
      iet    = iperm(lf) + n + 1
      jet    = iet + n + 1
      jat    = jet + nel*3
      lasti  = jat + 2*(nel + n + 5) + n
C
      ia(lf)    = iao
      ja(lf)    = jao
      idir(lf)  = idiro
C
      kal    = 1        
      ksol   = kal + 2*(nel + n + 5) + n
      krhs   = ksol + n
      kra    = krhs + n
      kx     = kra + 2*(nel + n + 5) + n
      ky     = kx + n
      krat   = ky + n
      lastr  = krat + 2*(nel + n + 5) + n
C
      kag(lf) = kra
C
      call rdmesh_perm(lformat,mt,m(ie),m(je),ned,m(inedg),
     >     m(idiro),m(iperm(lf)),r(kx),r(ky),nel,n,jdf)
C
      call iit(m(ie),m(je),nel,n,m(iet),m(jet))
C
      nnz = 0
      call smbasg(m(ie),m(je),m(iet),m(jet),m(iao),m(jao),
     >     n,nnz,m(idiro))
      nz(lf) = nnz
C
      iip = iet
C
      if (lmethod .le. 5) then
         istiff = lmethod
      else if (lmethod .eq. 6) then
         istiff = 3
      else
         write(*,*) ' Check the value of ICH(4).'
         stop
      end if
C
      call pde_choice(istiff,
     I     nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     O     m(iao),m(jao),r(kra),nnz,r(krhs),
     W     ned,m(inedg),m(idiro),m(iip))
C
C...  For power method, we need the Laplace matrix, A_0.
C
      call lapl_mass(1,nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     >     m(iao),m(jao),r(kal),nnz,m(idiro),m(iip))
C
cc      write(*,*) (r(krhs+ii-1),ii=1,n)
cc      write(*,*) (m(iao+ii-1),ii=1,n+1)
cc      write(*,*) (ij,'*',m(iord+ij-1),ij=1,n)
cc      call lprr(m(iao),m(jao),r(kal),n)
C
C********************************************************************
      if (iwr_a_rhs .eq. 1) then
C...     To print out the matrix A and the load vector RHS.
         call write_matrix_rhs(lf,n,m(iao),m(jao),r(kra),r(krhs),
     >        nnz,lformat)
         go to 950
      end if
C********************************************************************
C
cc      call outmat(m(iao),m(jao),n,200,r(kra),r(krhs))
cc      call outmat1(m(iao),m(jao),r(kra),n,nnz)

cc      if (lmethod .eq. 6) then
cc         call pde_choice(4,
cc     I        nel,n,jdf,m(ie),m(je),r(kx),r(ky),
cc     O        m(iao),m(jao),r(kra),nnz,r(krhs),
cc     W        ned,m(inedg),m(idiro),m(iip))
cc      end if
C
      if (iexact .eq. 2) then
         call soln_file_name(soln_file,lf)
         if (lformat .eq. 1) then
            open(ms,file=soln_file,status='unknown',form='formatted')
            read(ms,*) (soln(i),i=1,n)
         else 
            open(ms,file=soln_file,status='unknown',form='unformatted')
            read(ms) (soln(i),i=1,n)
         end if
      end if
C
      write(*,'(a)') mbell
C
      t_setup = dtime(timet)
C
cc      write(*,*) (r(krhs+ii), ii = 0, n-1)
C
      go to (100,200,300,400,500,600,700,800), lsolve
      write(*,*) ' Error in the choice of the linear system solver.'
      go to 990
C
 100  continue
C
C...  Compute the H^1 norm of the operator T = I - BA, i.e., ||T||_1.
C...  Here A is the EAFE scheme and B is the mg V-cycle on EAFE scheme.
C...  For this purpose, we use a power method.
C
      ifree = inedg
      kfree = krat
      maxit = 500
C
cc      call mg_s(m,iao,jao,idiro,ifree,r,ksol,krhs,kal,kfree,
cc     >     n,n1(lf),n2(lf),nnz,nel,lf,lc,tol) 
C
C...  It is needed to have zero rhs:
      call nullv(r(krhs),n)
C
      call power_new(m(1),m(iao),m(jao),m(idiro),
     >     r(1),r(ksol),r(krhs),r(kra),
     >     n,nnz,maxit,tol,iexact,
     >     nel,jdf,m(ie),m(je),r(kx),r(ky),
     >     iao,jao,idiro,kal,ksol,krhs,kra)
C
cc      call power_modify(m(1),m(iao),m(jao),m(idiro),
cc     >     r(1),r(ksol),r(krhs),r(kra),
cc     >     n,nnz,maxit,tol,iexact,
cc     >     nel,jdf,m(ie),m(je),r(kx),r(ky),
cc     >     iao,jao,idiro,kal,ksol,krhs,kra)
C
      t_iter = dtime(timet)
C
      go to 900
C
 200  continue
C
C...  Solve by Gaussian Elimination. (!!!!!!Check initial guess!!!!!!)
C
C...  Transposed system.
C
      iat = iet
      if(ich(1) .eq. 3 .and. ich(4) .eq. 2 .and. ich(6) .eq. -1) then
         write(*,*) ' transposed system::: '
         call aat(m(iao),m(jao),r(kra),n,n,m(iat),m(jat),r(krat))
         call icopyv(m(iat),m(iao),n+1)
         call icopyv(m(jat),m(jao),nnz)
         call copyv(r(krat),r(kra),nnz)
      endif
C
C...  Original system.
C
      kau = krat
      call copyv(r(krhs),r(ksol),n)
      call sgauss(m(iao),m(jao),r(kra),r(ksol),m(idiro),r(kau),n)
C
      t_iter = dtime(timet)
C
      go to 900
C
 300  continue
C
C...  Solve by Gauss-Seidel iteration. (!!!!!!Check initial guess!!!!!!)
C
      iwk = iet
      kwk = krat
      call gausei(m(iao),m(jao),m(idiro),m(iwk),r(ksol),r(krhs),r(kra),
     >     n,nnz)
C
      t_iter = dtime(timet)
C
      go to 900
C
 400  continue
C
C...  Solve by preconditioning with EAFE. (!!!!!!Check initial guess!!!!!!)
C
      ipre = inedg
      kpre = krat
C
      call eafe_iter(m(iao),m(jao),m(idiro),m(ipre),
     >     r(ksol),r(krhs),r(kra),r(kpre),n,nnz)
C
      t_iter = dtime(timet)
C
      go to 900
C
 500  continue
C
      ipre_normal = 0 !NO smoothings by G-S sweeps for the normal eq.
      go to 650
C
 600  continue
C     
      ipre_normal = 1 !DO smoothings by G-S sweeps for the normal eq.
      go to 650
C
 650  continue
C
C...  Solve by MG or MG_normal.
C
      ifree = inedg
      kfree = krat
      maxit = ich(39)
      mg_format = 1
C
C========================================================
cc      kasave = kfree
cc      kfree = kasave + nnz
ccC
cc      iasave = ifree
cc      jasave = iasave + n + 1
cc      ifree = jasave + nnz
ccC     
ccC...  Temporary save for iao, jao, ao
ccC
cc      call icopyv(m(iao),m(iasave),n+1)
cc      call icopyv(m(jao),m(jasave),nnz)
cc      call copyv(r(kra),r(kasave),nnz)
ccC
ccC...  Compute at = A^t.
cc      call aat(m(iasave),m(jasave),r(kasave),n,n,
cc     >     m(iao),m(jao),r(kra))
ccC
ccC...  Initial guess.
cc      call init_g0(m(iao),m(jao),r(kra),r(krhs),r(ksol),n)
ccC
cc      call mg_tr(m(1),m(iao),m(jao),m(idiro),
cc     >     r(1),r(ksol),r(krhs),r(kra),
cc     >     n,nnz,maxit,tol,iexact)
C=======================================================
C
C...  Initial guess.
      call init_g0(m(iao),m(jao),r(kra),r(krhs),r(ksol),n)
C
C...  The following two subroutines are identical.
C
      if (mg_format .eq. 0) then
         call mg(m(1),m(iao),m(jao),m(idiro),
     >        r(1),r(ksol),r(krhs),r(kra),
     >        n,nnz,maxit,tol,iexact)
      else
         call premg_normal(m(1),m(iao),m(jao),m(idiro),
     >        r(1),r(ksol),r(krhs),r(kra),
     >        n,nnz,maxit,tol,iexact,msolve)
      end if
C
      t_iter = dtime(timet)
C
      go to 900
C
 700  continue
C
      ipre_normal = 0 !NO smoothings by G-S sweeps for the normal eq.
      go to 850
C
 800  continue
C     
      ipre_normal = 1 !DO smoothings by G-S sweeps for the normal eq.
      go to 850
C
 850  continue
C
C...  Initial guess.
      call init_g0(m(iao),m(jao),r(kra),r(krhs),r(ksol),n)
cc      write(*,*) (r(ksol+ii), ii = 0, n-1)
C
C...  Solve by GMRES using the preconditioner MSOLVE.
C
      ifree = inedg
      kfree = krat
      maxit = ich(39)
C
C...  Get the information about preconditioning matrices.
      call precond_mat(m,m(iao),m(jao),m(idiro),r,r(kra),n,nnz) 
C
      isym  = 0
      itol  = 0 
      if (iexact .eq. 2)   itol = 11
      itmax = maxit
      iunit = 6
C
      ligw  = 100
      maxl  = 20  !80 !100
      kmp   = maxl
      jscal = 0
      jpre  = -1    !0 : no precond, + : right precond, - : left precond.
      nrmax = 10
C
      lrgw  = 2 + N*(MAXL+6) + MAXL*(MAXL+3)
C
      lasti = ifree
      lastr = kfree
C
      igwk  = lasti
      lasti = igwk + 100
      
      ksb   = lastr
      ksx   = ksb + 2
      krgwk = ksx + 2
      lastr = krgwk + lrgw
C
      ifree = lasti
      kfree = lastr
C
      m(igwk)     = maxl
      m(igwk + 1) = kmp
      m(igwk + 2) = jscal
      m(igwk + 3) = jpre
      m(igwk + 4) = nrmax
C
      call DGMRES(n, r(krhs), r(ksol), nnz,m(iao),m(jao),r(kra), 
     +     ISYM, MATVEC, MSOLVE,
     +     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, r(ksb), r(ksx), 
     +     r(krgwk), LRGW, m(IGWK), LIGW, r(1), m(1))
C
      write(*,*)
      write(*,*) '==================================================='
      if (ipre_normal .eq. 0) then
         write(*,*) ' METHOD : GMRES with MG precondition'
      else 
         write(*,*) ' METHOD : GMRES with MG + Normal precondition'
      end if
C
      write(*,*) '   Residual :        r := b-Au; '
      write(*,'(a,a,e12.4)')
     >     '    Convergence :     ||r||_0/||b||_0  and/or',
     >     '  ||u - sol||_0/||sol||_0 <', tol
      write(*,'(a,2(5X,i3))') '    Coarse and fine levels : ',lc,lf
      write(*,*)
      write(*,*) ' No. of iterations and error :', iter, err
      write(*,*) ' Return error flag : ', ierr
      write(*,*) ' Required minimum length of RGWK array : ',m(igwk+5)
      write(*,*) ' The total number of calls to MSOLVE : ', m(igwk+6)
      write(*,*) '==================================================='
      write(*,*)
C
      t_iter = dtime(timet)
C
      go to 900
C
 900  continue
C
C********************************************************************
C...  To print out the computed (exact) solution SOL.
C********************************************************************
      if (iexact .eq. 1 .and. n .le. 65*65) then
         call output(r(kx),r(ky),r(ksol),n)         
      else if (iexact .eq. 3) then
         call write_soln(lf,r(ksol),n,lformat)
      end if
C
      if (n .le. 65*65) then
cc         call output(r(kx),r(ky),r(ksol),n)
cc         call output_sol(m(ie),m(je),r(kx),r(ky),r(ksol),n,nel)
cc         call wrsol01(r(ksol),n,100)
         endfile(800)
         close(800)
         endfile(100)
         close(100)
         endfile(199)
         close(199)
      end if 
C
      smin = r(ksol)
      smax = smin
cc      if(smax .gt. 1d-05) write(*,*) 1
cc      write(100,'(e15.7$)') (r(ksol))
      do i = 2 , n
cc         if( i.le. 65*65) write(100,'(e15.7$)') (r(ksol+i-1))
         smin = dmin1(smin,r(ksol+i-1))
         smax = dmax1(smax,r(ksol+i-1))
cc         if(dabs(r(ksol+i-1)) .gt. 1d-05) then
cc             write(*,*) i,r(ksol+i-1)
cc             read(*,*)
cc         end if
      end do
C
cc      write(*,*) 
cc     >'=============================================================',
cc     >'================='
cc      write(*,*)
cc      write(*,'(a,f20.12)') ' MINIMUM of the solution: ', smin
cc      write(*,'(a,f20.12)') ' MAXIMUM of the solution: ', smax
cc      write(*,*)
cc      write(*,*) 
cc     >'=============================================================',
cc     >'================='

cc      write(*,'(10x,a)') ' CPU time:  '
cc      write(*,'(10x,a)') ' ======================== '
cc      write(*,'(10x,a,f10.3,a)') '   SETUP: ',t_setup, ' s '
cc      write(*,'(10x,a)') ' ------------------------'
cc      write(*,'(10x,a,f10.3,a)') '   SOLVE: ',t_iter, ' s '
cc      write(*,'(10x,a)') ' ------------------------'
cc      write(*,'(10x,a,f10.3,a)') '   TOTAL: ',t_setup+t_iter, ' s '
cc      write(*,'(10x,a)') ' ======================== '
C
cc      if (lsolve .le. 3) go to 990
C
 950  continue
C
cc      if ((lsolve .ge. 4) .and. lf .lt. 10) then
cc      if ((lsolve .ge. 4) .and. lf .lt. 9) then
      if (lf .lt. 9) then
         lf = lf + 1
         close(ms)
         close(mt)
cc         read(*,*)
         go to 50
      end if
      close(ms)
      close(mt)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 980  continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 990  close(mt)
      close(ms)
C
      stop
      end
C=====================================================================

