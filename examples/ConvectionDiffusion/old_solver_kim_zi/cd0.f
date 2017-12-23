C=====================================================================
      program ELLIPTIC
C=====================================================================
cc      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
cc      parameter (nlast = 16 000 000, nrlast = 8 000 000) !level = 9
cc      parameter (nlast = 60 000 000, nrlast = 30 000 000) !level = 10
      parameter (nlast = 128 000 000, nrlast = 64 000 000) !level = 11
      implicit real*8(a-h,o-z)
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      character*100 meshf(0:20)
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
      common /fname/ meshf
C---------------------------------------------------------------------
C...  This program solves the 2nd order linear elliptic PDEs. The 
C...  following finite element method can be applied: Standard Galerkin, 
C...  EAFE (Edge Average Finite Element), Streamline diffusion. 
C...  For the solver of linear systems the following method can be used:
C...  Gaussian Elimination, Gauss-Seidel, Multigrid. 
C...  Especially, if MGM is chosen, then we use the standard MG, which 
C...  means that we use A_{k-1} = PA_kP^t to get lower level matrices, 
C...  where P is the standard prolongation operator.
C...
C...  Parameter:
C...    istiff - choice of stiffness matrices in MG.
C...             1 = Standard Galerkin, 2 = EAFE, 
C...             3 = Streamline diffusion, 4 = Upwind, 5 = SUPG
C---------------------------------------------------------------------
      call input_data
      call quad_elt
      call quad_data
      call quad_data1
C 
      mt      = 16
      jdf     = ich(3)
      lmethod = ich(4)
      lread   = ich(5)
      lsolve  = ich(7)
      lfmg    = ich(36)
      lcmg    = ich(37)
      norder  = ich(41)
C
      write(*,1000) lread,ich(1),lmethod,ich(23),ich(11),lsolve
 1000 format(2x,' lread   = ', i7 / 2x,
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
                          do 980 ll = la,1
                          lcmg = ll
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      lf   = lfmg
      if (lcmg .ge. lf) lf = lcmg + 1
 50   continue 
C
      if(lread .eq. 1) then
         open(mt,file=meshf(lf),status='unknown',form='formatted')
cc         rewind(mt)
         read(mt,*) n1,n2
         read(mt,*) nel,n,ned
      else
         open(mt,file=meshf(lf),status='unknown',form='unformatted')
cc         rewind(mt)
         read(mt) n1,n2
         read(mt) nel,n,ned
      end if
      write(*,*) ' No. of elements, nodes, & bdry edges:', nel,n,ned
C
C...  Compute the addresses of each array in the arrays m and r.
C
      ie     = 1
      je     = ie + nel + 1
      ia     = je + nel*3
      ja     = ia + n + 1
      idir   = ja + 2*(nel + n + 5) + n
      inedg  = idir + n + 5
      iet    = inedg + ned * 2
      jet    = iet + n + 1
      jat    = jet + nel*3
      lasti  = jat + 2*(nel + n + 5) + n
C
      ksol   = 1
      krhs   = ksol + n
      kra    = krhs + n
      kx     = kra + 2*(nel + n + 5) + n
      ky     = kx + n
      krat   = ky + n
      lastr  = krat + 2*(nel + n + 5) + n
C
      call rdmesh(lread,mt,m(ie),m(je),ned,m(inedg),
     >     m(idir),r(kx),r(ky),nel,n,jdf)
C
      call iit(m(ie),m(je),nel,n,m(iet),m(jet))
C
      nnz = 0
      call smbasg(m(ie),m(je),m(iet),m(jet),m(ia),m(ja),n,nnz,m(idir))
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
     O     m(ia),m(ja),r(kra),nnz,r(krhs),
     W     ned,m(inedg),m(idir),m(iip))
C
cc      write(*,*) (ii,'*',r(krhs+ii-1),ii=1,n)
cc      write(*,*) (m(ia+ii-1),ii=1,n+1)
cc      write(*,*) (ij,'*',m(iord+ij-1),ij=1,n)
cc      call lprr(m(ia),m(ja),r(kra),n)
C
C********************************************************************
      iwr_a_rhs = 0
      if (iwr_a_rhs .eq. 1) then
C...     To print out the matrix A and the load vector RHS.
         call write_matrix_rhs(lf,n,m(ia),m(ja),r(kra),r(krhs),
     >        nnz,lread)
         go to 990
      end if
C********************************************************************
C
cc      call outmat(m(ia),m(ja),n,200,r(kra),r(krhs))
cc      call outmat1(m(ia),m(ja),r(kra),n,nnz)

cc      if (lmethod .eq. 6) then
cc         call pde_choice(4,
cc     I        nel,n,jdf,m(ie),m(je),r(kx),r(ky),
cc     O        m(ia),m(ja),r(kra),nnz,r(krhs),
cc     W        ned,m(inedg),m(idir),m(iip))
cc      end if
C         
      go to (200,300,400,500), lsolve
      write(*,*) ' Error in the choice of the linear system solver.'
      go to 990
C
 200  continue
C
C...  Solve by Gaussian Elimination.
C
C...  Transposed system.
C
      iat = iet
      if(ich(1) .eq. 3 .and. ich(4) .eq. 2 .and. ich(6) .eq. -1) then
         write(*,*) ' transposed system::: '
         call aat(m(ia),m(ja),r(kra),n,n,m(iat),m(jat),r(krat))
         call icopyv(m(iat),m(ia),n+1)
         call icopyv(m(jat),m(ja),nnz)
         call copyv(r(krat),r(kra),nnz)
      endif
C
C...  Original system.
C
      kau = krat
      call copyv(r(krhs),r(ksol),n)
      call sgauss(m(ia),m(ja),r(kra),r(ksol),m(idir),r(kau),n)
      go to 900
C
 300  continue
C
C...  Solve by Gauss-Seidel iteration.
C
      iwk = iet
      kwk = krat
      call gausei(m(ia),m(ja),m(idir),m(iwk),r(ksol),r(krhs),r(kra),
     >     n,nnz)
      go to 900
C
 400  continue
C
C...  Solve by MG.
C
      ifree = inedg
      kfree = krat
      ich_ord = 0
      tol  = 1.d-6
C
      if (ich_ord .eq. 0) then
         call mg_s(m(1),ia,ja,idir,ifree,r(1),ksol,krhs,kra,kfree,
     >        n,n1,n2,nnz,nel,lf,lcmg,tol) 
      else 
         call mg_s_ord(m(1),ia,ja,idir,ifree,r(1),ksol,krhs,kra,kfree,
     >        n,n1,n2,nnz,nel,lf,lcmg,tol)
      end if
C
      go to 900
C
 500  continue
C
C...  Solve by preconditioning with EAFE. 
C
      ipre = inedg
      kpre = krat
C
      call eafe_iter(m(ia),m(ja),m(idir),m(ipre),
     >     r(ksol),r(krhs),r(kra),r(kpre),n,nnz)
      go to 900
C
 900  continue
C
C...  Output of computed solution..
C
cc      call write_soln(lf,r(ksol),n,lread)

cc      if (n .lt. 8000 .and. lsolve .ne. 3) then
      if (n .le. 129*129) then
         call output(r(kx),r(ky),r(ksol),n)
cc         call output_sol(m(ie),m(je),r(kx),r(ky),r(ksol),n,nel)
         endfile(800)
         close(800)
         endfile(100)
         close(100)
         endfile(199)
         close(199)
      end if 
C
C********************************************************************
C...  To print out the computed solution SOL.
C********************************************************************
      iexact = 0
      if (iexact .eq. 1) then
         lvl = 0
         call write_soln(lvl,r(ksol),n,iformatted)
      end if
C
      if (n .le. 65*65) then
         call wrsol01(r(ksol),n,100)
         endfile(100)
         close(100)
      end if
C
      smin = r(ksol)
      smax = smin
C
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
      write(*,*) 
     >'=============================================================',
     >'================='
      write(*,*)
      write(*,'(a,f20.12)') ' MINIMUM of the solution: ', smin
      write(*,'(a,f20.12)') ' MAXIMUM of the solution: ', smax
      write(*,*)
      write(*,*) 
     >'=============================================================',
     >'================='

      write(*,'(10x,a)') ' CPU time:  '
      write(*,'(10x,a)') ' ======================== '
      write(*,'(10x,a,f10.3,a)') '   SETUP: ',t_setup, ' s '
      write(*,'(10x,a)') ' ------------------------'
      write(*,'(10x,a,f10.3,a)') '   SOLVE: ',t_iter, ' s '
      write(*,'(10x,a)') ' ------------------------'
      write(*,'(10x,a,f10.3,a)') '   TOTAL: ',t_setup+t_iter, ' s '
      write(*,'(10x,a)') ' ======================== '
C
C
      if (lsolve .ne. 3) go to 990
C
      if (lsolve .eq. 3 .and. lf .lt. 6) then
cc      if (lsolve .eq. 3 .and. lf .lt. lfmg) then
         lf = lf + 1
         close(mt)
cc         read(*,*)
         go to 50
      end if
      close(mt)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 980                      continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 990  close(mt)
      stop
      end
C=====================================================================

