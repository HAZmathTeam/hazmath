C=====================================================================
      program whatev
C=====================================================================
cc      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
cc      parameter (nlast = 16 000 000, nrlast = 8 000 000) !level = 9
      parameter (nlast = 50 000 000, nrlast = 50 000 000) !level = 10
C
      implicit real*8(a-h,o-z)
cc      common /top/ m(nlast)
cc      common /top1/ r(nrlast)
      dimension m(nlast), r(nrlast)
C
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
      character*100 meshf(0:20),soln_file
      common /fname/ meshf
C
      integer lmax
      parameter (lmax = 20)
C
      integer nd, n1, n2, nz,
     >     ia, ja, ka,
     >     ipp, jpp, kpp, nzp, 
     >     ku, kb, iord, 
     >     ifree, kfree
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax),
     >     ia(lmax), ja(lmax), ka(lmax), 
     >     ipp(lmax),jpp(lmax), kpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iord(lmax),
     >     ifree, kfree
C
      integer ipresm, ipossm, npresm, npossm, lcycle, lsolv_c,
     >     ich_ord, lf, lc
      common /mg_params/
     >     ipresm, ipossm, npresm, npossm, lcycle, lsolv_c, 
     >     ich_ord, lf, lc
C
      external msolve, matvec
C
      parameter (length = 1025*1025)
      double precision soln
      COMMON /DSLBLK/ SOLN(length)
C
      character*7 mbell
      common /ring/mbell
C
      real timet(2),t_setup,t_iter
      common /time0/t_setup,t_iter
C---------------------------------------------------------------------
C...  This program solves the 2nd order linear PDEs which are possibly 
C...  convection dominated. The following finite element method can be 
C...  applied: Standard Galerkin, EAFE (Edge Average Finite Element),
C...  Streamline diffusion. For the solver of linear systems the 
C...  following method can be used: Gaussian Elimination, Gauss-Seidel,
C...  Multigrid, GMRES, PCGM.
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
      tol = 2d0**(-10)
      tol = tol*tol*tol
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
      lcycle  = ich(31)
      ipresm  = ich(32)
      ipossm  = ich(33)
      npresm  = ich(34)
      npossm  = ich(35)
      lfmg    = ich(36)
      lcmg    = ich(37)
      lsolv_c = ich(38)
      ich_ord = ich(41)
C
      write(*,1000) lmethod,lsolve,lfmg,lcmg,lsolv_c,lformat,ich_ord
 1000 format(2x,' lmethod = ', i7 / 2x,
     >          ' lsolve  = ', i7 / 2x,
     >          ' lfmg    = ', i7 / 2x,
     >          ' lcmg    = ', i7 / 2x,
     >          ' lsolv_c = ', i7 / 2x,
     >          ' lformat = ', i7 / 2x,
     >          ' ich_ord = ', i7)
C
 10   continue
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      la = lcmg
      lc = lcmg
C  WE ALSO HAVE 
      mgend = lc
      call dtime(timet,t_setup)
C
      lf = lfmg
      if (lc .ge. lf) lf = lc + 1
 50   continue 
C
      write(*,*)
      write(*,*), ' LF, LC :  ', lf, lc, mgend
C
      if(lformat .eq. 1) then
         open(mt,file=meshf(lf),status='unknown',form='formatted')
cc         rewind(mt)
         read(mt,*) n1(lf),n2(lf)
         write(*,*) 'lf,n1,n2 = ',lf,n1(lf),n2(lf)
         read(mt,*) nel,n,ned
      else
         open(mt,file=meshf(lf),status='unknown',form='unformatted')
c     c        rewind(mt)
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
      iet    = inedg + ned * 2
      jet    = iet + n + 1
      jat    = jet + nel*3
      lasti  = jat + 2*(nel + n + 5) + n
C
      ia(lf)    = iao
      ja(lf)    = jao
C
      kao    = 1        
      ksol   = kao + 2*(nel + n + 5) + n
      krhs   = ksol + n
      kx     = krhs + n
      ky     = kx + n
      kuold = ky + n
      kaot   = kuold + n
      lastr  = kaot + 2*(nel + n + 5) + n
C
      call rdmesh(lformat,mt,m(ie),m(je),ned,m(inedg),
     >     m(idiro),r(kx),r(ky),nel,n,jdf)
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
      call pde_choice(lmethod,
     I     nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     O     m(iao),m(jao),r(kao),nnz,r(krhs),
     W     ned,m(inedg),m(idiro),m(iip))
C
      call dtime(timet,t_setup)
C
      ifree = inedg
      kfree = kaot
      maxit = ich(39)
C
      iwk   = iet
      kwk   = kaot
      nx = n1(lf)
      ny = n2(lf)
      h=1d0/dble(nx-1)
      call init_guess(m(iao),m(jao),r(kao),r(krhs),r(ksol),n,nx,ny,h)
cc      if (n .le. 257*257) then
         call output(18,r(kx),r(ky),r(ksol),n)
cc      end if
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(n .ne. nx*ny) stop 10
      lsub = (2**(mgend)+1)*(2**(mgend)+1)
      nnz0 = m(iao+n)-1
      ifree = jao+nnz0
      write(*,*)
      write(*,*), ' LF, LC :  ', lf, lc, mgend, lsub
c      write(*,*) ' nnz0 ', nnz0
c      read(*,*)
      iblock = iwk
      jblock = iwk + lsub + 1
      call blockit(lf,m(iblock),m(jblock),
     >     mgend,m(iao),m(jao),lsub,nx,ny,n)
C
      isize = m(iblock + lsub) - 1
      write(*,*) ' BLOCK SIZE : ', isize
      mask = jblock + isize
      ifree = mask  + n
      kfree = kwk
cc      do k = 1, n
cc         r(krhs+k-1) = 0d0
cc      end do
      maxsit=6
      mnp = 7
      do ijkl = 1 , maxsit
         call copyv(r(ksol),r(kuold),n)
         call subd(n,m(mask),m(iblock),m(jblock),lsub,
     >        m(iao),m(jao),r(kao),r(ksol),
     >        r(krhs),m(ifree),r(kfree),ijkl,mnp)
         call uuminv(r(kuold),r(ksol),n)
         smin = r(kuold)
         smax = smin
         do i = 2 , n
            smin = dmin1(smin,r(kuold+i-1))
            smax = dmax1(smax,r(kuold+i-1))
         end do
C     
         if (n .le. 257*257) then
            call output(mod(ijkl,6)+20,r(kx),r(ky),r(ksol),n)
         end if
         write(*,'(a,e16.8,a,e16.8,a,e16.8)') 
     >        ' MIN/MAX of the solution: ', 
     >        smin,';',smax,';', tol
         if(dabs(smin) .le. tol 
     >        .and. dabs(smax) .le. tol) go to 1234
      end do
 1234 continue
      stop
C
      igs   = 3
      iwr   = 1
C     
      call gauss_seidel(m(iao),m(jao),m(iwk),r(ksol),r(krhs),r(kao),
     >     r(kwk),iexact,n,nnz,igs,maxit,ich_ord,tol,iwr)
C
       call dtime(timet,t_iter)
C
 900  continue
C
C********************************************************************
C...  To print out the computed (exact) solution SOL.
C********************************************************************
c      if (n .le. 65*65) then
c         call output(r(kx),r(ky),r(ksol),n)         
c      end if
C
C...  Solve by Gauss-Seidel iteration. 
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
 950  continue
C
      stop
      end
