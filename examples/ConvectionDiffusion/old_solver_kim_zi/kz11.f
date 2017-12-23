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
      real dtime,timet(2),t_setup,t_iter
      external dtime
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
      write(*,1000) lmethod,lsolve,lfmg,lcmg,lsolv_c,ich_ord
 1000 format(2x,' lmethod = ', i7 / 2x,
     >          ' lsolve  = ', i7 / 2x,
     >          ' lfmg    = ', i7 / 2x,
     >          ' lcmg    = ', i7 / 2x,
     >          ' lsolv_c = ', i7 / 2x,
     >          ' ich_ord = ', i7)
C
 10   continue
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      la = lcmg
      lc = lcmg
      t_setup = dtime(timet)
C
      lf = lfmg
      if (lc .ge. lf) lf = lc + 1
 50   continue 
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
      iet    = inedg + ned * 2
      jet    = iet + n + 1
      jat    = jet + nel*3
      lasti  = jat + 2*(nel + n + 5) + n
C
      ia(lf)    = iao
      ja(lf)    = jao
C
      kal    = 1        
      ksol   = kal + 2*(nel + n + 5) + n
      krhs   = ksol + n
      kao    = krhs + n
      kx     = kao + 2*(nel + n + 5) + n
      ky     = kx + n
      kaot   = ky + n
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
      t_setup = dtime(timet)
C
      ifree = inedg
      kfree = kaot
      maxit = ich(39)
C
C...  Solve by Gauss-Seidel iteration. 
C
      iwk   = iet
      kwk   = kaot
      igs   = 3
      iwr   = 1
C     
      call init_guess(m(iao),m(jao),r(kal),r(krhs),r(ksol),n)
C      if (n .le. 65*65) then
C         call output(r(kx),r(ky),r(ksol),n)        
C         stop 
C      end if
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call gauss_seidel(m(iao),m(jao),m(iwk),r(ksol),r(krhs),r(kao),
     >     r(kwk),iexact,n,nnz,igs,maxit,ich_ord,tol,iwr)
C
      t_iter = dtime(timet)
C
 900  continue
C
C********************************************************************
C...  To print out the computed (exact) solution SOL.
C********************************************************************
      if (n .le. 65*65) then
         call output(r(kx),r(ky),r(ksol),n)         
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
C=====================================================================
      subroutine init_guess(ia,ja,a,b,sol,n)
C=====================================================================
      integer ia(1),ja(1),k,n,p1,p2, iseed(4)
      real*8 a(1),sol(1),b(1),
     >     solmax,solmin
      real time(2),t1,t2,t3
C--------------------------------------------------------------------
C... Initial guess:
C
      t1 = etime(time)
      call nullv(sol,n)
C
      do k = 1 , 4
         iseed(k) = 1
      end do
C... Get a uniformly distributed random numbers as initial guess
cccccccccccccccccc      call dlaruv( iseed, n, u )
      call dlarnv( 2, iseed, n, sol )
C
      do k = 1 , n
         p1 = ia(k)
         p2 = ia(k+1)-1
         if(p2-p1 .lt. 0) stop 10
         if(p2-p1 .eq. 0) then
            sol(k) = b(k)
         end if
      end do
 222  continue
      return
      end
