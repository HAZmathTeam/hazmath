C=====================================================================
      program Match_MG
C=====================================================================
      implicit none
C---------------------------------------------------------------------
cc      integer*4 nlast,nrlast
cc      integer m	
cc      double precision r	
cc      common /top/ m(nlast)
cc      common /top1/ r(nrlast)
C---------------------------------------------------------------------
      integer nlast,nrlast
cc      parameter (nlast = 140 000 000, nrlast = 100 000 000)  !limit
cc      parameter (nlast = 1 625 000, nrlast = 1 500 000)  !lvl = 8
cc      parameter (nlast = 6 500 000, nrlast = 6 000 000)  !lvl = 9
      parameter (nlast = 26 000 000, nrlast = 24 000 000)  !lvl = 10
cc      parameter (nlast = 105 000 000, nrlast = 96 000 000)  !lvl = 11
      integer m(nlast)
      double precision r(nrlast)
C      
      double precision zero, one
      character*100 matrix_file,soln_file
      character*7 mbell
      common /onezero/ zero,one
      common /namef/matrix_file,soln_file
      common /ring/mbell
C
      include 'paramsmg.h'
C
C...  memory addresses:::
C
      integer max_lvl
      integer nd,iord, 
     >     ia, ja, ka, 
     >     ipp,jpp, kpp,
     >     ku , kb, 
     >     ifree , kfree,lf,lc
C
      parameter (max_lvl = 100)
      common /point/ 
     >     nd(max_lvl),iord(max_lvl), 
     >     ia(max_lvl), ja(max_lvl), ka(max_lvl), 
     >     ipp(max_lvl),jpp(max_lvl), kpp(max_lvl),
     >     ku(max_lvl) , kb(max_lvl), 
     >     ifree , kfree ,lf,lc
C
      real dtime,timet(2),t_setup,t_iter
      external dtime
      common /time0/t_setup,t_iter
C
      integer lsolve,iformatted,maxitr,jj,kk,iexact,kex_sol
      integer iwk,igs,iwr,norder,max_sweeps,ms,lvl
      integer mt,ksol,krhs,ia0,ja0,ka0,n,nnz,i,kx,ky,ko,kc
      double precision smin,smax,tol
C
C--------------------------------------------------------------------
C...  Solve by MG method with matching.
C...  Paramaters:
C...    IEXACT - to determine whether reading or writing the computed
C...             exact solution is considered or not;
C...             0 - no reading and writing the computed exact soln;
C...             1 - writing the computed soln;
C...             2 - reading the computed exact soln;
C...             3 - writing the computed exact soln;
C--------------------------------------------------------------------
      iexact = 2
C
C********************************************************************
C...  Multiple running for different parameters.
cc      do 2000 jj = 7, 7
      do 2000 jj = 1, 10
      do 1000 kk = 4, 10
cc      do 1000 kk = 4, 8
C********************************************************************
C
      t_setup = dtime(timet)
C
C...  Initialize the pointers.
C
      call initial
C
      zero = .5d-13
      one = 1d0
      mt = 16
      ms = 17
      mbell(1:1) = char(7)
      mbell(2:7) = 'DONE.'
C
      lsolve = 0
      iformatted = 0
      maxitr = 0
C
cc      call read_param(lsolve,tol,maxitr,iformatted)
      call read_param_multi(lsolve,tol,maxitr,iformatted,kk,jj)
      if (iexact .eq. 3) tol = 1.d-13
C
 20   continue
C
      if(iformatted .eq. 2) then
         open(mt,file=matrix_file,status='unknown',form='unformatted')
         write(*,'(a$)') ' Reading... '
         read(mt) n
         nd(1) = n
         ia(1) = 1
         ja(1) = ia(1) + nd(1) + 1
         ia0 = ia(1)
         ja0 = ja(1)
C
         read(mt) (m(ia0+i-1),i=1,n+1)
C
         nnz   = m(ia0+n) - 1
C     
         ifree = ja0 + nnz
C     
cc         kx  = 1
cc         ky  = kx + n
C
         if (iexact .eq. 2) then
            kex_sol = 1
            ksol    = kex_sol + n
         else 
            kex_sol = 1
            ksol    = 1
         end if
         krhs  = ksol + n
         ka(1) = krhs + n
         ka0   = ka(1)
         kfree = ka0 + nnz
C     
         call read_u(mt,m(ja0),r(ka0),r(krhs),n,nnz)
         write(*,'(a)') mbell
C
         if (iexact .eq. 2) then
            open(ms,file=soln_file,status='unknown',form='unformatted')
            read(ms) (r(kex_sol+i-1),i=1,n)
         end if
      else
         open(mt,file=matrix_file,status='unknown',form='formatted')
         write(*,'(a$)') ' reading... '
         read(mt,*) n
         nd(1) = n
         ia(1) = 1
         ja(1) = ia(1) + nd(1) + 1
         ia0 = ia(1)
         ja0 = ja(1)
C
         read(mt,*) (m(ia0+i-1),i=1,n+1)
C
         nnz   = m(ia0+n) - 1
C     
         ifree = ja0 + nnz
C     
         if (iexact .eq. 2) then
            kex_sol = 1
            ksol    = kex_sol + n
         else 
            kex_sol = 1
            ksol    = 1
         end if
         krhs  = ksol + n
         ka(1) = krhs + n
         ka0   = ka(1)
         kfree = ka0 + nnz
C     
         call read_f(mt,m(ja0),r(ka0),r(krhs),n,nnz)
         write(*,'(a)') mbell
C
         if (iexact .eq. 2) then
            open(ms,file=soln_file,status='unknown',form='formatted')
            read(ms,*) (r(kex_sol+i-1),i=1,n)
         end if
      end if
C
      write(*,*) ' No. of nodes & no. of nonzero entries:',n,nnz
C
C...  Checking the matrix.
cc      ko = kfree
cc      kc = ko + n + 1
cc      call onesv(r(ko),n)
cc      call abyvg(m(ia0),m(ja0),r(ka0),r(ko),n,r(kc))
cc      write(*,*)  'Row sum:'
cc      write(*,*) (i,'*',r(kc+i-1),i=1,n)
cc      call vbya(m(ia0),m(ja0),r(ka0),r(ko),n,n,r(kc))
cc      write(*,*)  'Column sum:'
cc      write(*,*) (i,'*',r(kc+i-1),i=1,n)

cc      write(*,*) (r(krhs+i-1),i=1,n)
cc      write(*,*) (m(ia0+i-1),i=1,n+1)
cc      call lprr(m(ia0),m(ja0),r(ka0),n)
C
C...  INITIAL guess...
      lf = 1
      lc = 1
C
      t_setup = dtime(timet)
C
      call init_g0(m(ia0),m(ja0),r(ka0),r(krhs),r(ksol),n)
C
cc      write(*,*) (i,'*',r(ksol+i-1),i=1,n)
C
      go to (100,200,300) lsolve
C
 100  continue
C
C...  Solve by MG.
C
C...  SETUP: prolongations and coarse grid matrices definition...
C
      ku(1) = kfree
      kb(1) = ku(1) + n
      kfree = kb(1) + n
C
      call mg_match(r(ksol),r(krhs),n,maxitr,m(1),r(1),
     >     tol,r(kex_sol),iexact)
C
      t_iter = dtime(timet)
C
      go to 990
C
 200  continue
C
C...  Compute the spectral radius of (I-BA) for the mg_match.
C
C...  SETUP: prolongations and coarse grid matrices definition...
C
      ku(1) = kfree
      kb(1) = ku(1) + n
      kfree = kb(1) + n
C
cc      write(*,*) ifree, kfree
      call mg_match_power(r(ksol),r(krhs),n,maxitr,m(1),r(1),
     >     tol)
C
      go to 999
C
 300  continue
C
C...  Solve by Gauss-Siedel iteration.
C
      iwk = ifree
C
      igs = 1
      iwr = 1
      norder = ich_ord
      max_sweeps = 3*n
cc      tol = 1.d-9
C
      call gauss_seidel(m(ia0),m(ja0),m(iwk),
     >     r(kex_sol),r(ksol),r(krhs),r(ka0),iexact,
     >     n,nnz,igs,max_sweeps,norder,tol,iwr)
C
      t_iter = dtime(timet)
C
      go to 990
C
 990  continue
C
C********************************************************************
C...  To print out the computed (exact) solution SOL.
C********************************************************************
      if (iexact .eq. 1 .or. iexact .eq. 3) then
         lvl = kk
         call write_soln(lvl,r(ksol),n,iformatted)
      end if
C
      if (n .le. 65*65) then
cc         call output(r(kx),r(ky),r(ksol),n)
cc         call output_sol(m(ie),m(je),r(kx),r(ky),r(ksol),n,nel)
         call wrsol01(r(ksol),n,100)
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
 999  continue
C     
C********************************************************************
 1000 continue
 2000 continue
C********************************************************************
C     
      stop
      end
C====================================================================
      subroutine initial
C====================================================================
      include 'pointers.h'
      include 'paramsmg.h'
      integer k
C--------------------------------------------------------------------
C     Initialization of some pointers.
C--------------------------------------------------------------------
      do k = 1 , max_lvl
         ia(k) = 0
         ja(k) = 0
         ku(k) = 0
         kb(k) = 0
         ka(k) = 0
         nd(k) = 0
         ipp(k) = 0
         jpp(k) = 0
         kpp(k) = 0
         iord(k) = 0
      end do
C
      ifree = 0
      kfree = 0
C
      ipssm = 0
      iprsm = 0
      ipssm = 0
      nprsm = 0 
      npssm = 0
      icycle = 0
      isolv_coarse = 0
      ich_ord = 0
C
      lf = 1
      lc = 0
C
cc      write(*,'(a$)') 'press ENTER to start the solution process... '
cc      read(*,*)
cc      write(*,*)
C
      return
      end
C====================================================================
