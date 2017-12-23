C=====================================================================
      program CDFEM
C=====================================================================
      parameter (nlast = 5 000 000, nrlast = 3 000 000)
cc      parameter (nlast = 20 000 000, nrlast = 10 000 000)
cc      parameter (nlast = 60 000 000, nrlast = 30 000 000)
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
                          do 980 ll = la,5,2
                          lcmg = ll
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      lf   = lfmg
      if (lcmg .ge. lf) lf = lcmg + 1
 50   continue 
C
      if(lread .eq. 1) then
         if (lsolve .eq. 3) then
            open(mt,file=meshf(lf),status='unknown',form='formatted')
         else 
            open(mt,file=meshf(0),status='unknown',form='formatted')
         end if
cc         rewind(mt)
         read(mt,*) n1,n2
         read(mt,*) nel,n,ned
      else
         if (lsolve .eq. 3) then
            open(mt,file=meshf(lf),status='unknown',form='unformatted')
         else 
            open(mt,file=meshf(0),status='unknown',form='unformatted')
         end if
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
      iip   = iet
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
      if(ich(1) .eq. 3 .and. ich(4) .eq. 1 .and. ich(6) .eq. -1) then
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
      img = inedg
      kmg = krat
      tol   = 1.d-12
      maxit = ich(39)
C
      call mg(m(ia),m(ja),m(idir),m(img),r(ksol),r(krhs),r(kra),
     >        r(kmg),n,nnz,lf,lcmg,maxit,tol)
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
C...  OUTPUT.
C
      if (n .lt. 2000) then
         call output(r(kx),r(ky),r(ksol),n)
         call output_sol(m(ie),m(je),r(kx),r(ky),r(ksol),n,nel)
         endfile(800)
         close(800)
         endfile(100)
         close(100)
         endfile(199)
         close(199)
      end if 
C
      if (lsolve .ne. 3) go to 990
C
      if (lsolve .eq. 3 .and. lf .lt. 10) then
         lf = lf + 1
         close(mt)
cc         read(*,*)
         go to 50
      end if
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 980                      continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 990  close(mt)
      stop
      end
C=====================================================================
