C=====================================================================
      program BILMG
C=====================================================================
      implicit real*8(a-h,o-z)
      parameter(nsubmax = 10, max_dditer = 6)
cc      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
cc      parameter (nlast = 30 000 000, nrlast =15 000 000) !level = 9
       parameter (nlast = 52 000 000, nrlast = 55 000 000) !level = 10
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      dimension ns(nsubmax),idirs(nsubmax),n1s(nsubmax),n2s(nsubmax)
      dimension ias(nsubmax),jas(nsubmax),kas(nsubmax),ksols(nsubmax)
      dimension h1snorms(nsubmax),sl2norms(nsubmax),h1norms(nsubmax)
      character*100 meshf(0:20), solnf(1:20)
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
      common /fname/ meshf
      external exact_sol
C---------------------------------------------------------------------
C...  This program solves the 2nd order linear elliptic PDEs. 
C...  Multigrid method will be used on the bilinear element space..
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
      mtc     = 17
      mts     = 18
      jdf     = ich(3)
      lmethod = ich(4)
      lread   = ich(5)
      nsub    = ich(8)
      lc      = ich(9)
      lfmg    = ich(36)
      lcmg    = ich(37)
      lf      = lfmg
C
cc      write(*,1000) lread,ich(1),lmethod,nsub,lc,ich(23)
 1000 format(2x,' lread   = ', i7 / 2x,
     >          ' lpde    = ', i7 / 2x,
     >          ' lmethod = ', i7 / 2x,
     >          ' nsub    = ', i7 / 2x,
     >          ' lc      = ', i7 / 2x,
     >          ' lconv_f = ', i7)
C
         write(*,*) 
         write(*,*) '================================================'
         write(*,'(2(a,i3))') ' Fine Level:',lf,'   Coarse Level:',lc
         write(*,*) '================================================'
         write(*,*) 
C
      if(lread .eq. 1) then
         open(mt,file=meshf(lf),status='unknown',form='formatted')
         open(mtc,
     >        file=meshf(lc),status='unknown',form='formatted')
         read(mt,*) n1,n2
         read(mt,*) nel,n,ned
         read(mtc,*) n1c,n2c
         read(mtc,*) nelc,nc,nedc
      else
         open(mt,file=meshf(lf),status='unknown',form='unformatted')
         open(mtc,
     >        file=meshf(lc),status='unknown',form='unformatted')
c     c         rewind(mt)
         read(mt) n1,n2
         read(mt) nel,n,ned
         read(mtc) n1c,n2c
         read(mtc) nelc,nc,nedc
      end if
      write(*,*) ' No. of elements, nodes, & bdry edges:', nel,n,ned
C   
C.... Get the value of h
C
      hhhh = dsqrt(dble(2)/dble(nel))
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

C
      call pde_choice(lmethod,
     I     nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     O     m(ia),m(ja),r(kra),nnz,r(krhs),
     W     ned,m(inedg),m(idir),m(mfree))
C
C...  The choice of the right hand side and the exact FEM solution:
C...  1 - if the exact solution is known;
C...  0 - get the exact FEM solution from other source.
C
      ichoice_rhs = 0
      ku_exact = kfree
C
      if (ichoice_rhs .eq. 1) then
C
C...     Interpolation of the exact solution if it is known.
C
         call interp_func(exact_sol,r(ku_exact),n,r(kx),r(ky))
C
C...     Compute RHS = A*Exact_sol.
C
         call abyvg(m(ia),m(ja),r(kra),r(ku_exact),n,r(krhs))
      else 
         solnf(2) = 'sol2'
         solnf(3) = 'sol3'
         solnf(4) = 'sol4'
         solnf(5) = 'sol5'
         solnf(6) = 'sol6'
         solnf(7) = 'sol7'
         solnf(8) = 'sol8'
         solnf(9) = 'sol9'
         solnf(10)= 'sol10'
C
         if (lread .eq. 1) then
            open(mts,file=solnf(lf),status='unknown',form='formatted')
            read(mts,*) (r(ku_exact+l-1), l = 1, n)
         else
            open(mts,file=solnf(lf),status='unknown',form='unformatted')
            read(mts) (r(ku_exact+l-1), l = 1, n)
         end if
      end if
C
      kb     = ku_exact + n
      kfree  = kb + n
C

      
c
C...  Solve the SPD problem by MG V-cycle.
C     
      call mg_s(m(1),m(ia),m(ja),m(idir),mfree,
     >     r(1),r(ksol),r(krhs),r(kra),kfree,
     >     n,nx,ny,nnz,nel,lfmg,lcmg,tol) 


      subroutine mg_s(m,iao,jao,idiro,ifree,r,ksol,krhs,kra,kfree,
     >     n,nx,ny,nnz,nel,lf,lc) 

C     
C
         call wuminv(r(kerr),r(ku_exact),r(ksol),n)
cc         call l2norm(r(kerr),snorm,n,hhhh)
         call l2_h1s_norm(0,r(kerr),n,sl2norm,nel,jdf,m(ie),m(je),
     >        r(kx),r(ky))
         call l2_h1s_norm(1,r(kerr),n,h1snorm,nel,jdf,m(ie),m(je),
     >        r(kx),r(ky))
         h1norm = dsqrt(sl2norm*sl2norm + h1snorm*h1snorm)
C
cc         write(*,'(a,i6,4(e14.4))') 
cc     >        ' No.iter     L_2      H_1: ',kiter,sl2norm,h1norm
C
         write(*,'(i6,2(e14.4))') kiter, sl2norm, h1norm
C
      end do
C
C...  Write the exact solution after many iterations.
cc      call  write_soln(lf, r(ksol), n)
C
      if (n .lt. 18000) then
cc         call wuminv(r(kerr),r(ksol),r(ku_exact),n)     
cc         call abyvg(m(ia),m(ja),r(kra),r(kerr),n,r(krhs))
cc         call output(r(kx),r(ky),r(ksol),n)
      end if

      close(mt)
      close(mtc)
      close(mts)
C
      stop
      end
C=====================================================================

