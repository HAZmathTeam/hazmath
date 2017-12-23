C=====================================================================
      program TG_EIG
C=====================================================================
      implicit real*8(a-h,o-z)
cc      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
cc      parameter (nlast = 30 000 000, nrlast =15 000 000) !level = 9
       parameter (nlast = 52 000 000, nrlast = 55 000 000) !level = 10
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
C...  This program computes the smallest eigvalue of a variational 
C...  formulation of the Laplace equation using a two-grid algorithm:
C...
C...    (\nabla u_h, \nabla v) = \lamda_h (u_h, v) \forall v in V_h,
C...
C...  i.e., A_h*u_h = \lamda_h*M_h*u_h.
C...
C...  Standard Galerkin scheme is used. For the solver of linear systems 
C...  a V-cycle multigrid scheme is used.
C...
C...  Parameter:
C...    istiff - choice of stiffness matrices.
C...             1 = Standard Galerkin, (2 = EAFE.)
C---------------------------------------------------------------------
      call input_data
      call quad_elt
      call quad_data
      call quad_data1
C 
      mt      = 16
      mtc     = 17
      jdf     = ich(3)
      lread   = ich(5)
      lc      = ich(9)
      lfmg    = ich(36)
      lcmg    = ich(37)
      lf      = lfmg
      tol     = 1.d-6
C
cc      write(*,1000) lread,ich(1),lmethod,nsub,lc,ich(23)
 1000 format(2x,' lread   = ', i7 / 2x,
     >          ' lpde    = ', i7 / 2x,
     >          ' lmethod = ', i7 / 2x,
     >          ' nsub    = ', i7 / 2x,
     >          ' lc      = ', i7 / 2x,
     >          ' lconv_f = ', i7)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      lb = lf
      le = 8
 10   continue
      do 980 ll = lb, le, 2
         lf   = ll
         lfmg = lf
         lc   = lf/2
         write(*,*) 
         write(*,*) '================================================'
         write(*,'(2(a,i3))') ' Fine Level:',lf,'   Coarse Level:',lc
         write(*,*) '================================================'
         write(*,*) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
cc      hhhh = dsqrt(dble(2)/dble(nel))
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
C...  Compute the stiffness matrix on the fine grid.
C
      lm = 1              !1 - stiffness matrix,  2 -  mass matrix.
      call lapl_mass(lm,nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     >     m(ia),m(ja),r(kra),nnz,m(idir),m(iip))
C
C...  Compute the mass matrix on the fine grid.
C
      krm = krat
      lm = 2              !1 - stiffness matrix,  2 -  mass matrix.
      call lapl_mass(lm,nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     >     m(ia),m(ja),r(krm),nnz,m(idir),m(iip))
C
C...  Finding the stiffness matrix and the mass matrix on the coarse grid.
C
      ifree   = iet
      iac     = ifree
      jac     = iac + nc + 1
      iec     = jac + 2*(nelc + nc + 5) + nc
      jec     = iec + nelc + 1      
      idirc   = jec + nelc*3
      inedgc  = idirc + nc + 5
      ietc    = inedgc + nedc * 2
      jetc    = ietc + nc + 1
      lastic  = jetc + nelc*3
C
      kfree   = krm + nnz
      keigc   = kfree
      krac    = keigc + nc     !Use N instead of NC if updated later.
      krhsc   = krac + 2*(nelc + nc + 5) + nc
      ksolc   = krhsc + nc
      kxc     = ksolc + nc
      kyc     = kxc + nc
      lastrc  = kyc + nc
C
      call rdmesh(lread,mtc,m(iec),m(jec),nedc,m(inedgc),
     >     m(idirc),r(kxc),r(kyc),nelc,nc,jdf)
C
      call iit(m(iec),m(jec),nelc,nc,m(ietc),m(jetc))
C
      nnzc = 0
      call smbasg(m(iec),m(jec),m(ietc),
     >     m(jetc),m(iac),m(jac),nc,nnzc,m(idirc))
C
      iip = ietc
C
C...  Compute the stiffness matrix on the coarse grid.
C
      lm = 1              !1 - stiffness matrix,  2 -  mass matrix.
      call lapl_mass(lm,nelc,nc,jdf,m(iec),m(jec),r(kxc),r(kyc),
     >     m(iac),m(jac),r(krac),nnzc,m(idirc),m(iip))
C
      call wradj_r(m(iac),m(jac),r(krac),min0(nc,65*65),201)
C...  Compute the mass matrix on the coarse grid.
C
      krmc  = lastrc
C
      lm = 2              !1 - stiffness matrix,  2 -  mass matrix.
      call lapl_mass(lm,nelc,nc,jdf,m(iec),m(jec),r(kxc),r(kyc),
     >     m(iac),m(jac),r(krmc),nnzc,m(idirc),m(iip))
C
      call wradj_r(m(iac),m(jac),r(krmc),min0(nc,65*65),202)
C
cc      call lprr(m(iac),m(jac),r(krmc),nc)
cc      stop 1
C
C****************************************************************************
C...  STEP 1. Find the smallest eigenpair (eigvalc,eigc) on the coarse grid.*
C****************************************************************************
C...  Compute the smallest eigenpair (eigvalc,eigc) on the coarse grid.
C...  They are computed exactly.
C
      ifree = ietc
      kfree = krmc + nnzc
      tol   = 1.d-12
C
      call inv_pwr(eigvalc,m(1),iac,jac,idirc,ifree,
     >     r(1),ksolc,krhsc,krac,krmc,keigc,kfree,
     >     nc,n1c,n2c,nnzc,nelc,lc,lcmg,tol)
C
      eigvalc = 1.d0/eigvalc
C
C...  Compute the smallest eigenpair (eigvalf,eigf) on the fine grid.
C...  Thery are computed exactly and used for error estimates later.
C
      ifree = iet 
      keigf = krac
      kfree = keigf + n
      tol   = 1.d-12
C
      call inv_pwr(eigvalf,m(1),ia,ja,idir,ifree,
     >     r(1),ksol,krhs,kra,krm,keigf,kfree,
     >     n,n1,n2,nnz,nel,lfmg,lcmg,tol)
C
      eigvalf = 1.d0/eigvalf
C
C****************************************************************************
C... Step 2. Solve SPD system on the fine grid to get an eigenvector SOL.   *
C****************************************************************************
C...  Prolongate EIGC to a fine grid vector: sol = P*eigc.
C
      ipp = ifree
      jpp = ipp + n + 1
      lasti = jpp + 2*n
C     
      lftwo = lf
      lctwo = lc
C
      call copyv(r(keigc),r(ksol),nc)
C
      do i = lctwo+1,lftwo
         n1f = n1c*2 - 1
         n2f = n2c*2 - 1
         nf  = n1f*n2f
C     
         call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
         call pbyvg(m(ipp),m(jpp),r(ksol),nc,r(krhs),nf)
         call copyv(r(krhs),r(ksol),nf)
C     
         n1c = n1f
         n2c = n2f
         nc  = nf
      end do
C     
C...  Perform RHS = eigvalc*M*sol = eigvalc*M*(P*eigc).
      call abyvg(m(ia),m(ja),r(krm),r(ksol),n,r(krhs))
      call usmultu(r(krhs),n,eigvalc)
C
      tol  = 1.d-13
      call mg_s(m(1),ia,ja,idir,ifree,r(1),ksol,krhs,kra,kfree,
     >     n,n1,n2,nnz,nel,lfmg,lcmg,tol) 
C
C****************************************************************************
C...  Step 3. Rayleigh quotient to get an eigvalue: sol*A*sol/sol*M*sol     *
C****************************************************************************
C...  To compute H^1 seminorm squared, i.e., sol^t*A*sol
C
      call abyvg(m(ia),m(ja),r(kra),r(ksol),n,r(krhs))
      call scpro(r(ksol),r(krhs),eig_norm_h1,n)
C
C...  To compute L^2 norm squared, i.e., sol^t*M*sol
C
      call abyvg(m(ia),m(ja),r(krm),r(ksol),n,r(krhs))
      call scpro(r(ksol),r(krhs),eig_norm_l2,n)
C
      eigval = eig_norm_h1/eig_norm_l2
C
C...  Normalization if it is needed.
C
      eig_norm = dsqrt(eig_norm_h1)
      eig_norm_inv = 1.d0/eig_norm
      call usmultu(r(ksol),n,eig_norm_inv)
C
cc      print*, 'eigvector2 = ',(i,'*',r(ksol+i-1), i=1,n)
C
C****************************************************************************
C...  Error estimates: | eigval - eigvalf | and || sol - eigf ||_1          *
C****************************************************************************
C
      err_eigval = dabs(eigval - eigvalf)
C
C...  Compute sol = sol - eigf.
      call uuminv(r(ksol),r(keigf),n)
C
C...  Compute H^1 seminorm: err_h1s_sqrt = \sqrt(sol*A*sol)
      call abyvg(m(ia),m(ja),r(kra),r(ksol),n,r(krhs))
      call scpro(r(ksol),r(krhs),err_h1s,n)
      err_h1s_sqrt = dsqrt(err_h1s)
C
C...  Compute L^2 norm: err_l2_sqrt = \sqrt(sol*M*sol)
      call abyvg(m(ia),m(ja),r(krm),r(ksol),n,r(krhs))
      call scpro(r(ksol),r(krhs),err_l2,n)
      err_l2_sqrt = dsqrt(err_l2)
C
C...  Compute H^1 norm: err_h1 = \sqrt(err_h1s + err_l2)
      err_h1 = dsqrt(err_h1s + err_l2)
C
C...  Print out numerical results.
C
      write(*,*) 
      write(*,*) 
     >     ' ------------------------------------------------------',
     >     '--------------------'
      write(*,'(2(a))') 
     >     '  Coarse  Fine   Exact(coarse)    Exact(fine)',
     >     '       Computed        H^1 semi'
      write(*,*)
     >     ' ------------------------------------------------------',
     >     '--------------------'
      write(*, '(3x,i3,3X,i3,2X,3(3X,f13.10),3X,f11.8)')
     >     lc,lf,eigvalc,eigvalf,eigval,eig_norm
C
      write(*,*) 
      write(*,*) 
     >     ' ------------------------------------------------------',
     >     '--------------------'
      write(*,'(2(a))') 
     >     '  Coarse  Fine      L^2           H^1-semi          H^1',
     >     '          Eigenvalue'
      write(*,*)
     >     ' ------------------------------------------------------',
     >     '--------------------'
      write(*, '(3x,i3,3X,i3,4(2X,e14.7))')
     >     lc,lf,err_l2_sqrt,err_h1s_sqrt,err_h1,err_eigval
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 980  continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      stop
      end
C=====================================================================

