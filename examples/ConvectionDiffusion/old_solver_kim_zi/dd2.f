C=====================================================================
      program DDM
C=====================================================================
      implicit real*8(a-h,o-z)
      parameter(nsubmax = 10, max_dditer = 2, tol_dd = 1.d-06)
cc      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
cc        parameter (nlast = 30 000 000, nrlast =15 000 000) !level = 9
       parameter (nlast = 100 000 000, nrlast = 50 000 000) !level = 10
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
C...  This program solves the 2nd order linear elliptic PDEs. Domain
C...  decomposition is combined with two grid method and each subdomain
C...  is solved separately in a parallel manner.
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
      tol     = 1.d-6
C
      write(*,1000) lread,ich(1),lmethod,nsub,lc,ich(23)
 1000 format(2x,' lread   = ', i7 / 2x,
     >          ' lpde    = ', i7 / 2x,
     >          ' lmethod = ', i7 / 2x,
     >          ' nsub    = ', i7 / 2x,
     >          ' lc      = ', i7 / 2x,
     >          ' lconv_f = ', i7)
C
 10   continue
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
      iip = iet
C
C...  First we do the Laplace equation on the fine grid and then it will
C...  be used in the subroutine GET_BLOCKS.
C
      lm = 1              !1 - stiffness matrix,  2 -  mass matrix.
C
      call lapl_mass(lm,nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     >     m(ia),m(ja),r(kra),nnz,m(idir),m(iip))
C
C...  Finding stiffness matrix for the coarse grid.
C
      mfree   = iet
      iac     = mfree
      jac     = iac + nc + 1
      iec     = jac + 2*(nelc + nc + 5) + nc
      jec     = iec + nelc + 1      
      idirc   = jec + nelc*3
      inedgc  = idirc + nc + 5
      ietc    = inedgc + nedc * 2
      jetc    = ietc + nc + 1
      lastic  = jetc + nelc*3
C
      kfree   = krat
      krac    = kfree
      krhsc   = krac + 2*(nelc + nc + 5) + nc
      kxc     = krhsc + nc
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
C...  Finding the coarse grid stiffness matrix for the non-SPD problem 
C...  (SPD if solving the Laplace equation).
C
      call pde_choice(lmethod,
     I     nelc,nc,jdf,m(iec),m(jec),r(kxc),r(kyc),
     O     m(iac),m(jac),r(krac),nnzc,r(krhsc),
     W     nedc,m(inedgc),m(idirc),m(iip))
C
C...  Get the information about subdomains
C
      iwork   = idir
      mfree   = iec
      ielsub  = mfree
      mfree   = ielsub + nel
C
      kfree   = krhsc
C
      call formsub(lf,n,n1,n2,nsub,ns,n1s,n2s,m(mfree),mfree,
     >     isub_all,jsub_all,
     >     jsub_mid,idirs,m(ielsub),nel)
C
C...  Form the stiffness matrices for the Laplacian on each subdomain.
C
      call get_blocks(n,m(ia),m(ja),r(kra),
     >     ias,jas,kas,nsub,ns,
     >     m(isub_all),m(jsub_all),m(iwork),m(idirs(1)),
     >     m(mfree),r(kfree),mfree,kfree)
C
C...  Identify the overlapped nodes.
C
      node_overlap = mfree
      mfree        = node_overlap + n
C
      call inullv(m(node_overlap),n)
C
C...  ioverlap = 0 : non-overlapping domain (averaged on the 
C...                 common boundaries (midpoints) of subdomains only)
C...  ioverlap = 1 : overlapping domain (averaged on the 
C...                 overlapped regions of subdomains)
C
      ioverlap = 1
C
      if (ioverlap .eq. 0) then
         call overlap_mid(nsub,m(isub_all),m(jsub_all),ns,
     >        m(jsub_mid),m(node_overlap))
      else if (ioverlap .eq. 1) then
         call overlap_all(nsub,m(isub_all),m(jsub_all),ns,
     >        m(idirs(1)),m(node_overlap))
      end if
C
      iharmonic = 1                  !Harmonic extension or not.
C
      if (iharmonic .eq. 1) then
C...     !Information needed for Harmonic extension on the overlapping domain.
C
         ififth_inv = mfree
         ififth     = ififth_inv + n
C
         call formlast(m(node_overlap),m(ififth_inv),m(ififth),
     >        n,nsfifth)
C
         ia5 = ififth + nsfifth
         ja5 = ia5 + nsfifth + 1
         ka5 = kfree
C
         call get_fifth_block(m(ia),m(ja),r(kra),nsfifth,
     >        m(ififth_inv),m(ififth),m(ia5),m(ja5),r(ka5),nnz5)
C
         mfree = ja5+ + nnz5
         kfree = ka5+ + nnz5
      end if
C
C...  Form the fine grid stiffness matrix for non-SPD problem (SPD if
C...  solving the Laplace equation).
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
      kfree1 = kfree
C
      mfree1 = mfree
      nf = n
C
C****************** Iteration starts here. **************************
C
C...  Initial guess
C
      call nullv(r(ksol),n)
C
      kiter = 0
C
C...  The initial guess is 0, so the initial 'error' = ||u_exact||
C
cc      call l2norm(r(ku_exact),snorm,n,hhhh)
      call l2_h1s_norm(0,r(ku_exact),n,sl2norm,nel,jdf,m(ie),m(je),
     >     r(kx),r(ky))
      call l2_h1s_norm(1,r(ku_exact),n,h1snorm,nel,jdf,m(ie),m(je),
     >     r(kx),r(ky))
      h1norm = dsqrt(sl2norm*sl2norm + h1snorm*h1snorm)
C
      write(*,'(a,i6,4(e14.4))') 
     >     ' No.iter     L_2      H_1: ',  kiter, sl2norm, h1norm
C
      do while(kiter .lt. max_dditer)
         kiter = kiter + 1
C
C...     To compute the residual: b = rhs - A*sol.
C
         call abyvam(r(krhs),m(ia),m(ja),r(kra),r(ksol),n,r(kb))
C
C...     Restriction to the coarse grid: b = P^t*b.
C
         ipp = mfree
         jpp = ipp + nf + 1
         idirtwo = jpp + 2*nf
C     
         ktemp = kfree
C     
         nf = n
         n1f = n1
         n2f = n2
         lftwo = lf
         lctwo = lc
C     
         call icopyv(m(idir),m(idirtwo),nf)
C     
         do i = lftwo, lctwo+1, -1
            call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
            call rvbyp(m(ipp),m(jpp),r(kb),nf,r(ktemp),nc)
            call copyv(r(ktemp),r(kb),nc)
            call dir_two(r(kb),m(ipp),m(jpp),nf,nc,m(idirtwo))
C     
            n1f = n1c
            n2f = n2c
            nf = nc
         end do
C
C...     Solve the residual equation on the coarse grid by Gaussian 
C...     elimination: A*e_H = b.
C
         call sgauss(m(iac),m(jac),r(krac),r(kb),m(mfree1),r(kfree1),nc)
C
C...     Interpolation (prolongation) to the fine grid: u^h = P*e_H.
C
         do i = lctwo+1,lftwo
            n1f = n1c*2 - 1
            n2f = n2c*2 - 1
            nf  = n1f*n2f
C
            call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
            call pbyvg(m(ipp),m(jpp),r(kb),nc,r(ktemp),nf)
            call copyv(r(ktemp),r(kb),nf)
C     
            n1c = n1f
            n2c = n2f
            nc  = nf
         end do
C     
C...     Update the solution by adding the coarse grid correction: 
C...     sol = sol + u^h.
C
         call uupluv(r(ksol),r(kb),n)
C
C....    Compute the norm of the error after each half iteration.
C
         kerr  = kfree
         kfree = kerr + n
C
         call wuminv(r(kerr),r(ku_exact),r(ksol),n)
cc         call l2norm(r(kerr),snorm,n,hhhh)
         call l2_h1s_norm(0,r(kerr),n,sl2norm,nel,jdf,m(ie),m(je),
     >        r(kx),r(ky))
         call l2_h1s_norm(1,r(kerr),n,h1snorm,nel,jdf,m(ie),m(je),
     >        r(kx),r(ky))
         h1norm = dsqrt(sl2norm*sl2norm + h1snorm*h1snorm)
C
         write(*,'(a,f6.1,4(e14.4))') 
     >        ' No.iter     L_2      H_1: ',kiter-0.5,sl2norm,h1norm
C
C...     To compute the residual: b = rhs - A*sol for non-SPD.
C
         call abyvam(r(krhs),m(ia),m(ja),r(kra),r(ksol),n,r(kb))
C
C...     To solve the SPD system of the residual equation on each subdomain.
C
         call nullv(r(kerr),n)
C
         nelf  = nel
         mfr   = mfree
C
         ksols(1) = kfree
         do k = 2, nsub
            ksols(k) = ksols(k-1) + ns(k-1)
         end do
         kfree = ksols(nsub) + ns(nsub)
         kfr   = kfree 
C
         if (nsub .gt. 1) then
            isubb = 1
            isube = nsub
            do i = isubb,isube
               kbs   = kfree
               kfree = kbs + ns(i)
               iadr  = m(isub_all+i-1) - 1
C
C...           Distribute the residual to each subdomain.
C
               call scatter(m(jsub_all+iadr),r(kb),r(kbs),ns(i),
     >              m(idirs(i)))
C
C...           Solve the SPD problem on each subdomain by MG V-cycle.
C
               call mg_s(m(1),ias(i),jas(i),
     >              idirs(i),mfree,r(1),ksols(i),kbs,kas(i),kfree,
     >              ns(i),n1s(i),n2s(i),nnz,nel,lfmg,lcmg,tol) 
C
C...           Collect the coarse grid solutions to combine them into the
C...           global working array ERR.
C
               if (ioverlap .eq. 0) then
                  call gather_mid(m(jsub_all+iadr),r(ksols(i)),r(kerr),
     >                 ns(i),m(jsub_mid+iadr))
               else if (ioverlap .eq. 1) then
                  call gather_all(m(jsub_all+iadr),r(ksols(i)),r(kerr),
     >                 ns(i),m(idirs(i)))
               end if
C
               kfree = kfr
               mfree = mfr
               nel   = nelf
            end do
         else
            i     = 1
            kbs   = kfree
            kfree = kbs  +  ns(i)
            iadr  = m(isub_all+i-1) - 1
C
C..         Distribute the residual to one subdomain (i.e., domain itself).
C
            call scatter(m(jsub_all+iadr),r(kb),r(kbs),ns(i),
     >           m(idirs(i)))
C
C...        Solve the SPD problem on one subdomain by MG V-cycle.
C
            call mg_s(m(1),ias(i),jas(i),
     >           idirs(i),mfree,r(1),ksols(i),kbs,kas(i),kfree,
     >           ns(i),n1,n2,nnz,nel,lfmg,1,tol) 
C     
C...        Collect the coarse grid solution to combine it into the
C...        global working array ERR.
C
            call gather_all(m(jsub_all+iadr),r(ksols(i)),r(kerr),
     >           ns(i),m(idirs(i)))
C
            kfree = kfr
            mfree = mfr
            nel   = nelf
         end if
C     
C...     Define the errors within the subdomains if NSUB > 1.
C
         if (nsub .gt. 1) then
            globl2  = 0.0d00
            globh1s = 0.0d00
            globh1  = 0.0d00
            do i = isubb,isube
               iadr = m(isub_all+i-1) - 1
               call subd_norm(1,
     >              ns(i),r(ksols(i)),n,r(ksol),r(kb),r(ku_exact),
     >              h1snorms(i),i,m(ielsub),nel,
     >              m(ie),m(je),jdf,r(kx),r(ky),m(jsub_all+iadr))
C     
               call subd_norm(0,
     >              ns(i),r(ksols(i)),n,r(ksol),r(kb),r(ku_exact),
     >              sl2norms(i),i,m(ielsub),nel,
     >              m(ie),m(je),jdf,r(kx),r(ky),m(jsub_all+iadr))
C     
               h1norms(i) = sl2norms(i) + h1snorms(i)
               globl2     = globl2 + sl2norms(i)
               globh1s    = globh1s + h1snorms(i)
               globh1     = globh1 + h1norms(i)
C     
               h1snorms(i) = dsqrt(h1snorms(i))
               sl2norms(i) = dsqrt(sl2norms(i))
               h1norms(i)  = dsqrt(h1norms(i))
C
               write(*,'(a,2(i4),4(e13.4))') 
     >              ' No.iter  Subd  L^2  H^1s  H^1:',
     >              kiter,i,sl2norms(i),h1snorms(i),h1norms(i)
            end do
C
            write(*,'(a,i4,4(e13.4))') 
     >           ' Global  No.iter  L^2  H^1s  H^1:',
     >           kiter,dsqrt(globl2),dsqrt(globh1s),dsqrt(globh1)
         end if
C     
C...     Update the global solution by adding the solutions of the
C...     resudual equations of each subdomain.
C
         call update(r(ksol),r(kerr),n,m(node_overlap))
C
         if (iharmonic .eq. 1)  then
C...        Harmonic extension on the overlapping domain.
C
            kb5 = kfree
            kx5 = kb5 + nsfifth
            do kk = 1 , 1
               call abyvam(r(krhs),m(ia),m(ja),r(kra),r(ksol),n,r(kb))
               call scatter1(m(ififth),r(kb),r(kb5),nsfifth)
               call nullv(r(kx5),nsfifth)
               call sym_gs_ns(r(kx5),m(ia5),m(ja5),r(ka5),r(kb5),
     >              nsfifth,nsfifth*2,1)
C
               do k = 1 , nsfifth
                  i5 = m(ififth+k-1)
                  r(ksol+i5-1) = r(ksol+i5-1) + r(kx5+k-1)
               end do
            end do  
            kfree = kb5
         end if
C
C....    Compute the norm of the error after each iteration.
C
         call wuminv(r(kerr),r(ku_exact),r(ksol),n)
cc         call l2norm(r(kerr),snorm,n,hhhh)
         call l2_h1s_norm(0,r(kerr),n,sl2norm,nel,jdf,m(ie),m(je),
     >        r(kx),r(ky))
         call l2_h1s_norm(1,r(kerr),n,h1snorm,nel,jdf,m(ie),m(je),
     >        r(kx),r(ky))
         h1norm = dsqrt(sl2norm*sl2norm + h1snorm*h1snorm)
C
         write(*,'(a,i6,4(e14.4))') 
     >        ' No.iter     L_2      H_1: ',kiter,sl2norm,h1norm
C
         kfree = kfree1
         mfree = mfree1
      end do
C
      if (n .lt. 18000) then
         call wuminv(r(kerr),r(ksol),r(ku_exact),n)
C     
cc         call abyvg(m(ia),m(ja),r(kra),r(kerr),n,r(krhs))
C     
         call output(r(kx),r(ky),r(ksol),n)
      end if
 990  close(mt)
      close(mtc)
      close(mts)
      stop
      end
C=====================================================================

