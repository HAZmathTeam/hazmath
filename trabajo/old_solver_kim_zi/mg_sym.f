C=====================================================================
      subroutine mg_s(m,iao,jao,idiro,ifree,r,ksol,krhs,kra,kfree,
     >     n,nx,ny,nnz,nel,lf,lc,tol) 
C=====================================================================
      parameter (lmx = 20)
      implicit real*8(a-h,o-z)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx),n1(lmx),n2(lmx)
      character*100 meshf(0:lmx)
      common /menu/ ich(200)
      common /fname/ meshf
C---------------------------------------------------------------------
C...  Standard multigrid cycle for the SPD problems.
C...
C...  Parameter:
C...    N,NX,NY - dimension of the mesh, N = NX x NY
C...    ICYCLE  - choice of MG:  1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle,(4 = W-cylce, 5 = Full MG)
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    MAXIT   - maximum number of MG iterations
C...    TOL     - stopping criterion for iteration procedure.
C---------------------------------------------------------------------
      mt     = 15
      icycle = ich(31)
      maxit  = ich(39)
C
cc      tol    = 0.5d-6
      zero   = 1.d-13
C
C...  Generating the stiffness matrix A and the load vector B at each level.
C
      lasti = ifree
      lastr = kfree
      nelnew = nel
      nd(lf) = n
      n1(lf) = nx
      n2(lf) = ny
C
      do 50 i = lf, lc, -1

         if (i .eq. lf) then
            nz(i)   = nnz
            ia(i)   = iao
            ja(i)   = jao
            idir(i) = idiro
            ipp(i)  = lasti
         else
            nelnew  = nelnew/4
            ia(i)   = lasti
            ja(i)   = ia(i) + nd(i) + 1
            idir(i) = ja(i) + 2*(nelnew + nd(i) + 5) + nd(i)
            ipp(i)  = idir(i) + nd(i) + 5
         end if
C
         jpp(i)  = ipp(i) + nd(i) + 1
         lasti   = jpp(i) + 2*nd(i)
C
         ku(i)   = lastr
         kb(i)   = ku(i) + nd(i)
         if (i .eq. lf) then
            ka(i)   = kra
            lastr   = kb(i) + nd(i)
         else
            ka(i)   = kb(i) + nd(i)
            lastr   = ka(i) + 2*(nelnew + nd(i) + 5) + nd(i)
         end if
C
C...     To compute the symbolic prolongation matrix P from level i-1
C...     to level i.
C
         if (i .eq. lc) go to 40
         call prolng_symb(n1(i),n2(i),nd(i),n1(i-1),n2(i-1),nd(i-1),
     >        m(ipp(i)),m(jpp(i)),nzp(i))
 40      continue
C
C...     To compute stiffness matrices A at each level.
C
         if (i .eq. lf) then
            ia(i) = iao
            ja(i) = jao
            ka(i) = kra
            nz(i) = nnz
         else
            call formdir(m(idir(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),m(idir(i)),nd(i))
C
            ic = lasti
            ix = ic + nd(i+1) + 1
            jc = ix + nd(i) + 1
C
            call abybs(m(ia(i+1)),m(ja(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),nd(i+1),nd(i),m(ic),m(jc),m(ix))
C
            jce = m(ic+nd(i+1)) + 1
            kc  = lastr
            kx  = kc + jce
            call abyp(m(ia(i+1)),m(ja(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),m(ic),m(jc),r(ka(i+1)),r(kc),r(kx),nd(i+1))
C
            ict = jc + jce + 1
            jct = ict + nd(i+1) + 1
            kct = kx + nd(i+1) + 1
C
            call aat(m(ic),m(jc),r(kc),nd(i+1),nd(i),m(ict),m(jct),
     >           r(kct))
C
            call abybs(m(ict),m(jct),m(ipp(i+1)),m(jpp(i+1)),nd(i),
     >           nd(i+1),nd(i),m(ia(i)),m(ja(i)),m(ix))
C     
            call abyp(m(ict),m(jct),m(ipp(i+1)),m(jpp(i+1)),nd(i),
     >           m(ia(i)),m(ja(i)),r(kct),r(ka(i)),r(kx),nd(i+1))
C     
            nz(i) = m(ia(i)+nd(i))
C
            call adir_new(m(ia(i)),m(ja(i)),r(ka(i)),m(idir(i)),nd(i))
C
cc            write(*,*) 'printing ia, ja, ka'
cc            call lprr(m(ia(i)),m(ja(i)),r(ka(i)),nd(i))
cc            call outmat1(m(ia(i)),m(ja(i)),r(ka(i)),nd(i),nz(i))
         end if
 50   continue
C
C...  ------------------- MG starts here. ---------------------------
C
cc      write(*,*)
cc      write(*,*) '=============================', 
cc     >           '=================================================='
cc      if (icycle .eq. 1) then
cc         write(*,*) '     METHOD :          MG : V-cycle '
cc      else if (icycle .eq. 2) then
cc         write(*,*) '     METHOD :          MG : Backslash cycle '
cc      else if (icycle .eq. 3) then
cc         write(*,*) '     METHOD :          MG : Variable V-cycle '
cc      end if
ccC
cc      write(*,*) '     Residual :        r := b-Au; '
cc      write(*,'(a,a,e12.4)')
cc     >     '      Convergence :     ||r||_0/||r_0||_0  &',
cc     >     '  ||Br||_0/||u||_0 <', tol
cc      write(*,'(a,2(5X,i3))') '      Coarse and fine levels : ',lc,lf
cc      write(*,*) '-----------------------------',
cc     >           '--------------------------------------------------'
cc      write(*,*) 'Iteration  ||r||/||r_0||    ||r||  ',
cc     >                  '   ||Br||/||u||      ||Br||   ||Br||/||r_0||'
cc      write(*,*) '-----------------------------',
cc     >           '--------------------------------------------------'
C
      kbf   = kb(lf)
      kuf   = ku(lf)
C
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
      call init_g0(m(iao),m(jao),r(kra),r(krhs),r(ksol),n)
C
cc      write(*,*) (r(ksol+ii-1), ii = 1,n)
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
C...  To compute the initial residual: b = rhs - A*sol.
      call abyvam(r(krhs),m(iao),m(jao),r(kra),r(ksol),n,r(kbf))
C
      call scpro(r(kbf),r(kbf),err_res0,n)
      err_res0 = dsqrt(err_res0)
C
      if (err_res0 .lt. zero)  err_res0 = 1.d0
C
      kk = 0
cc      write(*,'(2X,i5,2X,5(3X,e11.4))') 
cc     >     kk,1.d0,err_res0,1.d0,1.d0,1.d0
C
C...  Iterate until convergence:  sol = sol + B*b = sol + B(rhs - A*sol).
C
      do 200 kk = 1, maxit
C
C...     Iterator B using an MG-cycle, i.e., u = B*b = B(rhs - A*sol).
C
         call premg_s(m(1),ia(1),ja(1),idir(1),
     >        ipp(1),jpp(1),nzp(1),nz(1),nd(1),
     >        r(1),ku(1),kb(1),ka(1),lmx,lf,lc,lasti,lastr)
C
C...     New solution:  sol = sol + u = sol + B(rhs - A*sol).
         call uupluv(r(ksol),r(kuf),n)
C
C...     To compute the relative error.
	 call scpro(r(kuf),r(kuf),xdiffn,n)
         call scpro(r(ksol),r(ksol),xnewn,n)
C
         xdiffn  = dsqrt(xdiffn)
         xnewn   = dsqrt(xnewn)
C
         if (xnewn .gt. xdiffn) then
            err_rel = xdiffn/xnewn
         else
            err_rel = xdiffn  
         end if
         err_brr = xdiffn/err_res0
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(r(krhs),m(iao),m(jao),r(kra),r(ksol),n,r(kbf))
C
C...     To compute the L_2 norm of the residual.
         call scpro(r(kbf),r(kbf),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
         nmod = mod(kk,3)
C
         if (nmod .eq. 1) then
cc            write(*,'(2X,i5,2X,5(3X,e11.4))')
cc     >           kk,err_relres,err_res,err_rel,xdiffn,err_brr
         end if
cc         if (err_relres .lt. tol .and. err_rel .lt. tol) go to 333
         if (err_rel .lt. tol) go to 333
 200  continue
C
cc      write(*,*), ' Iteration limit exceeded:', maxit
 333  continue
      niter = min(kk,maxit)
      write(*,'(2X,i5,2X,5(3X,e11.4))')
     >     niter,err_relres,err_res,err_rel,xdiffn,err_brr
C
cc      write(*,*) '=============================', 
cc     >           '=================================================='
C
cc      arfac  = (err_relres)**(1./dble(kk))
cc      arfac1 = (xdiffn/err_res0)**(1./dble(kk))
cc      write(*,*) 
cc      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
cc      write(*,*)
C
 500  return
      end
C=====================================================================
      subroutine premg_s(m,ia,ja,idir,ipp,jpp,nzp,nz,nd,
     >     r,ku,kb,ka,lmx,lf,lc,lasti,lastr)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  PREMG_S stands for an multigrid algorithm as a preconditioner, in 
C...  other words, it is an iterator B in an MG-cycle. For the solver
C...  of linear systems various multigrid method will be used: 
C...  V-cycle, \-cycle, W-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle, (4 = W-cycle, 5 = full MG)
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    ISOLV   - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C---------------------------------------------------------------------
      icycle = ich(31)
      iprsm = ich(32)
      ipssm = ich(33)
      nprsm = ich(34)
      npssm = ich(35)
      lsolv = ich(38)
C
      call nullv(r(ku(lf)),nd(lf))
      do 10 i = lc, lf - 1
         call nullv(r(ku(i)),nd(i))
         call nullv(r(kb(i)),nd(i))
 10   continue
C
C...  V(\)-cycle of the MG, i.e., the action of the iterator B.
C...  Going downward in the V-cycle.
C
      icount = -1
C
      do 100 k = lf, lc+1, -1
         iak = ia(k)
         jak = ja(k)
         idc = idir(k-1)
         kuk = ku(k)
         kbc = kb(k-1)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         nzk = nz(k)
         ipk = ipp(k)
         jpk = jpp(k)
         nzpk  = nzp(k)
         iend  = lasti
         kwk1  = lastr
         kend  = kwk1 + ndk*7

cc         write(*,*)'iend, kend in down-cycle', iend, kend
C
C...     Presmoothing by Gauss-seidel smoother noofsm times: u=u+R(b-Au)
C
         if (icycle .eq. 3) then !For variable V-cyle
            icount = icount + 1
            noofsm = nprsm*2**icount
         else
            noofsm = nprsm      !For V or \-cyle
         end if
C
	 if (iprsm .eq. 1) then 
            call fwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (iprsm .eq. 2) then 
            call bwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (iprsm .eq. 3) then 
            call sym_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
         end if
C
C...     To compute the residual wk = b - Au.
C
         call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
C
C...     Restriction to the lower level: b = P^t*wk = wk^t*P.
C
         call rvbyp(m(ipk),m(jpk),r(kwk1),ndk,r(kbc),ndc)
         call zero_dir(r(kbc),m(idc),ndc)
 100  continue
C
C...  Solving exactly at the coarsest level, i.e., u = A^(-1)*b.
C
      itmax = 1000
      kuc   = ku(lc)
      iac   = ia(lc)
      jac   = ja(lc)
      kac   = ka(lc)
C
      if (lsolv .eq. 1) then
         call fwd_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax,0)
      else if (lsolv .eq. 2) then
         call bwd_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax,0)
      else if (lsolv .eq. 3) then
         call sym_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax,0)
      else if (lsolv .ge. 4) then
         maxa = iend
         kau  = kwk1
         call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
         call copyv(r(kbc),r(kuc),ndc)
      end if
C
C...  Going upward in the V-cycle.
C
      do 200 k = lc+1, lf
         iak = ia(k)
         jak = ja(k)
         idc = idir(k-1)
         kuc = ku(k-1)
         kuk = ku(k)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         nzk = nz(k)
         ipk = ipp(k)
         jpk = jpp(k)
         nzpk  = nzp(k)
         iend  = lasti
         kend  = lastr

cc         write(*,*)'iend, kend in up-cycle', iend, kend
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C
         call pbyvcs(m(ipk),m(jpk),r(kuc),ndc,r(kuk),ndk)
C
C...     Postsmoothing by Gauss-Seidel smoother noofsm times: u=u+R(b-Au)
C
         if (icycle .eq. 2)  go to 200        !No post smoothing in \-cycle
C
         if (icycle .eq. 3) then !For variable V-cyle
            noofsm = npssm*2**icount
            icount = icount - 1
         else
            noofsm = npssm      !For V or \-cyle
         end if
C
	 if (ipssm .eq. 1) then 
            call fwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (ipssm .eq. 2) then 
            call bwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (ipssm .eq. 3) then 
            call sym_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
         end if
 200  continue 
C
      return
      end
C=====================================================================
      subroutine mg_s_ord(m,iao,jao,idiro,ifree,r,ksol,krhs,kra,kfree,
     >     n,nx,ny,nnz,nel,lf,lc,tol) 
C=====================================================================
      parameter (lmx = 20)
      implicit real*8(a-h,o-z)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx),iord(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx),n1(lmx),n2(lmx)
      character*100 meshf(0:lmx)
      common /menu/ ich(200)
      common /fname/ meshf
C---------------------------------------------------------------------
C...  Standard multigrid cycle for the elliptic problems. Tarjan's 
C...  ordering is considered. The following finite element method can 
C...  be applied: Standard Galerkin, EAFE (Edge Average Finite Element),
C...  Streamline diffusion. For the solver of linear systems various
C...  multigrid method will be used: V-cycle, \-cycle, etc
C...
C...  Parameter:
C...    N,NX,NY - dimension of the mesh, N = NX x NY
C...    ICYCLE  - choice of MG:  1 = V-cycle, 2 = \-cycle, 
C...              3 = Variable V-cycle, (4 = W-cylce, 5 = Full MG)
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    MAXIT   - maximum number of MG iterations
C...    TOL     - stopping criterion for iteration procedure.
C...    ISTIFF  - choice of the stiffness matrix (iterator B) in MG.
C...              1 = Standard Galerkin,    2 = EAFE, 
C...              3 = Streamline diffusion, 4 = Upwind, 5 = SUPG.
C---------------------------------------------------------------------
      mt     = 15
      jdf    = ich(3)
      ipde   = ich(4)
      lread  = ich(5)
      icycle = ich(31)
      maxit  = ich(39)
      norder = ich(41)
      istiff = ich(40)      
C
cc      tol    = 0.5d-6
      zero   = 1.d-13
C
C...  Generating the stiffness matrix A and the load vector B at each level.
C
      lasti = ifree
      lastr = kfree
      nelold = nel
C
      nd(lf)   = n
      n1(lf)   = nx
      n2(lf)   = ny
      nz(lf)   = nnz
      ia(lf)   = iao
      ja(lf)   = jao
      idir(lf) = idiro
      ipp(lf)  = lasti
      lasti    = ipp(lf) + nd(lf)
C
      ku(lf)   = lastr
      kb(lf)   = ku(lf) + nd(lf)
      lastr    = kb(lf) + nd(lf)
C
      if (istiff .eq. ipde) then 
         ka(lf) = kra
      else 
         if (lread .eq. 1) then 
            open(mt,file=meshf(lf),status='unknown',form='formatted')
cc            rewind(mt)
            read(mt,*) n1(lf),n2(lf)
            read(mt,*) nel,nd(lf),ned
         else
            open(mt,file=meshf(lf),status='unknown',form='unformatted')
cc            rewind(mt)
            read(mt) n1(lf),n2(lf)
            read(mt) nel,nd(lf),ned
         end if
         write(*,*)
cc     >        ' No. of elements, nodes, & bdry edges:',nel,nd(lf),ned
C
         ie      = lasti
         je      = ie + nel + 1
         inedg   = je + nel*3
         iet     = inedg + ned * 2
         jet     = iet + nd(lf) + 1
         jat     = jet + nel*3
         iend    = jat + 2*(nel + nd(lf) + 5) + nd(lf)
C
         ka(lf)  = lastr
         lastr   = ka(lf) + 2*(nel + nd(lf) + 5) + nd(lf)
         kx      = lastr
         ky      = kx + nd(lf)
         krat    = ky + nd(lf)
         kend    = krat + 2*(nel + nd(lf) + 5) + nd(lf)
C
         call rdmesh(lread,mt,m(ie),m(je),ned,m(inedg),
     >        m(idir(lf)),r(kx),r(ky),nel,nd(lf),jdf)
         close(mt)
C
         call iit(m(ie),m(je),nel,nd(lf),m(iet),m(jet))
C
         nz(i) = 0
         call smbasg(m(ie),m(je),m(iet),m(jet),m(ia(lf)),m(ja(lf)),
     >        nd(lf),nz(lf),m(idir(lf)))
C
         iip  = iet
         call pde_choice(istiff,
     I        nel,nd(lf),jdf,m(ie),m(je),r(kx),r(ky),
     O        m(ia(lf)),m(ja(lf)),r(ka(lf)),nz(lf),r(kb(lf)),
     W        ned,m(inedg),m(idir(lf)),m(iip))
      end if
C
      do 50 i = lf, lc, -1
         if (i .lt. lf) then
            nel     = nel/4
            ia(i)   = lasti
            ja(i)   = ia(i) + nd(i) + 1
            idir(i) = ja(i) + 2*(nel + nd(i) + 5) + nd(i)
            ipp(i)  = idir(i) + nd(i) + 5
C
            ku(i)   = lastr
            kb(i)   = ku(i) + nd(i)
            ka(i)   = kb(i) + nd(i)
            lastr   = ka(i) + 2*(nel + nd(i) + 5) + nd(i)
         end if
C
         jpp(i)  = ipp(i) + nd(i) + 1
         iord(i) = jpp(i) + 2*nd(i)
         lasti   = iord(i) + nd(i) + 1
C
C...     To compute the symbolic prolongation matrix P from level i-1
C...     to level i.
C
         if (i .eq. lc) go to 40
         call prolng_symb(n1(i),n2(i),nd(i),n1(i-1),n2(i-1),nd(i-1),
     >        m(ipp(i)),m(jpp(i)),nzp(i))
 40      continue
C
C...     To compute stiffness matrices A at the lower level.
C
         if (i .lt. lf) then
            call formdir(m(idir(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),m(idir(i)),nd(i))
C
            ic = lasti
            ix = ic + nd(i+1) + 1
            jc = ix + nd(i) + 1
C
            call abybs(m(ia(i+1)),m(ja(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),nd(i+1),nd(i),m(ic),m(jc),m(ix))
C
            jce = m(ic+nd(i+1)) + 1
            kc  = lastr
            kx  = kc + jce
C
            call abyp(m(ia(i+1)),m(ja(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),m(ic),m(jc),r(ka(i+1)),r(kc),r(kx),nd(i+1))
C
            ict = jc + jce + 1
            jct = ict + nd(i+1) + 1
            kct = kx + nd(i+1) + 1
C
            call aat(m(ic),m(jc),r(kc),nd(i+1),nd(i),m(ict),m(jct),
     >           r(kct))
c
            call abybs(m(ict),m(jct),m(ipp(i+1)),m(jpp(i+1)),nd(i),
     >           nd(i+1),nd(i),m(ia(i)),m(ja(i)),m(ix))
C     
            call abyp(m(ict),m(jct),m(ipp(i+1)),m(jpp(i+1)),nd(i),
     >           m(ia(i)),m(ja(i)),r(kct),r(ka(i)),r(kx),nd(i+1))
C     
            nz(i) = m(ia(i)+nd(i))
C
            call adir_new(m(ia(i)),m(ja(i)),r(ka(i)),m(idir(i)),nd(i))
C
cc            call lprr(m(ia(i)),m(ja(i)),r(ka(i)),nd(i))
cc            call outmat1(m(ia(i)),m(ja(i)),r(ka(i)),nd(i),nz(i))
         end if
C
C...     To determine whether Tarjan's ordering is considered or not.
C...     If ordering is considered, then set NORDER = 1.
C
         if (norder .eq. 1) then
            iaw    = lasti
            jaw    = iaw + nd(i) + 1 + 1
            isub   = jaw + nz(i)
            iwork1 = isub + nd(i) + 1
            iwork2 = iwork1 + nd(i) + 1
            iwork3 = iwork2 + nd(i) + 1
C     
            call cut_off(
     >           nd(i),m(ia(i)),m(ja(i)),r(ka(i)),m(iaw),m(jaw),
     >           m(iord(i)),m(isub),nblk,m(iwork1),m(iwork2),m(iwork3))
         else
            call iseqv(m(iord(i)),nd(i))
         end if
C
cc         write(*,*)(ij,'*',m(iord(i)+ij-1),ij=1,nd(i))
cc         read(*,*) 
 50   continue
C
      nel = nelold
C
C...  ------------------- MG starts here. ---------------------------
C
      write(*,*)
      write(*,*) '=============================', 
     >           '=================================================='
      if (icycle .eq. 1) then
         write(*,*) '     METHOD :          MG : V-cycle '
      else if (icycle .eq. 2) then
         write(*,*) '     METHOD :          MG : Backslash cycle '
      else if (icycle .eq. 3) then
         write(*,*) '     METHOD :          MG : Variable V-cycle '
      end if
C
      write(*,*) '     Residual :        r := b-Au; '
      write(*,'(a,a,e12.4)')
     >     '      Convergence :     ||r||_0/||r_0||_0  &',
     >     '  ||Br||_0/||u||_0 <', tol
      write(*,'(a,2(5X,i3))') '      Coarse and fine levels : ',lc,lf
      write(*,*) '-----------------------------',
     >           '--------------------------------------------------'
      write(*,*) 'Iteration  ||r||/||r_0||    ||r||  ',
     >                  '   ||Br||/||u||      ||Br||   ||Br||/||r_0||'
      write(*,*) '-----------------------------',
     >           '--------------------------------------------------'
C
      kbf   = kb(lf)
      kuf   = ku(lf)
C
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
      call init_g0(m(iao),m(jao),r(kra),r(krhs),r(ksol),n)
C
cc      write(*,*) (r(ksol+ii-1), ii = 1,n)
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
C...  To compute the initial residual: b = rhs - A*sol.
      call abyvam(r(krhs),m(iao),m(jao),r(kra),r(ksol),n,r(kbf))
C
      call scpro(r(kbf),r(kbf),err_res0,n)
      err_res0 = dsqrt(err_res0)
C
      if (err_res0 .lt. zero)  err_res0 = 1.d0
C
      kk = 0
      write(*,'(2X,i5,2X,5(3X,e11.4))') 
     >     kk,1.d0,err_res0,1.d0,1.d0,1.d0
C
C...  Iterate until convergence:  sol = sol + B*b = sol + B(rhs - A*sol).
C
      do 200 kk = 1, maxit
C
C...     Iterator B using an MG-cycle, i.e., u = B*b = B(rhs - A*sol).
C
         call premg_ord(m(1),ia(1),ja(1),idir(1),
     >        ipp(1),jpp(1),nzp(1),iord(1),nz(1),nd(1),
     >        r(1),ku(1),kb(1),ka(1),lmx,lf,lc,lasti,lastr)
C
C...     New solution:  sol = sol + u = sol + B(rhs - A*sol).
         call uupluv(r(ksol),r(kuf),n)
C
C...     To compute the relative error.
	 call scpro(r(kuf),r(kuf),xdiffn,n)
         call scpro(r(ksol),r(ksol),xnewn,n)
         xdiffn  = dsqrt(xdiffn)
         xnewn   = dsqrt(xnewn)
C
         if (xnewn .gt. xdiffn) then
            err_rel = xdiffn/xnewn
         else
            err_rel = xdiffn  
         end if
         err_brr = xdiffn/err_res0
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(r(krhs),m(iao),m(jao),r(kra),r(ksol),n,r(kbf))
C
C...     To compute the L_2 norm of the residual.
         call scpro(r(kbf),r(kbf),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
         nmod = mod(kk,3)
C
         if (nmod .eq. 1) then
            write(*,'(2X,i5,2X,5(3X,e11.4))')
     >           kk,err_relres,err_res,err_rel,xdiffn,err_brr
         end if
         if (err_relres .lt. tol .and. err_rel .lt. tol) go to 333
 200  continue
C
cc      write(*,*), ' Iteration limit exceeded:', maxit
 333  continue
      niter = min(kk,maxit)
      write(*,'(2X,i5,2X,5(3X,e11.4))')
     >     niter,err_relres,err_res,err_rel,xdiffn,err_brr
C
cc      write(*,*) '=============================', 
cc     >           '=================================================='
C
cc      arfac  = (err_relres)**(1./dble(kk))
cc      arfac1 = (xdiffn/err_res0)**(1./dble(kk))
cc      write(*,*) 
cc      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
cc      write(*,*)
C
 500  return
      end
C=====================================================================
      subroutine premg_ord(m,ia,ja,idir,ipp,jpp,nzp,iord,nz,nd,
     >     r,ku,kb,ka,lmx,lf,lc,lasti,lastr)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx),iord(lmx)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  PREMG stands for an multigrid algorithm as a preconditioner, in 
C...  other words, it is an iterator B in an MG-cycle. The following 
C...  finite element method can be applied: Standard Galerkin, EAFE 
C...  (Edge Average Finite Element), Streamline diffusion. For the 
C...  solver of linear systems various multigrid method will be used: 
C...  V-cycle, \-cycle, W-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle, (4 = W-cylce, 5 = Full MG)
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    ISOLV   - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...    ITERB   - choice of matrix A in the computation of an iterator
C...              B in MG iterations: 1 = EAFE, 2 = Standard Galerkin, 
C...              3 = Streamline diffusion
C---------------------------------------------------------------------
      icycle   = ich(31)
      iprsm = ich(32)
      ipssm = ich(33)
      nprsm = ich(34)
      npssm = ich(35)
      lsolv = ich(38)
      iterb = ich(40)
C
      call nullv(r(ku(lf)),nd(lf))
      do 10 i = lc, lf - 1
         call nullv(r(ku(i)),nd(i))
         call nullv(r(kb(i)),nd(i))
 10   continue
C
C...  V(\)-cycle of the MG, i.e., the action of the iterator B.
C...  Going downward in the V-cycle.
C
      icount = -1
C
      do 100 k = lf, lc+1, -1
         iak = ia(k)
         jak = ja(k)
         idc = idir(k-1)
         kuk = ku(k)
         kbc = kb(k-1)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         nzk = nz(k)
         ipk = ipp(k)
         jpk = jpp(k)
         nzpk  = nzp(k)
         iordk = iord(k)
         iend  = lasti
         kwk1  = lastr
         kend  = kwk1 + ndk*7

cc         write(*,*)'iend, kend in down-cycle', iend, kend
C
C...     Presmoothing by Gauss-seidel smoother noofsm times: u=u+R(b-Au)
C
         if (icycle .eq. 3) then !For variable V-cyle
            icount = icount + 1
            noofsm = nprsm*2**icount
         else
            noofsm = nprsm      !For V or \-cyle
         end if
C
	 if (iprsm .eq. 1) then 
            call fwd_gs_ns_ord_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (iprsm .eq. 2) then 
            call bwd_gs_ns_ord_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (iprsm .eq. 3) then 
            call sym_gs_ns_ord_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
         end if
C
C...     To compute the residual wk = b - Au.
C
         call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
C
C...     Restriction to the lower level: b = P^t*wk = wk^t*P.
C
         call rvbyp(m(ipk),m(jpk),r(kwk1),ndk,r(kbc),ndc)
         call zero_dir(r(kbc),m(idc),ndc)
C
 100  continue
C
C...  Solving exactly at the coarsest level, i.e., u = A^(-1)*b.
C
      itmax = 1000
      kuc   = ku(lc)
      iac   = ia(lc)
      jac   = ja(lc)
      kac   = ka(lc)
      iordc = iord(lc)
C
      if (lsolv .eq. 1) then
         call fwd_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 2) then
         call bwd_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 3) then
         call sym_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .ge. 4) then
         maxa = iend
         kau  = kwk1
         call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
         call copyv(r(kbc),r(kuc),ndc)
      end if
C
C...  Going upward in the V-cycle.
C
      do 200 k = lc+1, lf
         iak = ia(k)
         jak = ja(k)
         idc = idir(k-1)
         kuc = ku(k-1)
         kuk = ku(k)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         nzk = nz(k)
         ipk = ipp(k)
         jpk = jpp(k)
         nzpk  = nzp(k)
         iordk = iord(k)
         iend  = lasti
         kend  = lastr

cc         write(*,*)'iend, kend in up-cycle', iend, kend
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C
         call pbyvcs(m(ipk),m(jpk),r(kuc),ndc,r(kuk),ndk)
C
C...     Postsmoothing by Gauss-Seidel smoother noofsm times: u=u+R(b-Au)
C
         if (icycle .eq. 2)  go to 200        !No post smoothing in \-cycle
C
         if (icycle .eq. 3) then !For variable V-cyle
            noofsm = npssm*2**icount
            icount = icount - 1
         else
            noofsm = npssm      !For V or \-cyle
         end if
C
	 if (ipssm .eq. 1) then 
            call fwd_gs_ns_ord_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipssm .eq. 2) then 
            call bwd_gs_ns_ord_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipssm .eq. 3) then 
            call sym_gs_ns_ord_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
         end if
C
 200  continue 
C
      return
      end
C=====================================================================
