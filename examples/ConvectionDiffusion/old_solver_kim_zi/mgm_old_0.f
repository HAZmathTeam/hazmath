C=====================================================================
      subroutine mg(iao,jao,idiro,m,sol,rhs,ra,r,n,nnz,lf,lc) 
C=====================================================================
      parameter (lmx = 20)
      implicit real*8(a-h,o-z)
      dimension iao(1),jao(1),idiro(1),sol(1),rhs(1),ra(1)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx),iord(lmx),iperm(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx),n1(lmx),n2(lmx)
      character*100 meshf(0:lmx)
      common /menu/ ich(200)
      common /fname/ meshf
C---------------------------------------------------------------------
C...  This program solves the 2nd order linear PDEs which are possibly 
C...  convection dominated. The following finite element method can be 
C...  applied: Standard Galerkin, EAFE (Edge Average Finite Element),
C...  Streamline diffusion. For the solver of linear systems various
C...  multigrid method will be used: V-cycle, \-cycle, W-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE - choice of MG:  1 = V-cycle, 2 = \-cycle,
C...             3 = Variable V-cycle,(4 = W-cylce, 5 = Full MG)
C...    LF     - level of the finest mesh: 2,3,4, etc
C...    LC     - level of the coarsest mesh: 1,2,3, etc
C...    MAXIT  - maximum number of MG iterations
C...    TOL    - stopping iteration procedure.
C...    ISTIFF - choice of the stiffness matrix (iterator B) in MG.
C...             1 = Standard Galerkin,    2 = EAFE, 
C...             3 = Streamline diffusion, 4 = Upwind, 5 = SUPG.
C---------------------------------------------------------------------
      mt     = 15
      jdf    = ich(3)
      lread  = ich(5)
      icycle = ich(31)
      maxit  = ich(39)
      norder = ich(41)
      istiff = ich(40)
C
      tol    = 0.5d-6
      zero   = 1.d-13
C
C...  Generating the stiffness matrix A and the load vector B at each level.
C
      lasti = 1
      lastr = 1
cc      do 50 i = lc, lf
      do 50 i = lf, lc, -1
         if(lread .eq. 1) then
            open(mt,file=meshf(i),status='unknown',form='formatted')
cc            rewind(mt)
            read(mt,*) n1(i),n2(i)
            read(mt,*) nel,nd(i),ned
         else
            open(mt,file=meshf(i),status='unknown',form='unformatted')
cc            rewind(mt)
            read(mt) n1(i),n2(i)
            read(mt) nel,nd(i),ned
         end if
         write(*,*)
     >        ' No. of elements, nodes, & bdry edges:',nel,nd(i),ned
C
C...     Compute the addresses of each array in the arrays m and r.
C
         ia(i)   = lasti
         ja(i)   = ia(i) + nd(i) + 1
         idir(i) = ja(i) + 2*(nel + nd(i) + 5) + nd(i)
         ipp(i)  = idir(i) + nd(i) + 5
         jpp(i)  = ipp(i) + nd(i) + 1
         iord(i) = jpp(i) + 2*nd(i)
         iperm(i)= iord(i) + nd(i) + 1
         lasti   = iperm(i) + nd(i) + 1
         ie      = lasti
         je      = ie + nel + 1
         inedg   = je + nel*3
         iet     = inedg + ned * 2
         jet     = iet + nd(i) + 1
         jat     = jet + nel*3
         iend    = jat + 2*(nel + nd(i) + 5) + nd(i)
C
         ku(i)   = lastr
         kb(i)   = ku(i) + nd(i)
         ka(i)   = kb(i) + nd(i)
         lastr   = ka(i) + 2*(nel + nd(i) + 5) + nd(i)
         kx      = lastr
         ky      = kx + nd(i)
         krat    = ky + nd(i)
         kend    = krat + 2*(nel + nd(i) + 5) + nd(i)

cc         write(*,*) 'lasti,lastr,iend,kend', lasti,lastr,iend,kend
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
         call rdmesh_perm(lread,mt,m(ie),m(je),ned,m(inedg),
     >        m(idir(i)),m(iperm(i)),r(kx),r(ky),nel,nd(i),jdf)
         close(mt)
C
         call iit(m(ie),m(je),nel,nd(i),m(iet),m(jet))
C
         nz(i) = 0
         call smbasg(m(ie),m(je),m(iet),m(jet),m(ia(i)),m(ja(i)),
     >        nd(i),nz(i),m(idir(i)))
C
         iip    = iet
         call pde_choice(istiff,
     I        nel,nd(i),jdf,m(ie),m(je),r(kx),r(ky),
     O        m(ia(i)),m(ja(i)),r(ka(i)),nz(i),r(kb(i)),
     W        ned,m(inedg),m(idir(i)),m(iip))
         
cc         call outmat1(m(ia(i)),m(ja(i)),r(ka(i)),nd(i),nz(i))
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
C...  ------------------- MG starts here. ---------------------------
C
      write(*,*)
      write(*,'(a,a)') ' =============================',
     >      '=================================================='
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
      write(*,'(a,a)') ' -----------------------------',
     >           '--------------------------------------------------'
      write(*,*) 'Iteration  ||r||/||r_0||    ||r||  ',
     >                  '   ||Br||/||u||      ||Br||   ||Br||/||r_0||'
      write(*,'(a,a)') ' -----------------------------',
     >           '--------------------------------------------------'
C
      kbf   = kb(lf)
      kuf   = ku(lf)
C
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
      call init_g0(iao,jao,ra,rhs,sol,n)
C
cc      write(*,*) (r(ksol+ii-1), ii = 1,n)
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
C...  To compute the initial residual: b = rhs - A*sol.
      call abyvam(rhs,iao,jao,ra,sol,n,r(kbf))
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
         call premg(m(1),ia(1),ja(1),idir(1),iperm(1),
     >        ipp(1),jpp(1),nzp(1),iord(1),nz(1),nd(1),
     >        r(1),ku(1),kb(1),ka(1),lmx,lf,lc,lasti,lastr)
C
C...     New solution:  sol = sol + u = sol + B(rhs - A*sol).
         call uupluv(sol,r(kuf),n)
C
C...     To compute the relative error.
	 call scpro(r(kuf),r(kuf),xdiffn,n)
         call scpro(sol,sol,xnewn,n)
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
         call abyvam(rhs,iao,jao,ra,sol,n,r(kbf))
C
C...     To compute the residual error.
         call scpro(r(kbf),r(kbf),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
         if (lf .lt. 8) then
            nmod = mod(kk,100)
         else if (lf .lt. 10) then
            nmod = mod(kk,10)
         else 
            nmod = 1 
         end if
C
         if (nmod .eq. 1) then
            write(*,'(2X,i5,2X,5(3X,e11.4))')
     >           kk,err_relres,err_res,err_rel,xdiffn,err_brr
         end if
         if (err_relres .lt. tol .and. err_rel .lt. tol) go to 333
 200  continue
C
      write(*,*), ' Iteration limit exceeded:', maxit
 333  continue
      write(*,'(2X,i5,2X,5(3X,e11.4))')
     >     kk,err_relres,err_res,err_rel,xdiffn,err_brr
C
      write(*,'(a,a)') ' =============================',
     >      '=================================================='
C
      arfac  = (err_relres)**(1./dble(kk))
      arfac1 = (xdiffn/err_res0)**(1./dble(kk))
      write(*,*) 
      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
C
      niter = min(kk, maxit)
      write(*,*)
C
 500  return
      end
C=====================================================================
      subroutine premg(m,ia,ja,idir,iperm,ipp,jpp,nzp,iord,nz,nd,
     >     r,ku,kb,ka,lmx,lf,lc,lasti,lastr)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx),iord(lmx),iperm(lmx)
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
C...    ITERB   - choice of matrix A in the computation of an iterator
C...              B in MG iterations: 1 = EAFE, 2 = Standard Galerkin, 
C...              3 = Streamline diffusion
C---------------------------------------------------------------------
      icycle = ich(31)
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
         kwk2  = kwk1 + ndk + 1
         kend  = kwk2 + ndk*7

cc         write(*,*)'iend, kend in down-cycle', iend, kend
C
C...     Presmoothing by Gauss-seidel smoother nprsm times: u=u+R(b-Au)
C
         if (icycle .eq. 3) then !For variable V-cyle
            icount = icount + 1
            noofsm = nprsm*2**icount
         else
            noofsm = nprsm      !For V or \-cyle
         end if

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
C...     In order to perform the restriction on the uniform ordering,
C...     use IPERM to permute vectors back and forth.
C
         call pervec(m(iperm(k)),r(kwk1),r(kwk2),ndk)
         call rvbyp(m(ipk),m(jpk),r(kwk1),ndk,r(kbc),ndc)
         call perback(m(iperm(k-1)),r(kbc),r(kwk2),ndc)
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
      else if (lsolv .eq. 4) then
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
         kwk1  = lastr
         kwk2  = kwk1 + ndk + 1
         kend  = kwk2 + ndk + 1

cc         write(*,*)'iend, kend in up-cycle', iend, kend
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C...     In order to perform the prolongation on the uniform ordering,
C...     use IPERM to permute vectors back and forth.
C
C...     Either do this:
         call pervec(m(iperm(k-1)),r(kuc),r(kwk1),ndc)
         call pbyvg(m(ipk),m(jpk),r(kuc),ndc,r(kwk1),ndk)
         call perback(m(iperm(k)),r(kwk1),r(kwk2),ndk)
         call uupluv(r(kuk),r(kwk1),ndk)
C
C...     Or do the following:
cc         call pervec(m(iperm(k)),r(kuk),r(kwk1),ndk)
cc         call pervec(m(iperm(k-1)),r(kuc),r(kwk1),ndc)
cc         call pbyvcs(m(ipk),m(jpk),r(kuc),ndc,r(kuk),ndk)
cc         call perback(m(iperm(k)),r(kuk),r(kwk1),ndk)
C
C...     Postsmoothing by Gauss-Seidel smoother npssm times: u=u+R(b-Au)
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





