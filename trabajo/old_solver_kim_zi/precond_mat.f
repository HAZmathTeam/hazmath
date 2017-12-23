C=====================================================================
      subroutine mg(m,iao,jao,idiro,r,sol,rhs,ao,
     >     n,nnz,maxit,tol,ex_sol,iexact,ipre_normal) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1)
      dimension r(1)
C
      integer iao(1),jao(1),idiro(1)
      double precision sol(1),rhs(1),ao(1),ex_sol(1)
C
      parameter (lmax = 20)
      integer lmax
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ifree, kfree, lf, lc
C
      character*100 meshf(0:lmax)
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
C...             3 = Streamline diffusion, 4 = Upwind, 5 = SUPG,
C...             12 = double discretization: ich(4) & 2,
C...             13 = double discretization: ich(4) & 3,
C...             14 = double discretization: ich(4) & 4,
C...             15 = double discretization: ich(4) & 5.
C...    IEXACT - To determine the stoppping criterion;
C...             0,1,3 - compute the residual error: ||r||/||rhs||
C...                     where r = rhs - A*sol.
C...             2     - compute the relative error using the computed 
C...                     exact soln, i.e., ||sol - ex_sol|| / ||ex_sol||;
C...    IPRE_NORMAL - To detemine whether do smoothings by G-S sweeps 
C...                  for the normal equation: A^t*A*sol = A^t*rhs.
C...             0 - do not;
C...             1 - do;  in this case set the values: MAX_SWEEPS, IGS 
C---------------------------------------------------------------------
      mt     = 15
      jdf    = ich(3)
      lread  = ich(5)
      icycle = ich(31)
cc      maxit  = ich(39)
      norder = ich(41)
      istiff = ich(40)
C
cc      tol    = 0.5d-6
      zero   = 1.d-15
C
C...  Generating the stiffness matrix A and the load vector B at each level.
C
      lasti = ifree
      lastr = kfree
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
         if (i .eq. lf) then
            ipp(i)  = lasti
         else
            ia(i)   = lasti
            ja(i)   = ia(i) + nd(i) + 1
            idir(i) = ja(i) + 2*(nel + nd(i) + 5) + nd(i)
            ipp(i)  = idir(i) + nd(i) + 5
         end if
C
         jpp(i)  = ipp(i) + nd(i) + 1
         iord(i) = jpp(i) + 2*nd(i)
         iperm(i)= iord(i) + nd(i) + 1
         lasti   = iperm(i) + nd(i) + 1
C
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
C
         if (i .lt. lf .and. istiff .gt. 10) then
            kag(i)  = ka(i) + 2*(nel + nd(i) + 5) + nd(i)
            lastr   = kag(i) + 2*(nel + nd(i) + 5) + nd(i)
         else
            lastr   = ka(i) + 2*(nel + nd(i) + 5) + nd(i)
         end if
C
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
         if (i .lt. lf) then
            nz(i) = 0
            call smbasg(m(ie),m(je),m(iet),m(jet),m(ia(i)),m(ja(i)),
     >           nd(i),nz(i),m(idir(i)))
         end if
C
         iip  = iet
         ipde = mod(istiff, 10)
         call pde_choice(ipde,
     I        nel,nd(i),jdf,m(ie),m(je),r(kx),r(ky),
     O        m(ia(i)),m(ja(i)),r(ka(i)),nz(i),r(kb(i)),
     W        ned,m(inedg),m(idir(i)),m(iip))
         if (istiff .gt. 10 .and. i .lt. lf) then
            lmethod = ich(4)
            if (lmethod .le. 5) then
               ipde = lmethod
            else if (lmethod .eq. 6) then
               ipde = 3
            end if
            write(*,*) 'ipde', ipde
            call pde_choice(ipde,
     I           nel,nd(i),jdf,m(ie),m(je),r(kx),r(ky),
     O           m(ia(i)),m(ja(i)),r(kag(i)),nz(i),r(kb(i)),
     W           ned,m(inedg),m(idir(i)),m(iip))
         end if
C
cc         write(*,*) '============== Discretization 11111111 ==========='
cc         call outmat1(m(ia(i)),m(ja(i)),r(kag(i)),nd(i),nz(i))
C
cc         write(*,*) '============== Discretization 22222222 ==========='
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
cc      write(*,*) 'Check point 0000'
C...  G-S sweeps with normal equation if ipre_normal = 1.
C
      if (ipre_normal .eq. 0) go to 150
C
C...  Find C := A^t*A  &  bta := (rhs^t)*A = A^t*rhs so that G-S sweeps
C...  can be applied to the normal equation: A^t*A*sol = A^t*rhs.
C
      ic    = lasti
      iwk   = ic + n + 1
      iat   = iwk + n + 1 
      jat   = iat + n + 1
      jc    = jat + nnz + 1
C
      kbta  = lastr
      kwk   = kbta + n + 1
      kat   = kwk + n + 1
      kc    = kat + nnz + 1
C
cc      write(*,*) 'Check point 1111'
C...  Compute bta = (rhs^t)*A = A^t*rhs.
      call vbya(iao,jao,ao,rhs,n,n,r(kbta))
C
cc      write(*,*) 'Check point 2222'
C...  Compute at = A^t.
      call aat(iao,jao,ao,n,n,m(iat),m(jat),r(kat))
C
cc      write(*,*) 'Check point 3333'
C...  Compute C = A^t*A symbolically.
      call abybs(m(iat),m(jat),iao,jao,n,n,n,m(ic),m(jc),m(iwk))
C
cc      write(*,*) 'Check point 4444'
C...  Compute C = A^t*A numerically.
      call abyb(m(iat),m(jat),iao,jao,n,m(ic),m(jc),
     >     r(kat),ao,r(kc),r(kwk),n)
C
      nnz_c = m(ic+n) - 1
C
cc      write(*,*) 'Check point 5555'
C...  Move some arrays to the free spaces in order to save spaces.
      call icopyv(m(jc),m(iwk),nnz_c)
      call copyv(r(kc),r(kwk),nnz_c)
C
      jc = iwk
      kc = kwk
      lasti = jc + nnz_c + 1
      lastr = kc + nnz_c + 1
C
 150  continue
C
      arfac1 = (err_rel)**(1./dble(kk))
      write(*,*) 
      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
C
      write(*,*)
C
 500  return
      end
C=====================================================================
      subroutine premg(m,r,lasti,lastr)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1)
      dimension r(1)
C
      parameter (lmax = 20)
      integer lmax
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ifree, kfree, lf, lc
C
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
C...    ISTIFF  - choice of the stiffness matrix (iterator B) in MG.
C...              1 = Standard Galerkin,    2 = EAFE, 
C...              3 = Streamline diffusion, 4 = Upwind, 5 = SUPG,
C...              12 = double discretization: ich(4) & 2,
C...              13 = double discretization: ich(4) & 3,
C...              14 = double discretization: ich(4) & 4,
C...              15 = double discretization: ich(4) & 5.
C...    LACST   - choice of the stiffness matrix at the coarsest level
C...              1 = first order stable scheme, e.g., EAFE
C...              2 = second order (unstable) scheme, e.g., Galerkin
C---------------------------------------------------------------------
      icycle = ich(31)
      iprsm  = ich(32)
      ipssm  = ich(33)
      nprsm  = ich(34)
      npssm  = ich(35)
      lsolv  = ich(38)
      istiff = ich(40)
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
         iak  = ia(k)
         jak  = ja(k)
         idc  = idir(k-1)
         kuk  = ku(k)
         kbc  = kb(k-1)
         kbk  = kb(k)
         kak  = ka(k)
         kagk = kag(k)
         ndc  = nd(k-1)
         ndk  = nd(k)
         nzk  = nz(k)
         ipk  = ipp(k)
         jpk  = jpp(k)
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
         if (istiff .le. 10) then
            call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
         else
            call abyvam(r(kbk),m(iak),m(jak),r(kagk),r(kuk),ndk,r(kwk1))
         end if
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
      lacst = 1
      itmax = 1000
      kuc   = ku(lc)
      iac   = ia(lc)
      jac   = ja(lc)
      kac   = ka(lc)
      kagc  = kag(lc)
      iordc = iord(lc)
C
      if (lsolv .eq. 11) then
         call fwd_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 12) then
         call bwd_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 13) then
         call sym_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 14) then
         maxa = iend
         kau  = kwk1
         call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
         call copyv(r(kbc),r(kuc),ndc)
      else if (lsolv .eq. 21) then
         call fwd_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kagc),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 22) then
         call bwd_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kagc),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 23) then
         call sym_gs_ns_ord_wr(r(kuc),m(iac),m(jac),r(kagc),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv .eq. 24) then
         maxa = iend
         kau  = kwk1
         call sgauss(m(iac),m(jac),r(kagc),r(kbc),m(maxa),
     >        r(kau),ndc)
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





