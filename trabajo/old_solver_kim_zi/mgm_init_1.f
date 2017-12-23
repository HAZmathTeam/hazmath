C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     The only difference between MGM_INIT.F and MGM.F is that the 
C     subroutine MGM_INIT.F initializes the solution (SOL) inside the
C     the subroutine. On the other hand, the subroutine MGM.F does not 
C     initialize the solution (SOL) inside so that the solution, SOL,
C     should be initialized before calling the subroutine.
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

C=====================================================================
      subroutine mg(m,iao,jao,idiro,r,sol,rhs,ao,
     >     n,nnz,maxit,tol,iexact) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1)
      dimension r(1)
C
      integer iao(1),jao(1),idiro(1)
      double precision sol(1),rhs(1),ao(1)
C
      parameter (lmax = 20)
      integer lmax
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
C
      double precision soln(1)
      COMMON /DSLBLK/ SOLN
C
      character*100 meshf(0:lmax)
      common /menu/ ich(200)
      common /fname/ meshf
      common /prenormal/ ipre_normal
C---------------------------------------------------------------------
C...  This subroutine is identical to the subroutine PREMG_NORMAL.
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
C...                     exact soln, i.e., ||sol - soln|| / ||soln||;
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
C...  Find C := A^t*A so that G-S sweeps can be applied to the 
C...  normal equation: A^t*A*sol = A^t*rhs.
C
      ic    = lasti
      iwk   = ic + n + 1
      iat   = iwk + n + 1 
      jat   = iat + n + 1
      jc    = jat + nnz + 1
C
      kwk   = lastr
      kat   = kwk + n + 1
      kc    = kat + nnz + 1
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
C...  ------------------- MG_normal starts here. ---------------------------
C
      write(*,*)
      write(*,'(a,a)') ' =============================',
     >      '=================================================='
      if (icycle .eq. 1) then
         write(*,*) '   METHOD :  MG or MG_normal : V-cycle '
      else if (icycle .eq. 2) then
         write(*,*) '   METHOD :  MG or MG_normal : Backslash cycle '
      else if (icycle .eq. 3) then
         write(*,*) '   METHOD :  MG or MG_normal : Variable V-cycle '
      end if
C
      write(*,*) '     Residual :        r := b-Au; '
      write(*,'(a,a,e12.4)')
     >     '      Convergence :     ||r||_0/||b||_0  and/or',
     >     '  ||u - soln||_0/||soln||_0 <', tol
      write(*,'(a,2(5X,i3))') '      Coarse and fine levels : ',lc,lf
      write(*,'(a,a)') ' -----------------------------',
     >           '--------------------------------------------------'
      write(*,*) 'Iteration  ||r||/||b||      ||r||  ',
     >                  ' ||u-soln||/||soln||   ||Br||   |Br||/||Bb||'
      write(*,'(a,a)') ' -----------------------------',
     >           '--------------------------------------------------'
C
C...  Compute the L^2 norm of the computed exact solution if it is given.
C
      if (iexact .eq. 2) then
         call l2norm(soln,soln_norm,n,1d0)
         write(*,*) ' soln_norm: ', soln_norm
         if (soln_norm .le. zero) then
            soln_norm = 1.d0
            write(*,*) 'Zero exact solution is detected!!!!.'
         end if
      end if
C
      kbf   = kb(lf)
      kuf   = ku(lf)
C
cc      write(*,*)
cc      write(*,*) (rhs(ii), ii = 1,n)
C
      call init_g0(iao,jao,ao,rhs,sol,n)
C
cc      write(*,*) (sol(ii), ii = 1,n)
cc      write(*,*)
cc      write(*,*) (rhs(ii), ii = 1,n)
C
C...  To compute the initial residual: b = rhs - A*sol.
      call abyvam(rhs,iao,jao,ao,sol,n,r(kbf))
C
ccc      call scpro(rhs,rhs,err_res0,n)
      call scpro(r(kbf),r(kbf),err_res0,n)
      err_res0 = dsqrt(err_res0)
C
      if (err_res0 .le. zero) then
         err_res0 = 1.d0
         write(*,*) 'Zero right hand side is detected!!!!.'
      end if
C
      err_relres = err_res0
      err_rel = 0.d0
C
      kk = 0
      write(*,'(2X,i5,2X,3(3X,e11.4))') 
     >     kk,1.d0,err_res0,1.d0
C
      if (ipre_normal .eq. 0) go to 155
C
      kbta  = lastr
      lastr = kbta + n + 1
C
C...  Compute bta = (rhs^t)*A = A^t*rhs.
      call vbya(iao,jao,ao,rhs,n,n,r(kbta))
C
 155  continue
C
C...  Iterate until convergence:  sol = sol + B*b = sol + B(rhs - A*sol).
C
      do 200 kk = 1, maxit
C
         if (ipre_normal .eq. 0) go to 160
C
C...     Smoothing by G-S sweeps for the normal equation: A^t*A*sol = A^t*rhs.
C
         igs        = ich(42)
         max_sweeps = ich(43)
C
         if (igs .eq. 1) then
            call fwd_gs_ns_wr(sol,m(ic),m(jc),r(kc),r(kbta),n,
     >           max_sweeps,0)
         else if(igs .eq. 2) then
            call bwd_gs_ns_wr(sol,m(ic),m(jc),r(kc),r(kbta),n,
     >           max_sweeps,0) 
         else if(igs .eq. 3) then
            call sym_gs_ns_wr(sol,m(ic),m(jc),r(kc),r(kbta),n,
     >           max_sweeps,0)
         else
            write(*,*) ' Error in the choice of G-S! '
         end if
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(rhs,iao,jao,ao,sol,n,r(kbf))
C
 160     continue
C
C...     Iterator B using an MG-cycle, i.e., u = B*b = B(rhs - A*sol).
C
         call premg(m(1),r(1),lasti,lastr)
C
C...     New solution:  sol = sol + u = sol + B(rhs - A*sol).
         call uupluv(sol,r(kuf),n)
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(rhs,iao,jao,ao,sol,n,r(kbf))
C
C...     To compute the residual error.
         call scpro(r(kbf),r(kbf),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
C...     To compute the relative error if the exact solution is known.
         if (iexact .eq. 2) then
            call wuminv(r(kuf),soln,sol,n)
            call l2norm(r(kuf),xdiffn,n,1d0)
            err_rel = xdiffn/soln_norm  
         end if
C
         if (lf .lt. 3) then
            nmod = mod(kk,100)
         else if (lf .lt. 3) then
            nmod = mod(kk,10)
         else 
            nmod = 1 
         end if
C
         if (nmod .eq. 1) then
            write(*,'(2X,i5,2X,3(3X,e11.4))')
     >           kk,err_relres,err_res,err_rel
         end if
C
         if  (iexact .eq. 2) then
            if (err_rel .lt. tol) go to 333
         else
            if (err_relres .lt. tol .and. err_rel .lt. tol) go to 333
         end if
 200  continue
C
      write(*,*), ' Iteration limit exceeded:', maxit
 333  continue
      niter = min(kk, maxit)
      write(*,'(2X,i5,2X,3(3X,e11.4))')
     >     niter,err_relres,err_res,err_rel
C
      write(*,'(a,a)') ' =============================',
     >      '=================================================='
C
      arfac  = (err_relres)**(1./dble(kk))
      arfac1 = (err_rel)**(1./dble(kk))
      write(*,*) 
      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
C
      write(*,*)
C
 500  return
      end
C=====================================================================
      subroutine precond_mat(m,iao,jao,idiro,r,ao,n,nnz) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1)
      dimension r(1)
C
      integer iao(1),jao(1),idiro(1)
      double precision ao(1)
C
      parameter (lmax = 20)
      integer lmax
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord,
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
C
      character*100 meshf(0:lmax)
      common /menu/ ich(200)
      common /fname/ meshf
      common /prenormal/ ipre_normal
C---------------------------------------------------------------------
C...  Generate preconditioning matrices and other information. 
C...  The following finite element method can be applied: 
C...  Standard Galerkin, EAFE (Edge Average Finite Element).
C...  And the normal equation can be considered for smoothing purpose: 
C...  C*u = bta, where C = A^t*A and bta = A^t*rhs.
C...
C...  Parameter:
C...    LF     - level of the finest mesh: 2,3,4, etc
C...    LC     - level of the coarsest mesh: 1,2,3, etc
C...    ISTIFF - choice of the stiffness matrix (iterator B) in MG.
C...             1 = Standard Galerkin,    2 = EAFE, 
C...             3 = Streamline diffusion, 4 = Upwind, 5 = SUPG,
C...             12 = double discretization: ich(4) & 2,
C...             13 = double discretization: ich(4) & 3,
C...             14 = double discretization: ich(4) & 4,
C...             15 = double discretization: ich(4) & 5.
C...    IPRE_NORMAL - To detemine whether do smoothings by G-S sweeps 
C...                  for the normal equation: A^t*A*sol = A^t*rhs.
C...             0 - do not;
C...             1 - do;  in this case set the values: MAX_SWEEPS, IGS 
C---------------------------------------------------------------------
      mt     = 15
      jdf    = ich(3)
      lread  = ich(5)
      norder = ich(41)
      istiff = ich(40)
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
C...  Find C := A^t*A  so that G-S sweeps can be applied to the 
C...  normal equation: A^t*A*sol = A^t*rhs.
C
      ic    = lasti
      iwk   = ic + n + 1
      iat   = iwk + n + 1 
      jat   = iat + n + 1
      jc    = jat + nnz + 1
C
      kwk   = lastr
      kat   = kwk + n + 1
      kc    = kat + nnz + 1
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
      ifree = lasti
      kfree = lastr
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
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
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
      subroutine premg_normal(m,iao,jao,idiro,r,sol,rhs,ao,
     >     n,nnz,maxit,tol,iexact,msolve) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1)
      dimension r(1)
C
      integer iao(1),jao(1),idiro(1)
      double precision sol(1),rhs(1),ao(1)
C
      parameter (lmax = 20)
      integer lmax
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
C
      external msolve
C
      double precision soln(1)
      COMMON /DSLBLK/ SOLN
C
      character*100 meshf(0:lmax)
      common /menu/ ich(200)
      common /fname/ meshf
      common /prenormal/ ipre_normal
C---------------------------------------------------------------------
C...  This subroutine is identical to the subroutine MG.
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
C...                     exact soln, i.e., ||sol - soln|| / ||soln||;
C...    IPRE_NORMAL - To detemine whether do smoothings by G-S sweeps 
C...                  for the normal equation: A^t*A*sol = A^t*rhs.
C...             0 - do not;
C...             1 - do;  in this case set the values: MAX_SWEEPS, IGS 
C---------------------------------------------------------------------
      icycle = ich(31)
      zero   = 1.d-15
C
C...  Get the information about preconditioning matrices.
      call precond_mat(m,iao,jao,idiro,r,ao,n,nnz) 
C
      lasti = ifree
      lastr = kfree
C
C...  ------------------- MG_normal starts here. ---------------------------
C
      write(*,*)
      write(*,'(a,a)') ' =============================',
     >      '=================================================='
      if (icycle .eq. 1) then
         write(*,*) '   METHOD :  MG or MG_normal : V-cycle '
      else if (icycle .eq. 2) then
         write(*,*) '   METHOD :  MG or MG_normal : Backslash cycle '
      else if (icycle .eq. 3) then
         write(*,*) '   METHOD :  MG or MG_normal : Variable V-cycle '
      end if
C
      write(*,*) '     Residual :        r := b-Au; '
      write(*,'(a,a,e12.4)')
     >     '      Convergence :     ||r||_0/||b||_0  and/or',
     >     '  ||u - soln||_0/||soln||_0 <', tol
      write(*,'(a,2(5X,i3))') '      Coarse and fine levels : ',lc,lf
      write(*,'(a,a)') ' -----------------------------',
     >           '--------------------------------------------------'
      write(*,*) 'Iteration  ||r||/||b||      ||r||  ',
     >                  ' ||u-soln||/||soln||   ||Br||   |Br||/||Bb||'
      write(*,'(a,a)') ' -----------------------------',
     >           '--------------------------------------------------'
C
C...  Compute the L^2 norm of the computed exact solution if it is given.
C
      if (iexact .eq. 2) then
         call l2norm(soln,soln_norm,n,1d0)
         write(*,*) ' soln_norm: ', soln_norm
         if (soln_norm .le. zero) then
            soln_norm = 1.d0
            write(*,*) 'Zero exact solution is detected!!!!.'
         end if
      end if
C
cc      write(*,*)
cc      write(*,*) (rhs(ii), ii = 1,n)
C
      call init_g0(iao,jao,ao,rhs,sol,n)
C
cc      write(*,*) (sol(ii), ii = 1,n)
cc      write(*,*)
cc      write(*,*) (rhs(ii), ii = 1,n)
C
      krhs_tmp  = lastr
      ksol_tmp  = krhs_tmp + n + 1
      lastr = ksol_tmp + n + 1
C
C...  To compute the initial residual: b = rhs - A*sol.
      call abyvam(rhs,iao,jao,ao,sol,n,r(krhs_tmp))
C
ccc      call scpro(rhs,rhs,err_res0,n)
      call scpro(r(krhs_tmp),r(krhs_tmp),err_res0,n)
      err_res0 = dsqrt(err_res0)
C
      if (err_res0 .le. zero) then
         err_res0 = 1.d0
         write(*,*) 'Zero right hand side is detected!!!!.'
      end if
C
      err_relres = err_res0
      err_rel = 0.d0
C
      kk = 0
      write(*,'(2X,i5,2X,3(3X,e11.4))') 
     >     kk,1.d0,err_res0,1.d0
C
C...  Iterate until convergence:  sol = sol + B*b = sol + B(rhs - A*sol).
C
      do 200 kk = 1, maxit
C
         ifree = lasti
         kfree = lastr
C
C...     Iterator B using an MG + Precondition_Normal, 
C...     i.e., u = B*b = B(rhs - A*sol).
C
         call msolve(n, r(krhs_tmp), r(ksol_tmp), nnz, 
     >        iao, jao, ao, 0, r, m)
C     
C...     New solution:  sol = sol + u = sol + P(rhs - A*sol).
         call uupluv(sol,r(ksol_tmp),n)
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(rhs,iao,jao,ao,sol,n,r(krhs_tmp))
C
C...     To compute the residual error.
         call scpro(r(krhs_tmp),r(krhs_tmp),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
C...     To compute the relative error if the exact solution is known.
         if (iexact .eq. 2) then
            call wuminv(r(ksol_tmp),soln,sol,n)
            call l2norm(r(ksol_tmp),xdiffn,n,1d0)
            err_rel = xdiffn/soln_norm  
         end if
C
         if (lf .lt. 3) then
            nmod = mod(kk,100)
         else if (lf .lt. 3) then
            nmod = mod(kk,10)
         else 
            nmod = 1 
         end if
C
         if (nmod .eq. 1) then
            write(*,'(2X,i5,2X,3(3X,e11.4))')
     >           kk,err_relres,err_res,err_rel
         end if
C
         if  (iexact .eq. 2) then
            if (err_rel .lt. tol) go to 333
         else
            if (err_relres .lt. tol .and. err_rel .lt. tol) go to 333
         end if
 200  continue
C
      write(*,*), ' Iteration limit exceeded:', maxit
 333  continue
      niter = min(kk, maxit)
      write(*,'(2X,i5,2X,3(3X,e11.4))')
     >     niter,err_relres,err_res,err_rel
C
      write(*,'(a,a)') ' =============================',
     >      '=================================================='
C
      arfac  = (err_relres)**(1./dble(kk))
      arfac1 = (err_rel)**(1./dble(kk))
      write(*,*) 
      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
C
      write(*,*)
C
 500  return
      end
C=====================================================================
      subroutine msolve(n, rhs, sol, nnz, iao, jao, ao, isym, r, m)
C=====================================================================
      implicit real*8(a-h,o-z)
C
      integer n, nnz, iao(1), jao(1), isym, m(1)
      double precision rhs(n), sol(1), ao(1), r(1)
C
      parameter (lmax = 20)
      integer lmax
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
C
      common /menu/ ich(200)
      common /prenormal/ ipre_normal
C--------------------------------------------------------------------
C...  MSOLVE solves a linear system P*sol = rhs for sol given rhs with
C...  the preconditioning matrix P (P is supplied via R and M arrays). 
C...  (As an example, it is used as an external subroutine for DGMRES.)
C...
C...  Parameters:
C...    IAO,JAO,AO - contain the matrix data structure for A (It could 
C...                 take any form.)
C...    N       - the number of unknowns
C...    RHS     - the right-hand side vector
C...    SOL     - the solution upon return
C...    NNZ     - the number of non-zeros in the SLAP IAO, JAO, AO storage 
C...              for the matrix A. (It could take any form.)
C...    ISYM    - a flag which, if non-zero, denotes that A is symmetric
C...              and only the lower or upper triangle is stored. If it
C...              is zero, all non-zero entries of the matrix are stored.
C...    R       - can be used to pass necessary preconditioning 
C...              information and/or workspace 
C...    M       - the same purpose as RWORK.
C--------------------------------------------------------------------
C
      lasti = ifree
      lastr = kfree
C
      call nullv(sol,n)
C
      if (ipre_normal .eq. 0) then
         call copyv(rhs,r(kb(lf)),n)
         go to 160
      else 
C...     Compute b = (rhs^t)*A = A^t*rhs.
         call vbya(iao,jao,ao,rhs,n,n,r(kb(lf)))
      end if
C
C...  Precondition by G-S sweeps for the normal equation: 
C...  A^t*A*sol = A^t*rhs.
C
      igs        = ich(42)
      max_sweeps = ich(43)
C     
      if (igs .eq. 1) then
         call fwd_gs_ns_wr(sol,m(ic),m(jc),r(kc),r(kb(lf)),n,
     >        max_sweeps,0)
      else if(igs .eq. 2) then
         call bwd_gs_ns_wr(sol,m(ic),m(jc),r(kc),r(kb(lf)),n,
     >        max_sweeps,0) 
      else if(igs .eq. 3) then
         call sym_gs_ns_wr(sol,m(ic),m(jc),r(kc),r(kb(lf)),n,
     >        max_sweeps,0)
      else
         write(*,*) ' Error in the choice of G-S! '
      end if
C     
C...  To compute the residual: b = rhs - A*sol.
      call abyvam(rhs,iao,jao,ao,sol,n,r(kb(lf)))
C     
 160  continue
C     
C...  Iterator B using an MG-cycle, i.e., u = B*b = B(rhs - A*sol).
C     
      call premg(m(1),r(1),lasti,lastr)
C     
C...  New solution:  sol = sol + u = sol + B(rhs - A*sol).
      if (ipre_normal .eq. 0) then
         call copyv(r(ku(lf)),sol,n)
      else
         call uupluv(sol,r(ku(lf)),n)
      end if
C
      return
      end
C=====================================================================
      subroutine matvec(n,x,y,nelt,ia,ja,a,isym)
C=====================================================================
      integer n, nelt, ia(1), ja(1), isym
      double precision x(n), y(n), a(1)
cc      integer i, iaa, iab, k
cc      double precision u
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL SPARSE MATRIX BY VECTOR: y = A*x.
C...  MATVEC is used as an external subroutine for DGMRES subroutine.
C...
C...  Parameters:
C...    IA,JA,A - contain the matrix data structure for A (It could 
C...                take any form.)
C...    Y       - the product A*X
C...    X       - an input vector
C...    NELT    - the number of non-zeros in the SLAP IA, JA, A 
C...              storage for the matrix A. 
C...    ISYM    - a flag which, if non-zero, denotes that A is symmetric
C...              and only the lower or upper triangle is stored. If it
C...              is zero, all non-zero entries of the matrix are stored.
C--------------------------------------------------------------------
      call abyvg(ia,ja,a,x,n,y)
C
cc      do i = 1, n
cc        u = 0.0d00
cc         iaa = ia(i)
cc         iab = ia(i+1) - 1
cc         if(iab .ge. iaa) then
cc            do k = iaa, iab
cc               u = u + a(k) * x(ja(k))
cc            end do
cc         end if
cc         y(i) = u
cc      end do
C
      return
      end
C=====================================================================
