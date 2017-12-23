C=====================================================================
      subroutine form_submat(m,iao,jao,r,ao,n,nnz) 
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer m(1)
      double precision r(1)
C
      integer n, nnz, iao(1),jao(1)
      double precision ao(1)
C
      integer mt, jdf, lread, istiff
      integer lasti, lastr, iend, kend
      integer nel, ned, ie, je, inedg, iet, jet, i
      integer kx, ky, idir, iip, iaw, jaw, isub, nblk
      integer iwork1, iwork2, iwork3
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
      character*100 meshf(0:lmax)
      common /fname/ meshf
      integer ich
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  Generate preconditioning submatrices and other information. 
C...  The following finite element method can be applied: 
C...  Standard Galerkin, EAFE (Edge Average Finite Element).
C...
C...  Parameter:
C...    IAO,JAO,AO - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    NNZ     - the number of non-zeros in the IAO, JAO, AO storage 
C...              for the matrix A.
C...    M, R    - integer and real*8 spaces for storing all necessary
C...              information
C...    lmax    - maximal level for multi-level structure
C...    nd(k)   - the number of unknowns or the order of the Matrix for
C...              each level k on a multi-level structure
c...    n1(k),n2(k) - the dimensions in x/y direction of the uniform
C...              mesh for each level k on a multi-level structure
C...              They are needed for defining the prolongation and/or
C...              restriction matrices.
C...    IA(k),JA(k),KA(k) - pointers for controlling the data structure 
C...              of the matrix A_k in RRCU form for each level k.
C...              Numerical values of nonzeros of A_k in RRCU are stored
C...              in an array starting with r(ka(k)).
C...    NZ(k)   - the number of non-zeros in the matrix data structure,
C...              IA(k),JA(k),KA(k). More precisely, it will be the size
C...              of the storage related to the pointer JA(k) or KA(k)
C...    ipp(k),jpp(k),kpp(k)- pointers for controlling the prolongation
C...              matrix data structure, P_k. Because MG cycle is done 
C...              on the uniform mesh, it is not necessary to supply the
C...              numerical values of nonzeros of the prolongation matrix.
C...              That is, kpp(k) is unnecessary in this subroutine, instead
C...              only the actions are performed in the subroutines, 
C...              RVBYP and PBYVCS, using the information in ipp(k),jpp(k).
C...              They have the same data structure as IA(k),JA(k).
C...    nzp(k)  - the number of non-zeros in the matrix data structure,
C...              IPP(k), JPP(k). More precisely, it will be the size
C...              of the storage related to the pointer JPP(k).
C...    ku(k)   - working array (pointer) for the solution at level k
C...    kb(k)   - working array (pointer) for the load vector at level k
C...    iord(k) - a different ordering of unknowns is stored for G-S
C...              smoothing steps if ICH_ORD = 1. Please note that if 
C...              ICH_ORD is not 1, then simply assign m(iord(k)) = 0 
C...              with array size being 1 for each level k before calling
C...              this subroutine.
C...    ifree, kfree - pointers for free spaces for M and R 
C...    ICH_ORD - a different ordering of unknowns is used for G-S
C...              smoothing steps if it is 1. Please note that if it 
C...              is not 1, then simply assign m(iord(k)) = 0 with
C...              array size being 1 for each level k before calling
C...              this subroutine.
C...    LF      - level of the finest mesh: 4,5,6 etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    ISTIFF  - choice of the stiffness matrix (iterator B) in MG.
C...              1 = Standard Galerkin,    2 = EAFE, 
C...              3 = Streamline diffusion, 4 = Upwind, 5 = SUPG,
C...    JDF     - the choice of finite elements, i.e., degrees of freedom
C...              3 = linear element, 6 = quadratic element,
C...              4 = bilinear element (currently, it is 3 only)
C...    LREAD   - format of reading mesh file: 
C...              1 = formatted, 2 = unformatted
C---------------------------------------------------------------------
      mt     = 15
      jdf    = ich(3)
      lread  = ich(5)
      istiff = ich(40)
C
C...  Generating the stiffness matrix A and the load vector B at each level.
C
      lasti = ifree
      lastr = kfree
C
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
cc         write(*,*)
cc     >        ' No. of elements, nodes, & bdry edges:',nel,nd(i),ned
C
C...     Compute the addresses of each array in the arrays m and r.
C
         if (i .eq. lf) then
            ipp(i)  = lasti
         else
            ia(i)   = lasti
            ja(i)   = ia(i) + nd(i) + 1
            ipp(i)  = ja(i) + 2*(nel + nd(i) + 5) + nd(i)
         end if
C
         jpp(i)  = ipp(i) + nd(i) + 1
         iord(i) = jpp(i) + 2*nd(i)
         if (ich_ord .eq. 1) then
            lasti   = iord(i) + nd(i) + 1
         else
            lasti   = iord(i) + 1
         end if
C
         ie      = lasti
         je      = ie + nel + 1
         inedg   = je + nel*3
         iet     = inedg + ned * 2
         jet     = iet + nd(i) + 1
         iend    = jet + nel*3
C
         ku(i)   = lastr
         kb(i)   = ku(i) + nd(i)
         ka(i)   = kb(i) + nd(i)
         lastr   = ka(i) + 2*(nel + nd(i) + 5) + nd(i)
C
         kx      = lastr
         ky      = kx + nd(i)
         kend    = ky + nd(i)
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
         idir = iend
         call rdmesh(lread,mt,m(ie),m(je),ned,m(inedg),
     >        m(idir),r(kx),r(ky),nel,nd(i),jdf)
         close(mt)
C
         call iit(m(ie),m(je),nel,nd(i),m(iet),m(jet))
C
         if (i .lt. lf) then
            nz(i) = 0
            call smbasg(m(ie),m(je),m(iet),m(jet),m(ia(i)),m(ja(i)),
     >           nd(i),nz(i),m(idir))
         end if
C
         iip  = iet
         call pde_choice(istiff,
     I        nel,nd(i),jdf,m(ie),m(je),r(kx),r(ky),
     O        m(ia(i)),m(ja(i)),r(ka(i)),nz(i),r(kb(i)),
     W        ned,m(inedg),m(idir),m(iip))
C
cc         call outmat1(m(ia(i)),m(ja(i)),r(ka(i)),nd(i),nz(i))
C
C...     To determine whether Tarjan's ordering is considered or not.
C...     If ordering is considered, then set ICH_ORD = 1.
C
         if (ich_ord .eq. 1) then
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
            m(iord(i)) = 0
         end if
 50   continue
C
      ifree = lasti
      kfree = lastr
C
 500  return
      end
C=====================================================================
      subroutine mg(m,iao,jao,r,sol,rhs,ao,
     >     n,nnz,maxit,tol,iexact,msolve) 
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer m(1)
      double precision r(1)
C
      integer n,nnz,maxit,iexact,iao(1),jao(1)
      double precision sol(1),rhs(1),ao(1),tol
C
      integer kk, nmod, niter
      integer lasti, lastr
      double precision zero, soln_norm, xdiffn, arfac, arfac1
      double precision err_res0, err_relres, err_rel, err_res
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
      external msolve
C
      double precision soln(1)
      COMMON /DSLBLK/ SOLN
C
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*SOL = RHS by the
C...  standard linear iterative scheme
C...
C...       u_{i+1} = u_i + B(rhs - A*u_i)
C...
C...  The action of the iterator or preconditioner B is supplied by the 
C...  external subroutine MSOLVE (see below for its detailed description). 
C...  Various MG cycles (V-cycle, \-cycle, variable V-cycle) will be
C...  applied in the subroutine MSOLVE. This subroutine does not 
C...  initialize the solution (SOL), So the solution, SOL, should be 
C...  initialized before calling this subroutine.
C...
C...  Parameter:
C...    IAO,JAO,AO - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    RHS     - the right-hand side vector
C...    SOL     - the solution upon return
C...    NNZ     - the number of non-zeros in the IAO, JAO, AO storage 
C...              for the matrix A.
C...    ISYM    - a flag which, if non-zero, denotes that A is symmetric
C...              and only the lower or upper triangle is stored. If it
C...              is zero, all non-zero entries of the matrix are stored.
C...    MSOLVE  - External subroutine. Name of a routine which solves a 
C...              linear system MZ = S for Z given S with the 
C...              preconditioning matrix M (M is supplied via RWORK and
C...              IWORK arrays). The name of the MSOLVE routine must be
C...              declared external in the calling program. The calling
C...              sequence to MSOLVE is:
C...
C...              CALL MSOLVE(N, S, Z, NNZ, IAO, JAO, AO, ISYM, RWORK, IWORK)
C...
C...              Where N is the number of unknowns, S is the right-hand
C...              side vector and Z is the solution upon return. NNZ, IAO,
C...              JAO, AO and ISYM are defined as above. RWORK is a double 
C...              precision array that can be used to pass necessary 
C...              preconditioning information and/or workspace to MSOLVE.
C...              IWORK is an integer work array for the same purpose as 
C...              RWORK.
C...    MAXIT   - maximum number of MG iterations
C...    TOL     - stopping criterion for iteration procedure.
C...    IEXACT  - to determine whether reading or writing the computed
C...              exact solution is considered or not;
C...              0 - no reading and writing the computed (exact) soln;
C...              1 - writing the computed soln;
C...              2 - reading the computed exact soln;
C...              3 - writing the computed exact soln;
C...              If IEXACT = 0,1,3, iteration stops when 
C...                         ||r||/||rhs|| < TOL,
C...              where r = rhs - A*sol.
C...              IEXACT = 2 is often useful for checking and comparing 
C...              different routines. For this case, the user must supply
C...              the "exact" solution or a very accurate approximation 
C...              (one with an error much less than TOL) through a common 
C...              block,  
C...                         COMMON /DSLBLK/ SOLN( )
C...              If IEXACT = 2, iteration stops when 
C...                         ||sol-soln||/||soln|| < tol.
C...              Note that this requires the user to set up the 
C...              "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
C...              The routine with this declaration should be loaded before
C...              the stop test so that the correct length is used by the 
C...              loader.
C...    M, R    - integer and real*8 spaces for storing all necessary
C...              information
C...    lmax    - maximal level for multi-level structure
C...    nd(k)   - the number of unknowns or the order of the Matrix for
C...              each level k on a multi-level structure
c...    n1(k),n2(k) - the dimensions in x/y direction of the uniform
C...              mesh for each level k on a multi-level structure
C...              They are needed for defining the prolongation and/or
C...              restriction matrices.
C...    IA(k),JA(k),KA(k) - pointers for controlling the data structure 
C...              of the matrix A_k in RRCU form for each level k.
C...              Numerical values of nonzeros of A_k in RRCU are stored
C...              in an array starting with r(ka(k)).
C...    NZ(k)   - the number of non-zeros in the matrix data structure,
C...              IA(k),JA(k),KA(k). More precisely, it will be the size
C...              of the storage related to the pointer JA(k) or KA(k)
C...    ipp(k),jpp(k),kpp(k)- pointers for controlling the prolongation
C...              matrix data structure, P_k. Because MG cycle is done 
C...              on the uniform mesh, it is not necessary to supply the
C...              numerical values of nonzeros of the prolongation matrix.
C...              That is, kpp(k) is unnecessary in this subroutine, instead
C...              only the actions are performed in the subroutines, 
C...              RVBYP and PBYVCS, using the information in ipp(k),jpp(k).
C...              They have the same data structure as IA(k),JA(k).
C...    nzp(k)  - the number of non-zeros in the matrix data structure,
C...              IPP(k), JPP(k). More precisely, it will be the size
C...              of the storage related to the pointer JPP(k).
C...    ku(k)   - working array (pointer) for the solution at level k
C...    kb(k)   - working array (pointer) for the load vector at level k
C...    iord(k) - a different ordering of unknowns is stored for G-S
C...              smoothing steps if ICH_ORD = 1. Please note that if 
C...              ICH_ORD is not 1, then simply assign m(iord(k)) = 0 
C...              with array size being 1 for each level k before calling
C...              this subroutine.
C...    ifree, kfree - pointers for free spaces for M and R 
C...    IPRESM  - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPOSSM  - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRESM  - number of pre-smoothings: 1,2,3, etc
C...    NPOSSM  - number of post-smoothings: 1,2,3, etc
C...    LCYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle
C...    LSOLV_C - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...    ICH_ORD - a different ordering of unknowns is used for G-S
C...              smoothing steps if it is 1. Please note that if it 
C...              is not 1, then simply assign m(iord(k)) = 0 with
C...              array size being 1 for each level k before calling
C...              this subroutine.
C...    LF      - level of the finest mesh: 4,5,6 etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C---------------------------------------------------------------------
      zero   = 1.d-15
C
      lasti = ifree
      lastr = kfree
C
C...  ------------------- MG starts here. ---------------------------
C
      write(*,*)
      write(*,'(a,a)') ' =============================',
     >      '=================================================='
      if (lcycle .eq. 1) then
         write(*,*) '   METHOD :  MG  V-cycle '
      else if (lcycle .eq. 2) then
         write(*,*) '   METHOD :  MG  Backslash cycle '
      else if (lcycle .eq. 3) then
         write(*,*) '   METHOD :  MG  Variable V-cycle '
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
C...  To compute the initial residual: b = rhs - A*sol.
      call abyvam(rhs,iao,jao,ao,sol,n,r(kb(lf)))
C
      call scpro(r(kb(lf)),r(kb(lf)),err_res0,n)
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
C...     Iterator B using an MG cycle: either msolve or pre_mg. 
C...     i.e., u = B*b = B(rhs - A*sol).
C
         call msolve(n, r(kb(lf)), r(ku(lf)), nnz, 
     >        iao, jao, ao, 0, r, m)
C
cc         call pre_mg_standard(m,r)
C
C...     New solution:  sol = sol + u = sol + P(rhs - A*sol).
         call uupluv(sol,r(ku(lf)),n)
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(rhs,iao,jao,ao,sol,n,r(kb(lf)))
C
C...     To compute the residual error.
         call scpro(r(kb(lf)),r(kb(lf)),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
C...     To compute the relative error if the exact solution is known.
         if (iexact .eq. 2) then
            call wuminv(r(ku(lf)),soln,sol,n)
            call l2norm(r(ku(lf)),xdiffn,n,1d0)
            err_rel = xdiffn/soln_norm  
         end if
C
         nmod = 1 
C
         if (nmod .eq. 1) then
            write(*,'(2X,i5,2X,3(3X,e11.4))')
     >           kk,err_relres,err_res,err_rel
         end if
C
         if (iexact .eq. 2) then
            if (err_rel .lt. tol) go to 333
         else
            if (err_relres .lt. tol) go to 333
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
      subroutine pre_mg_standard(m,r)
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer m(1)
      double precision r(1)
C
      integer i, icount, k, kwk1, klast, noofsm, itmax, maxa, kau
      integer iak, jak, iac, kuk, kbc, kbk, kak, ndc, ndk
      integer ipk, jpk, nzpk, iordk, kuc, jac, kac, iordc
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
C---------------------------------------------------------------------
C...  PRE_MG_STANDARD stands for a standard multigrid cycle as a 
C...  preconditioner and it will be done on the structured uniform grid 
C...  with lexicographical ordering of unknowns. This action corresponds 
C...  to the iterator or preconditioner B in an MG-cycle. Various multigrid
C...  cycles can be used: V-cycle, \-cycle, variable V-cycle. 
C...
C...  Parameter:
C...    M, R    - integer and real*8 spaces for storing all necessary
C...              information
C...    lmax    - maximal level for multi-level structure
C...    nd(k)   - the number of unknowns or the order of the Matrix for
C...              each level k on a multi-level structure
c...    n1(k),n2(k) - the dimensions in x/y direction of the uniform
C...              mesh for each level k on a multi-level structure
C...              They are needed for defining the prolongation and/or
C...              restriction matrices.
C...    IA(k),JA(k),KA(k) - pointers for controlling the data structure 
C...              of the matrix A_k in RRCU form for each level k.
C...              Numerical values of nonzeros of A_k in RRCU are stored
C...              in an array starting with r(ka(k)).
C...    NZ(k)   - the number of non-zeros in the matrix data structure,
C...              IA(k),JA(k),KA(k). More precisely, it will be the size
C...              of the storage related to the pointer JA(k) or KA(k)
C...    ipp(k),jpp(k),kpp(k)- pointers for controlling the prolongation
C...              matrix data structure, P_k. Because MG cycle is done 
C...              on the uniform mesh, it is not necessary to supply the
C...              numerical values of nonzeros of the prolongation matrix.
C...              That is, kpp(k) is unnecessary in this subroutine, instead
C...              only the actions are performed in the subroutines, 
C...              RVBYP and PBYVCS, using the information in ipp(k),jpp(k).
C...              They have the same data structure as IA(k),JA(k).
C...    nzp(k)  - the number of non-zeros in the matrix data structure,
C...              IPP(k), JPP(k). More precisely, it will be the size
C...              of the storage related to the pointer JPP(k).
C...    ku(k)   - working array (pointer) for the solution at level k
C...    kb(k)   - working array (pointer) for the load vector at level k
C...    iord(k) - a different ordering of unknowns is stored for G-S
C...              smoothing steps if ICH_ORD = 1. Please note that if 
C...              ICH_ORD is not 1, then simply assign m(iord(k)) = 0 
C...              with array size being 1 for each level k before calling
C...              this subroutine.
C...    ifree, kfree - pointers for free spaces for M and R 
C...    IPRESM  - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPOSSM  - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRESM  - number of pre-smoothings: 1,2,3, etc
C...    NPOSSM  - number of post-smoothings: 1,2,3, etc
C...    LCYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle
C...    LSOLV_C - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...    ICH_ORD - a different ordering of unknowns is used for G-S
C...              smoothing steps if it is 1. Please note that if it 
C...              is not 1, then simply assign m(iord(k)) = 0 with
C...              array size being 1 for each level k before calling
C...              this subroutine.
C...    LF      - level of the finest mesh: 4,5,6 etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C---------------------------------------------------------------------
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
         iac  = ia(k-1)
         kuk  = ku(k)
         kbc  = kb(k-1)
         kbk  = kb(k)
         kak  = ka(k)
         ndc  = nd(k-1)
         ndk  = nd(k)
         ipk  = ipp(k)
         jpk  = jpp(k)
         nzpk  = nzp(k)
         iordk = iord(k)
C
         kwk1  = kfree
         klast  = kwk1 + ndk + 1
C
C...     Presmoothing by Gauss-seidel smoother npresm times: u=u+R(b-Au)
C
         if (lcycle .eq. 3) then !For variable V-cyle
            icount = icount + 1
            noofsm = npresm*2**icount
         else
            noofsm = npresm      !For V or \-cyle
         end if
C
	 if (ipresm .eq. 1) then 
            call fwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipresm .eq. 2) then 
            call bwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipresm .eq. 3) then 
            call sym_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
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
         call zero_bdry(r(kbc),m(iac),ndc)
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
      if (lsolv_c .eq. 1) then
         call fwd_gs_ord(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv_c .eq. 2) then
         call bwd_gs_ord(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv_c .eq. 3) then
         call sym_gs_ord(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv_c .eq. 4) then
         maxa = ifree
         kau  = kfree
         call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
         call copyv(r(kbc),r(kuc),ndc)
      end if
C     
C...  Going upward in the V-cycle.
C
      do 200 k = lc+1, lf
         iak = ia(k)
         jak = ja(k)
         iac = ia(k-1)
         kuc = ku(k-1)
         kuk = ku(k)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         ipk = ipp(k)
         jpk = jpp(k)
         nzpk  = nzp(k)
         iordk = iord(k)
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C
         call pbyvcs(m(ipk),m(jpk),r(kuc),ndc,r(kuk),ndk)
C
C...     Postsmoothing by Gauss-Seidel smoother npossm times: u=u+R(b-Au)
C
         if (lcycle .eq. 2)  go to 200        !No post smoothing in \-cycle
C
         if (lcycle .eq. 3) then !For variable V-cyle
            noofsm = npossm*2**icount
            icount = icount - 1
         else
            noofsm = npossm      !For V or \-cyle
         end if
C
	 if (ipossm .eq. 1) then 
            call fwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipossm .eq. 2) then 
            call bwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipossm .eq. 3) then 
            call sym_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
         end if
C
 200  continue 
C
      return
      end
C=====================================================================
      subroutine pre_mg_non_standard(m,r)
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer m(1)
      double precision r(1)
C
      integer i, icount, k, kwk1, klast, noofsm, itmax, maxa, kau
      integer iak, jak, iac, kuk, kbc, kbk, kak, ndc, ndk
      integer ipk, jpk, nzpk, iordk, kuc, jac, kac, iordc, kpk
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
C---------------------------------------------------------------------
C...  PRE_MG_NON_STANDARD stands for a multigrid cycle as a  
C...  preconditioner and it will be done on any grid with an arbitrary
C...  ordering of unknowns. This action corresponds to the iterator or 
C...  preconditioner B in an MG-cycle. Various multigrid cycles can be
C...  used: V-cycle, \-cycle, variable V-cycle. The information about 
C...  prolongation matrix should be stored using by the pointers, 
C...  ipp(k),jpp(k),kpp(k) for each level k.
C...
C...  Parameter:
C...    M, R    - integer and real*8 spaces for storing all necessary
C...              information
C...    lmax    - maximal level for multi-level structure
C...    nd(k)   - the number of unknowns or the order of the Matrix for
C...              each level k on a multi-level structure
c...    n1(k),n2(k) - the dimensions in x/y direction of the uniform
C...              mesh for each level k on a multi-level structure
C...              They are needed for defining the prolongation and/or
C...              restriction matrices.
C...    IA(k),JA(k),KA(k) - pointers for controlling the data structure 
C...              of the matrix A_k in RRCU form for each level k.
C...              Numerical values of nonzeros of A_k in RRCU are stored
C...              in an array starting with r(ka(k)).
C...    NZ(k)   - the number of non-zeros in the matrix data structure,
C...              IA(k),JA(k),KA(k). More precisely, it will be the size
C...              of the storage related to the pointer JA(k) or KA(k)
C...    ipp(k),jpp(k),kpp(k) - pointers for the data structure of the 
C...              prolongation matrix, P_k. They have the same data
C...              structure as IA(k),JA(k),KA(k). Because MG cycle is 
C...              done on the arbitrary mesh, we need all the information
C...              about the prolongation matrix.
C...    nzp(k)  - the number of non-zeros in the matrix data structure,
C...              IPP(k), JPP(k). More precisely, it will be the size
C...              of the storage related to the pointer JPP(k).
C...    ku(k)   - working array (pointer) for the solution at level k
C...    kb(k)   - working array (pointer) for the load vector at level k
C...    iord(k) - a different ordering of unknowns is stored for G-S
C...              smoothing steps if ICH_ORD = 1. Please note that if 
C...              ICH_ORD is not 1, then simply assign m(iord(k)) = 0 
C...              with array size being 1 for each level k before calling
C...              this subroutine.
C...    ifree, kfree - pointers for free spaces for M and R 
C...    IPRESM  - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPOSSM  - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRESM  - number of pre-smoothings: 1,2,3, etc
C...    NPOSSM  - number of post-smoothings: 1,2,3, etc
C...    LCYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle
C...    LSOLV_C - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...    ICH_ORD - a different ordering of unknowns is used for G-S
C...              smoothing steps if it is 1. Please note that if it 
C...              is not 1, then simply assign m(iord(k)) = 0 with
C...              array size being 1 for each level k before calling
C...              this subroutine.
C...    LF      - level of the finest mesh: 4,5,6 etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C---------------------------------------------------------------------
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
         iac  = ia(k-1)
         kuk  = ku(k)
         kbc  = kb(k-1)
         kbk  = kb(k)
         kak  = ka(k)
         ndc  = nd(k-1)
         ndk  = nd(k)
         ipk  = ipp(k)
         jpk  = jpp(k)
         kpk  = kpp(k)
         nzpk  = nzp(k)
         iordk = iord(k)
C
         kwk1  = kfree
         klast  = kwk1 + ndk + 1
C
C...     Presmoothing by Gauss-seidel smoother npresm times: u=u+R(b-Au)
C
         if (lcycle .eq. 3) then !For variable V-cyle
            icount = icount + 1
            noofsm = npresm*2**icount
         else
            noofsm = npresm      !For V or \-cyle
         end if
C
	 if (ipresm .eq. 1) then 
            call fwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipresm .eq. 2) then 
            call bwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipresm .eq. 3) then 
            call sym_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
         end if
C
C...     To compute the residual wk = b - Au.
C
         call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
C
C...     Restriction to the lower level: b = P^t*wk = wk^t*P.
C
         call vbya(m(ipk),m(jpk),r(kpk),r(kwk1),ndk,ndc,r(kbc))
         call zero_bdry(r(kbc),m(iac),ndc)
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
      if (lsolv_c .eq. 1) then
         call fwd_gs_ord(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv_c .eq. 2) then
         call bwd_gs_ord(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv_c .eq. 3) then
         call sym_gs_ord(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        m(iordc),ndc,itmax,0)
      else if (lsolv_c .eq. 4) then
         maxa = ifree
         kau  = kfree
         call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
         call copyv(r(kbc),r(kuc),ndc)
      end if
C     
C...  Going upward in the V-cycle.
C
      do 200 k = lc+1, lf
         iak = ia(k)
         jak = ja(k)
         iac = ia(k-1)
         kuc = ku(k-1)
         kuk = ku(k)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         ipk = ipp(k)
         jpk = jpp(k)
         kpk = kpp(k)
         nzpk  = nzp(k)
         iordk = iord(k)
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C
         call abyvcs(m(ipk),m(jpk),r(kpk),r(kuc),ndk,r(kuk))
cc         call zero_bdry(r(kuk),m(iac),ndc)
C
C...     Postsmoothing by Gauss-Seidel smoother npossm times: u=u+R(b-Au)
C
         if (lcycle .eq. 2)  go to 200        !No post smoothing in \-cycle
C
         if (lcycle .eq. 3) then !For variable V-cyle
            noofsm = npossm*2**icount
            icount = icount - 1
         else
            noofsm = npossm      !For V or \-cyle
         end if
C
	 if (ipossm .eq. 1) then 
            call fwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipossm .eq. 2) then 
            call bwd_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
	 else if (ipossm .eq. 3) then 
            call sym_gs_ord(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >                     m(iordk),ndk,noofsm,0) 
         end if
C
 200  continue 
C
      return
      end
C=====================================================================
      subroutine msolve(n, rhs, sol, nnz, iao, jao, ao, isym, r, m)
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer n, nnz, iao(1), jao(1), isym, m(1)
      double precision rhs(n), sol(1), ao(1), r(1)
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
C--------------------------------------------------------------------
C...  MSOLVE solves a linear system P*sol = rhs for sol given rhs with
C...  the preconditioning matrix P (P is supplied via R and M arrays). 
C...  (As an example, it is used as an external subroutine for DGMRES,
C...  DCG, and MG.) Various MG cycles can be used here for the call
C...  of PRE_MG_STANDARD.
C...
C...  Parameters:
C...    IAO,JAO,AO - contain the matrix data structure for A (It could 
C...                 take any form.)
C...    N       - the number of unknowns or the order of the Matrix
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
      call nullv(sol,n)
C
      call copyv(rhs,r(kb(lf)),n)
C     
C...  Iterator B using an MG-cycle, i.e., u = B*b = B(rhs - A*sol).
C     
      call pre_mg_standard(m(1),r(1))
C     
C...  New solution:  sol = sol + u = sol + B(rhs - A*sol).
      call copyv(r(ku(lf)),sol,n)
C
      return
      end
C=====================================================================
      subroutine matvec(n,x,y,nnz,ia,ja,a,isym)
C=====================================================================
      implicit none
C
      integer n, nnz, ia(1), ja(1), isym
      double precision x(n), y(n), a(1)
C
cc      integer i, iaa, iab, k
cc      double precision u
C--------------------------------------------------------------------
C...  PRODUCT - GENERAL SPARSE MATRIX BY VECTOR: y = A*x.
C...  MATVEC is used as an external subroutine for some subroutines,
C...  such as DGMRES and DCG.
C...
C...  Parameters:
C...    IA,JA,A - contain the matrix data structure for A (It could 
C...              take any form.)
C...    Y       - output : the product A*X
C...    X       - an input vector
C...    N       - the number of unknowns or the order of the Matrix
C...    NNZ     - the number of non-zeros in the SLAP IAO, JAO, AO storage 
C...              for the matrix A. (It could take any form.)
C...    ISYM    - a flag which, if non-zero, denotes that A is symmetric
C...              and only the lower or upper triangle is stored. If it
C...              is zero, all non-zero entries of the matrix are stored.
C--------------------------------------------------------------------
      call abyvg(ia,ja,a,x,n,y)
C
C...  The following is same as the subroutine call, ABYVG.
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

