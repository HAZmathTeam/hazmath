C=====================================================================
      subroutine msolve(n, rhs, sol, nnz, iao, jao, ao, isym, r, m)
C=====================================================================
      integer n, nnz, iao(1), jao(1), isym, m(1)
      double precision rhs(n), sol(1), ao(1), r(1)
      integer i, iaa, iab, k
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
C
C--------------------------------------------------------------------
C...  MSOLVE solves a linear system P*sol = rhs for sol given rhs with
C...  the preconditioning matrix P (P is supplied via R and M arrays). 
C...  It is used as an external subroutine for DGMRES 
C...  subroutine.
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
         copyv(rhs,r(kb(lf)),n)
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
      call uupluv(sol,r(ku(lf)),n)
C
      return
      end
C=====================================================================
